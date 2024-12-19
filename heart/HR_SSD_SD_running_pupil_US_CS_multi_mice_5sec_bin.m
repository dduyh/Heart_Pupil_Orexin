clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

folder = {'m483_L_VTA_R_NAc_oxLight\Apr_10_2024';
    'm485_L_VTA_R_NAc_oxLight\Apr_10_2024';
    'm486_L_NAc_nLightR\Apr_10_2024';
    'm487_L_VTA_nLightR\Apr_10_2024';
    'm1772\Apr_10_2024';
    'm1773\Apr_10_2024';
    'm1774\Apr_10_2024'};

Sample_Rate = 200;    % 200 scans per second.
FrameRate = 20;

trace_duration = 600;   % 600 seconds.

trials_num = size(folder,1);

run = NaN(trials_num,trace_duration/5);
HR = NaN(trials_num,trace_duration/5);    % 15 seconds per trial.
RMSSD = NaN(trials_num,trace_duration/5);
SD = NaN(trials_num,trace_duration/5);
% pupil = NaN(trials_num,trace_duration*FrameRate);
pupil = NaN(trials_num,trace_duration/5);

bin_xlims = 1:5:trace_duration;

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    % Load data

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    % Analyse ECG signals

    CS_start = ceil(step_timepoint(4)*Sample_Rate)+1;
    trace_start = CS_start-360*Sample_Rate;
    trace_end = trace_start+trace_duration*Sample_Rate-1;

    % Analyse running signals

    running = datas(trace_start:trace_end,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed = [Abs_speedDeg(1) Abs_speedDeg];

    for i = 1 : 5 : trace_duration

        speed_median(round(i/5)+1) = median(speed(((i-1)*Sample_Rate+1):(i+4)*Sample_Rate));

    end

    run(I,:) = speed_median;

    % Raw ECG signals

    ECG_raw = datas(trace_start:trace_end,2)';

    % Remove baseline wandering

    fpass=[11 13];

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.007;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm,1:trace_duration*Sample_Rate,'spline','extrap');

    for i = 1 : 5 : trace_duration

        heartRate_bpm_interp_median(round(i/5)+1) = median(heartRate_bpm_interp(((i-1)*Sample_Rate+1):(i+4)*Sample_Rate));

    end

    HR(I,:) = heartRate_bpm_interp_median;

    % Standard Deviation

    for i = 1 : 5 : trace_duration

        heartbeat_index = find(pksLocs(2:end)>=(i-1) & pksLocs(2:end)<i+4);
        RR_intervals_trace = RR_intervals(heartbeat_index);
        RR_intervals_std(round(i/5)+1) = std(RR_intervals_trace);

    end

    SD(I,:) = RR_intervals_std;

    % Square of Successive Differences

    SSD = abs(diff(RR_intervals));
    SSD_interp = interp1(pksLocs(3:end)*Sample_Rate,SSD,1:trace_duration*Sample_Rate,'spline','extrap');

    for i = 1 : 5 : trace_duration

        SSD_median(round(i/5)+1) = median(SSD_interp(((i-1)*Sample_Rate+1):(i+4)*Sample_Rate));

    end

    RMSSD(I,:) = SSD_median;

    % pupil size

    pupil_start = round((step_timepoint(4)-step_timepoint(1))*FrameRate)-360*FrameRate+1;
    pupil_end = pupil_start+trace_duration*FrameRate-1;

    raw_pupil = areas(pupil_start:pupil_end);

    for i = 1 : 5 : trace_duration

        medianArea(round(i/5)+1) = median(raw_pupil(((i-1)*FrameRate+1):(i+4)*FrameRate));

    end

    pupil_zscored = (medianArea - mean(medianArea)) / std(medianArea);
    pupil(I,:) = pupil_zscored;

    figure(1);
    subplot(5,1,1);
    hold on
    plot(bin_xlims,speed_median,'LineWidth',1)
    subplot(5,1,2)
    hold on
    plot(bin_xlims,pupil_zscored,'LineWidth',1)
    subplot(5,1,3)
    hold on
    plot(bin_xlims,heartRate_bpm_interp_median,'LineWidth',1)
    subplot(5,1,4)
    hold on
    plot(bin_xlims,SSD_median,'LineWidth',1)
    subplot(5,1,5)
    hold on
    plot(bin_xlims,RR_intervals_std,'LineWidth',1)

end

%%
figure(1);


subplot(5,1,1);
hold on
plot(bin_xlims,mean(run),'Color','k','LineWidth',2)
line([120,120],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([240,240],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([360,360],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
line([480,480],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
xlim([0 bin_xlims(end)])
axis off
title('Running Speed','FontSize',18,'FontWeight','bold')
hold off

subplot(5,1,2)
hold on
plot(bin_xlims,mean(pupil),'Color','k','LineWidth',2)
line([120,120],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([240,240],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([360,360],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
line([480,480],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
xlim([0 bin_xlims(end)])
% ylim([-2 4])
axis off
title('Pupil Size','FontSize',18,'FontWeight','bold','color',[229 114 190]./255)
hold off

subplot(5,1,3)
hold on
plot(bin_xlims,mean(HR),'Color','k','LineWidth',2)
line([120,120],[min(mean(HR)-std(HR)/sqrt(trials_num)),max(mean(HR)+std(HR)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([240,240],[min(mean(HR)-std(HR)/sqrt(trials_num)),max(mean(HR)+std(HR)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([360,360],[min(mean(HR)-std(HR)/sqrt(trials_num)),max(mean(HR)+std(HR)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
line([480,480],[min(mean(HR)-std(HR)/sqrt(trials_num)),max(mean(HR)+std(HR)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
xlim([0 bin_xlims(end)])
% ylim([11.3 12.7])
title('Heart Rate','FontSize',20,'FontWeight','bold','color',[255 128 128]./255)
ylabel('HR (bpm)','FontSize',15,'FontWeight','bold');
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

subplot(5,1,4)
hold on
plot(bin_xlims,mean(RMSSD),'Color','k','LineWidth',2)
line([120,120],[min(mean(RMSSD)-std(RMSSD)/sqrt(trials_num)),max(mean(RMSSD)+std(RMSSD)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([240,240],[min(mean(RMSSD)-std(RMSSD)/sqrt(trials_num)),max(mean(RMSSD)+std(RMSSD)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([360,360],[min(mean(RMSSD)-std(RMSSD)/sqrt(trials_num)),max(mean(RMSSD)+std(RMSSD)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
line([480,480],[min(mean(RMSSD)-std(RMSSD)/sqrt(trials_num)),max(mean(RMSSD)+std(RMSSD)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
xlim([0 bin_xlims(end)])
% ylim([0.002 0.006])
title('Root Mean Square of Successive Differences','FontSize',20,'FontWeight','bold','color','k')
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

subplot(5,1,5)
hold on
plot(bin_xlims,mean(SD),'Color','k','LineWidth',2)
line([120,120],[min(mean(SD)-std(SD)/sqrt(trials_num)),max(mean(SD)+std(SD)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([240,240],[min(mean(SD)-std(SD)/sqrt(trials_num)),max(mean(SD)+std(SD)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([360,360],[min(mean(SD)-std(SD)/sqrt(trials_num)),max(mean(SD)+std(SD)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
line([480,480],[min(mean(SD)-std(SD)/sqrt(trials_num)),max(mean(SD)+std(SD)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
xlim([0 bin_xlims(end)])
% ylim([0.002 0.008])
title('Standard Deviation','FontSize',20,'FontWeight','bold','color','k')
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

