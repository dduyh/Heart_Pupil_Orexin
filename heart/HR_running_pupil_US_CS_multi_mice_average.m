clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

folder = {'m485_L_VTA_R_NAc_oxLight\Apr_08_2024';
    'm486_L_NAc_nLightR\Apr_10_2024';
    'm487_L_VTA_nLightR\Apr_05_2024';
    'm487_L_VTA_nLightR\Apr_08_2024';
    'm1772\Apr_05_2024';
    'm1772\Apr_10_2024'};


Sample_Rate = 200;    % 200 scans per second.
FrameRate = 20;

trace_duration = 600;   % 600 seconds.

trials_num = size(folder,1);

run = NaN(trials_num,trace_duration*Sample_Rate);
trials = NaN(trials_num,trace_duration*Sample_Rate);    % 15 seconds per trial.
% pupil = NaN(trials_num,trace_duration*FrameRate);
pupil = NaN(trials_num,trace_duration);


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
    speed = movmean(Abs_speedDeg,100);
    run(I,:) = [speed(1) speed];

    % Raw ECG signals

    ECG_raw = datas(trace_start:trace_end,2)';

    % Remove baseline wandering

    fpass=[9 13];

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.01;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;

    heartRate_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate,1:trace_duration*Sample_Rate,'linear','extrap');
    heartRate_interp_smooth = movmean(heartRate_interp,5*Sample_Rate);
    trials(I,:) = heartRate_interp_smooth;

    % pupil size

    pupil_start = round((step_timepoint(4)-step_timepoint(1))*FrameRate)-360*FrameRate+1;
    pupil_end = pupil_start+trace_duration*FrameRate-1;

    raw_pupil = areas(pupil_start:pupil_end);

    for i = 1 : trace_duration
        medianArea(i) = median(raw_pupil(((i-1)*FrameRate+1):i*FrameRate));
    end

    pupil_zscored = (medianArea - mean(medianArea(241:360))) / std(medianArea(241:360));
    pupil(I,:) = pupil_zscored;

end

%%
figure;
xlims = (1:trace_duration*Sample_Rate)/Sample_Rate;
pupil_xlims = 1:trace_duration;

subplot(3,1,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(run)+std(run)/sqrt(trials_num) fliplr(mean(run)-std(run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',0.7);
plot(xlims,mean(run),'Color','k','LineWidth',2)
line([120,120],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([240,240],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([360,360],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
line([480,480],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
xlim([0 xlims(end)])
axis off
title('Abs. Running Speed','FontSize',20,'FontWeight','bold')
hold off

subplot(3,1,2)
hold on
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2)
line([120,120],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([240,240],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([360,360],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
line([480,480],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
xlim([0 pupil_xlims(end)])
ylim([-2 4])
axis off
title('Pupil Size','FontSize',20,'FontWeight','bold','color',[229 114 190]./255)
hold off

subplot(3,1,3)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(trials),'Color',[255 128 128]./255,'LineWidth',2)
line([120,120],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([240,240],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color',[33 113 181]/255,'linestyle','--','LineWidth',4);
line([360,360],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
line([480,480],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color',[227 26 28]/255,'linestyle','--','LineWidth',4);
xlim([0 xlims(end)])
ylim([11.3 12.7])
title('Heart Rate','FontSize',20,'FontWeight','bold','color',[255 128 128]./255)
ylabel('HR (Hz)','FontSize',15,'FontWeight','bold');

% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


