clc
clear
close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

folder = {'m1821\May_30_2024';
    'm1822\May_30_2024';
    'm1825\May_31_2024';
    'm1826\May_31_2024';
    'm2140\Mar_04_2025';
    'm2144\Mar_04_2025';
    'm2145\Mar_04_2025';
    'm2146\Mar_04_2025';
    'm2169\Mar_05_2025';
    'm2175\Mar_05_2025';
    'm2178\Mar_05_2025';
    'm2132\Mar_06_2025';
    'm2133\Mar_06_2025';
    'm2148\Mar_06_2025';
    'm2179\Mar_06_2025'
    'm839\Aug_17_2025';
    'm840\Aug_17_2025';
    'm843\Aug_22_2025';
    'm844\Aug_22_2025';
    'm833\Aug_22_2025';
    'm835\Aug_22_2025';
    'm836\Aug_22_2025'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

stimFreqs = 20;

pre_bout_duration = 100;   % 5 seconds.

trace_duration = 500;   % 90 seconds.

trials_num = 3*size(folder,1);

run = NaN(trials_num,trace_duration*FrameRate);
pupil = NaN(trials_num,trace_duration*FrameRate);
HR = NaN(trials_num,trace_duration*FrameRate);

xcorrwind = 5*FrameRate; % in samples

r_runpup = NaN(trials_num,2*xcorrwind+1);
r_runHR = NaN(trials_num,2*xcorrwind+1);
r_HRpup = NaN(trials_num,2*xcorrwind+1);

fpass_trials=[9 12;
    9 13;
    9.5 12;
    8 12;
    9 13;
    9 13;
    8 13;
    8 13;
    8 13;
    7 13;
    7 12;
    8 13;
    7 13;
    7 13;
    8 13;
    9 14;
    7 13.5;
    9 13;
    9.5 14;
    9 12.5;
    8.5 13;
    9.5 13];

outlier_pupil=[0 99;
    0 99;
    0 99;
    0 99;
    0 100;
    0 99.99;
    0 99.98;
    0 99.98;
    0 100;
    0 99.98;
    0 100;
    0 100;
    0 100;
    0 100;
    0 100;
    0 99.52;
    0 100;
    0 99.99;
    0 100;
    0 99.99;
    0 99.99;
    0 100];

smooth_window = 1;

%%

k = 1;

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    if I>15
        pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingJul16shuffle1_1000000_filtered.csv'],3,1);
    else
        pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);
    end

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs)-step_timepoint(1);

    data_start = ceil(step_timepoint(1)*FrameRate)+1;
    data_end = ceil(step_timepoint(1)*FrameRate+size(datas,1)*FrameRate/Sample_Rate);

    frame_time = (1/FrameRate):(1/FrameRate):size(datas,1)/Sample_Rate;

    %% Analyse running signals

    running = datas(:,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    Abs_speedDeg_outlier = filloutliers(Abs_speedDeg,"nearest","percentiles",[0 99.9]);
    speed_smooth = movmean(Abs_speedDeg_outlier,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);

    %% Analyse pupil size

    Pupil_up_x = pupil_data(data_start:data_end,1);
    Pupil_up_y = pupil_data(data_start:data_end,2);
    Pupil_left_x = pupil_data(data_start:data_end,4);
    Pupil_left_y = pupil_data(data_start:data_end,5);
    Pupil_down_x = pupil_data(data_start:data_end,7);
    Pupil_down_y = pupil_data(data_start:data_end,8);
    Pupil_right_x = pupil_data(data_start:data_end,10);
    Pupil_right_y = pupil_data(data_start:data_end,11);

    radii = zeros(size(Pupil_up_x,1),1);
    areas = zeros(size(Pupil_up_x,1),1);

    for i = 1:size(Pupil_up_x,1)
        X1(1) = Pupil_up_x(i,1);
        X1(2) = Pupil_up_y(i,1);

        X2(1) = Pupil_left_x(i,1);
        X2(2) = Pupil_left_y(i,1);

        Y1(1) = Pupil_down_x(i,1);
        Y1(2) = Pupil_down_y(i,1);

        Y2(1) = Pupil_right_x(i,1);
        Y2(2) = Pupil_right_y(i,1);

        MinorAxis = [X1(1), X1(2); Y1(1), Y1(2)];
        MinorAxisLength = pdist(MinorAxis, 'euclidean');

        MajorAxis = [X2(1), X2(2); Y2(1), Y2(2)];
        MajorAxisLength = pdist(MajorAxis, 'euclidean');

        diameters = mean([MajorAxisLength, MinorAxisLength]);
        radii(i,1) = diameters/2;
        areas(i,1) = pi*radii(i,1).^2;
    end

    pupil_areas = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

    pupil_smooth = movmean(pupil_areas,[smooth_window*FrameRate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    if I>4
        minPeakPromVal=0.07;
    else
        minPeakPromVal=0.007;
    end

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %%

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        raw_running = speed_smooth_resampled(t_start:t_end);
        run(k,:) = raw_running;

        raw_pupil = pupil_smooth(t_start:t_end);
        pupil(k,:) = raw_pupil;

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR(k,:) = raw_HR;

        [r_runpup(k,:),lags] = xcorr(raw_running-mean(raw_running),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');

        [r_runHR(k,:),lags] = xcorr(raw_running-mean(raw_running),raw_HR-mean(raw_HR),xcorrwind,'coeff');

        [r_HRpup(k,:),lags] = xcorr(raw_HR-mean(raw_HR),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');

        k = k+1;

    end

end

%%

xlims = (-5*FrameRate:5*FrameRate)/FrameRate;

fig = figure(1);
set(fig, 'Position', [7 773 2550 565]);

ax1 = subplot(1,3,1);
set(ax1, 'Position', [0.1300 0.1150 0.2134 0.7546]);
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runpup)+std(r_runpup)/sqrt(trials_num) fliplr(mean(r_runpup)-std(r_runpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',1);
plot(xlims,mean(r_runpup),'Color','k','LineWidth',2)
[maxcorr_runpup,l] = max(mean(r_runpup));
plot([xlims(l) xlims(l)],[0.04 maxcorr_runpup],'--','Color','k','LineWidth',2)
xlabel('Lag (seconds)','FontSize',18,'FontWeight','bold')
ylabel('Correlation coefficient','FontSize',18,'FontWeight','bold');
xlim([-5 5])
ylim([0.3 0.8])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

ax2 = subplot(1,3,2);
set(ax2, 'Position', [0.4108 0.1150 0.2134 0.7546]);
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_HRpup)+std(r_HRpup)/sqrt(trials_num) fliplr(mean(r_HRpup)-std(r_HRpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(xlims,mean(r_HRpup),'Color',[229 114 190]./255,'LineWidth',2)
[maxcorr_HRpup,l] = max(mean(r_HRpup));
plot([xlims(l) xlims(l)],[0.15 maxcorr_HRpup],'--','Color',[229 114 190]./255,'LineWidth',2)
xlabel('Lag (seconds)','FontSize',18,'FontWeight','bold')
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');
xlim([-5 5])
ylim([0.2 0.55])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

ax3 = subplot(1,3,3);
set(ax3, 'Position', [0.6916 0.1150 0.2134 0.7546]);
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runHR)+std(r_runHR)/sqrt(trials_num) fliplr(mean(r_runHR)-std(r_runHR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(r_runHR),'Color',[255 128 128]./255,'LineWidth',2)
[maxcorr_runHR,l] = max(mean(r_runHR));
plot([xlims(l) xlims(l)],[0.02 maxcorr_runHR],'--','Color',[255 128 128]./255,'LineWidth',2)
xlabel('Lag (seconds)','FontSize',18,'FontWeight','bold')
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');
xlim([-5 5])
ylim([0.2 0.5])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

%%

xlims = (-100*FrameRate+1:400*FrameRate)/FrameRate;

figure;

hold on
% patch('XData',[0, 0, 60, 60],'YData',[-5, 10, 10, -5],'EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR)+std(HR)/sqrt(trials_num) fliplr(mean(HR)-std(HR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(HR),'Color',[255 128 128]./255,'LineWidth',2)

title({'Heart Rate (Chrimson)',''},'FontSize',20,'FontWeight','bold')
ylabel('HR z-score','FontSize',15,'FontWeight','bold');

xlim([-100 400])
% ylim([-2 5])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off
