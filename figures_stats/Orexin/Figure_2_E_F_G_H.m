%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

Almo_folder = {
    'm1871\Oct_07_2024';
    'm1900\Oct_07_2024';
    'm1873\Oct_07_2024';
    'm1879\Oct_07_2024';
    'm1821\Oct_07_2024';
    'm1822\Oct_07_2024';
    'm1825\Oct_07_2024';
    'm1826\Oct_07_2024';
    'm2132\Apr_08_2025';
%     'm2133\Apr_08_2025';
    'm2146\Apr_08_2025';
    'm2169\Apr_08_2025';
    'm2140\Apr_07_2025';
    'm2148\Apr_07_2025';
    'm2144\Apr_07_2025';
    'm2145\Apr_07_2025';
    'm2178\Apr_08_2025';
    'm2179\Apr_08_2025';
    'm2175\Apr_08_2025'};

ctrl_folder = {
    'm1871\Oct_09_2024';
    'm1900\Oct_04_2024';
    'm1873\Oct_04_2024';
    'm1879\Oct_04_2024';
    'm1821\Oct_04_2024';
    'm1822\Oct_04_2024';
    'm1825\Oct_09_2024';
    'm1826\Oct_09_2024';
    'm2132\Apr_07_2025';
%     'm2133\Apr_07_2025';
    'm2146\Apr_07_2025';
    'm2169\Apr_07_2025';
    'm2140\Apr_14_2025';
    'm2148\Apr_14_2025';
    'm2144\Apr_14_2025';
    'm2145\Apr_14_2025';
    'm2178\Apr_14_2025';
    'm2179\Apr_14_2025';
    'm2175\Apr_14_2025'};

fpass_Almo=[
    10 13;
    10 13;
    9 12.5;
    8 13;
    8 13;
    10 14;
    9 13;
    9 12;
    8 13;
%     7 13;
    7 13;
    8 13;
    9 14;
    8 12.5;
    10 14;
    10 13;
    8.5 13;
    10 13;
    8 14];

fpass_ctrl=[
    10 13;
    11 13;
    10 12;
    9 13;
    8.5 13;
    10 14;
    8 13;
    8 13;
    8 13;
%     9.5 13;
    8.5 13;
    8 13;
    8 12.5;
    7 12;
    8 13;
    9 13;
    6 11;
    8 13;
    9.5 13.5];

outlier_pupil_Almo=[
    0 100;
    0 100;
    0 99.2;
    0 100;
    0 99.9;
    0 100;
    0 100;
    0 22;
    0 99.97;
%     0 93.1;
    0 100;
    0 99.95;
    0 100;
    0 99.99;
    0 99.99;
    0 100;
    0 99.99;
    0 99.99;
    0 100];

outlier_pupil_ctrl=[
    0 100;
    0 99.9;
    0 99.9;
    0 99.9;
    0 100;
    0 99.9;
    0 100;
    0 95;
    0 100;
%     0 100;
    0 93.4;
    0 96.9;
    0 100;
    0 100;
    0 100;
    0 99.99;
    0 100;
    0 100;
    0 99.99];

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

stimFreqs = 20;

pre_bout_duration = 100;   % 5 seconds.

trace_duration = 500;   % 30 seconds.

wheel_radius = 10;     % cm

unit = (pi/180)*wheel_radius*Sample_Rate;

smooth_window = 3;

Run_Almo = [];
Run_ctrl = [];

Pupil_Almo = [];
Pupil_ctrl = [];

HR_Almo = [];
HR_ctrl = [];

HRV_Almo = [];
HRV_ctrl = [];

RMSSD_on_Almo = [];
RMSSD_off_Almo = [];
RMSSD_on_ctrl = [];
RMSSD_off_ctrl = [];

for I=1:size(Almo_folder,1)

    Data_Folder = [Directory Almo_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

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
    Abs_speedDeg_outlier = filloutliers(Abs_speedDeg,"nearest","percentiles",[0 99.99]);
    speed_smooth = movmean(Abs_speedDeg_outlier,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);
    speed_smooth_resampled_cm = speed_smooth_resampled*unit;

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

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil_Almo(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_Almo(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %% HRV

    SD = abs(diff(RR_intervals)).*1000;
    SD_interp = interp1(pksLocs(3:end),SD,frame_time,'nearest','extrap');
    SD_interp_smooth = movmean(SD_interp,[smooth_window*FrameRate 0]);
    
    %% Optogenetic stimulation bouts

    run_trials = [];
    pupil_trials = [];
    HR_trials = [];
    HRV_trials = [];
    laser_on_RMSSD_trials = [];
    laser_off_RMSSD_trials = [];

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_trials = [run_trials; speed_smooth_resampled_cm(t_start:t_end)];

        raw_pupil = pupil_smooth(t_start:t_end)';
        pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
        pupil_trials = [pupil_trials; raw_pupil];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR_trials = [HR_trials; raw_HR];

        raw_HRV = SD_interp_smooth(t_start:t_end);
        HRV_trials = [HRV_trials; raw_HRV];

        baseline_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)-pre_bout_duration) & pksLocs(2:end)<laser_onsets(II));
        baseline_RR_intervals = RR_intervals(baseline_heartbeat_index);
        baseline_RMSSD = sqrt(mean((diff(baseline_RR_intervals)).^2));

        laser_on_heartbeat_index = find(pksLocs(2:end)>=laser_onsets(II) & pksLocs(2:end)<(laser_onsets(II)+60));
        laser_on_RR_intervals = RR_intervals(laser_on_heartbeat_index);
        laser_on_RMSSD_trials = [laser_on_RMSSD_trials; sqrt(mean((diff(laser_on_RR_intervals)).^2))-baseline_RMSSD];

        laser_off_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)+60) & pksLocs(2:end)<(laser_onsets(II)+400));
        laser_off_RR_intervals = RR_intervals(laser_off_heartbeat_index);
        laser_off_RMSSD_trials = [laser_off_RMSSD_trials; sqrt(mean((diff(laser_off_RR_intervals)).^2))-baseline_RMSSD];

    end

    Run_Almo = [Run_Almo; mean(run_trials,1)];

    Pupil_Almo = [Pupil_Almo; mean(pupil_trials,1)];

    HR_Almo = [HR_Almo; mean(HR_trials,1)];

    HRV_Almo = [HRV_Almo; mean(HRV_trials)];

    RMSSD_on_Almo = [RMSSD_on_Almo; mean(laser_on_RMSSD_trials)];
    RMSSD_off_Almo = [RMSSD_off_Almo; mean(laser_off_RMSSD_trials)];

end

for I=1:size(ctrl_folder,1)

    Data_Folder = [Directory ctrl_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

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
    Abs_speedDeg_outlier = filloutliers(Abs_speedDeg,"nearest","percentiles",[0 99.99]);
    speed_smooth = movmean(Abs_speedDeg_outlier,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);
    speed_smooth_resampled_cm = speed_smooth_resampled*unit;

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

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil_ctrl(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_ctrl(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %% HRV

    SD = abs(diff(RR_intervals)).*1000;
    SD_interp = interp1(pksLocs(3:end),SD,frame_time,'nearest','extrap');
    SD_interp_smooth = movmean(SD_interp,[smooth_window*FrameRate 0]);

    %% Optogenetic stimulation bouts

    run_trials = [];
    pupil_trials = [];
    HR_trials = [];
    HRV_trials = [];
    laser_on_RMSSD_trials = [];
    laser_off_RMSSD_trials = [];

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_trials = [run_trials; speed_smooth_resampled_cm(t_start:t_end)];

        raw_pupil = pupil_smooth(t_start:t_end)';
        pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
        pupil_trials = [pupil_trials; raw_pupil];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR_trials = [HR_trials; raw_HR];

        raw_HRV = SD_interp_smooth(t_start:t_end);
        HRV_trials = [HRV_trials; raw_HRV];

        baseline_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)-pre_bout_duration) & pksLocs(2:end)<laser_onsets(II));
        baseline_RR_intervals = RR_intervals(baseline_heartbeat_index);
        baseline_RMSSD = sqrt(mean((diff(baseline_RR_intervals)).^2));

        laser_on_heartbeat_index = find(pksLocs(2:end)>=laser_onsets(II) & pksLocs(2:end)<(laser_onsets(II)+60));
        laser_on_RR_intervals = RR_intervals(laser_on_heartbeat_index);
        laser_on_RMSSD_trials = [laser_on_RMSSD_trials; sqrt(mean((diff(laser_on_RR_intervals)).^2))-baseline_RMSSD];

        laser_off_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)+60) & pksLocs(2:end)<(laser_onsets(II)+400));
        laser_off_RR_intervals = RR_intervals(laser_off_heartbeat_index);
        laser_off_RMSSD_trials = [laser_off_RMSSD_trials; sqrt(mean((diff(laser_off_RR_intervals)).^2))-baseline_RMSSD];

    end

    Run_ctrl = [Run_ctrl; mean(run_trials,1)];

    Pupil_ctrl = [Pupil_ctrl; mean(pupil_trials,1)];

    HR_ctrl = [HR_ctrl; mean(HR_trials,1)];

    HRV_ctrl = [HRV_ctrl; mean(HRV_trials)];

    RMSSD_on_ctrl = [RMSSD_on_ctrl; mean(laser_on_RMSSD_trials)];
    RMSSD_off_ctrl = [RMSSD_off_ctrl; mean(laser_off_RMSSD_trials)];

end

%%

xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

trials_num = size(Run_Almo,1);

fig = figure(1);
set(fig, 'Position', [2561 49 922 1315]);

subplot(4,2,1);
hold on
patch('XData',[0, 0, 60, 60],'YData',[0, 8, 8, 0],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_ctrl)+std(Run_ctrl)/sqrt(trials_num) fliplr(mean(Run_ctrl)-std(Run_ctrl)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.7);
plot1 = plot(xlims,mean(Run_ctrl),'Color',[0.5 0.5 0.5],'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_Almo)+std(Run_Almo)/sqrt(trials_num) fliplr(mean(Run_Almo)-std(Run_Almo)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.7);
plot2 = plot(xlims,mean(Run_Almo),'Color','k','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([0 6])
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Running Speed','FontSize',15,'FontWeight','bold','color','k')
ylabel('cm/s','FontSize',12,'FontWeight','bold');
legend([plot2, plot1], {'Almorexant', 'Vehicle'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
hold off

subplot(4,2,3)
hold on
patch('XData',[0, 0, 60, 60],'YData',[200, 900, 900, 200],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_ctrl)+std(Pupil_ctrl)/sqrt(trials_num) fliplr(mean(Pupil_ctrl)-std(Pupil_ctrl)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.7);
plot1 = plot(xlims,mean(Pupil_ctrl),'Color',[139 92 158]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_Almo)+std(Pupil_Almo)/sqrt(trials_num) fliplr(mean(Pupil_Almo)-std(Pupil_Almo)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot2 = plot(xlims,mean(Pupil_Almo),'Color',[229 114 190]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([200 600])
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Pupil Size','FontSize',15,'FontWeight','bold','color','k')
ylabel('pixels','FontSize',12,'FontWeight','bold');
legend([plot2, plot1], {'Almorexant', 'Vehicle'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
hold off

subplot(4,2,5)
hold on
patch('XData',[0, 0, 60, 60],'YData',[580, 720, 720, 580],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_ctrl)+std(HR_ctrl)/sqrt(trials_num) fliplr(mean(HR_ctrl)-std(HR_ctrl)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255, 214, 170]./255,'FaceAlpha',0.7);
plot1 = plot(xlims,mean(HR_ctrl),'Color',[250, 157, 86]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_Almo)+std(HR_Almo)/sqrt(trials_num) fliplr(mean(HR_Almo)-std(HR_Almo)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot2 = plot(xlims,mean(HR_Almo),'Color',[255 128 128]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([600 720])
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Heart Rate','FontSize',15,'FontWeight','bold')
ylabel('bpm','FontSize',12,'FontWeight','bold');
legend([plot2, plot1], {'Almorexant', 'Vehicle'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
hold off

subplot(4,2,7);
hold on
patch('XData',[0, 0, 60, 60],'YData',[0, 20, 20, 0],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HRV_ctrl)+std(HRV_ctrl)/sqrt(trials_num) fliplr(mean(HRV_ctrl)-std(HRV_ctrl)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255, 214, 170]./255,'FaceAlpha',0.7);
plot1 = plot(xlims,mean(HRV_ctrl),'Color',[250, 157, 86]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HRV_Almo)+std(HRV_Almo)/sqrt(trials_num) fliplr(mean(HRV_Almo)-std(HRV_Almo)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot2 = plot(xlims,mean(HRV_Almo),'Color',[255 128 128]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([2 18])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Heart Rate Variability','FontSize',15,'FontWeight','bold')
ylabel('SD (ms)','FontSize',12,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold')
legend([plot2, plot1], {'Almorexant', 'Vehicle'},'FontSize',12,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

%%

peak_Run_Almo = mean(Run_Almo(:,100*FrameRate+1:160*FrameRate),2) - mean(Run_Almo(:,40*FrameRate+1:100*FrameRate),2);

peak_Run_ctrl = mean(Run_ctrl(:,100*FrameRate+1:160*FrameRate),2) - mean(Run_ctrl(:,40*FrameRate+1:100*FrameRate),2);

peak_Run = [peak_Run_Almo peak_Run_ctrl];

figure(1);

axes('Position', [0.5703 0.7673 0.15 0.1561]);
hold on

for k = 1:size(peak_Run,1)
    plot([1.3 1.7],peak_Run(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_Run_Almo,1)
    plot(1.3,peak_Run_Almo(k),'marker','o','markersize',4,...
        'markeredgecolor','k','markerfacecolor','k',...
        'linestyle','none');
end

for k = 1:size(peak_Run_ctrl,1)
    plot(1.7,peak_Run_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_Run_Almo,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_Run_Almo,1));

[S_chow,M_chow] = std(peak_Run_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_Run_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",'k','LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[0.5 0.5 0.5],'LineWidth',2);
plot(1.1, M_PB,'marker','o','color','k','linewidth',3,'markeredgecolor','k','markerfacecolor','k','markersize',5);
plot(1.9, M_chow,'marker','o','color',[0.5 0.5 0.5],'linewidth',3,'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5);

line([1.3 1.7], [7, 7], 'Color', 'k', 'LineWidth', 2);
text(1.5, 7.3, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, 3.5, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -4.4, 'Almo', 'Color', 'k', 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -4.4, 'Veh', 'Color', [0.5 0.5 0.5], 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-4 8])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Locomotion (cm/s)','FontSize',12,'FontWeight','bold');

hold off

[h_Run, p_Run, ~, stats_Run] = ttest(peak_Run_Almo,peak_Run_ctrl,'Tail','left')
[h_Run_Almo, p_Run_Almo, ~, stats_Run_Almo] = ttest(peak_Run_Almo,0,'Tail','right')
[h_Run_ctrl, p_Run_ctrl, ~, stats_Run_ctrl] = ttest(peak_Run_ctrl,0,'Tail','right')

%%

peak_Pupil_Almo = mean(Pupil_Almo(:,100*FrameRate+1:160*FrameRate),2) - mean(Pupil_Almo(:,40*FrameRate+1:100*FrameRate),2);

peak_Pupil_ctrl = mean(Pupil_ctrl(:,100*FrameRate+1:160*FrameRate),2) - mean(Pupil_ctrl(:,40*FrameRate+1:100*FrameRate),2);

peak_Pupil = [peak_Pupil_Almo peak_Pupil_ctrl];

figure(1);

axes('Position', [0.5703 0.5482 0.15 0.1561]);
hold on

for k = 1:size(peak_Pupil,1)
    plot([1.3 1.7],peak_Pupil(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_Pupil_Almo,1)
    plot(1.3,peak_Pupil_Almo(k),'marker','o','markersize',4,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Pupil_ctrl,1)
    plot(1.7,peak_Pupil_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_Pupil_Almo,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_Pupil_Almo,1));

[S_chow,M_chow] = std(peak_Pupil_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_Pupil_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[229 114 190]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[139 92 158]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[229 114 190]./255,'linewidth',3,'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[139 92 158]./255,'linewidth',3,'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,'markersize',5);

line([1.3 1.7], [600, 600], 'Color', 'k', 'LineWidth', 2);
text(1.5, 612, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 210, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, 290, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -75, 'Almo', 'Color', [229 114 190]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -75, 'Veh', 'Color', [139 92 158]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-50 650])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Pupil (pixels)','FontSize',12,'FontWeight','bold');

hold off

[h_Pupil, p_Pupil, ~, stats_Pupil] = ttest(peak_Pupil_Almo,peak_Pupil_ctrl,'Tail','left')
[h_Pupil_Almo, p_Pupil_Almo, ~, stats_Pupil_Almo] = ttest(peak_Pupil_Almo,0,'Tail','right')
[h_Pupil_ctrl, p_Pupil_ctrl, ~, stats_Pupil_ctrl] = ttest(peak_Pupil_ctrl,0,'Tail','right')

%%

peak_HR_Almo = mean(HR_Almo(:,100*FrameRate+1:160*FrameRate),2) - mean(HR_Almo(:,40*FrameRate+1:100*FrameRate),2);

peak_HR_ctrl = mean(HR_ctrl(:,100*FrameRate+1:160*FrameRate),2) - mean(HR_ctrl(:,40*FrameRate+1:100*FrameRate),2);

peak_HR = [peak_HR_Almo peak_HR_ctrl];

figure(1);

axes('Position', [0.5703 0.3291 0.15 0.1561]);
hold on

for k = 1:size(peak_HR,1)
    plot([1.3 1.7],peak_HR(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_HR_Almo,1)
    plot(1.3,peak_HR_Almo(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HR_ctrl,1)
    plot(1.7,peak_HR_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HR_Almo,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HR_Almo,1));

[S_chow,M_chow] = std(peak_HR_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HR_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[250, 157, 86]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[250, 157, 86]./255,'linewidth',3,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',5);

line([1.3 1.7], [85, 85], 'Color', 'k', 'LineWidth', 2);
text(1.5, 93, 'p=0.0647', 'FontSize', 10,'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 27, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, 40, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -46, 'Almo', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -46, 'Veh', 'Color', [250, 157, 86]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-40 100])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Heart Rate (bpm)','FontSize',12,'FontWeight','bold');

hold off

[h_HR, p_HR, ~, stats_HR] = ttest(peak_HR_Almo,peak_HR_ctrl,'Tail','left')
[h_HR_Almo, p_HR_Almo, ~, stats_HR_Almo] = ttest(peak_HR_Almo,0,'Tail','right')
[h_HR_ctrl, p_HR_ctrl, ~, stats_HR_ctrl] = ttest(peak_HR_ctrl,0,'Tail','right')

%%

% peak_HR_20Hz = mean(HR_20Hz(:,160*FrameRate+1:460*FrameRate),2) - mean(HR_20Hz(:,40*FrameRate+1:100*FrameRate),2);
% 
% peak_HR_5Hz = mean(HR_5Hz(:,160*FrameRate+1:460*FrameRate),2) - mean(HR_5Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_HR_Almo = mean(HR_Almo(:,160*FrameRate+1:end),2) - mean(HR_Almo(:,40*FrameRate+1:100*FrameRate),2);

peak_HR_ctrl = mean(HR_ctrl(:,160*FrameRate+1:end),2) - mean(HR_ctrl(:,40*FrameRate+1:100*FrameRate),2);

peak_HR = [peak_HR_Almo peak_HR_ctrl];

figure(1);

axes('Position', [0.83 0.3291 0.15 0.1561]);
hold on

for k = 1:size(peak_HR,1)
    plot([1.3 1.7],peak_HR(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_HR_Almo,1)
    plot(1.3,peak_HR_Almo(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HR_ctrl,1)
    plot(1.7,peak_HR_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HR_Almo,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HR_Almo,1));

[S_chow,M_chow] = std(peak_HR_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HR_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[250, 157, 86]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[250, 157, 86]./255,'linewidth',3,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',5);

line([1.3 1.7], [30, 30], 'Color', 'k', 'LineWidth', 2);
text(1.5, 36, 'p=0.0520', 'FontSize', 10,'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, -33, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -65, 'Almo', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -65, 'Veh', 'Color', [250, 157, 86]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-60 40])
xlim([0.9 2.1])

title({'Laser OFF'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Heart Rate (bpm)','FontSize',12,'FontWeight','bold');

hold off

[h_HR, p_HR, ~, stats_HR] = ttest(peak_HR_Almo,peak_HR_ctrl,'Tail','right')
[h_HR_Almo, p_HR_Almo, ~, stats_HR_Almo] = ttest(peak_HR_Almo,0,'Tail','left')
[h_HR_ctrl, p_HR_ctrl, ~, stats_HR_ctrl] = ttest(peak_HR_ctrl,0,'Tail','left')

%%

peak_HRV_Almo = RMSSD_on_Almo.*1000;

peak_HRV_ctrl = RMSSD_on_ctrl.*1000;

peak_HRV = [peak_HRV_Almo peak_HRV_ctrl];

figure(1);

axes('Position', [0.5703 0.1100 0.15 0.1561]);
hold on

for k = 1:size(peak_HRV,1)
    plot([1.3 1.7],peak_HRV(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_HRV_Almo,1)
    plot(1.3,peak_HRV_Almo(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HRV_ctrl,1)
    plot(1.7,peak_HRV_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HRV_Almo,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HRV_Almo,1));

[S_chow,M_chow] = std(peak_HRV_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HRV_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[250, 157, 86]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[250, 157, 86]./255,'linewidth',3,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',5);

text(1.1, -5, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -13, 'Almo', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -13, 'Veh', 'Color', [250, 157, 86]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-12 12])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ RMSSD (ms)','FontSize',12,'FontWeight','bold');

hold off

[h_HRV, p_HRV, ~, stats_HRV] = ttest(peak_HRV_Almo,peak_HRV_ctrl,'Tail','left')
[h_HRV_Almo, p_HRV_Almo, ~, stats_HRV_Almo] = ttest(peak_HRV_Almo,0,'Tail','left')
[h_HRV_ctrl, p_HRV_ctrl, ~, stats_HRV_ctrl] = ttest(peak_HRV_ctrl,0,'Tail','left')

%%

peak_HRV_Almo = RMSSD_off_Almo.*1000;

peak_HRV_ctrl = RMSSD_off_ctrl.*1000;

peak_HRV = [peak_HRV_Almo peak_HRV_ctrl];

figure(1);

axes('Position', [0.83 0.1100 0.15 0.1561]);
hold on

for k = 1:size(peak_HRV,1)
    plot([1.3 1.7],peak_HRV(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_HRV_Almo,1)
    plot(1.3,peak_HRV_Almo(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HRV_ctrl,1)
    plot(1.7,peak_HRV_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HRV_Almo,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HRV_Almo,1));

[S_chow,M_chow] = std(peak_HRV_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HRV_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[250, 157, 86]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[250, 157, 86]./255,'linewidth',3,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',5);

line([1.3 1.7], [13, 13], 'Color', 'k', 'LineWidth', 2);
text(1.5, 13.5, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 4, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, 7, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -6, 'Almo', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -6, 'Veh', 'Color', [250, 157, 86]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-5 15])
xlim([0.9 2.1])

title({'Laser OFF'},'FontSize',15,'FontWeight','bold')
ylabel('Δ RMSSD (ms)','FontSize',12,'FontWeight','bold');

hold off

[h_HRV, p_HRV, ~, stats_HRV] = ttest(peak_HRV_Almo,peak_HRV_ctrl,'Tail','left')
[h_HRV_Almo, p_HRV_Almo, ~, stats_HRV_Almo] = ttest(peak_HRV_Almo,0,'Tail','right')
[h_HRV_ctrl, p_HRV_ctrl, ~, stats_HRV_ctrl] = ttest(peak_HRV_ctrl,0,'Tail','right')

%%

saveas(gcf, 'Orexin_Almo_OptoStim_Run_Pupil_HR_HRV.svg')

saveas(gcf, 'Orexin_Almo_OptoStim_Run_Pupil_HR_HRV.png')

