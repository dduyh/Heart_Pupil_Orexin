%%

clear; close all; clc;

%% Set the path for output data

Directory = 'P:\Yihui\data\';    % Main directory\

saline_folder = {'m1900\Sep_05_2024';
    'm1879\Sep_05_2024';
    'm1822\Sep_06_2024';
    'm1826\Sep_06_2024';
    'm1871\Sep_07_2024';
    'm1873\Sep_08_2024';
    'm2132\Mar_10_2025';
%     'm2133\Mar_10_2025';
    'm2140\Mar_10_2025';
    'm2148\Mar_13_2025';
    'm2144\Mar_13_2025';
    'm2178\Mar_13_2025';
    'm2179\Mar_13_2025';
    'm2145\Mar_13_2025';
    'm2146\Mar_13_2025';
    'm2169\Mar_13_2025';
    'm2175\Mar_13_2025'};

propranolol_folder = {'m1871\Sep_05_2024';
    'm1873\Sep_05_2024';
    'm1821\Sep_06_2024';
    'm1825\Sep_06_2024';
    'm1900\Sep_07_2024';
    'm1879\Sep_07_2024';
    'm2148\Mar_10_2025';
    'm2144\Mar_10_2025';
    'm2145\Mar_10_2025';
    'm2146\Mar_10_2025';
    'm2169\Mar_10_2025';
    'm2132\Mar_17_2025';
%     'm2133\Mar_17_2025';
    'm2140\Mar_17_2025';
    'm2175\Mar_17_2025';
    'm2178\Mar_17_2025';
    'm2179\Mar_17_2025'};

atropinum_folder = {'m1871\Oct_14_2024';
    'm1873\Oct_14_2024';
    'm1879\Oct_14_2024';
    'm1900\Oct_14_2024';
    'm1821\Oct_14_2024';
    'm1822\Oct_14_2024';
    'm1825\Oct_14_2024';
    'm1826\Oct_14_2024'};

atenolol_folder = {'m1871\Sep_28_2024';
    'm1873\Sep_28_2024';
    'm1879\Sep_28_2024';
    'm1900\Sep_28_2024';
    'm1821\Sep_30_2024';
    'm1822\Sep_30_2024';
    'm1825\Sep_30_2024';
    'm1826\Sep_30_2024'};

fpass_saline=[10 14;
    9 15;
    10 13;
    8 12;
    9 13;
    9 13;
    9 13;
%     8.5 13.5;
    9 13;
    7 13;
    9 13.5;
    7 12.5;
    7 13;
    7 13;
    9 14;
    9 13;
    9 14];

fpass_propranolol=[7 12;
    3 9;
    6 9;
    5 9;
    7 12;
    6 11;
    4 9;
    7 11;
    5 10;
    4 8;
    6 10;
    4 8;
%     8 13;
    3 8;
    5 9;
    7 10;
    6 10];

fpass_atropinum=[11 13;
    10 12.5;
    8 13;
    11 13;
    10 12.5;
    11 13;
    9 12;
    9 12];

fpass_atenolol=[7 10;
    5 10;
    7 10;
    8 11;
    6 10;
    8 11;
    6 10;
    6 9.5];

outlier_pupil_saline=[0.1 97.5;
    0 98.7;
    0 99.98;
    0 95.7;
    0 99.4;
    0 99.7;
    0 98.5;
%     0 98.2;
    0 93.8;
    0 99.1;
    0 98.9;
    0 98.8;
    0 99.1;
    0 95.0;
    0 99.1;
    0 99.1;
    0 99.1];

outlier_pupil_propranolol=[0 96.7;
    0.1 96.7;
    0 99.8;
    0 99.4;
    0.5 96.5;
    0 99.8;
    0 99.91;
    0 98.6;
    0 98.9;
    0 99.1;
    0 98.6;
    0 99.1;
%     0 99.0;
    0 99.2;
    0 99.1;
    0 99.3;
    0 99.2];

outlier_pupil_atropinum=[0 99.52;
    0 98.8;
    0 99.64;
    0 98.6;
    0 99.05;
    0 98.79;
    0 99.82;
    0 98.21];

outlier_pupil_atenolol=[0 99.2;
    0 98.2;
    0 99.2;
    0.2 96.8;
    0 98.6;
    0 99.3;
    0 99.3;
    0 99.4];

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

stimFreqs = 20;

pre_bout_duration = 100;   % 5 seconds.

trace_duration = 500;   % 30 seconds.

wheel_radius = 10;     % cm

unit = (pi/180)*wheel_radius*Sample_Rate;

smooth_window = 3;

Run_saline = [];
Run_propranolol = [];
Run_atropinum = [];
Run_atenolol = [];

HR_saline = [];
HR_saline_delta = [];
HR_propranolol = [];
HR_propranolol_delta = [];
HR_atropinum = [];
HR_atropinum_delta = [];
HR_atenolol = [];
HR_atenolol_delta = [];

Pupil_saline = [];
Pupil_saline_zscore = [];
Pupil_propranolol = [];
Pupil_propranolol_zscore = [];
Pupil_atropinum = [];
Pupil_atropinum_zscore = [];
Pupil_atenolol = [];
Pupil_atenolol_zscore = [];

%%

for I=1:size(saline_folder,1)

    Data_Folder = [Directory saline_folder{I} '\'];

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

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil_saline(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    pupil_smooth_zscore = zscore(pupil_smooth);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_saline(I,:);

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

    %% Optogenetic stimulation bouts

    run_trials = [];
    pupil_trials = [];
    pupil_trials_zscore = [];
    HR_trials = [];
    HR_trials_delta = [];

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_trials = [run_trials; speed_smooth_resampled_cm(t_start:t_end)];

        pupil_trials = [pupil_trials; pupil_smooth(t_start:t_end)'];
        pupil_trials_zscore = [pupil_trials_zscore; pupil_smooth_zscore(t_start:t_end)'];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        delta_HR = raw_HR - mean(raw_HR(1:pre_bout_duration*FrameRate));
        HR_trials = [HR_trials; raw_HR];
        HR_trials_delta = [HR_trials_delta; delta_HR];

    end

    Run_saline = [Run_saline; mean(run_trials,1)];

    Pupil_saline = [Pupil_saline; mean(pupil_trials,1)];
    Pupil_saline_zscore = [Pupil_saline_zscore; mean(pupil_trials_zscore,1)];

    HR_saline = [HR_saline; mean(HR_trials,1)];
    HR_saline_delta = [HR_saline_delta; mean(HR_trials_delta,1)];

end

%%

for I=1:size(propranolol_folder,1)

    Data_Folder = [Directory propranolol_folder{I} '\'];

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

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil_propranolol(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    pupil_smooth_zscore = zscore(pupil_smooth);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_propranolol(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.3,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %% Optogenetic stimulation bouts

    run_trials = [];
    pupil_trials = [];
    pupil_trials_zscore = [];
    HR_trials = [];
    HR_trials_delta = [];

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_trials = [run_trials; speed_smooth_resampled_cm(t_start:t_end)];

        pupil_trials = [pupil_trials; pupil_smooth(t_start:t_end)'];
        pupil_trials_zscore = [pupil_trials_zscore; pupil_smooth_zscore(t_start:t_end)'];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        delta_HR = raw_HR - mean(raw_HR(1:pre_bout_duration*FrameRate));
        HR_trials = [HR_trials; raw_HR];
        HR_trials_delta = [HR_trials_delta; delta_HR];

    end

    Run_propranolol = [Run_propranolol; mean(run_trials,1)];

    Pupil_propranolol = [Pupil_propranolol; mean(pupil_trials,1)];
    Pupil_propranolol_zscore = [Pupil_propranolol_zscore; mean(pupil_trials_zscore,1)];

    HR_propranolol = [HR_propranolol; mean(HR_trials,1)];
    HR_propranolol_delta = [HR_propranolol_delta; mean(HR_trials_delta,1)];

end

%%

for I=1:size(atropinum_folder,1)

    Data_Folder = [Directory atropinum_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder atropinum_folder{I}(1:5) '_' atropinum_folder{I}(7:end) 'DLC_resnet50_Pupil_trackingDec18shuffle1_1000000_filtered.csv'],3,1);

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

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil_atropinum(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    pupil_smooth_zscore = zscore(pupil_smooth);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_atropinum(I,:);

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

    %% Optogenetic stimulation bouts

    run_trials = [];
    pupil_trials = [];
    pupil_trials_zscore = [];
    HR_trials = [];
    HR_trials_delta = [];

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_trials = [run_trials; speed_smooth_resampled_cm(t_start:t_end)];

        pupil_trials = [pupil_trials; pupil_smooth(t_start:t_end)'];
        pupil_trials_zscore = [pupil_trials_zscore; pupil_smooth_zscore(t_start:t_end)'];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        delta_HR = raw_HR - mean(raw_HR(1:pre_bout_duration*FrameRate));
        HR_trials = [HR_trials; raw_HR];
        HR_trials_delta = [HR_trials_delta; delta_HR];

    end

    Run_atropinum = [Run_atropinum; mean(run_trials,1)];

    Pupil_atropinum = [Pupil_atropinum; mean(pupil_trials,1)];
    Pupil_atropinum_zscore = [Pupil_atropinum_zscore; mean(pupil_trials_zscore,1)];

    HR_atropinum = [HR_atropinum; mean(HR_trials,1)];
    HR_atropinum_delta = [HR_atropinum_delta; mean(HR_trials_delta,1)];

end

%%

for I=1:size(atenolol_folder,1)

    Data_Folder = [Directory atenolol_folder{I} '\'];

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

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil_atenolol(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    pupil_smooth_zscore = zscore(pupil_smooth);
    
    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_atenolol(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.3,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %% Optogenetic stimulation bouts

    run_trials = [];
    pupil_trials = [];
    pupil_trials_zscore = [];
    HR_trials = [];
    HR_trials_delta = [];

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_trials = [run_trials; speed_smooth_resampled_cm(t_start:t_end)];

        pupil_trials = [pupil_trials; pupil_smooth(t_start:t_end)'];
        pupil_trials_zscore = [pupil_trials_zscore; pupil_smooth_zscore(t_start:t_end)'];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        delta_HR = raw_HR - mean(raw_HR(1:pre_bout_duration*FrameRate));
        HR_trials = [HR_trials; raw_HR];
        HR_trials_delta = [HR_trials_delta; delta_HR];

    end

    Run_atenolol = [Run_atenolol; mean(run_trials,1)];

    Pupil_atenolol = [Pupil_atenolol; mean(pupil_trials,1)];
    Pupil_atenolol_zscore = [Pupil_atenolol_zscore; mean(pupil_trials_zscore,1)];

    HR_atenolol = [HR_atenolol; mean(HR_trials,1)];
    HR_atenolol_delta = [HR_atenolol_delta; mean(HR_trials_delta,1)];

end

%%

xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

trials_num_saline = size(Run_saline,1);
trials_num_propranolol = size(Run_propranolol,1);
trials_num_atropinum = size(Run_atropinum,1);
trials_num_atenolol = size(Run_atenolol,1);

fig = figure(1);
set(fig, 'Position', [2561 49 922 1315]);

subplot(4,2,3);
hold on
patch('XData',[0, 0, 60, 60],'YData',[0, 12, 12, 0],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_saline)+std(Run_saline)/sqrt(trials_num_saline) fliplr(mean(Run_saline)-std(Run_saline)/sqrt(trials_num_saline))],'EdgeColor','none','FaceColor',[255 209 181]./255,'FaceAlpha',0.3);
plot1 = plot(xlims,mean(Run_saline),'Color',[250 153 94]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_atropinum)+std(Run_atropinum)/sqrt(trials_num_atropinum) fliplr(mean(Run_atropinum)-std(Run_atropinum)/sqrt(trials_num_atropinum))],'EdgeColor','none','FaceColor',[232 150 158]./255,'FaceAlpha',0.3);
plot2 = plot(xlims,mean(Run_atropinum),'Color',[248 82 101]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_atenolol)+std(Run_atenolol)/sqrt(trials_num_atenolol) fliplr(mean(Run_atenolol)-std(Run_atenolol)/sqrt(trials_num_atenolol))],'EdgeColor','none','FaceColor',[189 106 170]./255,'FaceAlpha',0.3);
plot3 = plot(xlims,mean(Run_atenolol),'Color',[125 49 136]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_propranolol)+std(Run_propranolol)/sqrt(trials_num_propranolol) fliplr(mean(Run_propranolol)-std(Run_propranolol)/sqrt(trials_num_propranolol))],'EdgeColor','none','FaceColor',[148 141 173]./255,'FaceAlpha',0.3);
plot4 = plot(xlims,mean(Run_propranolol),'Color',[58 56 105]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([0 12])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Running Speed','FontSize',15,'FontWeight','bold')
ylabel('cm/s','FontSize',12,'FontWeight','bold');
legend([plot4, plot3, plot2, plot1], {'PRO', 'ATE', 'ATR', 'SAL'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
hold off

subplot(4,2,5)
hold on
patch('XData',[0, 0, 60, 60],'YData',[-1, 3, 3, -1],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_saline_zscore)+std(Pupil_saline_zscore)/sqrt(trials_num_saline) fliplr(mean(Pupil_saline_zscore)-std(Pupil_saline_zscore)/sqrt(trials_num_saline))],'EdgeColor','none','FaceColor',[255 209 181]./255,'FaceAlpha',0.3);
plot1 = plot(xlims,mean(Pupil_saline_zscore),'Color',[250 153 94]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_atropinum_zscore)+std(Pupil_atropinum_zscore)/sqrt(trials_num_atropinum) fliplr(mean(Pupil_atropinum_zscore)-std(Pupil_atropinum_zscore)/sqrt(trials_num_atropinum))],'EdgeColor','none','FaceColor',[232 150 158]./255,'FaceAlpha',0.3);
plot2 = plot(xlims,mean(Pupil_atropinum_zscore),'Color',[248 82 101]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_atenolol_zscore)+std(Pupil_atenolol_zscore)/sqrt(trials_num_atenolol) fliplr(mean(Pupil_atenolol_zscore)-std(Pupil_atenolol_zscore)/sqrt(trials_num_atenolol))],'EdgeColor','none','FaceColor',[189 106 170]./255,'FaceAlpha',0.3);
plot3 = plot(xlims,mean(Pupil_atenolol_zscore),'Color',[125 49 136]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_propranolol_zscore)+std(Pupil_propranolol_zscore)/sqrt(trials_num_propranolol) fliplr(mean(Pupil_propranolol_zscore)-std(Pupil_propranolol_zscore)/sqrt(trials_num_propranolol))],'EdgeColor','none','FaceColor',[148 141 173]./255,'FaceAlpha',0.3);
plot4 = plot(xlims,mean(Pupil_propranolol_zscore),'Color',[58 56 105]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([-1 3])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Pupil Size','FontSize',15,'FontWeight','bold')
ylabel('z-score (s.d.)','FontSize',12,'FontWeight','bold');
legend([plot4, plot3, plot2, plot1], {'PRO', 'ATE', 'ATR', 'SAL'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
hold off

subplot(4,2,7)
hold on
patch('XData',[0, 0, 60, 60],'YData',[-100, 100, 100, -100],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_saline_delta)+std(HR_saline_delta)/sqrt(trials_num_saline) fliplr(mean(HR_saline_delta)-std(HR_saline_delta)/sqrt(trials_num_saline))],'EdgeColor','none','FaceColor',[255 209 181]./255,'FaceAlpha',0.3);
plot1 = plot(xlims,mean(HR_saline_delta),'Color',[250 153 94]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_atropinum_delta)+std(HR_atropinum_delta)/sqrt(trials_num_atropinum) fliplr(mean(HR_atropinum_delta)-std(HR_atropinum_delta)/sqrt(trials_num_atropinum))],'EdgeColor','none','FaceColor',[232 150 158]./255,'FaceAlpha',0.3);
plot2 = plot(xlims,mean(HR_atropinum_delta),'Color',[248 82 101]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_atenolol_delta)+std(HR_atenolol_delta)/sqrt(trials_num_atenolol) fliplr(mean(HR_atenolol_delta)-std(HR_atenolol_delta)/sqrt(trials_num_atenolol))],'EdgeColor','none','FaceColor',[189 106 170]./255,'FaceAlpha',0.3);
plot3 = plot(xlims,mean(HR_atenolol_delta),'Color',[125 49 136]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_propranolol_delta)+std(HR_propranolol_delta)/sqrt(trials_num_propranolol) fliplr(mean(HR_propranolol_delta)-std(HR_propranolol_delta)/sqrt(trials_num_propranolol))],'EdgeColor','none','FaceColor',[148 141 173]./255,'FaceAlpha',0.3);
plot4 = plot(xlims,mean(HR_propranolol_delta),'Color',[58 56 105]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([-80 70])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Heart Rate','FontSize',15,'FontWeight','bold')
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold')
ylabel('Δ bpm','FontSize',12,'FontWeight','bold');
legend([plot4, plot3, plot2, plot1], {'PRO', 'ATE', 'ATR', 'SAL'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
hold off

%%

peak_Run_saline = mean(Run_saline(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Run_saline(:,1:pre_bout_duration*FrameRate),2);

peak_Run_atropinum = mean(Run_atropinum(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Run_atropinum(:,1:pre_bout_duration*FrameRate),2);

peak_Run_atenolol = mean(Run_atenolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Run_atenolol(:,1:pre_bout_duration*FrameRate),2);

peak_Run_propranolol = mean(Run_propranolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Run_propranolol(:,1:pre_bout_duration*FrameRate),2);

Run_mean = [mean(peak_Run_saline) mean(peak_Run_atropinum) mean(peak_Run_atenolol) mean(peak_Run_propranolol)];
Run_sem = [std(peak_Run_saline)/sqrt(trials_num_saline) std(peak_Run_atropinum)/sqrt(trials_num_atropinum) std(peak_Run_atenolol)/sqrt(trials_num_atenolol) std(peak_Run_propranolol)/sqrt(trials_num_propranolol)];

figure(1);

axes('Position', [0.5703 0.5482 0.15 0.1561]);
hold on

b = bar(Run_mean,'EdgeColor','none');
b.FaceColor = 'flat';
b.CData(1,:) = [255 209 181]./255;
b.CData(2,:) = [232 150 158]./255;
b.CData(3,:) = [189 106 170]./255;
b.CData(4,:) = [148 141 173]./255;

for k = 1:size(peak_Run_saline,1)
    plot(1,peak_Run_saline(k),'marker','o','markersize',5,...
        'markeredgecolor',[250 153 94]./255,'markerfacecolor',[250 153 94]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Run_atropinum,1)
    plot(2,peak_Run_atropinum(k),'marker','o','markersize',5,...
        'markeredgecolor',[248 82 101]./255,'markerfacecolor',[248 82 101]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Run_atenolol,1)
    plot(3,peak_Run_atenolol(k),'marker','o','markersize',5,...
        'markeredgecolor',[125 49 136]./255,'markerfacecolor',[125 49 136]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Run_propranolol,1)
    plot(4,peak_Run_propranolol(k),'marker','o','markersize',5,...
        'markeredgecolor',[58 56 105]./255,'markerfacecolor',[58 56 105]./255,...
        'linestyle','none');
end

errorbar(1:4,Run_mean,Run_sem,'k','linestyle','none','linewidth',2,'CapSize',15);

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2 3 4])
xticklabels({'SAL','ATR','ATE','PRO'})

ylim([0 15])
xlim([0.4 4.6])

title({'Δ Locomotion'},'FontSize',15,'FontWeight','bold')
ylabel('Δ cm/s','FontSize',12,'FontWeight','bold');

hold off

[h_Run_sal_atr, p_Run_sal_atr, ~, stats_Run_sal_atr] = ttest2(peak_Run_saline,peak_Run_atropinum,'Tail','left')
[h_Run_sal_ate, p_Run_sal_ate, ~, stats_Run_sal_ate] = ttest2(peak_Run_saline,peak_Run_atenolol,'Tail','right')
[h_Run_sal_pro, p_Run_sal_pro, ~, stats_Run_sal_pro] = ttest2(peak_Run_saline,peak_Run_propranolol,'Tail','right')

%%

peak_Pupil_saline = mean(Pupil_saline(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Pupil_saline(:,1:pre_bout_duration*FrameRate),2);

peak_Pupil_atropinum = mean(Pupil_atropinum(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Pupil_atropinum(:,1:pre_bout_duration*FrameRate),2);

peak_Pupil_atenolol = mean(Pupil_atenolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Pupil_atenolol(:,1:pre_bout_duration*FrameRate),2);

peak_Pupil_propranolol = mean(Pupil_propranolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Pupil_propranolol(:,1:pre_bout_duration*FrameRate),2);

Pupil_mean = [mean(peak_Pupil_saline) mean(peak_Pupil_atropinum) mean(peak_Pupil_atenolol) mean(peak_Pupil_propranolol)];
Pupil_sem = [std(peak_Pupil_saline)/sqrt(trials_num_saline) std(peak_Pupil_atropinum)/sqrt(trials_num_atropinum) std(peak_Pupil_atenolol)/sqrt(trials_num_atenolol) std(peak_Pupil_propranolol)/sqrt(trials_num_propranolol)];

figure(1);

axes('Position', [0.5703 0.3291 0.15 0.1561]);
hold on

b = bar(Pupil_mean,'EdgeColor','none');
b.FaceColor = 'flat';
b.CData(1,:) = [255 209 181]./255;
b.CData(2,:) = [232 150 158]./255;
b.CData(3,:) = [189 106 170]./255;
b.CData(4,:) = [148 141 173]./255;

for k = 1:size(peak_Pupil_saline,1)
    plot(1,peak_Pupil_saline(k),'marker','o','markersize',5,...
        'markeredgecolor',[250 153 94]./255,'markerfacecolor',[250 153 94]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Pupil_atropinum,1)
    plot(2,peak_Pupil_atropinum(k),'marker','o','markersize',5,...
        'markeredgecolor',[248 82 101]./255,'markerfacecolor',[248 82 101]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Pupil_atenolol,1)
    plot(3,peak_Pupil_atenolol(k),'marker','o','markersize',5,...
        'markeredgecolor',[125 49 136]./255,'markerfacecolor',[125 49 136]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Pupil_propranolol,1)
    plot(4,peak_Pupil_propranolol(k),'marker','o','markersize',5,...
        'markeredgecolor',[58 56 105]./255,'markerfacecolor',[58 56 105]./255,...
        'linestyle','none');
end

errorbar(1:4,Pupil_mean,Pupil_sem,'k','linestyle','none','linewidth',2,'CapSize',15);

text(2, 810, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2 3 4])
xticklabels({'SAL','ATR','ATE','PRO'})

ylim([0 900])
xlim([0.4 4.6])

title({'Δ Pupil Size'},'FontSize',15,'FontWeight','bold')
ylabel('Δ pixels','FontSize',12,'FontWeight','bold');

hold off

[h_Pupil_sal_atr, p_Pupil_sal_atr, ~, stats_Pupil_sal_atr] = ttest2(peak_Pupil_saline,peak_Pupil_atropinum,'Tail','left')
[h_Pupil_sal_ate, p_Pupil_sal_ate, ~, stats_Pupil_sal_ate] = ttest2(peak_Pupil_saline,peak_Pupil_atenolol,'Tail','right')
[h_Pupil_sal_pro, p_Pupil_sal_pro, ~, stats_Pupil_sal_pro] = ttest2(peak_Pupil_saline,peak_Pupil_propranolol,'Tail','right')

%%

peak_HR_saline = mean(HR_saline(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(HR_saline(:,1:pre_bout_duration*FrameRate),2);

peak_HR_atropinum = mean(HR_atropinum(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(HR_atropinum(:,1:pre_bout_duration*FrameRate),2);

peak_HR_atenolol = mean(HR_atenolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(HR_atenolol(:,1:pre_bout_duration*FrameRate),2);

peak_HR_propranolol = mean(HR_propranolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(HR_propranolol(:,1:pre_bout_duration*FrameRate),2);

HR_mean = [mean(peak_HR_saline) mean(peak_HR_atropinum) mean(peak_HR_atenolol) mean(peak_HR_propranolol)];
HR_sem = [std(peak_HR_saline)/sqrt(trials_num_saline) std(peak_HR_atropinum)/sqrt(trials_num_atropinum) std(peak_HR_atenolol)/sqrt(trials_num_atenolol) std(peak_HR_propranolol)/sqrt(trials_num_propranolol)];

figure(1);

axes('Position', [0.5703 0.1100 0.15 0.1561]);
hold on

b = bar(HR_mean,'EdgeColor','none');
b.FaceColor = 'flat';
b.CData(1,:) = [255 209 181]./255;
b.CData(2,:) = [232 150 158]./255;
b.CData(3,:) = [189 106 170]./255;
b.CData(4,:) = [148 141 173]./255;

for k = 1:size(peak_HR_saline,1)
    plot(1,peak_HR_saline(k),'marker','o','markersize',5,...
        'markeredgecolor',[250 153 94]./255,'markerfacecolor',[250 153 94]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HR_atropinum,1)
    plot(2,peak_HR_atropinum(k),'marker','o','markersize',5,...
        'markeredgecolor',[248 82 101]./255,'markerfacecolor',[248 82 101]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HR_atenolol,1)
    plot(3,peak_HR_atenolol(k),'marker','o','markersize',5,...
        'markeredgecolor',[125 49 136]./255,'markerfacecolor',[125 49 136]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HR_propranolol,1)
    plot(4,peak_HR_propranolol(k),'marker','o','markersize',5,...
        'markeredgecolor',[58 56 105]./255,'markerfacecolor',[58 56 105]./255,...
        'linestyle','none');
end

errorbar(1:4,HR_mean,HR_sem,'k','linestyle','none','linewidth',2,'CapSize',15);

text(3, 30, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(4, 40, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2 3 4])
xticklabels({'SAL','ATR','ATE','PRO'})

ylim([-80 60])
xlim([0.4 4.6])

title({'Δ Heart Rate'},'FontSize',15,'FontWeight','bold')
ylabel('Δ bpm','FontSize',12,'FontWeight','bold');

hold off

[h_HR_sal_atr, p_HR_sal_atr, ~, stats_HR_sal_atr] = ttest2(peak_HR_saline,peak_HR_atropinum,'Tail','right')
[h_HR_sal_ate, p_HR_sal_ate, ~, stats_HR_sal_ate] = ttest2(peak_HR_saline,peak_HR_atenolol,'Tail','right')
[h_HR_sal_pro, p_HR_sal_pro, ~, stats_HR_sal_pro] = ttest2(peak_HR_saline,peak_HR_propranolol,'Tail','right')

