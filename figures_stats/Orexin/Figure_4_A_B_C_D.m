%%

clear; close all; clc;

%% Set the path for output data

Directory = 'P:\Yihui\data\';    % Main directory\

saline_folder = {'m2070\Nov_24_2024';
    'm2071\Nov_24_2024';
    'm2072\Nov_24_2024';
    'm2151\Jan_06_2025';
    'm2152\Jan_06_2025';
    'm2154\Jan_06_2025'};

propranolol_folder = {'m2070\Nov_25_2024';
    'm2071\Nov_25_2024';
    'm2072\Nov_25_2024';
    'm2151\Feb_25_2025';
    'm2152\Feb_25_2025';
    'm2154\Feb_25_2025'};

atropinum_folder = {'m2070\Dec_09_2024';
    'm2071\Dec_09_2024';
    'm2072\Dec_09_2024';
    'm2151\Feb_28_2025';
    'm2152\Feb_28_2025';
    'm2154\Feb_28_2025'};

atenolol_folder = {'m2070\Nov_29_2024';
    'm2071\Nov_29_2024';
    'm2072\Nov_29_2024';
    'm2151\Apr_03_2025';
    'm2152\Apr_03_2025';
    'm2154\Apr_03_2025'};

fpass_saline=[8 13;
    10 13;
    9.5 13.5;
    7 13;
    8 13;
    7 13];

fpass_propranolol=[6.5 10;
    6 10;
    6.5 9;
    6 8;
    6 9;
    5 8];

fpass_atropinum=[8 12.5;
    9 13;
    10 13;
    9.5 13;
    9 13;
    9 13];

fpass_atenolol=[7.5 11;
    6 11;
    6.5 11;
    5 9;
    5.5 9.5
    4 9];

outlier_pupil_saline=[0 99.5;
    0 99.2;
    0 99.6;
    0 98.7;
    0 99.2;
    0 99.75];

outlier_pupil_propranolol=[0 99.5;
    0 99.2;
    0 98.5;
    0 94.2;
    0 97.6;
    0 98.5];

outlier_pupil_atropinum=[0 98.3;
    0 99.4;
    0 98.1;
    0 98.9;
    0 99.2;
    0 99.4];

outlier_pupil_atenolol=[0 99.7;
    0 99.1;
    0 99.0;
    0 99.0;
    0 98.1;
    0 99.1];

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

pre_bout_duration = 20;   % 5 seconds.

trace_duration = 100;   % 30 seconds.

wheel_radius = 10;     % cm

unit = (pi/180)*wheel_radius*Sample_Rate;

speed_thresh = 2;     % cm/s

smooth_window = 1;

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

Cell_saline = [];
Cell_propranolol = [];
Cell_atropinum = [];
Cell_atenolol = [];

Cell_mouse_saline = [];
Cell_mouse_propranolol = [];
Cell_mouse_atropinum = [];
Cell_mouse_atenolol = [];

%%

cells_num_saline = 0;
running_bouts_num_saline = 0;

for I=1:size(saline_folder,1)

    Data_Folder = [Directory saline_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    Cell_Status_data = readtable([Data_Folder saline_folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder saline_folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);
    cell_times = cell_traces_data(:,1);
    cell_traces = cell_traces_data(:,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times = cell_times(~nanRows, :);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
    vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

    timepoint = times(data_onset:data_offset,1)';
    time = timepoint(1,:)-timepoint(1,1);

    frame_time = 0:(1/FrameRate):time(end);

    %% Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num_saline = cells_num_saline + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
%         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
    end

    cell_traces_interpolated_smooth = movmean(cell_traces_interpolated,[smooth_window*FrameRate 0],2);

    %% Analyse running signals

    running = datas(data_onset:data_offset,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed_smooth = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);
    speed_smooth_resampled_cm = speed_smooth_resampled*unit;

    %% Analyse pupil size

    Pupil_up_x = pupil_data(vid_onset:min(end,vid_offset),1);
    Pupil_up_y = pupil_data(vid_onset:min(end,vid_offset),2);
    Pupil_left_x = pupil_data(vid_onset:min(end,vid_offset),4);
    Pupil_left_y = pupil_data(vid_onset:min(end,vid_offset),5);
    Pupil_down_x = pupil_data(vid_onset:min(end,vid_offset),7);
    Pupil_down_y = pupil_data(vid_onset:min(end,vid_offset),8);
    Pupil_right_x = pupil_data(vid_onset:min(end,vid_offset),10);
    Pupil_right_y = pupil_data(vid_onset:min(end,vid_offset),11);

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

    ECG_raw = datas(data_onset:data_offset,2)';

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

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %% Running bouts

    binarized_speed = speed_smooth_resampled_cm > speed_thresh;

    running_onsets = find(diff(binarized_speed)==1);

    Injection_onset = find(datas(data_onset:data_offset,4),1,'last');
    running_onsets(running_onsets<Injection_onset*FrameRate/Sample_Rate)=[];
    running_onsets(running_onsets>(length(speed_smooth_resampled_cm)-trace_duration*FrameRate+pre_bout_duration*FrameRate))=[];

    good_running_onsets = [];

    for i = 1:length(running_onsets)

        ii = sum(binarized_speed(running_onsets(i)-pre_bout_duration*FrameRate+1:running_onsets(i)));

        if ii == 0
            good_running_onsets = [good_running_onsets,running_onsets(i)];
        end
    end

    running_bouts_num_saline = running_bouts_num_saline + length(good_running_onsets);

    run_all_bouts = [];
    pupil_all_bouts = [];
    pupil_all_bouts_zscore = [];
    hr_all_bouts = [];
    hr_all_bouts_delta = [];
    cell_all_bouts = [];

    for i = 1:length(good_running_onsets)

        t_start = good_running_onsets(i) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_all_bouts = [run_all_bouts; speed_smooth_resampled_cm(t_start:t_end)];

        pupil_all_bouts = [pupil_all_bouts; pupil_smooth(t_start:t_end)'];
        pupil_all_bouts_zscore = [pupil_all_bouts_zscore; pupil_smooth_zscore(t_start:t_end)'];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        delta_HR = raw_HR - mean(raw_HR(1:pre_bout_duration*FrameRate));
        hr_all_bouts = [hr_all_bouts; raw_HR];
        hr_all_bouts_delta = [hr_all_bouts_delta; delta_HR];

        cell_traces_bouts = cell_traces_interpolated_smooth(:,t_start:t_end);
        cell_all_bouts = cat(3, cell_all_bouts, cell_traces_bouts);

    end

    Run_saline = [Run_saline; mean(run_all_bouts,1)];

    Pupil_saline = [Pupil_saline; mean(pupil_all_bouts,1)];
    Pupil_saline_zscore = [Pupil_saline_zscore; mean(pupil_all_bouts_zscore,1)];

    HR_saline = [HR_saline; mean(hr_all_bouts,1)];
    HR_saline_delta = [HR_saline_delta; mean(hr_all_bouts_delta,1)];

    Cell_saline = [Cell_saline; mean(cell_all_bouts,3)];

    Cell_mouse_saline  = [Cell_mouse_saline; mean(mean(cell_all_bouts,3),1)];

end

%%

cells_num_propranolol = 0;
running_bouts_num_propranolol = 0;

for I=1:size(propranolol_folder,1)

    Data_Folder = [Directory propranolol_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    Cell_Status_data = readtable([Data_Folder propranolol_folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder propranolol_folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);
    cell_times = cell_traces_data(:,1);
    cell_traces = cell_traces_data(:,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times = cell_times(~nanRows, :);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
    vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

    timepoint = times(data_onset:data_offset,1)';
    time = timepoint(1,:)-timepoint(1,1);

    frame_time = 0:(1/FrameRate):time(end);

    %% Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num_propranolol = cells_num_propranolol + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
%         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
    end

    cell_traces_interpolated_smooth = movmean(cell_traces_interpolated,[smooth_window*FrameRate 0],2);

    %% Analyse running signals

    running = datas(data_onset:data_offset,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed_smooth = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);
    speed_smooth_resampled_cm = speed_smooth_resampled*unit;

    %% Analyse pupil size

    Pupil_up_x = pupil_data(vid_onset:min(end,vid_offset),1);
    Pupil_up_y = pupil_data(vid_onset:min(end,vid_offset),2);
    Pupil_left_x = pupil_data(vid_onset:min(end,vid_offset),4);
    Pupil_left_y = pupil_data(vid_onset:min(end,vid_offset),5);
    Pupil_down_x = pupil_data(vid_onset:min(end,vid_offset),7);
    Pupil_down_y = pupil_data(vid_onset:min(end,vid_offset),8);
    Pupil_right_x = pupil_data(vid_onset:min(end,vid_offset),10);
    Pupil_right_y = pupil_data(vid_onset:min(end,vid_offset),11);

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

    ECG_raw = datas(data_onset:data_offset,2)';

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

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    heartRate_bpm_interp_smooth_zscore = zscore(heartRate_bpm_interp_smooth);

    %% Running bouts

    binarized_speed = speed_smooth_resampled_cm > speed_thresh;

    running_onsets = find(diff(binarized_speed)==1);

    Injection_onset = find(datas(data_onset:data_offset,4),1,'last');
    running_onsets(running_onsets<Injection_onset*FrameRate/Sample_Rate)=[];
    running_onsets(running_onsets>(length(speed_smooth_resampled_cm)-trace_duration*FrameRate+pre_bout_duration*FrameRate))=[];

    good_running_onsets = [];

    for i = 1:length(running_onsets)

        ii = sum(binarized_speed(running_onsets(i)-pre_bout_duration*FrameRate+1:running_onsets(i)));

        if ii == 0
            good_running_onsets = [good_running_onsets,running_onsets(i)];
        end
    end

    running_bouts_num_propranolol = running_bouts_num_propranolol + length(good_running_onsets);

    run_all_bouts   = [];
    pupil_all_bouts = [];
    pupil_all_bouts_zscore = [];
    hr_all_bouts    = [];
    hr_all_bouts_delta    = [];
    cell_all_bouts  = [];

    for i = 1:length(good_running_onsets)

        t_start = good_running_onsets(i) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_all_bouts = [run_all_bouts; speed_smooth_resampled_cm(t_start:t_end)];

        pupil_all_bouts = [pupil_all_bouts; pupil_smooth(t_start:t_end)'];
        pupil_all_bouts_zscore = [pupil_all_bouts_zscore; pupil_smooth_zscore(t_start:t_end)'];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        delta_HR = raw_HR - mean(raw_HR(1:pre_bout_duration*FrameRate));
        hr_all_bouts = [hr_all_bouts; raw_HR];
        hr_all_bouts_delta = [hr_all_bouts_delta; delta_HR];

        cell_traces_bouts = cell_traces_interpolated_smooth(:,t_start:t_end);
        cell_all_bouts = cat(3, cell_all_bouts, cell_traces_bouts);

    end

    Run_propranolol = [Run_propranolol; mean(run_all_bouts,1)];

    Pupil_propranolol = [Pupil_propranolol; mean(pupil_all_bouts,1)];
    Pupil_propranolol_zscore = [Pupil_propranolol_zscore; mean(pupil_all_bouts_zscore,1)];

    HR_propranolol = [HR_propranolol; mean(hr_all_bouts,1)];
    HR_propranolol_delta = [HR_propranolol_delta; mean(hr_all_bouts_delta,1)];

    Cell_propranolol = [Cell_propranolol; mean(cell_all_bouts,3)];

    Cell_mouse_propranolol  = [Cell_mouse_propranolol; mean(mean(cell_all_bouts,3),1)];

end

%%

cells_num_atropinum = 0;
running_bouts_num_atropinum = 0;

for I=1:size(atropinum_folder,1)

    Data_Folder = [Directory atropinum_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec18shuffle1_1000000_filtered.csv'],3,1);

    Cell_Status_data = readtable([Data_Folder atropinum_folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder atropinum_folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);
    cell_times = cell_traces_data(:,1);
    cell_traces = cell_traces_data(:,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times = cell_times(~nanRows, :);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
    vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

    timepoint = times(data_onset:data_offset,1)';
    time = timepoint(1,:)-timepoint(1,1);

    frame_time = 0:(1/FrameRate):time(end);

    %% Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num_atropinum = cells_num_atropinum + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
%         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
    end

    cell_traces_interpolated_smooth = movmean(cell_traces_interpolated,[smooth_window*FrameRate 0],2);

    %% Analyse running signals

    running = datas(data_onset:data_offset,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed_smooth = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);
    speed_smooth_resampled_cm = speed_smooth_resampled*unit;

    %% Analyse pupil size

    Pupil_up_x = pupil_data(vid_onset:min(end,vid_offset),1);
    Pupil_up_y = pupil_data(vid_onset:min(end,vid_offset),2);
    Pupil_left_x = pupil_data(vid_onset:min(end,vid_offset),4);
    Pupil_left_y = pupil_data(vid_onset:min(end,vid_offset),5);
    Pupil_down_x = pupil_data(vid_onset:min(end,vid_offset),7);
    Pupil_down_y = pupil_data(vid_onset:min(end,vid_offset),8);
    Pupil_right_x = pupil_data(vid_onset:min(end,vid_offset),10);
    Pupil_right_y = pupil_data(vid_onset:min(end,vid_offset),11);

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

    ECG_raw = datas(data_onset:data_offset,2)';

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

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    heartRate_bpm_interp_smooth_zscore = zscore(heartRate_bpm_interp_smooth);

    %% Running bouts

    binarized_speed = speed_smooth_resampled_cm > speed_thresh;

    running_onsets = find(diff(binarized_speed)==1);

    Injection_onset = find(datas(data_onset:data_offset,4),1,'last');
    running_onsets(running_onsets<Injection_onset*FrameRate/Sample_Rate)=[];
    running_onsets(running_onsets>(length(speed_smooth_resampled_cm)-trace_duration*FrameRate+pre_bout_duration*FrameRate))=[];

    good_running_onsets = [];

    for i = 1:length(running_onsets)

        ii = sum(binarized_speed(running_onsets(i)-pre_bout_duration*FrameRate+1:running_onsets(i)));

        if ii == 0
            good_running_onsets = [good_running_onsets,running_onsets(i)];
        end
    end

    running_bouts_num_atropinum = running_bouts_num_atropinum + length(good_running_onsets);

    run_all_bouts   = [];
    pupil_all_bouts = [];
    pupil_all_bouts_zscore = [];
    hr_all_bouts    = [];
    hr_all_bouts_delta    = [];
    cell_all_bouts  = [];

    for i = 1:length(good_running_onsets)

        t_start = good_running_onsets(i) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_all_bouts = [run_all_bouts; speed_smooth_resampled_cm(t_start:t_end)];

        pupil_all_bouts = [pupil_all_bouts; pupil_smooth(t_start:t_end)'];
        pupil_all_bouts_zscore = [pupil_all_bouts_zscore; pupil_smooth_zscore(t_start:t_end)'];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        delta_HR = raw_HR - mean(raw_HR(1:pre_bout_duration*FrameRate));
        hr_all_bouts = [hr_all_bouts; raw_HR];
        hr_all_bouts_delta = [hr_all_bouts_delta; delta_HR];

        cell_traces_bouts = cell_traces_interpolated_smooth(:,t_start:t_end);
        cell_all_bouts = cat(3, cell_all_bouts, cell_traces_bouts);

    end

    Run_atropinum = [Run_atropinum; mean(run_all_bouts,1)];

    Pupil_atropinum = [Pupil_atropinum; mean(pupil_all_bouts,1)];
    Pupil_atropinum_zscore = [Pupil_atropinum_zscore; mean(pupil_all_bouts_zscore,1)];

    HR_atropinum = [HR_atropinum; mean(hr_all_bouts,1)];
    HR_atropinum_delta = [HR_atropinum_delta; mean(hr_all_bouts_delta,1)];

    Cell_atropinum = [Cell_atropinum; mean(cell_all_bouts,3)];

    Cell_mouse_atropinum  = [Cell_mouse_atropinum; mean(mean(cell_all_bouts,3),1)];

end

%%

cells_num_atenolol = 0;
running_bouts_num_atenolol = 0;

for I=1:size(atenolol_folder,1)

    Data_Folder = [Directory atenolol_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    Cell_Status_data = readtable([Data_Folder atenolol_folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder atenolol_folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);
    cell_times = cell_traces_data(:,1);
    cell_traces = cell_traces_data(:,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times = cell_times(~nanRows, :);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
    vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

    timepoint = times(data_onset:data_offset,1)';
    time = timepoint(1,:)-timepoint(1,1);

    frame_time = 0:(1/FrameRate):time(end);

    %% Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num_atenolol = cells_num_atenolol + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
%         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
    end

    cell_traces_interpolated_smooth = movmean(cell_traces_interpolated,[smooth_window*FrameRate 0],2);

    %% Analyse running signals

    running = datas(data_onset:data_offset,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed_smooth = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);
    speed_smooth_resampled_cm = speed_smooth_resampled*unit;

    %% Analyse pupil size

    Pupil_up_x = pupil_data(vid_onset:min(end,vid_offset),1);
    Pupil_up_y = pupil_data(vid_onset:min(end,vid_offset),2);
    Pupil_left_x = pupil_data(vid_onset:min(end,vid_offset),4);
    Pupil_left_y = pupil_data(vid_onset:min(end,vid_offset),5);
    Pupil_down_x = pupil_data(vid_onset:min(end,vid_offset),7);
    Pupil_down_y = pupil_data(vid_onset:min(end,vid_offset),8);
    Pupil_right_x = pupil_data(vid_onset:min(end,vid_offset),10);
    Pupil_right_y = pupil_data(vid_onset:min(end,vid_offset),11);

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

    ECG_raw = datas(data_onset:data_offset,2)';

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

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    heartRate_bpm_interp_smooth_zscore = zscore(heartRate_bpm_interp_smooth);

    %% Running bouts

    binarized_speed = speed_smooth_resampled_cm > speed_thresh;

    running_onsets = find(diff(binarized_speed)==1);

    Injection_onset = find(datas(data_onset:data_offset,4),1,'last');
    running_onsets(running_onsets<Injection_onset*FrameRate/Sample_Rate)=[];
    running_onsets(running_onsets>(length(speed_smooth_resampled_cm)-trace_duration*FrameRate+pre_bout_duration*FrameRate))=[];

    good_running_onsets = [];

    for i = 1:length(running_onsets)

        ii = sum(binarized_speed(running_onsets(i)-pre_bout_duration*FrameRate+1:running_onsets(i)));

        if ii == 0
            good_running_onsets = [good_running_onsets,running_onsets(i)];
        end
    end

    running_bouts_num_atenolol = running_bouts_num_atenolol + length(good_running_onsets);

    run_all_bouts   = [];
    pupil_all_bouts = [];
    pupil_all_bouts_zscore = [];
    hr_all_bouts    = [];
    hr_all_bouts_delta    = [];
    cell_all_bouts  = [];

    for i = 1:length(good_running_onsets)

        t_start = good_running_onsets(i) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_all_bouts = [run_all_bouts; speed_smooth_resampled_cm(t_start:t_end)];

        pupil_all_bouts = [pupil_all_bouts; pupil_smooth(t_start:t_end)'];
        pupil_all_bouts_zscore = [pupil_all_bouts_zscore; pupil_smooth_zscore(t_start:t_end)'];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        delta_HR = raw_HR - mean(raw_HR(1:pre_bout_duration*FrameRate));
        hr_all_bouts = [hr_all_bouts; raw_HR];
        hr_all_bouts_delta = [hr_all_bouts_delta; delta_HR];

        cell_traces_bouts = cell_traces_interpolated_smooth(:,t_start:t_end);
        cell_all_bouts = cat(3, cell_all_bouts, cell_traces_bouts);

    end

    Run_atenolol = [Run_atenolol; mean(run_all_bouts,1)];

    Pupil_atenolol = [Pupil_atenolol; mean(pupil_all_bouts,1)];
    Pupil_atenolol_zscore = [Pupil_atenolol_zscore; mean(pupil_all_bouts_zscore,1)];

    HR_atenolol = [HR_atenolol; mean(hr_all_bouts,1)];
    HR_atenolol_delta = [HR_atenolol_delta; mean(hr_all_bouts_delta,1)];

    Cell_atenolol = [Cell_atenolol; mean(cell_all_bouts,3)];

    Cell_mouse_atenolol  = [Cell_mouse_atenolol; mean(mean(cell_all_bouts,3),1)];

end

%%

xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

trials_num_saline = size(Run_saline,1);
trials_num_propranolol = size(Run_propranolol,1);
trials_num_atropinum = size(Run_atropinum,1);
trials_num_atenolol = size(Run_atenolol,1);

fig = figure(1);
set(fig, 'Position', [2561 49 922 1315]);

subplot(4,2,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Cell_mouse_saline)+std(Cell_mouse_saline)/sqrt(trials_num_saline) fliplr(mean(Cell_mouse_saline)-std(Cell_mouse_saline)/sqrt(trials_num_saline))],'EdgeColor','none','FaceColor',[255 209 181]./255,'FaceAlpha',0.3);
plot1 = plot(xlims,mean(Cell_mouse_saline),'Color',[250 153 94]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Cell_mouse_atropinum)+std(Cell_mouse_atropinum)/sqrt(trials_num_atropinum) fliplr(mean(Cell_mouse_atropinum)-std(Cell_mouse_atropinum)/sqrt(trials_num_atropinum))],'EdgeColor','none','FaceColor',[232 150 158]./255,'FaceAlpha',0.3);
plot2 = plot(xlims,mean(Cell_mouse_atropinum),'Color',[248 82 101]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Cell_mouse_atenolol)+std(Cell_mouse_atenolol)/sqrt(trials_num_atenolol) fliplr(mean(Cell_mouse_atenolol)-std(Cell_mouse_atenolol)/sqrt(trials_num_atenolol))],'EdgeColor','none','FaceColor',[189 106 170]./255,'FaceAlpha',0.3);
plot3 = plot(xlims,mean(Cell_mouse_atenolol),'Color',[125 49 136]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Cell_mouse_propranolol)+std(Cell_mouse_propranolol)/sqrt(trials_num_propranolol) fliplr(mean(Cell_mouse_propranolol)-std(Cell_mouse_propranolol)/sqrt(trials_num_propranolol))],'EdgeColor','none','FaceColor',[148 141 173]./255,'FaceAlpha',0.3);
plot4 = plot(xlims,mean(Cell_mouse_propranolol),'Color',[58 56 105]./255,'LineWidth',2);
line([0,0],[-0.5,4.5],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([-0.5 1.5])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Orx-GCaMP6s','FontSize',15,'FontWeight','bold')
ylabel('z-score (s.d.)','FontSize',12,'FontWeight','bold');
legend([plot4, plot3, plot2, plot1], {'PRO', 'ATE', 'ATR', 'SAL'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
hold off

subplot(4,2,3);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_saline)+std(Run_saline)/sqrt(trials_num_saline) fliplr(mean(Run_saline)-std(Run_saline)/sqrt(trials_num_saline))],'EdgeColor','none','FaceColor',[255 209 181]./255,'FaceAlpha',0.3);
plot1 = plot(xlims,mean(Run_saline),'Color',[250 153 94]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_atropinum)+std(Run_atropinum)/sqrt(trials_num_atropinum) fliplr(mean(Run_atropinum)-std(Run_atropinum)/sqrt(trials_num_atropinum))],'EdgeColor','none','FaceColor',[232 150 158]./255,'FaceAlpha',0.3);
plot2 = plot(xlims,mean(Run_atropinum),'Color',[248 82 101]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_atenolol)+std(Run_atenolol)/sqrt(trials_num_atenolol) fliplr(mean(Run_atenolol)-std(Run_atenolol)/sqrt(trials_num_atenolol))],'EdgeColor','none','FaceColor',[189 106 170]./255,'FaceAlpha',0.3);
plot3 = plot(xlims,mean(Run_atenolol),'Color',[125 49 136]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_propranolol)+std(Run_propranolol)/sqrt(trials_num_propranolol) fliplr(mean(Run_propranolol)-std(Run_propranolol)/sqrt(trials_num_propranolol))],'EdgeColor','none','FaceColor',[148 141 173]./255,'FaceAlpha',0.3);
plot4 = plot(xlims,mean(Run_propranolol),'Color',[58 56 105]./255,'LineWidth',2);
line([0,0],[0,14],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([0 14])
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
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_saline_zscore)+std(Pupil_saline_zscore)/sqrt(trials_num_saline) fliplr(mean(Pupil_saline_zscore)-std(Pupil_saline_zscore)/sqrt(trials_num_saline))],'EdgeColor','none','FaceColor',[255 209 181]./255,'FaceAlpha',0.3);
plot1 = plot(xlims,mean(Pupil_saline_zscore),'Color',[250 153 94]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_atropinum_zscore)+std(Pupil_atropinum_zscore)/sqrt(trials_num_atropinum) fliplr(mean(Pupil_atropinum_zscore)-std(Pupil_atropinum_zscore)/sqrt(trials_num_atropinum))],'EdgeColor','none','FaceColor',[232 150 158]./255,'FaceAlpha',0.3);
plot2 = plot(xlims,mean(Pupil_atropinum_zscore),'Color',[248 82 101]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_atenolol_zscore)+std(Pupil_atenolol_zscore)/sqrt(trials_num_atenolol) fliplr(mean(Pupil_atenolol_zscore)-std(Pupil_atenolol_zscore)/sqrt(trials_num_atenolol))],'EdgeColor','none','FaceColor',[189 106 170]./255,'FaceAlpha',0.3);
plot3 = plot(xlims,mean(Pupil_atenolol_zscore),'Color',[125 49 136]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_propranolol_zscore)+std(Pupil_propranolol_zscore)/sqrt(trials_num_propranolol) fliplr(mean(Pupil_propranolol_zscore)-std(Pupil_propranolol_zscore)/sqrt(trials_num_propranolol))],'EdgeColor','none','FaceColor',[148 141 173]./255,'FaceAlpha',0.3);
plot4 = plot(xlims,mean(Pupil_propranolol_zscore),'Color',[58 56 105]./255,'LineWidth',2);
line([0,0],[-1,3.2],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([-1 3.2])
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
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_saline_delta)+std(HR_saline_delta)/sqrt(trials_num_saline) fliplr(mean(HR_saline_delta)-std(HR_saline_delta)/sqrt(trials_num_saline))],'EdgeColor','none','FaceColor',[255 209 181]./255,'FaceAlpha',0.3);
plot1 = plot(xlims,mean(HR_saline_delta),'Color',[250 153 94]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_atropinum_delta)+std(HR_atropinum_delta)/sqrt(trials_num_atropinum) fliplr(mean(HR_atropinum_delta)-std(HR_atropinum_delta)/sqrt(trials_num_atropinum))],'EdgeColor','none','FaceColor',[232 150 158]./255,'FaceAlpha',0.3);
plot2 = plot(xlims,mean(HR_atropinum_delta),'Color',[248 82 101]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_atenolol_delta)+std(HR_atenolol_delta)/sqrt(trials_num_atenolol) fliplr(mean(HR_atenolol_delta)-std(HR_atenolol_delta)/sqrt(trials_num_atenolol))],'EdgeColor','none','FaceColor',[189 106 170]./255,'FaceAlpha',0.3);
plot3 = plot(xlims,mean(HR_atenolol_delta),'Color',[125 49 136]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_propranolol_delta)+std(HR_propranolol_delta)/sqrt(trials_num_propranolol) fliplr(mean(HR_propranolol_delta)-std(HR_propranolol_delta)/sqrt(trials_num_propranolol))],'EdgeColor','none','FaceColor',[148 141 173]./255,'FaceAlpha',0.3);
plot4 = plot(xlims,mean(HR_propranolol_delta),'Color',[58 56 105]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
line([0,0],[-100,200],'Color','k','linestyle','--','LineWidth',2);
ylim([-50 100])
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

peak_Cell_saline = mean(Cell_mouse_saline(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Cell_mouse_saline(:,1:pre_bout_duration*FrameRate),2);

peak_Cell_atropinum = mean(Cell_mouse_atropinum(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Cell_mouse_atropinum(:,1:pre_bout_duration*FrameRate),2);

peak_Cell_atenolol = mean(Cell_mouse_atenolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Cell_mouse_atenolol(:,1:pre_bout_duration*FrameRate),2);

peak_Cell_propranolol = mean(Cell_mouse_propranolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Cell_mouse_propranolol(:,1:pre_bout_duration*FrameRate),2);

peak_Cell = [peak_Cell_saline peak_Cell_atropinum peak_Cell_atenolol peak_Cell_propranolol];

Cell_mean = [mean(peak_Cell_saline) mean(peak_Cell_atropinum) mean(peak_Cell_atenolol) mean(peak_Cell_propranolol)];
Cell_sem = [std(peak_Cell_saline)/sqrt(trials_num_saline) std(peak_Cell_atropinum)/sqrt(trials_num_atropinum) std(peak_Cell_atenolol)/sqrt(trials_num_atenolol) std(peak_Cell_propranolol)/sqrt(trials_num_propranolol)];

figure(1);

axes('Position', [0.5703 0.7673 0.15 0.1561]);
hold on

b = bar(Cell_mean,'EdgeColor','none');
b.FaceColor = 'flat';
b.CData(1,:) = [255 209 181]./255;
b.CData(2,:) = [232 150 158]./255;
b.CData(3,:) = [189 106 170]./255;
b.CData(4,:) = [148 141 173]./255;

for k = 1:size(peak_Cell,1)
    plot(1:4,peak_Cell(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_Cell_saline,1)
    plot(1,peak_Cell_saline(k),'marker','o','markersize',5,...
        'markeredgecolor',[250 153 94]./255,'markerfacecolor',[250 153 94]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Cell_atropinum,1)
    plot(2,peak_Cell_atropinum(k),'marker','o','markersize',5,...
        'markeredgecolor',[248 82 101]./255,'markerfacecolor',[248 82 101]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Cell_atenolol,1)
    plot(3,peak_Cell_atenolol(k),'marker','o','markersize',5,...
        'markeredgecolor',[125 49 136]./255,'markerfacecolor',[125 49 136]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Cell_propranolol,1)
    plot(4,peak_Cell_propranolol(k),'marker','o','markersize',5,...
        'markeredgecolor',[58 56 105]./255,'markerfacecolor',[58 56 105]./255,...
        'linestyle','none');
end

errorbar(1:4,Cell_mean,Cell_sem,'k','linestyle','none','linewidth',2,'CapSize',15);

text(3, 1.25, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2 3 4])
xticklabels({'SAL','ATR','ATE','PRO'})

ylim([0 1.4])
xlim([0.4 4.6])

title({'Δ Orexin'},'FontSize',15,'FontWeight','bold')
ylabel('Δ z-score (s.d.)','FontSize',12,'FontWeight','bold');

hold off

[h_Cell_sal_atr, p_Cell_sal_atr, ~, stats_Cell_sal_atr] = ttest(peak_Cell_saline,peak_Cell_atropinum,'Tail','right')
[h_Cell_sal_ate, p_Cell_sal_ate, ~, stats_Cell_sal_ate] = ttest(peak_Cell_saline,peak_Cell_atenolol,'Tail','left')
[h_Cell_sal_pro, p_Cell_sal_pro, ~, stats_Cell_sal_pro] = ttest(peak_Cell_saline,peak_Cell_propranolol,'Tail','right')

%%

peak_Run_saline = mean(Run_saline(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Run_saline(:,1:pre_bout_duration*FrameRate),2);

peak_Run_atropinum = mean(Run_atropinum(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Run_atropinum(:,1:pre_bout_duration*FrameRate),2);

peak_Run_atenolol = mean(Run_atenolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Run_atenolol(:,1:pre_bout_duration*FrameRate),2);

peak_Run_propranolol = mean(Run_propranolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Run_propranolol(:,1:pre_bout_duration*FrameRate),2);

peak_Run = [peak_Run_saline peak_Run_atropinum peak_Run_atenolol peak_Run_propranolol];

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

for k = 1:size(peak_Run,1)
    plot(1:4,peak_Run(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

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

text(4, 4, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2 3 4])
xticklabels({'SAL','ATR','ATE','PRO'})

ylim([0 12])
xlim([0.4 4.6])

title({'Δ Locomotion'},'FontSize',15,'FontWeight','bold')
ylabel('Δ cm/s','FontSize',12,'FontWeight','bold');

hold off

[h_Run_sal_atr, p_Run_sal_atr, ~, stats_Run_sal_atr] = ttest(peak_Run_saline,peak_Run_atropinum,'Tail','left')
[h_Run_sal_ate, p_Run_sal_ate, ~, stats_Run_sal_ate] = ttest(peak_Run_saline,peak_Run_atenolol,'Tail','right')
[h_Run_sal_pro, p_Run_sal_pro, ~, stats_Run_sal_pro] = ttest(peak_Run_saline,peak_Run_propranolol,'Tail','right')


%%

peak_Pupil_saline = mean(Pupil_saline(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Pupil_saline(:,1:pre_bout_duration*FrameRate),2);

peak_Pupil_atropinum = mean(Pupil_atropinum(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Pupil_atropinum(:,1:pre_bout_duration*FrameRate),2);

peak_Pupil_atenolol = mean(Pupil_atenolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Pupil_atenolol(:,1:pre_bout_duration*FrameRate),2);

peak_Pupil_propranolol = mean(Pupil_propranolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(Pupil_propranolol(:,1:pre_bout_duration*FrameRate),2);

peak_Pupil = [peak_Pupil_saline peak_Pupil_atropinum peak_Pupil_atenolol peak_Pupil_propranolol];

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

for k = 1:size(peak_Pupil,1)
    plot(1:4,peak_Pupil(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

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

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2 3 4])
xticklabels({'SAL','ATR','ATE','PRO'})

ylim([0 600])
xlim([0.4 4.6])

title({'Δ Pupil Size'},'FontSize',15,'FontWeight','bold')
ylabel('Δ pixels','FontSize',12,'FontWeight','bold');

hold off

[h_Pupil_sal_atr, p_Pupil_sal_atr, ~, stats_Pupil_sal_atr] = ttest(peak_Pupil_saline,peak_Pupil_atropinum,'Tail','right')
[h_Pupil_sal_ate, p_Pupil_sal_ate, ~, stats_Pupil_sal_ate] = ttest(peak_Pupil_saline,peak_Pupil_atenolol,'Tail','right')
[h_Pupil_sal_pro, p_Pupil_sal_pro, ~, stats_Pupil_sal_pro] = ttest(peak_Pupil_saline,peak_Pupil_propranolol,'Tail','right')

%%

peak_HR_saline = mean(HR_saline(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(HR_saline(:,1:pre_bout_duration*FrameRate),2);

peak_HR_atropinum = mean(HR_atropinum(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(HR_atropinum(:,1:pre_bout_duration*FrameRate),2);

peak_HR_atenolol = mean(HR_atenolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(HR_atenolol(:,1:pre_bout_duration*FrameRate),2);

peak_HR_propranolol = mean(HR_propranolol(:,pre_bout_duration*FrameRate+1:2*pre_bout_duration*FrameRate),2) - mean(HR_propranolol(:,1:pre_bout_duration*FrameRate),2);

peak_HR = [peak_HR_saline peak_HR_atropinum peak_HR_atenolol peak_HR_propranolol];

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

for k = 1:size(peak_HR,1)
    plot(1:4,peak_HR(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

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

text(2, 60, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(3, 70, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(4, 60, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2 3 4])
xticklabels({'SAL','ATR','ATE','PRO'})

ylim([0 110])
xlim([0.4 4.6])

title({'Δ Heart Rate'},'FontSize',15,'FontWeight','bold')
ylabel('Δ bpm','FontSize',12,'FontWeight','bold');

hold off

[h_HR_sal_atr, p_HR_sal_atr, ~, stats_HR_sal_atr] = ttest(peak_HR_saline,peak_HR_atropinum,'Tail','right')
[h_HR_sal_ate, p_HR_sal_ate, ~, stats_HR_sal_ate] = ttest(peak_HR_saline,peak_HR_atenolol,'Tail','right')
[h_HR_sal_pro, p_HR_sal_pro, ~, stats_HR_sal_pro] = ttest(peak_HR_saline,peak_HR_propranolol,'Tail','right')
