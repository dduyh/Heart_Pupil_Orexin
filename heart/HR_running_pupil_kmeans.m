%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';                     % Main directory\
mouse_name = 'm1026';            % Mouse name\
date = 'Jul_07_2024';                             % Date\
stim = 'Fear_Conditioning'; % Sucrose Quinine

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);
load([Data_Folder 'Pupil.mat'], 'areas');

%% Analyse ECG signals

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

vid_start = ceil(step_timepoint(1)*Sample_Rate)+1;

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);

smooth_window = 1;

%% Analyse running signals

running = datas(vid_start:end,1)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
Abs_speedDeg = abs(speedDeg);
speed = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
speed = [speed(1) speed];

%% k-means clustering running bouts

pre_bout_duration = 5;   % 5 seconds.

trace_duration = 30;   % 30 seconds.

speed_thresh = 0.02;

binarized_speed = speed > speed_thresh;

running_onsets = find(diff(binarized_speed)==1);

running_onsets(running_onsets<pre_bout_duration*Sample_Rate)=[];
running_onsets(running_onsets>(length(speed)-trace_duration*Sample_Rate+pre_bout_duration*Sample_Rate))=[];

good_running_onsets = [];

for i = 1:length(running_onsets)

    ii = sum(binarized_speed(running_onsets(i)-pre_bout_duration*Sample_Rate+1:running_onsets(i)));

    if ii == 0
        good_running_onsets = [good_running_onsets,running_onsets(i)];
    end
end

run = NaN(length(good_running_onsets),trace_duration*Sample_Rate);


for i = 1:length(good_running_onsets)

    t_start = good_running_onsets(i) - pre_bout_duration*Sample_Rate +1;
    t_end = t_start + trace_duration*Sample_Rate -1;

    run(i,:) = speed(t_start:t_end);

end

run_input_zscore = zscore(run);

km = kmeans(run_input_zscore,2);

figure;
plot(run(km==1,:)')
figure;
plot(run(km==2,:)')


%% Analyse pupil size

pupil = filloutliers(areas,"nearest","percentiles",[0 100]);

%% Analyse ECG signals

% Raw ECG signals

ECG_raw = datas(vid_start:end,2)';

% Remove baseline wandering

fpass=[11 12.5];

ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

% find peaks

minPeakPromVal=0.007;

[pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

% heart rate in time

RR_intervals = diff(pksLocs);
heartRate=1./RR_intervals;
heartRate_bpm=heartRate*60;

heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",[0 100]);

heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm_outlier,0:time(end)*Sample_Rate,'nearest','extrap');


