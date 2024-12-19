%% ECG data processing

clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm1746';            % Mouse name\
date = 'Mar_05_2024';                             % Date\
session = 'session_1';

Data_Folder = [Directory mouse_name '\' date '\' session '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

%% Analyse ECG signals

Sample_Rate = 1000;    % 1000 scans per second.

trace_duration = 360;   % 360 seconds.

vid_start = ceil(step_timepoint(1))*Sample_Rate+1;

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

trace_end = vid_start+trace_duration*Sample_Rate-1;

%% Raw ECG signals

ECG_raw = datas(vid_start:trace_end,2)';

%% Remove baseline wandering

fpass=[7 14];

ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

%% find peaks

minPeakPromVal=0.01;
figure;
findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

[pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

%% heart rate in time

RR_intervals = diff(pksLocs);
heartRate=1./RR_intervals;
heartRate_smooth = movmean(heartRate,100);


%% heart rate in time

% Analyse running signals
running = datas(vid_start:trace_end,1)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
Abs_speedDeg = abs(speedDeg);
speed = movmean(Abs_speedDeg,500);

speedThreshold = 0.002;
running_index = zeros(1,trace_duration*Sample_Rate);
running_index(speed > speedThreshold) = 1;

heartbeat_index = round(pksLocs * 1000);
movement_state = running_index(heartbeat_index);
movement_state(1) = [];

running_RR_intervals = RR_intervals(movement_state==1);
resting_RR_intervals = RR_intervals(movement_state==0);

running_heartRate = mean(1./running_RR_intervals)
resting_heartRate = mean(1./resting_RR_intervals)

%% Poincare Plots and HRV


running_RMSSD = sqrt(mean((diff(running_RR_intervals)).^2))
resting_RMSSD = sqrt(mean((diff(resting_RR_intervals)).^2))

%% Plot ECG signals
figure
subplot(2,1,1);
plot(pksLocs(2:end),heartRate_smooth)
ylim([0 30])
xlim([0 trace_duration])

subplot(2,1,2);
imagesc(movement_state)

%
% figure
% subplot(2,1,1);
% plot(time(300000:310000),ECG_raw(300000:310000),'color',[34 75 160]./255)
% ylim([1.7 2.3])
% title('ECG Raw Signal','FontSize',15,'FontWeight','bold','color',[34 75 160]./255)
% % axis off
%
% subplot(2,1,2);
% plot(time(300000:310000),ECG_Bandpass(300000:310000),'color',[34 75 160]./255)
% ylim([-0.3 0.3])
% title('ECG Raw Signal','FontSize',15,'FontWeight','bold','color',[34 75 160]./255)
% xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
% % axis off
%
% % xlim([0 time(end)])


