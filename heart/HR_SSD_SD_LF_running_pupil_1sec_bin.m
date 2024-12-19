%% ECG data processing

clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm483_L_VTA_R_NAc_oxLight';            % Mouse name\
date = 'Apr_05_2024';                             % Date\
stim = 'Fear_Conditioning'; % Sucrose Quinine

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

%% Analyse ECG signals

Sample_Rate = 200;    % 200 scans per second.

vid_start = ceil(step_timepoint(1))*Sample_Rate+1;

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

%% Analyse running signals

running = datas(vid_start:end,1)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
Abs_speedDeg = abs(speedDeg);
speed = movmean(Abs_speedDeg,100);

%% Analyse pupil size

Pupil_up_x = pupil_data(:,1);
Pupil_up_y = pupil_data(:,2);
Pupil_left_x = pupil_data(:,4);
Pupil_left_y = pupil_data(:,5);
Pupil_down_x = pupil_data(:,7);
Pupil_down_y = pupil_data(:,8);
Pupil_right_x = pupil_data(:,10);
Pupil_right_y = pupil_data(:,11);

center_x = zeros(size(pupil_data,1),1);
center_y = zeros(size(pupil_data,1),1);
radii = zeros(size(pupil_data,1),1);
areas = zeros(size(pupil_data,1),1);

for i = 1:size(pupil_data,1)
    X1(1) = Pupil_up_x(i,1);
    X1(2) = Pupil_up_y(i,1);

    X2(1) = Pupil_left_x(i,1);
    X2(2) = Pupil_left_y(i,1);

    Y1(1) = Pupil_down_x(i,1);
    Y1(2) = Pupil_down_y(i,1);

    Y2(1) = Pupil_right_x(i,1);
    Y2(2) = Pupil_right_y(i,1);

    [center_x(i,1), center_y(i,1)]= node(X1,Y1,X2,Y2);

    MinorAxis = [X1(1), X1(2); Y1(1), Y1(2)];
    MinorAxisLength = pdist(MinorAxis, 'euclidean');

    MajorAxis = [X2(1), X2(2); Y2(1), Y2(2)];
    MajorAxisLength = pdist(MajorAxis, 'euclidean');

    diameters = mean([MajorAxisLength, MinorAxisLength]);
    radii(i,1) = diameters/2;
    areas(i,1) = pi*radii(i,1).^2;
end

save([Data_Folder 'Pupil.mat'], 'radii', 'areas', 'center_x', 'center_y');

FrameRate = 20;

medianArea = [];

for i = 1 : total_time
    medianArea(i) = median(areas(((i-1)*FrameRate+1):i*FrameRate));
end

%% Raw ECG signals

ECG_raw = datas(vid_start:end,2)';

%% Remove baseline wandering

fpass=[11 13];

ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

%% find peaks

minPeakPromVal=0.007;
figure;
findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

[pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

%% heart rate in time

RR_intervals = diff(pksLocs);
% heartRate=1./RR_intervals;
% heartRate_bpm=heartRate*60;

RR_intervals_median = NaN(1,total_time);
RR_intervals_std = NaN(1,total_time);

for i = 1 : total_time
    heartbeat_index = find(pksLocs(2:end)>=(i-1) & pksLocs(2:end)<i);
    RR_intervals_trace = RR_intervals(heartbeat_index);
    RR_intervals_median(i) = median(RR_intervals_trace);
    RR_intervals_std(i) = std(RR_intervals_trace);
end

heartRate_median = 1./RR_intervals_median;
heartRate_median_bpm = heartRate_median*60;
heartRate_median_bpm_smooth = movmean(heartRate_median_bpm,5);

% heartRate_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm,0:time(end)*Sample_Rate,'spline','extrap');
% heartRate_interp_smooth = movmean(heartRate_interp,2*Sample_Rate);

%% HRV

% SSD = (diff(RR_intervals)).^2;
% 
% LF = bandpass(heartRate_interp,[0.4 0.8],Sample_Rate);

% RR_intervals_mean = mean(RR_intervals_median);
% SD = abs(RR_intervals_median-RR_intervals_mean);

RR_intervals_std_smooth = movmean(RR_intervals_std,5);

SSD = abs(diff(RR_intervals_median));
SSD_smooth = movmean(SSD,5);

LF = bandpass(heartRate_median,[0.4 0.8],1);

%% Plot running, heartRate, HRV and pupil size

figure
subplot(6,1,1);
plot(time,[speed speed(end)],'k')
xlim([0 time(end)])
% ylim([0 0.15])
title('Running Speed','FontSize',15,'FontWeight','bold')
ylabel('Speed','FontSize',15,'FontWeight','bold')
axis off

subplot(6,1,2)
plot(1:total_time, medianArea,'color',[229 114 190]./255,'LineWidth',2)
xlim([0 time(end)])
% ylim([0 700])
ylabel('pupil area (pixels)','FontSize',15,'FontWeight','bold')
title('Pupil Size','FontSize',15,'FontWeight','bold','color',[229 114 190]./255)
axis off

subplot(6,1,3)
% plot(time,heartRate_interp,'LineWidth',2)
plot(1:total_time, heartRate_median_bpm_smooth, 'LineWidth',2)
xlim([0 time(end)])
ylim([600 800])
% yticks(5:2:14)
title('Heart Rate','FontSize',15,'FontWeight','bold','color',[34 75 160]./255)
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(6,1,4)
plot(2:total_time,SSD_smooth,'color','k','LineWidth',2)
% plot(pksLocs(3:end),SSD,'LineWidth',2)
title('Root Mean Square of Successive Differences','FontSize',15,'FontWeight','bold','color','k')
xlim([0 time(end)])
% ylim([0 0.00001])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(6,1,5)
plot(1:total_time,RR_intervals_std_smooth,'color','k','LineWidth',2)
% plot(pksLocs(3:end),SSD,'LineWidth',2)
title('Standard Deviation','FontSize',15,'FontWeight','bold','color','k')
xlim([0 time(end)])
ylim([0 0.005])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(6,1,6)
plot(1:total_time,LF,'color','k','LineWidth',2)
% plot(time,LF,'color','k','LineWidth',2)
title('LF Signal (0.4-0.8Hz)','FontSize',15,'FontWeight','bold','color','k')
xlim([0 time(end)])
ylim([-0.40 0.40])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')



