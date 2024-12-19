%% ECG data processing

clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm1772';            % Mouse name\
date = 'Mar_07_2024\session_1';                             % Date\
stim = 'heart'; % Sucrose Quinine

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

%% Analyse ECG signals

Sample_Rate = 1000;    % 200 scans per second.

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
stdArea = [];

for i = 1 : total_time
medianArea(i) = median(areas(((i-1)*FrameRate+1):i*FrameRate));
stdArea(i) = std(areas(((i-1)*FrameRate+1):i*FrameRate));
end

%% Plot running, ECG signals and pupil size 

ECG_raw = datas(vid_start:end,2)';

figure(1)
subplot(3,1,1);
plot(time,[speed speed(end)],'k')
xlim([0 time(end)])
% ylim([0 0.15])
title('Running Speed','FontSize',15,'FontWeight','bold')
ylabel('Speed','FontSize',15,'FontWeight','bold')
axis off

subplot(3,1,2)
plot(1:total_time, medianArea,'color',[229 114 190]./255,'LineWidth',2)
xlim([0 time(end)])
ylim([0 900])
ylabel('pupil area (pixels)','FontSize',15,'FontWeight','bold')
title('Pupil Size','FontSize',15,'FontWeight','bold','color',[229 114 190]./255)
axis off

subplot(3,1,3)
plot(time,ECG_raw,'color',[34 75 160]./255)
xlim([0 time(end)])
% ylim([1.7 2.3])
title('ECG Raw Signal','FontSize',15,'FontWeight','bold','color',[34 75 160]./255)
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
axis off



