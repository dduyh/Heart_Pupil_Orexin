%% Fiber photometry data processing

clear; close all; clc;

%% set the path for output data

Directory = 'E:\data\';                     % Main directory\
mouse_name = 'm487_R_NAc_nLightR';            % Mouse name\
date = 'Dec_09_2023';                             % Date\
session = 'session_2';
stim = 'Quinine'; % Sucrose Quinine

single_plots = 0;

Data_Folder = [Directory mouse_name '\' date '\' session '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

pupil_data = csvread([Data_Folder mouse_name '_' date '_' session '_' stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

%% Analyse photometry signals

Sample_Rate = 200; % 200 scans per second.
vid_start = ceil(step_timepoint(1))*Sample_Rate+1;

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

raw = datas(vid_start:end,4)';

if single_plots
    % Plot raw data
    figure(1);
    subplot(4,1,1)
    plot(time,raw,'b')
    title('Raw signals (m487, Dec 09 2023, session 2, Quinine)','FontSize',10,'FontWeight','bold')
end

% Smooth data

smooth_win = 10;
smooth = movmean(raw,smooth_win);

% Remove bleaching slope

F = @(x,t)x(1)*exp(-t/x(2)) + x(3)*exp(-t/(x(4)*x(2))) + x(5);

x0 = [max(smooth)/4, 3600, max(smooth)/4, 0.1, max(smooth)/2] ;
lb = [0, 600, 0, 0, 0];
ub = [max(smooth), 36000, max(smooth), 1, max(smooth)];
x = lsqcurvefit(F, x0, time, smooth, lb, ub);
base = F(x,time);
ch1 = smooth - base;

if single_plots
    % Plot smoothed signals with bleaching slopes
    figure(1)
    subplot(4,1,2)
    hold on
    title('Smoothed signals with double exponential bleaching fits','FontSize',10,'FontWeight','bold')
    plot(time, smooth, 'color', [231 208 84]./255)
    plot(time, base, 'k-')
    hold off
    
    % Plot bleaching baseline substracted signals
    figure(1)
    subplot(4,1,3)
    plot(time,ch1,'color',[231 208 84]./255)
    legend('560 nm','FontSize',10,'TextColor',[231 208 84]./255)
    legend('boxoff')
    title('Bleaching Correction by Double Exponential Fit','FontSize',10,'FontWeight','bold')
end

% Standardize signals

ch1_zscored = (ch1 - mean(ch1)) / std(ch1);

if single_plots
    % Plot normalized signals
    figure(1)
    subplot(4,1,4)
    plot(time,ch1_zscored,'color',[255 128 128]./255)
    title('nLightR NAc R','FontSize',10,'FontWeight','bold')
    ylabel('z-score','FontSize',10,'FontWeight','bold')
    xlabel('Time (seconds)','FontSize',10,'FontWeight','bold')
end

%% Analyse licking and running signals

licking = datas(vid_start:end,1)';

running = datas(vid_start:end,2)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
speed = movmean(speedDeg,100);

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

%% Plot running, licking, photometry signals and pupil size 

figure
subplot(4,1,1);
plot(time,licking,'b')
xlim([0 time(end)])
ylim([-0.5 1.5])
title('Licking','FontSize',15,'FontWeight','bold')

subplot(4,1,2);
plot(time,[speed speed(end)],'k')
xlim([0 time(end)])
% ylim([-0.3 0.2])
title('Running','FontSize',15,'FontWeight','bold')
ylabel('Speed','FontSize',15,'FontWeight','bold')

subplot(4,1,3)
plot(time,ch1_zscored,'color',[255 128 128]./255)
xlim([0 time(end)])
title('nLightR NAc R','FontSize',15,'FontWeight','bold')
ylabel('z-score','FontSize',15,'FontWeight','bold')

subplot(4,1,4)
fill([1:total_time fliplr(1:total_time)],[medianArea+stdArea fliplr(medianArea-stdArea)], [0.8 0.8 0.8])
hold on
plot(1:total_time, medianArea, 'LineWidth',2)
xlim([0 time(end)])
% ylim([1500 4000])
ylabel('pupil area (pixels)','FontSize',15,'FontWeight','bold')
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
title('Pupil Size','FontSize',15,'FontWeight','bold')



