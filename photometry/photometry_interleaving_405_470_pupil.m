%% Fiber photometry data processing

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';                     % Main directory\
mouse_name = 'm758_L_MCH_R_Orx_GCaMP';            % Mouse name\
date = 'Dec_10_2023';                             % Date\
session = 'session_3';
stim = 'Quinine'; % Sucrose Quinine

Data_Folder = [Directory mouse_name '\' date '\' session '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

pupil_data = csvread([Data_Folder mouse_name '_' date '_' session '_' stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

%% Analyse photometry signals

Sample_Rate = 200; % 200 scans per second.
vid_start = ceil(step_timepoint(1))*Sample_Rate+2;

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

raw_ch1 = datas(vid_start:end,3)';
raw_ch2 = datas(vid_start:end,4)';

% Seperate two channels
index_470 = find(~repmat([ones(1,10) zeros(1,10)],1,10*total_time));
index_405 = find(repmat([ones(1,10) zeros(1,10)],1,10*total_time));

raw_470_ch1 = raw_ch1(index_470);
raw_405_ch1 = raw_ch1(index_405);
raw_470_ch2 = raw_ch2(index_470);
raw_405_ch2 = raw_ch2(index_405);

time_470 = time(index_470);
time_405 = time(index_405);

% Smooth data
smooth_win = 10;
smooth_470_ch1 = movmean(raw_470_ch1,smooth_win);
smooth_405_ch1 = movmean(raw_405_ch1,smooth_win);
smooth_470_ch2 = movmean(raw_470_ch2,smooth_win);
smooth_405_ch2 = movmean(raw_405_ch2,smooth_win);

% Remove bleaching slope
F = @(x,t)x(1)*exp(-t/x(2)) + x(3)*exp(-t/(x(4)*x(2))) + x(5);

x0_470_ch1 = [max(smooth_470_ch1)/4, 3600, max(smooth_470_ch1)/4, 0.1, max(smooth_470_ch1)/2] ;
lb_470_ch1 = [0, 600, 0, 0, 0];
ub_470_ch1 = [max(smooth_470_ch1), 36000, max(smooth_470_ch1), 1, max(smooth_470_ch1)];
x_470_ch1 = lsqcurvefit(F, x0_470_ch1, time_470, smooth_470_ch1, lb_470_ch1, ub_470_ch1);
base_470_ch1 = F(x_470_ch1,time_470);
ch1_470 = smooth_470_ch1 - base_470_ch1;

x0_405_ch1 = [max(smooth_405_ch1)/4, 3600, max(smooth_405_ch1)/4, 0.1, max(smooth_405_ch1)/2] ;
lb_405_ch1 = [0, 600, 0, 0, 0];
ub_405_ch1 = [max(smooth_405_ch1), 36000, max(smooth_405_ch1), 1, max(smooth_405_ch1)];
x_405_ch1 = lsqcurvefit(F, x0_405_ch1, time_405, smooth_405_ch1, lb_405_ch1, ub_405_ch1);
base_405_ch1 = F(x_405_ch1,time_405);
ch1_405 = smooth_405_ch1 - base_405_ch1;

x0_470_ch2 = [max(smooth_470_ch2)/4, 3600, max(smooth_470_ch2)/4, 0.1, max(smooth_470_ch2)/2] ;
lb_470_ch2 = [0, 600, 0, 0, 0];
ub_470_ch2 = [max(smooth_470_ch2), 36000, max(smooth_470_ch2), 1, max(smooth_470_ch2)];
x_470_ch2 = lsqcurvefit(F, x0_470_ch2, time_470, smooth_470_ch2, lb_470_ch2, ub_470_ch2);
base_470_ch2 = F(x_470_ch2,time_470);
ch2_470 = smooth_470_ch2 - base_470_ch2;

x0_405_ch2 = [max(smooth_405_ch2)/4, 3600, max(smooth_405_ch2)/4, 0.1, max(smooth_405_ch2)/2] ;
lb_405_ch2 = [0, 600, 0, 0, 0];
ub_405_ch2 = [max(smooth_405_ch2), 36000, max(smooth_405_ch2), 1, max(smooth_405_ch2)];
x_405_ch2 = lsqcurvefit(F, x0_405_ch2, time_405, smooth_405_ch2, lb_405_ch2, ub_405_ch2);
base_405_ch2 = F(x_405_ch2,time_405);
ch2_405 = smooth_405_ch2 - base_405_ch2;

% Calculate difference between 470nm and 405nm signals
fitdata_ch1 = fit(ch1_405',ch1_470',fittype('poly1'),'Robust','on');
fitdata_ch2 = fit(ch2_405',ch2_470',fittype('poly1'),'Robust','on');

if fitdata_ch1.p1 > 0
    ch1_405 = fitdata_ch1(ch1_405)';
    ch1 = ch1_470 - ch1_405;
else
    ch1 = ch1_470;
end

if fitdata_ch2.p1 > 0
    ch2_405 = fitdata_ch2(ch2_405)';
    ch2 = ch2_470 - ch2_405;
else
    ch2 = ch2_470;
end

% Standardize signals
ch1_zscored = (ch1 - mean(ch1)) / std(ch1);
ch2_zscored = (ch2 - mean(ch2)) / std(ch2);

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

figure(1)
subplot(5,1,1);
plot(time,licking,'b')
xlim([0 time(end)])
ylim([-0.5 1.5])
title('Licking','FontSize',15,'FontWeight','bold')

subplot(5,1,2);
plot(time,[speed speed(end)],'k')
xlim([0 time(end)])
% ylim([-0.4 0.4])
title('Running','FontSize',15,'FontWeight','bold')
ylabel('Speed','FontSize',15,'FontWeight','bold')

subplot(5,1,3)
plot(time_470,ch1_zscored,'color',[56 181 72]./255)
xlim([0 time(end)])
title('MCH-GCaMP6s L','FontSize',15,'FontWeight','bold')
ylabel('z-score','FontSize',15,'FontWeight','bold')

subplot(5,1,4)
plot(time_470,ch2_zscored,'color',[56 181 72]./255)
xlim([0 time(end)])
title('Orx-GCaMP6s R','FontSize',15,'FontWeight','bold')
ylabel('z-score','FontSize',15,'FontWeight','bold')

subplot(5,1,5)
fill([1:total_time fliplr(1:total_time)],[medianArea+stdArea fliplr(medianArea-stdArea)], [0.8 0.8 0.8])
hold on
plot(1:total_time, medianArea, 'LineWidth',2)
xlim([0 time(end)])
% ylim([0 1500])
ylabel('pupil area (pixels)','FontSize',15,'FontWeight','bold')
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
title('Pupil Size','FontSize',15,'FontWeight','bold')



