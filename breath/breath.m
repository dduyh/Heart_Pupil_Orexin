%% 

clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm1822';            % Mouse name\
date = 'May_24_2024';                             % Date\

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% 

load([Data_Folder 'step_timepoint.mat']);

data = readmatrix([Data_Folder 'Values.csv']);

FrameRate = 20;

data_start = ceil(step_timepoint(1)*FrameRate)+1;

breath_data = data(data_start:end,2);

pixDiff = diff(breath_data);

% plot((2:size(A,1))/20,pixDiff')

breath_Bandpass = lowpass(pixDiff,1.5,FrameRate);

% findpeaks(breath_Bandpass,FrameRate,'MaxPeakWidth',0.5,'MinPeakProminence',0.2);
[pksVal, pksLocs]=findpeaks(breath_Bandpass,FrameRate,'MaxPeakWidth',0.5,'MinPeakProminence',0.2);

RR_intervals = diff(pksLocs);
breathRate=1./RR_intervals;

breath_outlier = filloutliers(breathRate,"nearest","percentiles",[5 95]);

breathRate_interp = interp1(pksLocs(2:end)*FrameRate,breath_outlier,1:size(breath_data,1),'spline','extrap');
breathRate_interp_smooth = movmean(breathRate_interp,FrameRate*30);

%%
figure

time = (1/FrameRate):(1/FrameRate):(size(breath_data,1)/FrameRate);
plot(time,breathRate_interp_smooth)
xlim([0 time(end)])

% %%
% figure
% plot((2:size(A,1))/20,movsum(movmean(abs(pixDiff)>0.2,10)>0.05,[1200 0]))


