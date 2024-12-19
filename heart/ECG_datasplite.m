%% ECG data splite

clear; close all; clc;

%% set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm1746';            % Mouse name\
date = 'Jun_04_2024';                             % Date\

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'step_timepoint.mat']);

%%

Sample_Rate = 1000;    % 200 scans per second.

vid_start = ceil(step_timepoint(1))*Sample_Rate+1;
% vid_start = 1;

EEG = datas(vid_start:end,2);

running = datas(vid_start:end,1);
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
speed = movmean(speedDeg,100);

EMG = [speed; speed(end)];

%% save data files

save([Data_Folder 'ECG.mat'],'EEG','EMG');

%%

fig1 = figure(1);
fig2 = figure(2);
copyobj(fig1.Children, fig2)
