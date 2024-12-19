%%

clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm2072';            % Mouse name\
date = 'Nov_09_2024';                             % Date\
stim = 'Fear_Conditioning'; % Sucrose Quinine

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

Cell_Status_data = readtable([Data_Folder mouse_name '_cell_traces.csv'],'Format','auto');
Cell_Status = Cell_Status_data(1,2:end);
Cell_Status_transposed = rows2vars(Cell_Status);
Cell_Status_column = Cell_Status_transposed{:, 2};
cell_status = strcmp(Cell_Status_column, 'accepted');

cell_traces_data = readmatrix([Data_Folder mouse_name '_cell_traces.csv'],'Range',[3 1]);
cell_times = cell_traces_data(:,1);
cell_traces = cell_traces_data(:,2:end);
cell_traces_accepted = cell_traces(:,cell_status);

nanRows = any(isnan(cell_traces_accepted), 2);
cell_times = cell_times(~nanRows, :);
cell_traces_accepted = cell_traces_accepted(~nanRows, :);

cell_traces_accepted_zscore = zscore(cell_traces_accepted);

%% Analyse ECG signals

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

data_onset = find(datas(:,3),1)-50;
data_offset = find(datas(:,3),1,'last')-10;

vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

timepoint = times(data_onset:data_offset,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

smooth_window = 1;

%% Analyse running signals

running = datas(data_onset:data_offset,1)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
Abs_speedDeg = abs(speedDeg);
speed = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);

%% Analyse pupil size

Pupil_up_x = pupil_data(vid_onset:vid_offset,1);
Pupil_up_y = pupil_data(vid_onset:vid_offset,2);
Pupil_left_x = pupil_data(vid_onset:vid_offset,4);
Pupil_left_y = pupil_data(vid_onset:vid_offset,5);
Pupil_down_x = pupil_data(vid_onset:vid_offset,7);
Pupil_down_y = pupil_data(vid_onset:vid_offset,8);
Pupil_right_x = pupil_data(vid_onset:vid_offset,10);
Pupil_right_y = pupil_data(vid_onset:vid_offset,11);

center_x = zeros(size(Pupil_up_x,1),1);
center_y = zeros(size(Pupil_up_x,1),1);
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

    [center_x(i,1), center_y(i,1)]= node(X1,Y1,X2,Y2);

    MinorAxis = [X1(1), X1(2); Y1(1), Y1(2)];
    MinorAxisLength = pdist(MinorAxis, 'euclidean');

    MajorAxis = [X2(1), X2(2); Y2(1), Y2(2)];
    MajorAxisLength = pdist(MajorAxis, 'euclidean');

    diameters = mean([MajorAxisLength, MinorAxisLength]);
    radii(i,1) = diameters/2;
    areas(i,1) = pi*radii(i,1).^2;
end

pupil = filloutliers(areas,"nearest","percentiles",[0 99.99]);

pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

%% Raw ECG signals

ECG_raw = datas(data_onset:data_offset,2)';

%% Remove baseline wandering

fpass=[9.5 14.5];

ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

%% find peaks

minPeakPromVal=0.07;

% figure;
% findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

[pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

%% heart rate in time

RR_intervals = diff(pksLocs);
heartRate=1./RR_intervals;
heartRate_bpm=heartRate*60;

heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",[0 100]);

heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm_outlier,1:size(ECG_raw,2),'nearest','extrap');

heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[5*Sample_Rate 0]);

%%

figure

numCells = size(cell_traces_accepted,2);

for i = 1:numCells

    subplot(numCells, 1, i);
    plot(cell_times,cell_traces_accepted_zscore(:,i),'Color',[56 181 72]./255,'LineWidth',2.5)
    axis off

end

%% Plot running, pupil, heartRate and cell traces

% time = (1/Sample_Rate):(1/Sample_Rate):(size(ECG_raw,2)/Sample_Rate);
video_time = (1/FrameRate):(1/FrameRate):(size(pupil,1)/FrameRate);

figure

subplot(4,1,1);
plot(time,[speed speed(end)],'k')
xlim([0 time(end)])
% ylim([0 0.15])
title('Running Speed','FontSize',15,'FontWeight','bold')
% ylabel('Speed','FontSize',15,'FontWeight','bold')
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(4,1,2)
imagesc(cell_traces_accepted_zscore'); 
colormap('turbo'); 
title('Orx-GCaMP6s','FontSize',15,'FontWeight','bold','color',[64 145 48]./255)
hbar = colorbar('east');
set(hbar, 'Position', [0.1172    0.5572    0.0083    0.1374])
hbar.FontSize = 15;
hbar.Label.String = 'z-score';
hbar.FontWeight = 'bold';
axis off

subplot(4,1,3)
plot(video_time, pupil_smooth,'color',[229 114 190]./255,'LineWidth',2)
xlim([0 time(end)])
% ylim([100 600])
ylabel('pixels','FontSize',15,'FontWeight','bold')
title('Pupil Size','FontSize',15,'FontWeight','bold','color',[229 114 190]./255)
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(4,1,4)
plot(time,heartRate_bpm_interp_smooth,'LineWidth',2)
% plot(1:total_time, heartRate_median_bpm_smooth, 'LineWidth',2)
xlim([0 time(end)])
% ylim([480 600])
% yticks(5:2:14)
ylabel('bpm','FontSize',15,'FontWeight','bold')
title('Heart Rate','FontSize',15,'FontWeight','bold','color',[34 75 160]./255)
% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
