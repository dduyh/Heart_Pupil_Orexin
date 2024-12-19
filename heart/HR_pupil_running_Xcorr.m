%% 

clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm1826';            % Mouse name\
date = 'May_31_2024';                             % Date\
stim = 'Fear_Conditioning'; % Sucrose Quinine

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'step_timepoint.mat']);

pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

%% 

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

data_start = round(step_timepoint(1)*FrameRate)+1;
data_end = round(step_timepoint(1)*FrameRate + size(datas,1)/Sample_Rate*FrameRate);

smooth_window = 1;

%% Analyse running signals

running = datas(:,1)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
Abs_speedDeg = abs(speedDeg);
Abs_speedDeg_outlier = filloutliers(Abs_speedDeg,"nearest","percentiles",[0 99.9]);
speed = [Abs_speedDeg_outlier(1) Abs_speedDeg_outlier];
speed_smooth = movmean(speed,[smooth_window*Sample_Rate 0]);

%% Analyse pupil size

Pupil_up_x = pupil_data(data_start:data_end,1);
Pupil_up_y = pupil_data(data_start:data_end,2);
Pupil_left_x = pupil_data(data_start:data_end,4);
Pupil_left_y = pupil_data(data_start:data_end,5);
Pupil_down_x = pupil_data(data_start:data_end,7);
Pupil_down_y = pupil_data(data_start:data_end,8);
Pupil_right_x = pupil_data(data_start:data_end,10);
Pupil_right_y = pupil_data(data_start:data_end,11);

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

    MinorAxis = [X1(1), X1(2); Y1(1), Y1(2)];
    MinorAxisLength = pdist(MinorAxis, 'euclidean');

    MajorAxis = [X2(1), X2(2); Y2(1), Y2(2)];
    MajorAxisLength = pdist(MajorAxis, 'euclidean');

    diameters = mean([MajorAxisLength, MinorAxisLength]);
    radii(i,1) = diameters/2;
    areas(i,1) = pi*radii(i,1).^2;
end

pupil = filloutliers(areas,"nearest","percentiles",[0.1 99.9]);

video_time = (1/FrameRate):(1/FrameRate):(size(pupil,1)/FrameRate);

pupil_interp = interp1(video_time*Sample_Rate,pupil,1:size(datas,1),'nearest','extrap');

pupil_interp_smooth = movmean(pupil_interp,[smooth_window*Sample_Rate 0]);

%% Analyse Heart Rate

% Raw ECG signals

ECG_raw = datas(:,2)';

% Remove baseline wandering

fpass=[8 12];

ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

% find peaks

minPeakPromVal=0.007;

[pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

% heart rate in time

RR_intervals = diff(pksLocs);
heartRate=1./RR_intervals;
heartRate_bpm=heartRate*60;

heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",[0 100]);

heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm_outlier,1:size(datas,1),'nearest','extrap');

heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*Sample_Rate 0]);

%% Plot pulse, running, pupil size and heartRate

time = (1/Sample_Rate):(1/Sample_Rate):(size(datas,1)/Sample_Rate);

pulse = datas(:,3)';

figure

subplot(4,1,1);
plot(time,pulse,'k')
xlim([0 time(end)])
% ylim([0 0.2])
title('Laser Pulse','FontSize',15,'FontWeight','bold')
ylabel('Speed','FontSize',15,'FontWeight','bold')
axis off

subplot(4,1,2);
plot(time,speed_smooth,'k')
xlim([0 time(end)])
% ylim([0 0.12])
title('Running Speed','FontSize',15,'FontWeight','bold')
ylabel('Speed','FontSize',15,'FontWeight','bold')
axis off

subplot(4,1,3)
plot(time, pupil_interp_smooth,'color',[229 114 190]./255,'LineWidth',2)
xlim([0 time(end)])
% ylim([0 1500])
ylabel('pupil area (pixels)','FontSize',15,'FontWeight','bold')
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
% ylim([520 600])
% yticks(5:2:14)
title('Heart Rate','FontSize',15,'FontWeight','bold','color',[34 75 160]./255)
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

%%

xcorrwind = 5*Sample_Rate; % in samples

[r,lags] = xcorr(speed-mean(speed),heartRate_bpm_interp-mean(heartRate_bpm_interp),xcorrwind,'coeff');
figure
stem(lags,r)

[r,lags] = xcorr(speed-mean(speed),pupil_interp-mean(pupil_interp),xcorrwind,'coeff');
figure
stem(lags,r)


% % do pupil-orexin xcorrs for all sessions
% [c_orxpup(a,:),lags] = xcorr(currOrx,currPup,xcorrwind,'coeff');
% [maxcorr_puporxbest(a),l] = max(c_orxpup(a,:));
% maxlag_puporxbest(a) = lags(l);
% nMice = size(maxcorr_puporxbest,1);
% 
% % xcorr plot orx-pup mean + sds
% mean_orxpup_xcorr = mean(c_orxpup,1);
% sd_orxpup_xcorr = std(c_orxpup,0,1)/sqrt(nMice);
% 
% close all
% f = figure(1)
% subplot(1,3,1)
% s = shadedErrorBar((-xcorrwind:xcorrwind)/20,mean_orxpup_xcorr, [sd_orxpup_xcorr],'PatchSaturation',0.7,'lineprops', 'w')
% s.patch.FaceColor = (orxColor + pupColor)/2;hold on
% plot((-xcorrwind:xcorrwind)/20,mean_orxpup_xcorr,'Color',(orxColor + pupColor)/2)
% hold on; xline(0), title('Orexin-Pupil X-corr')
% ylim([0,1])
