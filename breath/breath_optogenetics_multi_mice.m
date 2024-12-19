clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

breath_folder = {'m1825\May_23_2024';
    'm1826\May_23_2024';
    'm1821\May_24_2024_2';
    'm1822\May_24_2024'};

ctrl_folder = {'m1840\May_23_2024';
    'm1821\May_24_2024'};

FrameRate = 20;

stimFreqs = 20;

trace_duration = 500;   % 90 seconds.

trials_num = size(breath_folder,1);

ctrl_trials_num = size(ctrl_folder,1);

breath = NaN(trials_num,trace_duration*FrameRate);

breath_ctrl = NaN(ctrl_trials_num,trace_duration*FrameRate);

MinPeakHeights=[1 1.5 1.4 1.5];

MinPeakHeights_ctrl=[0.5 0.5];

outlier_breath=[0.5 99.5;
    0.5 99.5;
    0.5 99.5;
    0 99];

outlier_ctrl=[0.5 99.5;
    0.5 99.5];

xlims = (-100*FrameRate+1:400*FrameRate)/FrameRate;

for I=1:size(breath_folder,1)

    Data_Folder = [Directory breath_folder{I} '\'];

    breath_data = csvread([Data_Folder 'Values.csv'],1,1);

    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs);

    %% Analyse breathing rate

    pixDiff = diff(breath_data);

    MinPeakHeight=MinPeakHeights(I);

    [pksVal, pksLocs]=findpeaks(pixDiff,FrameRate,'MaxPeakWidth',0.2,'MinPeakHeight',MinPeakHeight,'Threshold',1e-4,'WidthReference','halfheight');

    BB_intervals = diff(pksLocs);
    breathRate=1./BB_intervals;

    breath_outlier = filloutliers(breathRate,"nearest","percentiles",outlier_breath(I,:));

    breathRate_interp = interp1(pksLocs(2:end)*FrameRate,breath_outlier,1:size(breath_data,1),'spline','extrap');
    breathRate_interp_smooth = movmean(breathRate_interp,FrameRate*5);

    %%

    breath_trials = NaN(length(idxFreqs),trace_duration*FrameRate);

    for II = 1:length(laser_onsets)

        breath_start = round(laser_onsets(II)*FrameRate) - 100*FrameRate +1;
        breath_end = round(laser_onsets(II)*FrameRate) + 400*FrameRate;

        raw_breath = breathRate_interp_smooth(breath_start:breath_end);

%         breath_zscored = (raw_breath - mean(raw_breath(1:100*FrameRate))) / std(raw_breath(1:100*FrameRate));

        breath_trials(II,:) = raw_breath;

    end

    figure(1);
    hold on
    plot(xlims,mean(breath_trials),'LineWidth',1)

    breath(I,:) = mean(breath_trials);

end

for I=1:size(ctrl_folder,1)

    Data_Folder = [Directory ctrl_folder{I} '\'];

    breath_data = csvread([Data_Folder 'Values.csv'],1,1);

    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs);

    %% Analyse breathing rate

    pixDiff = diff(breath_data);

    MinPeakHeight=MinPeakHeights_ctrl(I);

    [pksVal, pksLocs]=findpeaks(pixDiff,FrameRate,'MaxPeakWidth',0.2,'MinPeakHeight',MinPeakHeight,'Threshold',1e-4,'WidthReference','halfheight');

    BB_intervals = diff(pksLocs);
    breathRate=1./BB_intervals;

    breath_outlier = filloutliers(breathRate,"nearest","percentiles",outlier_ctrl(I,:));

    breathRate_interp = interp1(pksLocs(2:end)*FrameRate,breath_outlier,1:size(breath_data,1),'spline','extrap');
    breathRate_interp_smooth = movmean(breathRate_interp,FrameRate*5);

    %%

    breath_trials = NaN(length(idxFreqs),trace_duration*FrameRate);

    for II = 1:length(laser_onsets)

        breath_start = round(laser_onsets(II)*FrameRate) - 100*FrameRate +1;
        breath_end = round(laser_onsets(II)*FrameRate) + 400*FrameRate;

        raw_breath = breathRate_interp_smooth(breath_start:breath_end);

%         breath_zscored = (raw_breath - mean(raw_breath(1:100*FrameRate))) / std(raw_breath(1:100*FrameRate));

        breath_trials(II,:) = raw_breath;

    end

    figure(2);
    hold on
    plot(xlims,mean(breath_trials),'LineWidth',1)

    breath_ctrl(I,:) = mean(breath_trials);

end

%%

figure(1);

patch('XData',[0, 0, 60, 60],'YData',[0.6, 1.6, 1.6, 0.6],'EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

plot(xlims,mean(breath),'Color','k','LineWidth',2)
title('Respiratory Rate (Chrimson)','FontSize',20,'FontWeight','bold')
ylabel('Hz','FontSize',15,'FontWeight','bold');

xlim([-100 400])
ylim([0.6 1.6])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


figure(2);

patch('XData',[0, 0, 60, 60],'YData',[0.6, 1.6, 1.6, 0.6],'EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

plot(xlims,mean(breath_ctrl),'Color','k','LineWidth',2)
title('Respiratory Rate (Control)','FontSize',20,'FontWeight','bold')
ylabel('Hz','FontSize',15,'FontWeight','bold');

xlim([-100 400])
ylim([0.6 1.6])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

