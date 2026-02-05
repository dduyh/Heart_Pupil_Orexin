%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

folder = {'m1821\May_30_2024';
    'm1822\May_30_2024';
    'm1825\May_31_2024';
    'm1826\May_31_2024';
    'm2140\Mar_04_2025';
    'm2144\Mar_04_2025';
    'm2145\Mar_04_2025';
    'm2146\Mar_04_2025';
    'm2169\Mar_05_2025';
    'm2175\Mar_05_2025';
    'm2178\Mar_05_2025';
    'm2132\Mar_06_2025';
    %     'm2133\Mar_06_2025';
    'm2148\Mar_06_2025';
    'm2179\Mar_06_2025';
    'm839\Aug_17_2025';
    'm840\Aug_17_2025';
    'm843\Aug_22_2025';
    'm844\Aug_22_2025';
    'm833\Aug_22_2025';
    'm835\Aug_22_2025';
    'm836\Aug_22_2025'};

fpass_trials=[9 12;
    9 13;
    9.5 12;
    8 12;
    9 13;
    9 13;
    8 13;
    8 13;
    8 13;
    7 13;
    7 12;
    8 13;
    %     7 13;
    7 13;
    8 13;
    9 14;
    7 13.5;
    9 13;
    9.5 14;
    9 12.5;
    8.5 13;
    9.5 13];

outlier_pupil=[0 99;
    0 99;
    0 99;
    0 99;
    0 100;
    0 99.99;
    0 99.98;
    0 99.98;
    0 100;
    0 99.98;
    0 100;
    0 100;
    %     0 100;
    0 100;
    0 100;
    0 99.52;
    0 100;
    0 99.99;
    0 100;
    0 99.99;
    0 99.99;
    0 100];

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

pre_bout_duration = 100;   % 5 seconds.

trace_duration = 500;   % 30 seconds.

wheel_radius = 10;     % cm

unit = (pi/180)*wheel_radius*Sample_Rate;

smooth_window = 3;

Run_20Hz = [];
Run_5Hz = [];

Pupil_20Hz = [];
Pupil_5Hz = [];

HR_20Hz = [];
HR_5Hz = [];

HRV_20Hz = [];
HRV_5Hz = [];

RMSSD_on_20Hz = [];
RMSSD_off_20Hz = [];
RMSSD_on_5Hz = [];
RMSSD_off_5Hz = [];

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    if I>14
        pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingJul16shuffle1_1000000_filtered.csv'],3,1);
    else
        pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);
    end

    data_start = ceil(step_timepoint(1)*FrameRate)+1;
    data_end = ceil(step_timepoint(1)*FrameRate+size(datas,1)*FrameRate/Sample_Rate);

    frame_time = (1/FrameRate):(1/FrameRate):size(datas,1)/Sample_Rate;

    %% Analyse running signals

    running = datas(:,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    Abs_speedDeg_outlier = filloutliers(Abs_speedDeg,"nearest","percentiles",[0 99.99]);
    speed_smooth = movmean(Abs_speedDeg_outlier,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);
    speed_smooth_resampled_cm = speed_smooth_resampled*unit;

    %% Analyse pupil size

    Pupil_up_x = pupil_data(data_start:data_end,1);
    Pupil_up_y = pupil_data(data_start:data_end,2);
    Pupil_left_x = pupil_data(data_start:data_end,4);
    Pupil_left_y = pupil_data(data_start:data_end,5);
    Pupil_down_x = pupil_data(data_start:data_end,7);
    Pupil_down_y = pupil_data(data_start:data_end,8);
    Pupil_right_x = pupil_data(data_start:data_end,10);
    Pupil_right_y = pupil_data(data_start:data_end,11);

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

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    if I>4
        minPeakPromVal=0.07;
    else
        minPeakPromVal=0.007;
    end

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %% HRV

    SD = abs(diff(RR_intervals)).*1000;
    SD_interp = interp1(pksLocs(3:end),SD,frame_time,'nearest','extrap');
    SD_interp_smooth = movmean(SD_interp,[smooth_window*FrameRate 0]);

    %% Optogenetic stimulation 20 Hz bouts

    run_trials = [];
    pupil_trials = [];
    HR_trials = [];
    HRV_trials = [];
    laser_on_RMSSD_trials = [];
    laser_off_RMSSD_trials = [];

    idxFreqs = find(stimSchedule == 20);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs)-step_timepoint(1);

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_trials = [run_trials; speed_smooth_resampled_cm(t_start:t_end)];

        raw_pupil = pupil_smooth(t_start:t_end)';
        pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
        pupil_trials = [pupil_trials; raw_pupil];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR_trials = [HR_trials; raw_HR];

        raw_HRV = SD_interp_smooth(t_start:t_end);
        HRV_trials = [HRV_trials; raw_HRV];

        baseline_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)-pre_bout_duration) & pksLocs(2:end)<laser_onsets(II));
        baseline_RR_intervals = RR_intervals(baseline_heartbeat_index);
        baseline_RMSSD = sqrt(mean((diff(baseline_RR_intervals)).^2));

        laser_on_heartbeat_index = find(pksLocs(2:end)>=laser_onsets(II) & pksLocs(2:end)<(laser_onsets(II)+60));
        laser_on_RR_intervals = RR_intervals(laser_on_heartbeat_index);
        laser_on_RMSSD_trials = [laser_on_RMSSD_trials; sqrt(mean((diff(laser_on_RR_intervals)).^2))-baseline_RMSSD];

        laser_off_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)+60) & pksLocs(2:end)<(laser_onsets(II)+400));
        laser_off_RR_intervals = RR_intervals(laser_off_heartbeat_index);
        laser_off_RMSSD_trials = [laser_off_RMSSD_trials; sqrt(mean((diff(laser_off_RR_intervals)).^2))-baseline_RMSSD];

    end

    Run_20Hz = [Run_20Hz; mean(run_trials)];

    Pupil_20Hz = [Pupil_20Hz; mean(pupil_trials)];

    HR_20Hz = [HR_20Hz; mean(HR_trials)];

    HRV_20Hz = [HRV_20Hz; mean(HRV_trials)];

    RMSSD_on_20Hz = [RMSSD_on_20Hz; mean(laser_on_RMSSD_trials)];
    RMSSD_off_20Hz = [RMSSD_off_20Hz; mean(laser_off_RMSSD_trials)];

    %% Optogenetic stimulation 5 Hz bouts

    run_trials = [];
    pupil_trials = [];
    HR_trials = [];
    HRV_trials = [];
    laser_on_RMSSD_trials = [];
    laser_off_RMSSD_trials = [];

    idxFreqs = find(stimSchedule == 5);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs)-step_timepoint(1);

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        run_trials = [run_trials; speed_smooth_resampled_cm(t_start:t_end)];

        raw_pupil = pupil_smooth(t_start:t_end)';
        pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
        pupil_trials = [pupil_trials; raw_pupil];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR_trials = [HR_trials; raw_HR];

        raw_HRV = SD_interp_smooth(t_start:t_end);
        HRV_trials = [HRV_trials; raw_HRV];

        baseline_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)-pre_bout_duration) & pksLocs(2:end)<laser_onsets(II));
        baseline_RR_intervals = RR_intervals(baseline_heartbeat_index);
        baseline_RMSSD = sqrt(mean((diff(baseline_RR_intervals)).^2));

        laser_on_heartbeat_index = find(pksLocs(2:end)>=laser_onsets(II) & pksLocs(2:end)<(laser_onsets(II)+60));
        laser_on_RR_intervals = RR_intervals(laser_on_heartbeat_index);
        laser_on_RMSSD_trials = [laser_on_RMSSD_trials; sqrt(mean((diff(laser_on_RR_intervals)).^2))-baseline_RMSSD];

        laser_off_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)+60) & pksLocs(2:end)<(laser_onsets(II)+400));
        laser_off_RR_intervals = RR_intervals(laser_off_heartbeat_index);
        laser_off_RMSSD_trials = [laser_off_RMSSD_trials; sqrt(mean((diff(laser_off_RR_intervals)).^2))-baseline_RMSSD];

    end

    Run_5Hz = [Run_5Hz; mean(run_trials)];

    Pupil_5Hz = [Pupil_5Hz; mean(pupil_trials)];

    HR_5Hz = [HR_5Hz; mean(HR_trials)];

    HRV_5Hz = [HRV_5Hz; mean(HRV_trials)];

    RMSSD_on_5Hz = [RMSSD_on_5Hz; mean(laser_on_RMSSD_trials)];
    RMSSD_off_5Hz = [RMSSD_off_5Hz; mean(laser_off_RMSSD_trials)];

end

%%

xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

trials_num = size(Run_20Hz,1);

fig = figure(1);
set(fig, 'Position', [2561 49 922 1315]);

subplot(4,2,1);
hold on
patch('XData',[0, 0, 60, 60],'YData',[0, 8, 8, 0],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_5Hz)+std(Run_5Hz)/sqrt(trials_num) fliplr(mean(Run_5Hz)-std(Run_5Hz)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.7);
plot1 = plot(xlims,mean(Run_5Hz),'Color',[0.5 0.5 0.5],'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run_20Hz)+std(Run_20Hz)/sqrt(trials_num) fliplr(mean(Run_20Hz)-std(Run_20Hz)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.7);
plot2 = plot(xlims,mean(Run_20Hz),'Color','k','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([0 5])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Running Speed','FontSize',15,'FontWeight','bold')
ylabel('cm/s','FontSize',12,'FontWeight','bold');
legend([plot2, plot1], {'20Hz', '5Hz'},'FontSize',12,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

subplot(4,2,3)
hold on
patch('XData',[0, 0, 60, 60],'YData',[300, 900, 900, 300],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_5Hz)+std(Pupil_5Hz)/sqrt(trials_num) fliplr(mean(Pupil_5Hz)-std(Pupil_5Hz)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.7);
plot1 = plot(xlims,mean(Pupil_5Hz),'Color',[139 92 158]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil_20Hz)+std(Pupil_20Hz)/sqrt(trials_num) fliplr(mean(Pupil_20Hz)-std(Pupil_20Hz)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot2 = plot(xlims,mean(Pupil_20Hz),'Color',[229 114 190]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([300 900])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Pupil Size','FontSize',15,'FontWeight','bold')
ylabel('pixels','FontSize',12,'FontWeight','bold');
legend([plot2, plot1], {'20Hz', '5Hz'},'FontSize',12,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

subplot(4,2,5)
hold on
patch('XData',[0, 0, 60, 60],'YData',[580, 720, 720, 580],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_5Hz)+std(HR_5Hz)/sqrt(trials_num) fliplr(mean(HR_5Hz)-std(HR_5Hz)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255, 214, 170]./255,'FaceAlpha',0.7);
plot1 = plot(xlims,mean(HR_5Hz),'Color',[250, 157, 86]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR_20Hz)+std(HR_20Hz)/sqrt(trials_num) fliplr(mean(HR_20Hz)-std(HR_20Hz)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot2 = plot(xlims,mean(HR_20Hz),'Color',[255 128 128]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([600 720])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Heart Rate','FontSize',15,'FontWeight','bold')
ylabel('bpm','FontSize',12,'FontWeight','bold');
legend([plot2, plot1], {'20Hz', '5Hz'},'FontSize',12,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

subplot(4,2,7);
hold on
patch('XData',[0, 0, 60, 60],'YData',[0, 20, 20, 0],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HRV_5Hz)+std(HRV_5Hz)/sqrt(trials_num) fliplr(mean(HRV_5Hz)-std(HRV_5Hz)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255, 214, 170]./255,'FaceAlpha',0.7);
plot1 = plot(xlims,mean(HRV_5Hz),'Color',[250, 157, 86]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HRV_20Hz)+std(HRV_20Hz)/sqrt(trials_num) fliplr(mean(HRV_20Hz)-std(HRV_20Hz)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot2 = plot(xlims,mean(HRV_20Hz),'Color',[255 128 128]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([2 15])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Heart Rate Variability','FontSize',15,'FontWeight','bold')
ylabel('SD (ms)','FontSize',12,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold')
legend([plot2, plot1], {'20Hz', '5Hz'},'FontSize',12,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

%%

peak_Run_20Hz = mean(Run_20Hz(:,100*FrameRate+1:160*FrameRate),2) - mean(Run_20Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_Run_5Hz = mean(Run_5Hz(:,100*FrameRate+1:160*FrameRate),2) - mean(Run_5Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_Run = [peak_Run_20Hz peak_Run_5Hz];

figure(1);

axes('Position', [0.5703 0.7673 0.15 0.1561]);
hold on

for k = 1:size(peak_Run,1)
    plot([1.3 1.7],peak_Run(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_Run_20Hz,1)
    plot(1.3,peak_Run_20Hz(k),'marker','o','markersize',4,...
        'markeredgecolor','k','markerfacecolor','k',...
        'linestyle','none');
end

for k = 1:size(peak_Run_5Hz,1)
    plot(1.7,peak_Run_5Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_Run_20Hz,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_Run_20Hz,1));

[S_chow,M_chow] = std(peak_Run_5Hz,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_Run_5Hz,1));

errorbar(1.1, M_PB, SEM_PB, "Color",'k','LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[0.5 0.5 0.5],'LineWidth',2);
plot(1.1, M_PB,'marker','o','color','k','linewidth',3,'markeredgecolor','k','markerfacecolor','k','markersize',5);
plot(1.9, M_chow,'marker','o','color',[0.5 0.5 0.5],'linewidth',3,'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5);

line([1.3 1.7], [5.3, 5.3], 'Color', 'k', 'LineWidth', 2);
text(1.5, 5.45, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 3, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -2.3, '20Hz', 'Color', 'k', 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -2.3, '5Hz', 'Color', [0.5 0.5 0.5], 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-2 6])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Locomotion (cm/s)','FontSize',12,'FontWeight','bold');

hold off

[h_Run_20Hz_5Hz, p_Run_20Hz_5Hz, ~, stats_Run_20Hz_5Hz] = ttest(peak_Run_20Hz,peak_Run_5Hz,'Tail','right')
[h_Run_20Hz, p_Run_20Hz, ~, stats_Run_20Hz] = ttest(peak_Run_20Hz,0,'Tail','right')
[h_Run_5Hz, p_Run_5Hz, ~, stats_Run_5Hz] = ttest(peak_Run_5Hz,0,'Tail','right')

%%

peak_Pupil_20Hz = mean(Pupil_20Hz(:,100*FrameRate+1:160*FrameRate),2) - mean(Pupil_20Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_Pupil_5Hz = mean(Pupil_5Hz(:,100*FrameRate+1:160*FrameRate),2) - mean(Pupil_5Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_Pupil = [peak_Pupil_20Hz peak_Pupil_5Hz];

figure(1);

axes('Position', [0.5703 0.5482 0.15 0.1561]);
hold on

for k = 1:size(peak_Pupil,1)
    plot([1.3 1.7],peak_Pupil(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_Pupil_20Hz,1)
    plot(1.3,peak_Pupil_20Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

for k = 1:size(peak_Pupil_5Hz,1)
    plot(1.7,peak_Pupil_5Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_Pupil_20Hz,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_Pupil_20Hz,1));

[S_chow,M_chow] = std(peak_Pupil_5Hz,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_Pupil_5Hz,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[229 114 190]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[139 92 158]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[229 114 190]./255,'linewidth',3,'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[139 92 158]./255,'linewidth',3,'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,'markersize',5);

line([1.3 1.7], [630, 630], 'Color', 'k', 'LineWidth', 2);
text(1.5, 645, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 450, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, 150, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -80, '20Hz', 'Color', [229 114 190]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -80, '5Hz', 'Color', [139 92 158]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-50 700])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Pupil (pixels)','FontSize',12,'FontWeight','bold');

hold off

[h_Pupil, p_Pupil, ~, stats_Pupil] = ttest(peak_Pupil_20Hz,peak_Pupil_5Hz,'Tail','right')
[h_Pupil_20Hz, p_Pupil_20Hz, ~, stats_Pupil_20Hz] = ttest(peak_Pupil_20Hz,0,'Tail','right')
[h_Pupil_5Hz, p_Pupil_5Hz, ~, stats_Pupil_5Hz] = ttest(peak_Pupil_5Hz,0,'Tail','right')

%%

peak_HR_20Hz = mean(HR_20Hz(:,100*FrameRate+1:160*FrameRate),2) - mean(HR_20Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_HR_5Hz = mean(HR_5Hz(:,100*FrameRate+1:160*FrameRate),2) - mean(HR_5Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_HR = [peak_HR_20Hz peak_HR_5Hz];

figure(1);

axes('Position', [0.5703 0.3291 0.15 0.1561]);
hold on

for k = 1:size(peak_HR,1)
    plot([1.3 1.7],peak_HR(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_HR_20Hz,1)
    plot(1.3,peak_HR_20Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HR_5Hz,1)
    plot(1.7,peak_HR_5Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HR_20Hz,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HR_20Hz,1));

[S_chow,M_chow] = std(peak_HR_5Hz,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HR_5Hz,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[250, 157, 86]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[250, 157, 86]./255,'linewidth',3,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',5);

line([1.3 1.7], [91, 91], 'Color', 'k', 'LineWidth', 2);
text(1.5, 94, '**', 'FontSize', 20,'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 40, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, 20, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -47, '20Hz', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -47, '5Hz', 'Color', [250, 157, 86]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-40 100])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Heart Rate (bpm)','FontSize',12,'FontWeight','bold');

hold off

[h_HR, p_HR, ~, stats_HR] = ttest(peak_HR_20Hz,peak_HR_5Hz,'Tail','right')
[h_HR_20Hz, p_HR_20Hz, ~, stats_HR_20Hz] = ttest(peak_HR_20Hz,0,'Tail','right')
[h_HR_5Hz, p_HR_5Hz, ~, stats_HR_5Hz] = ttest(peak_HR_5Hz,0,'Tail','right')

%%

% peak_HR_20Hz = mean(HR_20Hz(:,160*FrameRate+1:460*FrameRate),2) - mean(HR_20Hz(:,40*FrameRate+1:100*FrameRate),2);
%
% peak_HR_5Hz = mean(HR_5Hz(:,160*FrameRate+1:460*FrameRate),2) - mean(HR_5Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_HR_20Hz = mean(HR_20Hz(:,160*FrameRate+1:end),2) - mean(HR_20Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_HR_5Hz = mean(HR_5Hz(:,160*FrameRate+1:end),2) - mean(HR_5Hz(:,40*FrameRate+1:100*FrameRate),2);

peak_HR = [peak_HR_20Hz peak_HR_5Hz];

figure(1);

axes('Position', [0.83 0.3291 0.15 0.1561]);
hold on

for k = 1:size(peak_HR,1)
    plot([1.3 1.7],peak_HR(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_HR_20Hz,1)
    plot(1.3,peak_HR_20Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HR_5Hz,1)
    plot(1.7,peak_HR_5Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HR_20Hz,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HR_20Hz,1));

[S_chow,M_chow] = std(peak_HR_5Hz,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HR_5Hz,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[250, 157, 86]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[250, 157, 86]./255,'linewidth',3,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',5);

line([1.3 1.7], [36, 36], 'Color', 'k', 'LineWidth', 2);
text(1.5, 45, 'p=0.0506', 'FontSize', 10,'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, -45, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, -30, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -128, '20Hz', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -128, '5Hz', 'Color', [250, 157, 86]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-120 50])
xlim([0.9 2.1])

title({'Laser OFF'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Heart Rate (bpm)','FontSize',12,'FontWeight','bold');

hold off

[h_HR, p_HR, ~, stats_HR] = ttest(peak_HR_20Hz,peak_HR_5Hz,'Tail','left')
[h_HR_20Hz, p_HR_20Hz, ~, stats_HR_20Hz] = ttest(peak_HR_20Hz,0,'Tail','left')
[h_HR_5Hz, p_HR_5Hz, ~, stats_HR_5Hz] = ttest(peak_HR_5Hz,0,'Tail','left')

%%

peak_HRV_20Hz = RMSSD_on_20Hz.*1000;

peak_HRV_5Hz = RMSSD_on_5Hz.*1000;

peak_HRV = [peak_HRV_20Hz peak_HRV_5Hz];

figure(1);

axes('Position', [0.5703 0.1100 0.15 0.1561]);
hold on

for k = 1:size(peak_HRV,1)
    plot([1.3 1.7],peak_HRV(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_HRV_20Hz,1)
    plot(1.3,peak_HRV_20Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HRV_5Hz,1)
    plot(1.7,peak_HRV_5Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HRV_20Hz,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HRV_20Hz,1));

[S_chow,M_chow] = std(peak_HRV_5Hz,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HRV_5Hz,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[250, 157, 86]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[250, 157, 86]./255,'linewidth',3,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',5);

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -16, '20Hz', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -16, '5Hz', 'Color', [250, 157, 86]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-15 10])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ RMSSD (ms)','FontSize',12,'FontWeight','bold');

hold off

[h_HRV, p_HRV, ~, stats_HRV] = ttest(peak_HRV_20Hz,peak_HRV_5Hz,'Tail','left')
[h_HRV_20Hz, p_HRV_20Hz, ~, stats_HRV_20Hz] = ttest(peak_HRV_20Hz,0,'Tail','left')
[h_HRV_5Hz, p_HRV_5Hz, ~, stats_HRV_5Hz] = ttest(peak_HRV_5Hz,0,'Tail','left')

%%

peak_HRV_20Hz = RMSSD_off_20Hz.*1000;

peak_HRV_5Hz = RMSSD_off_5Hz.*1000;

peak_HRV = [peak_HRV_20Hz peak_HRV_5Hz];

figure(1);

axes('Position', [0.83 0.1100 0.15 0.1561]);
hold on

for k = 1:size(peak_HRV,1)
    plot([1.3 1.7],peak_HRV(k,:),'marker','none','markersize',3,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',2);
end

for k = 1:size(peak_HRV_20Hz,1)
    plot(1.3,peak_HRV_20Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HRV_5Hz,1)
    plot(1.7,peak_HRV_5Hz(k),'marker','o','markersize',4,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HRV_20Hz,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HRV_20Hz,1));

[S_chow,M_chow] = std(peak_HRV_5Hz,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HRV_5Hz,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[250, 157, 86]./255,'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[250, 157, 86]./255,'linewidth',3,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',5);

text(1.1, 7, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.9, 6, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -11, '20Hz', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -11, '5Hz', 'Color', [250, 157, 86]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-10 15])
xlim([0.9 2.1])

title({'Laser OFF'},'FontSize',15,'FontWeight','bold')
ylabel('Δ RMSSD (ms)','FontSize',12,'FontWeight','bold');

hold off

[h_HRV, p_HRV, ~, stats_HRV] = ttest(peak_HRV_20Hz,peak_HRV_5Hz,'Tail','right')
[h_HRV_20Hz, p_HRV_20Hz, ~, stats_HRV_20Hz] = ttest(peak_HRV_20Hz,0,'Tail','right')
[h_HRV_5Hz, p_HRV_5Hz, ~, stats_HRV_5Hz] = ttest(peak_HRV_5Hz,0,'Tail','right')

%%

saveas(gcf, 'Orexin_OptoStim_Run_Pupil_HR_HRV.svg')

saveas(gcf, 'Orexin_OptoStim_Run_Pupil_HR_HRV.png')

