%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

folder = {'m2070\Feb_19_2025';
    'm2071\Feb_19_2025';
    'm2072\Feb_19_2025';
    'm2151\Feb_19_2025';
    'm2152\Feb_19_2025';
    'm2154\Feb_19_2025'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

fpass_trials=[5 10.5;
    6 11;
    6 12;
    6 12;
    7 12;
    5 11];

outlier_pupil=[0 99.97;
    0 99.99;
    0 99.99;
    0 99.95;
    0 99.99;
    0 100];

PB_onsets = [1541 1252 1315 1440 1277 1263];

pre_bout_duration = 20;   % 5 seconds.

trace_duration = 120;   % 30 seconds.

smooth_window = 1;

PB_Run = [];
startle_Run = [];

PB_HR = [];
startle_HR = [];

PB_Pupil = [];
startle_Pupil = [];

PB_Cell = [];
startle_Cell = [];

cells_num = 0;

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    Cell_Status_data = readtable([Data_Folder folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);
    cell_times = cell_traces_data(:,1);
    cell_traces = cell_traces_data(:,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times = cell_times(~nanRows, :);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    %     cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
    vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

    tone_onset = step_timepoint(3)-data_onset/Sample_Rate;

    timepoint = times(data_onset:data_offset,1)';
    time = timepoint(1,:)-timepoint(1,1);

    frame_time = 0:(1/FrameRate):time(end);

    %% Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num = cells_num + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        %         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
        cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
    end

    cell_traces_interpolated_smooth = movmean(cell_traces_interpolated,[smooth_window*FrameRate 0],2);

    %% Analyse running signals

    running = datas(data_onset:data_offset,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed_smooth = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);

    %% Analyse pupil size

    Pupil_up_x = pupil_data(vid_onset:min(end,vid_offset),1);
    Pupil_up_y = pupil_data(vid_onset:min(end,vid_offset),2);
    Pupil_left_x = pupil_data(vid_onset:min(end,vid_offset),4);
    Pupil_left_y = pupil_data(vid_onset:min(end,vid_offset),5);
    Pupil_down_x = pupil_data(vid_onset:min(end,vid_offset),7);
    Pupil_down_y = pupil_data(vid_onset:min(end,vid_offset),8);
    Pupil_right_x = pupil_data(vid_onset:min(end,vid_offset),10);
    Pupil_right_y = pupil_data(vid_onset:min(end,vid_offset),11);

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

    ECG_raw = datas(data_onset:data_offset,2)';

    % Remove baseline wandering

    fpass=fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %% Running bouts

    PB_t_start = round(PB_onsets(I)*FrameRate) - pre_bout_duration*FrameRate +1;
    PB_t_end = PB_t_start + trace_duration*FrameRate -1;

    PB_Run = [PB_Run; speed_smooth_resampled(PB_t_start:PB_t_end)];

    raw_pupil = pupil_smooth(PB_t_start:PB_t_end)';
    %     pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
    PB_Pupil = [PB_Pupil; raw_pupil];

    PB_HR = [PB_HR; heartRate_bpm_interp_smooth(PB_t_start:PB_t_end)];

    cell_traces_bouts = cell_traces_interpolated_smooth(:,PB_t_start:PB_t_end);
    cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),2))./std(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),0,2);
    %         cell_traces_bouts_zscored = zscore(cell_traces_bouts,0,2);
    PB_Cell = [PB_Cell; mean(cell_traces_bouts_zscored)];
    %     PB_Cell = [PB_Cell; cell_traces_bouts_zscored];


    startle_t_start = round(tone_onset*FrameRate) - pre_bout_duration*FrameRate +1;
    startle_t_end = startle_t_start + trace_duration*FrameRate -1;

    startle_Run = [startle_Run; speed_smooth_resampled(startle_t_start:startle_t_end)];

    raw_pupil = pupil_smooth(startle_t_start:startle_t_end)';
    %     pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
    startle_Pupil = [startle_Pupil; raw_pupil];

    startle_HR = [startle_HR; heartRate_bpm_interp_smooth(startle_t_start:startle_t_end)];

    cell_traces_bouts = cell_traces_interpolated_smooth(:,startle_t_start:startle_t_end);
    cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),2))./std(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),0,2);
    %         cell_traces_bouts_zscored = zscore(cell_traces_bouts,0,2);
    startle_Cell = [startle_Cell; mean(cell_traces_bouts_zscored)];
    %     startle_Cell = [startle_Cell; cell_traces_bouts_zscored];

end

%%

xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

trials_num = size(PB_Run,1);

traces_num = size(PB_Cell,1);

fig = figure(2);
set(fig, 'Position', [2561 49 2560 1315]);

ax1 = subplot(2,4,1);
set(ax1, 'Position', [0.1300 0.5838 0.2 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(startle_Cell)+std(startle_Cell)/sqrt(traces_num) fliplr(mean(startle_Cell)-std(startle_Cell)/sqrt(traces_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',1);
plot1 = plot(xlims,mean(startle_Cell),'Color','#547DB1','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(PB_Cell)+std(PB_Cell)/sqrt(traces_num) fliplr(mean(PB_Cell)-std(PB_Cell)/sqrt(traces_num))],'EdgeColor','none','FaceColor',[211 229 209]./255,'FaceAlpha',0.7);
plot2 = plot(xlims,mean(PB_Cell),'Color',[64 145 48]./255,'LineWidth',2);
line([0 0],[-5 25],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
% ylim([-1.5 3])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Orexin Neural Activities','FontSize',20,'FontWeight','bold')
ylabel('z-score (s.d.)','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'PB', 'Startle'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

ax3 = subplot(2,4,3);
set(ax3, 'Position', [0.5422 0.5838 0.2 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(startle_Run)+std(startle_Run)/sqrt(trials_num) fliplr(mean(startle_Run)-std(startle_Run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.7);
plot1 = plot(xlims,mean(startle_Run),'Color',[0.5 0.5 0.5],'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(PB_Run)+std(PB_Run)/sqrt(trials_num) fliplr(mean(PB_Run)-std(PB_Run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.7);
plot2 = plot(xlims,mean(PB_Run),'Color','k','LineWidth',2);
line([0,0],[0,0.07],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([0 0.07])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Running Speed','FontSize',20,'FontWeight','bold')
ylabel('Running Speed','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'PB', 'Startle'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

ax5 = subplot(2,4,5);
set(ax5, 'Position', [0.1300 0.1100 0.2 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(startle_Pupil)+std(startle_Pupil)/sqrt(trials_num) fliplr(mean(startle_Pupil)-std(startle_Pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',1);
plot1 = plot(xlims,mean(startle_Pupil),'Color',[139 92 158]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(PB_Pupil)+std(PB_Pupil)/sqrt(trials_num) fliplr(mean(PB_Pupil)-std(PB_Pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',1);
plot2 = plot(xlims,mean(PB_Pupil),'Color',[229 114 190]./255,'LineWidth',2);
line([0 0],[100 700],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([100 700])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Pupil Size Dynamics','FontSize',20,'FontWeight','bold')
ylabel('pixels','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'PB', 'Startle'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

ax7 = subplot(2,4,7);
set(ax7, 'Position', [0.5422 0.1100 0.2 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(startle_HR)+std(startle_HR)/sqrt(trials_num) fliplr(mean(startle_HR)-std(startle_HR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',1);
plot1 = plot(xlims,mean(startle_HR),'Color',[255 128 128]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(PB_HR)+std(PB_HR)/sqrt(trials_num) fliplr(mean(PB_HR)-std(PB_HR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255, 214, 170]./255,'FaceAlpha',1);
plot2 = plot(xlims,mean(PB_HR),'Color',[250, 157, 86]./255,'LineWidth',2);
line([0,0],[400,700],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([400 700])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Heart Rate','FontSize',20,'FontWeight','bold')
ylabel('HR (bpm)','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'PB', 'Startle'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

%%

% peak_PB_Cell = max(PB_Cell,[],2);
%
% peak_chow_Cell = max(chow_Cell,[],2);

% peak_PB_Cell = trapz(PB_Cell,2);
%
% peak_startle_Cell = trapz(startle_Cell,2);

% peak_PB_Cell = trapz(PB_Cell(:,25*FrameRate+1:50*FrameRate),2);
%
% peak_startle_Cell = trapz(startle_Cell(:,25*FrameRate+1:50*FrameRate),2);

peak_PB_Cell = mean(PB_Cell(:,20*FrameRate+1:end),2);

peak_startle_Cell = mean(startle_Cell(:,20*FrameRate+1:end),2);

peak_Cell = [peak_PB_Cell peak_startle_Cell];

figure(2);

ax2 = subplot(2,4,2);
set(ax2, 'Position', [0.38 0.5838 0.1 0.3412]);
hold on

for k = 1:size(peak_Cell,1)
    plot([1.3 1.7],peak_Cell(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(peak_PB_Cell,1)
    plot(1.3,peak_PB_Cell(k),'marker','o','markersize',5,...
        'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,...
        'linestyle','none');
end

for k = 1:size(peak_startle_Cell,1)
    plot(1.7,peak_startle_Cell(k),'marker','o','markersize',5,...
        'markeredgecolor','#547DB1','markerfacecolor','#547DB1',...
        'linestyle','none');
end

[S_PB,M_PB] = std(peak_PB_Cell,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_PB_Cell,1));

[S_startle,M_startle] = std(peak_startle_Cell,'omitnan');
SEM_startle = S_startle/sqrt(size(peak_startle_Cell,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[64 145 48]./255,'LineWidth',3);
errorbar(1.9, M_startle, SEM_startle, "Color",'#547DB1','LineWidth',3);
plot(1.1, M_PB,'marker','o','color',[64 145 48]./255,'linewidth',4,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,'markersize',6);
plot(1.9, M_startle,'marker','o','color','#547DB1','linewidth',4,'markeredgecolor','#547DB1','markerfacecolor','#547DB1','markersize',6);

line([1.3 1.7], [23.5, 23.5], 'Color', 'k', 'LineWidth', 2);
text(1.5, 24, '*', 'FontSize', 25, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};
% ax.YTickLabel = {};

xticks([1.2 1.8])

text(1.2, -1, 'PB', 'Color', [64 145 48]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -1, 'Startle', 'Color', '#547DB1', 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([0 25])
xlim([0.9 2.1])

ylabel('Mean Orexin Response','FontSize',18,'FontWeight','bold');

hold off

% [h_Cell_PB, p_Cell_PB, ~, stats_Cell_PB] = ttest(peak_PB_Cell,0,'Tail','right')
% [h_Cell_startle, p_Cell_startle, ~, stats_Cell_startle] = ttest(peak_startle_Cell,0,'Tail','right')
[h_Cell, p_Cell, ~, stats_Cell] = ttest(peak_PB_Cell,peak_startle_Cell,'Tail','right')

%%

% peak_PB_Run = max(PB_Run,[],2);
%
% peak_startle_Run = max(startle_Run,[],2);

% peak_PB_Run = trapz(PB_Run,2);
%
% peak_startle_Run = trapz(startle_Run,2);

peak_PB_Run = mean(PB_Run(:,20*FrameRate+1:end),2);

peak_startle_Run = mean(startle_Run(:,20*FrameRate+1:end),2);

peak_Run = [peak_PB_Run peak_startle_Run];

figure(2);

ax4 = subplot(2,4,4);
set(ax4, 'Position', [0.7923 0.5838 0.1 0.3412]);
hold on

for k = 1:size(peak_Run,1)
    plot([1.3 1.7],peak_Run(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(peak_PB_Run,1)
    plot(1.3,peak_PB_Run(k),'marker','o','markersize',5,...
        'markeredgecolor','k','markerfacecolor','k',...
        'linestyle','none');
end

for k = 1:size(peak_startle_Run,1)
    plot(1.7,peak_startle_Run(k),'marker','o','markersize',5,...
        'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],...
        'linestyle','none');
end

[S_PB,M_PB] = std(peak_PB_Run,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_PB_Run,1));

[S_startle,M_startle] = std(peak_startle_Run,'omitnan');
SEM_startle = S_startle/sqrt(size(peak_startle_Run,1));

errorbar(1.1, M_PB, SEM_PB, "Color",'k','LineWidth',3);
errorbar(1.9, M_startle, SEM_startle, "Color",[0.5 0.5 0.5],'LineWidth',3);
plot(1.1, M_PB,'marker','o','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',6);
plot(1.9, M_startle,'marker','o','color',[0.5 0.5 0.5],'linewidth',4,'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',6);

line([1.3 1.7], [0.028, 0.028], 'Color', 'k', 'LineWidth', 2);
text(1.5, 0.0286, '**', 'FontSize', 25,'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};
% ax.YTickLabel = {};

xticks([1.2 1.8])

text(1.2, -0.001, 'PB', 'Color', 'k', 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -0.001, 'Startle', 'Color', [0.5 0.5 0.5], 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([0 0.030])
xlim([0.9 2.1])

ylabel('Mean Locomotion Response','FontSize',18,'FontWeight','bold');

hold off

[h_Run, p_Run, ~, stats_Run] = ttest(peak_PB_Run,peak_startle_Run,'Tail','right')

%%

% peak_PB_Pupil = max(PB_Pupil,[],2);
%
% peak_chow_Pupil = max(chow_Pupil,[],2);

% peak_PB_Pupil = trapz(PB_Pupil,2);
%
% peak_startle_Pupil = trapz(startle_Pupil,2);

peak_PB_Pupil = mean(PB_Pupil(:,20*FrameRate+1:end),2);

peak_startle_Pupil = mean(startle_Pupil(:,20*FrameRate+1:end),2);

peak_Pupil = [peak_PB_Pupil peak_startle_Pupil];

figure(2);

ax6 = subplot(2,4,6);
set(ax6, 'Position', [0.38 0.1100 0.1 0.3412]);
hold on

for k = 1:size(peak_Pupil,1)
    plot([1.3 1.7],peak_Pupil(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(peak_PB_Pupil,1)
    plot(1.3,peak_PB_Pupil(k),'marker','o','markersize',5,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

for k = 1:size(peak_startle_Pupil,1)
    plot(1.7,peak_startle_Pupil(k),'marker','o','markersize',5,...
        'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,...
        'linestyle','none');
end

[S_PB,M_PB] = std(peak_PB_Pupil,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_PB_Pupil,1));

[S_startle,M_startle] = std(peak_startle_Pupil,'omitnan');
SEM_startle = S_startle/sqrt(size(peak_startle_Pupil,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[229 114 190]./255,'LineWidth',3);
errorbar(1.9, M_startle, SEM_startle, "Color",[139 92 158]./255,'LineWidth',3);
plot(1.1, M_PB,'marker','o','color',[229 114 190]./255,'linewidth',4,'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,'markersize',6);
plot(1.9, M_startle,'marker','o','color',[139 92 158]./255,'linewidth',4,'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,'markersize',6);

line([1.3 1.7], [760, 760], 'Color', 'k', 'LineWidth', 2);
text(1.5, 770, '*', 'FontSize', 25,'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};
% ax.YTickLabel = {};

xticks([1.2 1.8])

text(1.2, 80, 'PB', 'Color', [229 114 190]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, 80, 'Startle', 'Color', [139 92 158]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([100 800])
xlim([0.9 2.1])

ylabel('Mean Pupil Response','FontSize',18,'FontWeight','bold');

hold off

[h_Pupil, p_Pupil, ~, stats_Pupil] = ttest(peak_PB_Pupil,peak_startle_Pupil,'Tail','right')

%%

% peak_PB_HR = max(PB_HR,[],2);
%
% peak_startle_HR = max(startle_HR,[],2);

peak_PB_HR = mean(PB_HR(:,20*FrameRate+1:end),2);

peak_startle_HR = mean(startle_HR(:,20*FrameRate+1:end),2);

peak_HR = [peak_PB_HR peak_startle_HR];

figure(2);

ax8 = subplot(2,4,8);
set(ax8, 'Position', [0.7923 0.1100 0.1 0.3412]);
hold on

for k = 1:size(peak_HR,1)
    plot([1.3 1.7],peak_HR(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(peak_PB_HR,1)
    plot(1.3,peak_PB_HR(k),'marker','o','markersize',5,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

for k = 1:size(peak_startle_HR,1)
    plot(1.7,peak_startle_HR(k),'marker','o','markersize',5,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

[S_PB,M_PB] = std(peak_PB_HR,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_PB_HR,1));

[S_startle,M_startle] = std(peak_startle_HR,'omitnan');
SEM_startle = S_startle/sqrt(size(peak_startle_HR,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[250, 157, 86]./255,'LineWidth',3);
errorbar(1.9, M_startle, SEM_startle, "Color",[255 128 128]./255,'LineWidth',3);
plot(1.1, M_PB,'marker','o','color',[250, 157, 86]./255,'linewidth',4,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',6);
plot(1.9, M_startle,'marker','o','color',[255 128 128]./255,'linewidth',4,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',6);

line([1.3 1.7], [660, 660], 'Color', 'k', 'LineWidth', 2);
text(1.5, 665, 'NS', 'FontSize', 15,'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};
% ax.YTickLabel = {};

xticks([1.2 1.8])

text(1.2, 495, 'PB', 'Color', [250, 157, 86]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, 495, 'Startle', 'Color', [255 128 128]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([500 680])
xlim([0.9 2.1])

ylabel('Mean Cardiac Response','FontSize',18,'FontWeight','bold');

hold off

[h_HR, p_HR, ~, stats_HR] = ttest(peak_PB_HR,peak_startle_HR,'Tail','right')

%%

% saveas(gcf, 'PB_Startle_Orexin_Pupil_Run_HR_integral.svg')
%
% saveas(gcf, 'PB_Startle_Orexin_Pupil_Run_HR_integral.png')
