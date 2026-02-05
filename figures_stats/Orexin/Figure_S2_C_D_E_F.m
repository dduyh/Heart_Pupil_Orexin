%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

novel_folder = {'m2070\Mar_20_2025';
    'm2072\Mar_20_2025';
    'm2151\Mar_20_2025';
    'm2152\Mar_20_2025'};

familiar_folder = {'m2070\Mar_24_2025';
    'm2072\Mar_24_2025';
    'm2151\Mar_24_2025';
    'm2152\Mar_24_2025'};

novel_fpass_trials=[7 12;
    6 11;
    6 11;
    7 11];

familiar_fpass_trials=[6 10.5;
    7 12;
    6.5 11;
    6.5 11];

novel_outlier_pupil=[0 100;
    0 99.99;
    0 99.99;
    0 99.97];

familiar_outlier_pupil=[0 99.8;
    0 99.9;
    0 99.9;
    0 99.92];

novel_Raisins_onsets = [1379 1196 1238 1140];

familiar_Raisins_onsets = [1262 1378 1311 1275];

Novel_Run = [];
Familiar_Run = [];

Novel_HR = [];
Familiar_HR = [];

Novel_Pupil = [];
Familiar_Pupil = [];

Novel_Cell = [];
Familiar_Cell = [];

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

pre_bout_duration = 20;   % 5 seconds.

trace_duration = 120;   % 30 seconds.

smooth_window = 1;

cells_num = 0;

for I=1:size(novel_folder,1)

    Data_Folder = [Directory novel_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    Cell_Status_data = readtable([Data_Folder novel_folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder novel_folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);

    cell_times = cell_traces_data(:,1);
    dt = diff(cell_times);
    [~, jump_frame] = max(dt);

    cell_times_half = cell_times(1:jump_frame,1);
    cell_traces = cell_traces_data(1:jump_frame,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times_half = cell_times_half(~nanRows, :);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    %             cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
    vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

    timepoint = times(data_onset:data_offset,1)';
    time = timepoint(1,:)-timepoint(1,1);

    frame_time = 0:(1/FrameRate):time(end);

    %% Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num = cells_num + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        %         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
        cell_traces_interpolated(i,:) = interp1(cell_times_half,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
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

    pupil = filloutliers(areas,"nearest","percentiles",novel_outlier_pupil(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(data_onset:data_offset,2)';

    % Remove baseline wandering

    fpass=novel_fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[5*FrameRate 0]);

    %% Running bouts

    Raisins_t_start = round(novel_Raisins_onsets(I)*FrameRate) - pre_bout_duration*FrameRate +1;
    Raisins_t_end = Raisins_t_start + trace_duration*FrameRate -1;

    Novel_Run = [Novel_Run; speed_smooth_resampled(Raisins_t_start:Raisins_t_end)];

    raw_pupil = pupil_smooth(Raisins_t_start:Raisins_t_end)';
    pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
    Novel_Pupil = [Novel_Pupil; pupil_zscored];

    Novel_HR = [Novel_HR; heartRate_bpm_interp_smooth(Raisins_t_start:Raisins_t_end)];

    cell_traces_bouts = cell_traces_interpolated_smooth(:,Raisins_t_start:Raisins_t_end);
    cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),2))./std(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),0,2);
    %         cell_traces_bouts_zscored = zscore(cell_traces_bouts,0,2);
    Novel_Cell = [Novel_Cell; mean(cell_traces_bouts_zscored)];
    %     Novel_Cell = [Novel_Cell; cell_traces_bouts_zscored];

end

for I=1:size(familiar_folder,1)

    Data_Folder = [Directory familiar_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    Cell_Status_data = readtable([Data_Folder familiar_folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder familiar_folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);

    cell_times = cell_traces_data(:,1);
    dt = diff(cell_times);
    [~, jump_frame] = max(dt);

    cell_times_half = cell_times(jump_frame+1:end,1);
    cell_times_half = cell_times_half(:,1)-cell_times_half(1,1);
    cell_traces = cell_traces_data(jump_frame+1:end,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times_half = cell_times_half(~nanRows, :);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    %             cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
    vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

    timepoint = times(data_onset:data_offset,1)';
    time = timepoint(1,:)-timepoint(1,1);

    frame_time = 0:(1/FrameRate):time(end);

    %% Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num = cells_num + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        %         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
        cell_traces_interpolated(i,:) = interp1(cell_times_half,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
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

    pupil = filloutliers(areas,"nearest","percentiles",familiar_outlier_pupil(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(data_onset:data_offset,2)';

    % Remove baseline wandering

    fpass=familiar_fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[5*FrameRate 0]);

    %% Running bouts

    Raisins_t_start = round(familiar_Raisins_onsets(I)*FrameRate) - pre_bout_duration*FrameRate +1;
    Raisins_t_end = Raisins_t_start + trace_duration*FrameRate -1;

    Familiar_Run = [Familiar_Run; speed_smooth_resampled(Raisins_t_start:Raisins_t_end)];

    raw_pupil = pupil_smooth(Raisins_t_start:Raisins_t_end)';
    pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
    Familiar_Pupil = [Familiar_Pupil; pupil_zscored];

    Familiar_HR = [Familiar_HR; heartRate_bpm_interp_smooth(Raisins_t_start:Raisins_t_end)];

    cell_traces_bouts = cell_traces_interpolated_smooth(:,Raisins_t_start:Raisins_t_end);
    cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),2))./std(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),0,2);
    %         cell_traces_bouts_zscored = zscore(cell_traces_bouts,0,2);
    Familiar_Cell = [Familiar_Cell; mean(cell_traces_bouts_zscored)];
    %     Familiar_Cell = [Familiar_Cell; cell_traces_bouts_zscored];

end

%% rank sorting

traces_num = size(Familiar_Cell,1);

mean_response_Novel = mean(Novel_Cell(:, 20*FrameRate:40*FrameRate-1), 2);
mean_response_Familiar = mean(Familiar_Cell(:, 20*FrameRate:40*FrameRate-1), 2);

[~, novel_sorted_indices] = sort(mean_response_Novel);
[~, rank_Novel] = sort(novel_sorted_indices);

[~, familiar_sorted_indices] = sort(mean_response_Familiar);
[~, rank_Familiar] = sort(familiar_sorted_indices);

rank_Cell = [rank_Novel rank_Familiar];

figure(3);

hold on

for k = 1:size(rank_Cell,1)
    plot([1.3 1.7],rank_Cell(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(rank_Novel,1)
    plot(1.3,rank_Novel(k),'marker','o','markersize',5,...
        'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,...
        'linestyle','none');
end

for k = 1:size(rank_Familiar,1)
    plot(1.7,rank_Familiar(k),'marker','o','markersize',5,...
        'markeredgecolor','#547DB1','markerfacecolor','#547DB1',...
        'linestyle','none');
end

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.3 1.7])

text(1.3, -9, 'Novel', 'Color', [64 145 48]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.7, -9, 'Familiar', 'Color', '#547DB1', 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-5 70])
xlim([1.1 1.9])

ylabel('Ranks of Orexin Response','FontSize',18,'FontWeight','bold');

hold off

% [h_rank, p_rank, ~, stats_rank] = ttest(rank_Raisins,rank_Startle)

% saveas(gcf, 'Raisins_Startle_Orexin_rank.svg')
%
% saveas(gcf, 'Raisins_Startle_Orexin_rank.png')

%%

novel_sorted_traces = Novel_Cell(novel_sorted_indices, :);

fig1 = figure(1);
set(fig1, 'Position', [1 49 2560 1315]);

ax1 = subplot(1,2,1);
set(ax1,'YDir','reverse');

hold on

% imagesc(sorted_traces);
imagesc(novel_sorted_traces,[-10 40]);
colormap('turbo');

line([pre_bout_duration*FrameRate,pre_bout_duration*FrameRate],[1,traces_num],'Color','white','linestyle','--','LineWidth',2);
ylim([1,traces_num])

xlim([1,trace_duration*FrameRate])
xticks([1:20*FrameRate:trace_duration*FrameRate trace_duration*FrameRate])
xticklabels({'-20','0','20','40','60','80','100'})

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 18;
ax.FontWeight = 'bold';

title('Novel Raisins','FontSize',20,'FontWeight','bold')
ylabel('sorted cell traces','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')

hold off

%%

familiar_sorted_traces = Familiar_Cell(familiar_sorted_indices, :);

figure(1);

ax2 = subplot(1,2,2);
set(ax2,'YDir','reverse');

hold on

% imagesc(sorted_traces);
imagesc(familiar_sorted_traces,[-10 40]);
colormap('turbo');

hbar = colorbar('north');
set(hbar, 'Position', [0.8238 0.0494 0.0814 0.0238])
hbar.AxisLocation = 'out';
hbar.FontSize = 15;
hbar.Label.String = 'z-score';
hbar.FontWeight = 'bold';

line([pre_bout_duration*FrameRate,pre_bout_duration*FrameRate],[1,traces_num],'Color','white','linestyle','--','LineWidth',2);
ylim([1,traces_num])

xlim([1,trace_duration*FrameRate])
xticks([1:20*FrameRate:trace_duration*FrameRate trace_duration*FrameRate])
xticklabels({'-20','0','20','40','60','80','100'})

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 18;
ax.FontWeight = 'bold';

title('Familiar Raisins','FontSize',20,'FontWeight','bold')
ylabel('sorted cell traces','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')

hold off

% saveas(gcf, 'Raisins_Startle_Orexin_heatmap.svg')
%
% saveas(gcf, 'Raisins_Startle_Orexin_heatmap.png')

%%

xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

trials_num = size(Novel_Run,1);

traces_num = size(Novel_Cell,1);

fig2 = figure(2);
set(fig2, 'Position', [2561 49 2560 1315]);

ax1 = subplot(2,4,1);
set(ax1, 'Position', [0.1300 0.5838 0.2 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Familiar_Cell)+std(Familiar_Cell)/sqrt(traces_num) fliplr(mean(Familiar_Cell)-std(Familiar_Cell)/sqrt(traces_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',1);
plot1 = plot(xlims,mean(Familiar_Cell),'Color','#547DB1','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Novel_Cell)+std(Novel_Cell)/sqrt(traces_num) fliplr(mean(Novel_Cell)-std(Novel_Cell)/sqrt(traces_num))],'EdgeColor','none','FaceColor',[211 229 209]./255,'FaceAlpha',0.7);
plot2 = plot(xlims,mean(Novel_Cell),'Color',[64 145 48]./255,'LineWidth',2);
line([0 0],[-5 35],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([-5 35])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Orexin Neural Activities','FontSize',20,'FontWeight','bold')
ylabel('z-score (s.d.)','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

ax3 = subplot(2,4,3);
set(ax3, 'Position', [0.5422 0.5838 0.2 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Familiar_Run)+std(Familiar_Run)/sqrt(trials_num) fliplr(mean(Familiar_Run)-std(Familiar_Run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.7);
plot1 = plot(xlims,mean(Familiar_Run),'Color',[0.5 0.5 0.5],'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Novel_Run)+std(Novel_Run)/sqrt(trials_num) fliplr(mean(Novel_Run)-std(Novel_Run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.7);
plot2 = plot(xlims,mean(Novel_Run),'Color','k','LineWidth',2);
line([0,0],[-0.02,0.14],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([0 0.14])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Running Speed','FontSize',20,'FontWeight','bold')
ylabel('Running Speed (AU)','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

ax5 = subplot(2,4,5);
set(ax5, 'Position', [0.1300 0.1100 0.2 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Familiar_Pupil)+std(Familiar_Pupil)/sqrt(trials_num) fliplr(mean(Familiar_Pupil)-std(Familiar_Pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.5);
plot1 = plot(xlims,mean(Familiar_Pupil),'Color',[139 92 158]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Novel_Pupil)+std(Novel_Pupil)/sqrt(trials_num) fliplr(mean(Novel_Pupil)-std(Novel_Pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.5);
plot2 = plot(xlims,mean(Novel_Pupil),'Color',[229 114 190]./255,'LineWidth',2);
line([0 0],[-10 100],'Color','k','linestyle','--','LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([-10 60])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Pupil Size Dynamics','FontSize',20,'FontWeight','bold')
ylabel('z-score (s.d.)','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

ax7 = subplot(2,4,7);
set(ax7, 'Position', [0.5422 0.1100 0.2 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Familiar_HR)+std(Familiar_HR)/sqrt(trials_num) fliplr(mean(Familiar_HR)-std(Familiar_HR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.5);
plot1 = plot(xlims,mean(Familiar_HR),'Color',[255 128 128]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Novel_HR)+std(Novel_HR)/sqrt(trials_num) fliplr(mean(Novel_HR)-std(Novel_HR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255, 214, 170]./255,'FaceAlpha',0.5);
plot2 = plot(xlims,mean(Novel_HR),'Color',[250, 157, 86]./255,'LineWidth',2);
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
legend([plot2, plot1], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

%%

peak_novel_Cell = mean(Novel_Cell(:,20*FrameRate+1:end),2);

peak_familiar_Cell = mean(Familiar_Cell(:,20*FrameRate+1:end),2);

peak_Cell = [peak_novel_Cell peak_familiar_Cell];

figure(2);

ax2 = subplot(2,4,2);
set(ax2, 'Position', [0.38 0.5838 0.1 0.3412]);
hold on

for k = 1:size(peak_Cell,1)
    plot([1.3 1.7],peak_Cell(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(peak_novel_Cell,1)
    plot(1.3,peak_novel_Cell(k),'marker','o','markersize',5,...
        'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,...
        'linestyle','none');
end

for k = 1:size(peak_familiar_Cell,1)
    plot(1.7,peak_familiar_Cell(k),'marker','o','markersize',5,...
        'markeredgecolor','#547DB1','markerfacecolor','#547DB1',...
        'linestyle','none');
end

[S_novel,M_novel] = std(peak_novel_Cell,'omitnan');
SEM_novel = S_novel/sqrt(size(peak_novel_Cell,1));

[S_familiar,M_familiar] = std(peak_familiar_Cell,'omitnan');
SEM_familiar = S_familiar/sqrt(size(peak_familiar_Cell,1));

errorbar(1.1, M_novel, SEM_novel, "Color",[64 145 48]./255,'LineWidth',3);
errorbar(1.9, M_familiar, SEM_familiar, "Color",'#547DB1','LineWidth',3);
plot(1.1, M_novel,'marker','o','color',[64 145 48]./255,'linewidth',4,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,'markersize',6);
plot(1.9, M_familiar,'marker','o','color','#547DB1','linewidth',4,'markeredgecolor','#547DB1','markerfacecolor','#547DB1','markersize',6);

line([1.3 1.7], [33, 33], 'Color', 'k', 'LineWidth', 2);
text(1.5, 34, 'NS', 'FontSize', 15, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};
% ax.YTickLabel = {};

xticks([1.2 1.8])

text(1.2, -6, 'Novel', 'Color', [64 145 48]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -6, 'Familiar', 'Color', '#547DB1', 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-5 35])
xlim([0.9 2.1])

ylabel('Mean Orexin Response','FontSize',18,'FontWeight','bold');

hold off

% [h_Cell_Raisins, p_Cell_Raisins, ~, stats_Cell_Raisins] = ttest(peak_Raisins_Cell,0,'Tail','right')
% [h_Cell_startle, p_Cell_startle, ~, stats_Cell_startle] = ttest(peak_startle_Cell,0,'Tail','right')
[h_Cell, p_Cell, ~, stats_Cell] = ttest(peak_novel_Cell,peak_familiar_Cell,'Tail','left')

%%

peak_novel_Run = mean(Novel_Run(:,20*FrameRate+1:end),2);

peak_familiar_Run = mean(Familiar_Run(:,20*FrameRate+1:end),2);

peak_Run = [peak_novel_Run peak_familiar_Run];

figure(2);

ax4 = subplot(2,4,4);
set(ax4, 'Position', [0.7923 0.5838 0.1 0.3412]);
hold on

for k = 1:size(peak_Run,1)
    plot([1.3 1.7],peak_Run(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(peak_novel_Run,1)
    plot(1.3,peak_novel_Run(k),'marker','o','markersize',5,...
        'markeredgecolor','k','markerfacecolor','k',...
        'linestyle','none');
end

for k = 1:size(peak_familiar_Run,1)
    plot(1.7,peak_familiar_Run(k),'marker','o','markersize',5,...
        'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],...
        'linestyle','none');
end

[S_novel,M_novel] = std(peak_novel_Run,'omitnan');
SEM_novel = S_novel/sqrt(size(peak_novel_Run,1));

[S_familiar,M_familiar] = std(peak_familiar_Run,'omitnan');
SEM_familiar = S_familiar/sqrt(size(peak_familiar_Run,1));

errorbar(1.1, M_novel, SEM_novel, "Color",'k','LineWidth',3);
errorbar(1.9, M_familiar, SEM_familiar, "Color",[0.5 0.5 0.5],'LineWidth',3);
plot(1.1, M_novel,'marker','o','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',6);
plot(1.9, M_familiar,'marker','o','color',[0.5 0.5 0.5],'linewidth',4,'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',6);

line([1.3 1.7], [0.038, 0.038], 'Color', 'k', 'LineWidth', 2);
text(1.5, 0.0388, 'NS', 'FontSize', 15,'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};
% ax.YTickLabel = {};

xticks([1.2 1.8])

text(1.2, 0.0143, 'Novel', 'Color', 'k', 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, 0.0143, 'Familiar', 'Color', [0.5 0.5 0.5], 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([0.015 0.040])
xlim([0.9 2.1])

ylabel('Mean Locomotion Response','FontSize',18,'FontWeight','bold');

hold off

[h_Run, p_Run, ~, stats_Run] = ttest(peak_novel_Run,peak_familiar_Run,'Tail','left')

%%

peak_novel_Pupil = mean(Novel_Pupil(:,20*FrameRate+1:end),2);

peak_familiar_Pupil = mean(Familiar_Pupil(:,20*FrameRate+1:end),2);

peak_Pupil = [peak_novel_Pupil peak_familiar_Pupil];

figure(2);

ax6 = subplot(2,4,6);
set(ax6, 'Position', [0.38 0.1100 0.1 0.3412]);
hold on

for k = 1:size(peak_Pupil,1)
    plot([1.3 1.7],peak_Pupil(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(peak_novel_Pupil,1)
    plot(1.3,peak_novel_Pupil(k),'marker','o','markersize',5,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

for k = 1:size(peak_familiar_Pupil,1)
    plot(1.7,peak_familiar_Pupil(k),'marker','o','markersize',5,...
        'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,...
        'linestyle','none');
end

[S_novel,M_novel] = std(peak_novel_Pupil,'omitnan');
SEM_novel = S_novel/sqrt(size(peak_novel_Pupil,1));

[S_familiar,M_familiar] = std(peak_familiar_Pupil,'omitnan');
SEM_familiar = S_familiar/sqrt(size(peak_familiar_Pupil,1));

errorbar(1.1, M_novel, SEM_novel, "Color",[229 114 190]./255,'LineWidth',3);
errorbar(1.9, M_familiar, SEM_familiar, "Color",[139 92 158]./255,'LineWidth',3);
plot(1.1, M_novel,'marker','o','color',[229 114 190]./255,'linewidth',4,'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,'markersize',6);
plot(1.9, M_familiar,'marker','o','color',[139 92 158]./255,'linewidth',4,'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,'markersize',6);

line([1.3 1.7], [37, 37], 'Color', 'k', 'LineWidth', 2);
text(1.5, 38.5, 'NS', 'FontSize', 15,'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};
% ax.YTickLabel = {};

xticks([1.2 1.8])

text(1.2, -6, 'Novel', 'Color', [229 114 190]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -6, 'Familiar', 'Color', [139 92 158]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-5 40])
xlim([0.9 2.1])

ylabel('Mean Pupil Response','FontSize',18,'FontWeight','bold');

hold off

[h_Pupil, p_Pupil, ~, stats_Pupil] = ttest(peak_novel_Pupil,peak_familiar_Pupil,'Tail','left')

%%

peak_novel_HR = mean(Novel_HR(:,20*FrameRate+1:end),2);

peak_familiar_HR = mean(Familiar_HR(:,20*FrameRate+1:end),2);

peak_HR = [peak_novel_HR peak_familiar_HR];

figure(2);

ax8 = subplot(2,4,8);
set(ax8, 'Position', [0.7923 0.1100 0.1 0.3412]);
hold on

for k = 1:size(peak_HR,1)
    plot([1.3 1.7],peak_HR(k,:),'marker','none','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 1],'linewidth',3);
end

for k = 1:size(peak_novel_HR,1)
    plot(1.3,peak_novel_HR(k),'marker','o','markersize',5,...
        'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,...
        'linestyle','none');
end

for k = 1:size(peak_familiar_HR,1)
    plot(1.7,peak_familiar_HR(k),'marker','o','markersize',5,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

[S_novel,M_novel] = std(peak_novel_HR,'omitnan');
SEM_novel = S_novel/sqrt(size(peak_novel_HR,1));

[S_familiar,M_familiar] = std(peak_familiar_HR,'omitnan');
SEM_familiar = S_familiar/sqrt(size(peak_familiar_HR,1));

errorbar(1.1, M_novel, SEM_novel, "Color",[250, 157, 86]./255,'LineWidth',3);
errorbar(1.9, M_familiar, SEM_familiar, "Color",[255 128 128]./255,'LineWidth',3);
plot(1.1, M_novel,'marker','o','color',[250, 157, 86]./255,'linewidth',4,'markeredgecolor',[250, 157, 86]./255,'markerfacecolor',[250, 157, 86]./255,'markersize',6);
plot(1.9, M_familiar,'marker','o','color',[255 128 128]./255,'linewidth',4,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',6);

line([1.3 1.7], [610, 610], 'Color', 'k', 'LineWidth', 2);
text(1.5, 615, 'NS', 'FontSize', 15,'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};
% ax.YTickLabel = {};

xticks([1.2 1.8])

text(1.2, 455, 'Novel', 'Color', [250, 157, 86]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, 455, 'Familiar', 'Color', [255 128 128]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([460 620])
xlim([0.9 2.1])

ylabel('Mean Cardiac Response','FontSize',18,'FontWeight','bold');

hold off

[h_HR, p_HR, ~, stats_HR] = ttest(peak_novel_HR,peak_familiar_HR,'Tail','left')

%%

% saveas(gcf, 'Raisins_Startle_Orexin_Pupil_Run_HR_mean.svg')
%
% saveas(gcf, 'Raisins_Startle_Orexin_Pupil_Run_HR_mean.png')
