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

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

novel_Raisins_onsets = [1379 1196 1238 1140];

familiar_Raisins_onsets = [1262 1378 1311 1275];

pre_bout_duration = 20;   % 5 seconds.

trace_duration = 120;   % 30 seconds.

smooth_window = 1;

Novel_Cell = [];
Familiar_Cell = [];

cells_num = 0;

for I=1:size(novel_folder,1)

    Data_Folder = [Directory novel_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

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

    %% Stimulation bouts

    Raisins_t_start = round(novel_Raisins_onsets(I)*FrameRate) - pre_bout_duration*FrameRate +1;
    Raisins_t_end = Raisins_t_start + trace_duration*FrameRate -1;

    cell_traces_bouts = cell_traces_interpolated_smooth(:,Raisins_t_start:Raisins_t_end);
    cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),2))./std(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),0,2);
    %         cell_traces_bouts_zscored = zscore(cell_traces_bouts,0,2);
    %     Novel_Cell = [Novel_Cell; mean(cell_traces_bouts_zscored)];
    Novel_Cell = [Novel_Cell; cell_traces_bouts_zscored];

end

for I=1:size(familiar_folder,1)

    Data_Folder = [Directory familiar_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

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

    %% Stimulation bouts

    Raisins_t_start = round(familiar_Raisins_onsets(I)*FrameRate) - pre_bout_duration*FrameRate +1;
    Raisins_t_end = Raisins_t_start + trace_duration*FrameRate -1;

    cell_traces_bouts = cell_traces_interpolated_smooth(:,Raisins_t_start:Raisins_t_end);
    cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),2))./std(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),0,2);
    %         cell_traces_bouts_zscored = zscore(cell_traces_bouts,0,2);
    %     Familiar_Cell = [Familiar_Cell; mean(cell_traces_bouts_zscored)];
    Familiar_Cell = [Familiar_Cell; cell_traces_bouts_zscored];

end


%% rank sorting

traces_num = size(Familiar_Cell,1);

mean_response_Novel = mean(Novel_Cell(:, 20*FrameRate:40*FrameRate-1), 2);
mean_response_Familiar = mean(Familiar_Cell(:, 20*FrameRate:40*FrameRate-1), 2);

[~, novel_sorted_indices] = sort(mean_response_Novel);
[~, familiar_sorted_indices] = sort(mean_response_Familiar);

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

%% PCA

concatenated_traces = [novel_sorted_traces familiar_sorted_traces];

% concatenated_traces = [Novel_Cell Familiar_Cell];

concatenated_traces_zscore = zscore(concatenated_traces, 0, 2);

[coeff, score, ~, ~, explained] = pca(concatenated_traces_zscore');

%% PCA 2D plots

trajectory_novel = score(1:trace_duration*FrameRate, 1:2);
trajectory_familiar = score(trace_duration*FrameRate+1:end, 1:2);

fig2 = figure(2);
set(fig2, 'Position', [1 49 680 630]);

plot1 = plot(trajectory_novel(:,1), trajectory_novel(:,2),'-o','color','#547DB1','linewidth',3,'markersize',4,'markeredgecolor','#547DB1','markerfacecolor','#547DB1');
hold on;
plot2 = plot(trajectory_familiar(:,1), trajectory_familiar(:,2),'-o','color',[64 145 48]./255,'linewidth',3,'markersize',4,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255);

xlim([-10 20])
ylim([-8 12])
xlabel('PC1','FontSize',18,'FontWeight','bold');
ylabel('PC2','FontSize',18,'FontWeight','bold');
axis square
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
grid off;
title('Neuronal Trajectories in PCA Space','FontSize',20,'FontWeight','bold');
legend([plot1, plot2], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
legend('boxoff')
hold off

%% PCA 2D plots

trajectory_novel = score(1:trace_duration*FrameRate, [1 3]);
trajectory_familiar = score(trace_duration*FrameRate+1:end, [1 3]);

fig3 = figure(3);
set(fig3, 'Position', [1 49 680 630]);

plot1 = plot(trajectory_novel(:,1), trajectory_novel(:,2),'-o','color','#547DB1','linewidth',3,'markersize',4,'markeredgecolor','#547DB1','markerfacecolor','#547DB1');
hold on;
plot2 = plot(trajectory_familiar(:,1), trajectory_familiar(:,2),'-o','color',[64 145 48]./255,'linewidth',3,'markersize',4,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255);

xlim([-10 20])
ylim([-10 10])
xlabel('PC1','FontSize',18,'FontWeight','bold');
ylabel('PC3','FontSize',18,'FontWeight','bold');
axis square
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
grid off;
title('Neuronal Trajectories in PCA Space','FontSize',20,'FontWeight','bold');
legend([plot1, plot2], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
legend('boxoff')
hold off

%% PCA 3D plots

trajectory_novel_3D = score(1:trace_duration*FrameRate, 1:3);
trajectory_familiar_3D   = score(trace_duration*FrameRate+1:end, 1:3);

fig4 = figure(4);
set(fig4, 'Position', [1 49 1080 748]);

plot1 = plot3(trajectory_novel_3D(:,1), trajectory_novel_3D(:,2), trajectory_novel_3D(:,3),'-o','color','#547DB1','linewidth',3,'markersize',4,'markeredgecolor','#547DB1','markerfacecolor','#547DB1');
hold on;
plot2 = plot3(trajectory_familiar_3D(:,1), trajectory_familiar_3D(:,2), trajectory_familiar_3D(:,3),'-o','color',[64 145 48]./255,'linewidth',3,'markersize',4,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255);

xlim([-10 20])
ylim([-10 10])
zlim([-10 10])
xlabel('PC1','FontSize',18,'FontWeight','bold');
ylabel('PC2','FontSize',18,'FontWeight','bold');
zlabel('PC3','FontSize',18,'FontWeight','bold');
axis square
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
grid on;
title('Neural Trajectories in PCA Space (3D)','FontSize',20,'FontWeight','bold');
lgd = legend([plot1, plot2], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
set(lgd, 'Position', [0.6960 0.7733 0.0974 0.0755]);
legend('boxoff')
view(-30, 60);
hold off

%% PCA 3D dynamic plots MP4

video_filename = 'neural_trajectory.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = FrameRate;
open(v);

fig5 = figure(5);
set(fig5, 'Position', [1 49 1080 748]);

hold on;
xlim([-10 20])
ylim([-10 10])
zlim([-10 10])
xlabel('PC1','FontSize',18,'FontWeight','bold');
ylabel('PC2','FontSize',18,'FontWeight','bold');
zlabel('PC3','FontSize',18,'FontWeight','bold');
grid on;
axis square
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

title('Neural Trajectories','FontSize',20,'FontWeight','bold');

view(-30, 60);

line_novel = plot3(NaN, NaN, NaN,'-o','color','#547DB1','linewidth',3,'markersize',4,'markeredgecolor','#547DB1','markerfacecolor','#547DB1');
line_familiar   = plot3(NaN, NaN, NaN,'-o','color',[64 145 48]./255,'linewidth',3,'markersize',4,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255);

time_label = text(-3, 16, 10, '', 'FontSize',18,'Color','k','FontWeight','bold','HorizontalAlignment','center');

lgd = legend([line_novel, line_familiar], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
set(lgd, 'Position', [0.6960 0.7733 0.0974 0.0755]);
legend('boxoff')

for t = 1:trace_duration*FrameRate/10
    set(line_novel, 'XData', trajectory_novel_3D(1:10*t,1), 'YData', trajectory_novel_3D(1:10*t,2), 'ZData', trajectory_novel_3D(1:10*t,3));

    set(line_familiar,   'XData', trajectory_familiar_3D(1:10*t,1), 'YData', trajectory_familiar_3D(1:10*t,2), 'ZData', trajectory_familiar_3D(1:10*t,3));

    set(time_label, 'String', sprintf('Time: %d s', round(t/2)));

    drawnow;

    frame = getframe(gcf);
    writeVideo(v, frame);

end

close(v);

%% PCA 3D dynamic plots GIF

filename = 'neural_trajectory.gif';

fig6 = figure(6);
set(fig6, 'Position', [1 49 1080 748]);

hold on;
xlim([-10 20])
ylim([-10 10])
zlim([-10 10])
xlabel('PC1','FontSize',18,'FontWeight','bold');
ylabel('PC2','FontSize',18,'FontWeight','bold');
zlabel('PC3','FontSize',18,'FontWeight','bold');
grid on;
axis square
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

title('Neural Trajectories','FontSize',20,'FontWeight','bold');

view(-30, 60);

line_novel = plot3(NaN, NaN, NaN,'-o','color','#547DB1','linewidth',3,'markersize',4,'markeredgecolor','#547DB1','markerfacecolor','#547DB1');
line_familiar   = plot3(NaN, NaN, NaN,'-o','color',[64 145 48]./255,'linewidth',3,'markersize',4,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255);

time_label = text(-3, 16, 10, '', 'FontSize',18,'Color','k','FontWeight','bold','HorizontalAlignment','center');

lgd = legend([line_novel, line_familiar], {'Novel', 'Familiar'},'FontSize',15,'FontWeight','bold');
set(lgd, 'Position', [0.6960 0.7733 0.0974 0.0755]);
legend('boxoff')

for t = 1:trace_duration*FrameRate/10
    set(line_novel, 'XData', trajectory_novel_3D(1:10*t,1), 'YData', trajectory_novel_3D(1:10*t,2), 'ZData', trajectory_novel_3D(1:10*t,3));

    set(line_familiar,   'XData', trajectory_familiar_3D(1:10*t,1), 'YData', trajectory_familiar_3D(1:10*t,2), 'ZData', trajectory_familiar_3D(1:10*t,3));

    set(time_label, 'String', sprintf('Time: %d s', round(t/2)));

    drawnow;

    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if t == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1/FrameRate);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/FrameRate);
    end

end


