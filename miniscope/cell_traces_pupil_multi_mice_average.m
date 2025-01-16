%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

folder = {'m2070\Jan_11_2025';
    'm2071\Jan_11_2025';
    'm2072\Jan_11_2025';
    'm2151\Jan_11_2025';
    'm2152\Jan_11_2025';
    'm2154\Jan_11_2025'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

outlier_pupil=[0 100;
    0 99.99;
    0 99.99;
    0 99.97;
    0 99.98;
    0 99.99];

PB_onsets = [960 1590 1053 1844 1495 1952];

chow_onsets = [1846 948 1790 1108 2411 1116];

pre_bout_duration = 20;   % 10 seconds.

trace_duration = 120;   % 60 seconds.

smooth_window = 1;

PB_Pupil = [];
chow_Pupil = [];

PB_Cell = [];
chow_Cell = [];

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

    %% Analyse pupil size

    Pupil_up_x = pupil_data(vid_onset:min(end,vid_offset),1);
    Pupil_up_y = pupil_data(vid_onset:min(end,vid_offset),2);
    Pupil_left_x = pupil_data(vid_onset:min(end,vid_offset),4);
    Pupil_left_y = pupil_data(vid_onset:min(end,vid_offset),5);
    Pupil_down_x = pupil_data(vid_onset:min(end,vid_offset),7);
    Pupil_down_y = pupil_data(vid_onset:min(end,vid_offset),8);
    Pupil_right_x = pupil_data(vid_onset:min(end,vid_offset),10);
    Pupil_right_y = pupil_data(vid_onset:min(end,vid_offset),11);

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

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    %% Running bouts

    PB_t_start = round(PB_onsets(I)*FrameRate) - pre_bout_duration*FrameRate +1;
    PB_t_end = PB_t_start + trace_duration*FrameRate -1;

    raw_pupil = pupil_smooth(PB_t_start:PB_t_end)';
    pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
    PB_Pupil = [PB_Pupil; pupil_zscored];

    cell_traces_bouts = cell_traces_interpolated_smooth(:,PB_t_start:PB_t_end);
            cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),2))./std(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),0,2);
    %         cell_traces_bouts_zscored = zscore(cell_traces_bouts,0,2);
    PB_Cell = [PB_Cell; mean(cell_traces_bouts_zscored)];

    chow_t_start = round(chow_onsets(I)*FrameRate) - pre_bout_duration*FrameRate +1;
    chow_t_end = chow_t_start + trace_duration*FrameRate -1;

    raw_pupil = pupil_smooth(chow_t_start:chow_t_end)';
    pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
    chow_Pupil = [chow_Pupil; pupil_zscored];

    cell_traces_bouts = cell_traces_interpolated_smooth(:,chow_t_start:chow_t_end);
            cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),2))./std(cell_traces_bouts(:,1:pre_bout_duration*FrameRate),0,2);
    %         cell_traces_bouts_zscored = zscore(cell_traces_bouts,0,2);
    chow_Cell = [chow_Cell; mean(cell_traces_bouts_zscored)];

end

%% plot figure

xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

trials_num = size(PB_Pupil,1);

traces_num = size(PB_Cell,1);

fig = figure(1);
set(fig, 'Position', [1 1 1400 1300]);

ax1 = subplot(2,2,1);
set(ax1, 'Position', [0.1300 0.5838 0.3347 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(chow_Cell)+std(chow_Cell)/sqrt(traces_num) fliplr(mean(chow_Cell)-std(chow_Cell)/sqrt(traces_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',1);
plot1 = plot(xlims,mean(chow_Cell),'Color','#547DB1','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(PB_Cell)+std(PB_Cell)/sqrt(traces_num) fliplr(mean(PB_Cell)-std(PB_Cell)/sqrt(traces_num))],'EdgeColor','none','FaceColor',[211 229 209]./255,'FaceAlpha',1);
plot2 = plot(xlims,mean(PB_Cell),'Color',[64 145 48]./255,'LineWidth',2);
line([0 0],[-5 20],'Color','k','linestyle','--','LineWidth',2);
ylim([-5 20])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Orexin Neural Activities','FontSize',20,'FontWeight','bold')
ylabel('z-score (s.d.)','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'PB', 'Chow'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

ax3 = subplot(2,2,3);
set(ax3, 'Position', [0.1300 0.1100 0.3347 0.3412]);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(chow_Pupil)+std(chow_Pupil)/sqrt(trials_num) fliplr(mean(chow_Pupil)-std(chow_Pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',1);
plot1 = plot(xlims,mean(chow_Pupil),'Color',[139 92 158]./255,'LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(PB_Pupil)+std(PB_Pupil)/sqrt(trials_num) fliplr(mean(PB_Pupil)-std(PB_Pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',1);
plot2 = plot(xlims,mean(PB_Pupil),'Color',[229 114 190]./255,'LineWidth',2);
line([0 0],[-5 30],'Color','k','linestyle','--','LineWidth',2);
ylim([-5 30])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Pupil Size Dynamics','FontSize',20,'FontWeight','bold')
ylabel('z-score (s.d.)','FontSize',18,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',18,'FontWeight','bold')
legend([plot2, plot1], {'PB', 'Chow'},'FontSize',15,'FontWeight','bold');
% legend('Location','bestoutside')
legend('boxoff')
hold off

%%

peak_PB_Cell = max(PB_Cell,[],2);

peak_chow_Cell = max(chow_Cell,[],2);

peak_Cell = [peak_PB_Cell peak_chow_Cell];

figure(1);

ax2 = subplot(2,2,2);
set(ax2, 'Position', [0.5703 0.5838 0.16 0.3412]);
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

for k = 1:size(peak_chow_Cell,1)
    plot(1.7,peak_chow_Cell(k),'marker','o','markersize',5,...
        'markeredgecolor','#547DB1','markerfacecolor','#547DB1',...
        'linestyle','none');
end

[S_PB,M_PB] = std(peak_PB_Cell,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_PB_Cell,1));

[S_chow,M_chow] = std(peak_chow_Cell,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_chow_Cell,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[64 145 48]./255,'LineWidth',3);
errorbar(1.9, M_chow, SEM_chow, "Color",'#547DB1','LineWidth',3);
plot(1.1, M_PB,'marker','o','color',[64 145 48]./255,'linewidth',4,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,'markersize',6);
plot(1.9, M_chow,'marker','o','color','#547DB1','linewidth',4,'markeredgecolor','#547DB1','markerfacecolor','#547DB1','markersize',6);

line([1.3 1.7], [32, 32], 'Color', 'k', 'LineWidth', 2);
text(1.5, 32.5, '***', 'FontSize', 25, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -2, 'PB', 'Color', [64 145 48]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -2, 'Chow', 'Color', '#547DB1', 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([0 35])
xlim([0.9 2.1])

ylabel('Max z-score (s.d.)','FontSize',18,'FontWeight','bold');

hold off

[h_Cell_PB, p_Cell_PB, ~, stats_Cell_PB] = ttest(peak_PB_Cell,0,'Tail','right')
[h_Cell_chow, p_Cell_chow, ~, stats_Cell_chow] = ttest(peak_chow_Cell,0,'Tail','right')
[h_Cell, p_Cell, ~, stats_Cell] = ttest(peak_PB_Cell,peak_chow_Cell,'Tail','right')


%%

peak_PB_Pupil = max(PB_Pupil,[],2);

peak_chow_Pupil = max(chow_Pupil,[],2);

peak_Pupil = [peak_PB_Pupil peak_chow_Pupil];

figure(1);

ax4 = subplot(2,2,4);
set(ax4, 'Position', [0.5703 0.1100 0.16 0.3412]);
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

for k = 1:size(peak_chow_Pupil,1)
    plot(1.7,peak_chow_Pupil(k),'marker','o','markersize',5,...
        'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,...
        'linestyle','none');
end

[S_PB,M_PB] = std(peak_PB_Pupil,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_PB_Pupil,1));

[S_chow,M_chow] = std(peak_chow_Pupil,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_chow_Pupil,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[229 114 190]./255,'LineWidth',3);
errorbar(1.9, M_chow, SEM_chow, "Color",[139 92 158]./255,'LineWidth',3);
plot(1.1, M_PB,'marker','o','color',[229 114 190]./255,'linewidth',4,'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,'markersize',6);
plot(1.9, M_chow,'marker','o','color',[139 92 158]./255,'linewidth',4,'markeredgecolor',[139 92 158]./255,'markerfacecolor',[139 92 158]./255,'markersize',6);

line([1.3 1.7], [45, 45], 'Color', 'k', 'LineWidth', 2);
text(1.5, 46.3, 'NS', 'FontSize', 15,'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -2, 'PB', 'Color', [229 114 190]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -2, 'Chow', 'Color', [139 92 158]./255, 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([0 50])
xlim([0.9 2.1])

ylabel('Max z-score (s.d.)','FontSize',18,'FontWeight','bold');

hold off

[h_Pupil, p_Pupil, ~, stats_Pupil] = ttest(peak_PB_Pupil,peak_chow_Pupil,'Tail','right')

%%

saveas(gcf, 'PB_Chow_Orexin_Pupil.svg')

saveas(gcf, 'PB_Chow_Orexin_Pupil.png')
