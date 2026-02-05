%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

sema_folder = {'m2070\Apr_15_2025';
    'm2072\Apr_15_2025';
    'm2151\Apr_15_2025';
    'm2152\Apr_15_2025';
    'm2154\Apr_15_2025'};

glu_folder = {'m2070\May_01_2025';
    'm2072\May_01_2025';
    'm2151\May_01_2025';
    'm2152\May_01_2025'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate_sema = 20;

FrameRate_glu = 15;

trace_duration_sema = 4500;

trace_duration_glu = 3600;

pre_duration = 1300;

smooth_window = 1;

Cell_sema = [];

Cell_glu = [];

%%

cells_num_sema = 0;

for I=1:size(sema_folder,1)

    Data_Folder = [Directory sema_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    Injection_onset = find(datas(data_onset:data_offset,4),1);

    if I==2
        trace_start = round(Injection_onset*10/Sample_Rate)-pre_duration*10;
        trace_end = trace_start+trace_duration_sema*10-1;
    else
        trace_start = round(Injection_onset*20/Sample_Rate)-pre_duration*20;
        trace_end = trace_start+trace_duration_sema*20-1;
    end

    Cell_Status_data = readtable([Data_Folder sema_folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder sema_folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);
    cell_times = cell_traces_data(trace_start:trace_end,1);
    cell_traces = cell_traces_data(trace_start:trace_end,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times = cell_times(~nanRows, :);
    cell_times = cell_times(:,1)-cell_times(1,1);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    frame_time = 0:(1/FrameRate_sema):trace_duration_sema-1/FrameRate_sema;

    % Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num_sema = cells_num_sema + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
        %         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
    end

    cell_traces_interpolated_smooth = movmean(cell_traces_interpolated,[smooth_window*FrameRate_sema 0],2);

    % Remove bleaching slope
    photometry_detrend_all = zeros(size(cell_traces_interpolated_smooth));

    lastN = 60000; % number of last points to use for tail correction

    for i = 1:numCells
        trace = cell_traces_interpolated_smooth(i, :);

        % Only use first pre_duration seconds for exponential fitting
        idx = frame_time <= pre_duration;
        t_fit = frame_time(idx);
        y_fit = trace(idx);

        % Exponential model: A*exp(-t/tau) + C
        F = @(x,t) x(1)*exp(-t./x(2)) + x(3);
        x0 = [max(y_fit)-min(y_fit), 1000, median(y_fit)];
        lb = [0, 1, -Inf];
        ub = [Inf, Inf, Inf];
        opts = optimoptions('lsqcurvefit', 'Display', 'off');

        x = lsqcurvefit(F, x0, t_fit, y_fit, lb, ub, opts);

        % Compute baseline
        baseline = F(x, frame_time);

        % --- Fix: tail anchoring (no 'end' misuse) ---
        N = length(frame_time);
        if N >= lastN
            tail_idx = N - lastN + 1 : N;
            shift = median(trace(tail_idx)) - median(baseline(tail_idx));
            baseline = baseline + shift;
        end

        % Detrend
        photometry_detrend_all(i, :) = trace - baseline;
    end

    % Average (or keep individual detrended traces)
    photometry_detrend = mean(photometry_detrend_all);
    Cell_sema = [Cell_sema; photometry_detrend];

end

%%

cells_num_glu = 0;

for I=1:size(glu_folder,1)

    Data_Folder = [Directory glu_folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    Injection_onset = find(datas(data_onset:data_offset,4),1);

    trace_start = round(Injection_onset*FrameRate_glu/Sample_Rate)-pre_duration*FrameRate_glu;
    trace_end = trace_start+trace_duration_glu*FrameRate_glu-1;

    Cell_Status_data = readtable([Data_Folder glu_folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder glu_folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);
    cell_times = cell_traces_data(trace_start:trace_end,1);
    cell_traces = cell_traces_data(trace_start:trace_end,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times = cell_times(~nanRows, :);
    cell_times = cell_times(:,1)-cell_times(1,1);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    cell_traces_accepted_zscore = zscore(cell_traces_accepted);

    frame_time = 0:(1/FrameRate_glu):trace_duration_glu-1/FrameRate_glu;

    % Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cells_num_glu = cells_num_glu + numCells;

    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted_zscore(:, i),frame_time,'nearest','extrap');
        %         cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
    end

    cell_traces_interpolated_smooth = movmean(cell_traces_interpolated,[smooth_window*FrameRate_glu 0],2);

    % Remove bleaching slope
    photometry_detrend_all = zeros(size(cell_traces_interpolated_smooth));

    lastN = 60000; % number of last points to use for tail correction

    for i = 1:numCells
        trace = cell_traces_interpolated_smooth(i, :);

        % Only use first pre_duration seconds for exponential fitting
        idx = frame_time <= pre_duration;
        t_fit = frame_time(idx);
        y_fit = trace(idx);

        % Exponential model: A*exp(-t/tau) + C
        F = @(x,t) x(1)*exp(-t./x(2)) + x(3);
        x0 = [max(y_fit)-min(y_fit), 1000, median(y_fit)];
        lb = [0, 1, -Inf];
        ub = [Inf, Inf, Inf];
        opts = optimoptions('lsqcurvefit', 'Display', 'off');

        x = lsqcurvefit(F, x0, t_fit, y_fit, lb, ub, opts);

        % Compute baseline
        baseline = F(x, frame_time);

        % --- Fix: tail anchoring (no 'end' misuse) ---
        N = length(frame_time);
        if N >= lastN
            tail_idx = N - lastN + 1 : N;
            shift = median(trace(tail_idx)) - median(baseline(tail_idx));
            baseline = baseline + shift;
        end

        % Detrend
        photometry_detrend_all(i, :) = trace - baseline;
    end

    % Average (or keep individual detrended traces)
    photometry_detrend = mean(photometry_detrend_all);
    Cell_glu = [Cell_glu; photometry_detrend];

end

%%

xlims = (-pre_duration*FrameRate_glu+1:(trace_duration_glu-pre_duration)*FrameRate_glu)/FrameRate_glu;

traces_num_glu = size(Cell_glu,1);

fig = figure(1);
set(fig, 'Position', [2655 665 2017 593]);

axes('Position', [0.1300 0.1100 0.2134 0.8099]);
hold on
patch('XData',[0, 0, 60, 60],'YData',[-2 4 4 -2],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',0.7);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Cell_glu)+std(Cell_glu)/sqrt(traces_num_glu) fliplr(mean(Cell_glu)-std(Cell_glu)/sqrt(traces_num_glu))],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',1);
plot(xlims,mean(Cell_glu),'Color',[64 145 48]./255,'LineWidth',2)

xlim([-pre_duration trace_duration_glu-pre_duration])
ylim([-2 2])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Glucose','FontSize',18,'FontWeight','bold')
ylabel('Orx-GCaMP6s (z-score)','FontSize',15,'FontWeight','bold');
hold off


%%

xlims = (-pre_duration*FrameRate_sema+1:(trace_duration_sema-pre_duration)*FrameRate_sema)/FrameRate_sema;

traces_num_sema = size(Cell_sema,1);

figure(1);

axes('Position', [0.4108 0.1100 0.2134 0.8099]);
hold on
patch('XData',[0, 0, 60, 60],'YData',[-2 4 4 -2],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',0.7);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Cell_sema)+std(Cell_sema)/sqrt(traces_num_sema) fliplr(mean(Cell_sema)-std(Cell_sema)/sqrt(traces_num_sema))],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',1);
plot(xlims,mean(Cell_sema),'Color',[64 145 48]./255,'LineWidth',2)

xlim([-pre_duration trace_duration_sema-pre_duration])
ylim([-1 3])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Semaglutide','FontSize',18,'FontWeight','bold')
ylabel('Orx-GCaMP6s (z-score)','FontSize',15,'FontWeight','bold');
hold off

%%

% peak_Cell_glu = mean(Cell_glu(:,1360*FrameRate_glu+1:2360*FrameRate_glu),2) - mean(Cell_glu(:,pre_duration*FrameRate_glu-1000*FrameRate_glu+1:pre_duration*FrameRate_glu),2);
% 
% peak_Cell_sema = mean(Cell_sema(:,end-1000*FrameRate_sema+1:end),2) - mean(Cell_sema(:,pre_duration*FrameRate_sema-1000*FrameRate_sema+1:pre_duration*FrameRate_sema),2);

% peak_Cell_glu = mean(Cell_glu(:,1360*FrameRate_glu+1:end),2) - mean(Cell_glu(:,1:pre_duration*FrameRate_glu),2);
% 
% peak_Cell_sema = mean(Cell_sema(:,1360*FrameRate_sema+1:end),2) - mean(Cell_sema(:,1:pre_duration*FrameRate_sema),2);

peak_Cell_glu = min(Cell_glu(:,1360*FrameRate_glu+1:2360*FrameRate_glu),[],2) - min(Cell_glu(:,pre_duration*FrameRate_glu-1000*FrameRate_glu+1:pre_duration*FrameRate_glu),[],2);

peak_Cell_sema = min(Cell_sema(:,end-1000*FrameRate_sema+1:end),[],2) - min(Cell_sema(:,pre_duration*FrameRate_sema-1000*FrameRate_sema+1:pre_duration*FrameRate_sema),[],2);

Cell_mean = [mean(peak_Cell_glu) mean(peak_Cell_sema)];
Cell_sem = [std(peak_Cell_glu)/sqrt(traces_num_glu) std(peak_Cell_sema)/sqrt(traces_num_sema)];

figure(1);

axes('Position', [0.6916 0.1100 0.1 0.8099]);
hold on

b = bar(Cell_mean,'EdgeColor','none');
b.FaceColor = 'flat';
b.CData(1,:) = [222, 202, 224] / 255;
b.CData(2,:) = [198, 214, 234] / 255;

for k = 1:size(peak_Cell_glu,1)
    plot(1,peak_Cell_glu(k),'marker','o','markersize',5,...
        'markeredgecolor','#8B5C9E','markerfacecolor','#8B5C9E',...
        'linestyle','none');
end

for k = 1:size(peak_Cell_sema,1)
    plot(2,peak_Cell_sema(k),'marker','o','markersize',5,...
        'markeredgecolor','#719DC9','markerfacecolor','#719DC9',...
        'linestyle','none');
end

errorbar(1:2,Cell_mean,Cell_sem,'k','linestyle','none','linewidth',2,'CapSize',15);

% text(3, 1.25, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xticks([1 2])
xticklabels({'GLUC','SEMA'})

% ylim([0 1.4])
xlim([0.3 2.7])

title({'Δ Orexin'},'FontSize',18,'FontWeight','bold')
ylabel('Δ z-score (s.d.)','FontSize',15,'FontWeight','bold');

hold off

[h_Cell_glu_sema, p_Cell_glu_sema, ~, stats_Cell_glu_sema] = ttest2(peak_Cell_glu,peak_Cell_sema,'Tail','left')

