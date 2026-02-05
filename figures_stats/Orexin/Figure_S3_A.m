clc
clear
close all

%% set the path for data.

Directory = 'P:\Yihui\data\';    % Main directory\

breath_folder = {'m1821\Jun_13_2024';
    'm1822\Jun_14_2024';
    'm1825\Jun_15_2024';
    'm1826\Jun_15_2024';
    'm2178\Jun_19_2025';
    'm2179\Jun_19_2025';
    'm2144\Jun_19_2025';
    'm2145\Jun_19_2025';
    'm2146\Jun_19_2025';
    'm2140\Jun_19_2025';
    'm2148\Jun_19_2025';
    'm2169\Jun_19_2025'};

ctrl_folder = {'m1840\Jun_13_2024';
    'm1841\Jun_13_2024';
    'm2155\Jun_20_2025';
    'm2156\Jun_20_2025'};

FrameRate = 20;

stimFreqs = 20;

pre_bout_duration = 100;

trace_duration = 500;   % 90 seconds.

trials_num = size(breath_folder,1);

ctrl_trials_num = size(ctrl_folder,1);

breath = NaN(trials_num,trace_duration*FrameRate);

breath_ctrl = NaN(ctrl_trials_num,trace_duration*FrameRate);

MinPeakProminences=[0.3 0.05 0.2 0.06 0.3 0.25 0.3 0.3 0.2 0.2 0.2 0.08];

MinPeakProminences_ctrl=[0.1 0.04 0.12 0.05];

smooth_window = 5;

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

    MinPeakProminence=MinPeakProminences(I);

    [pksVal, pksLocs]=findpeaks([pixDiff(1); pixDiff],FrameRate,'MaxPeakWidth',0.2,'MinPeakProminence',MinPeakProminence,'Threshold',1e-4,'WidthReference','halfheight');

    BB_intervals = diff(pksLocs);
    breathRate=1./BB_intervals;

    breath_outlier = filloutliers(breathRate,"nearest","percentiles",[0 100]);

    breathRate_interp = interp1(pksLocs(2:end)*FrameRate,breath_outlier,1:size(breath_data,1),'nearest','extrap');

    breathRate_interp_smooth = movmean(breathRate_interp,[smooth_window*FrameRate 0]);

    %%

    if I>4
        stim_num = length(idxFreqs);
    else
        stim_num = 2;
    end
    
    breath_trials = NaN(stim_num,trace_duration*FrameRate);

    for II = 1:stim_num

        breath_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        breath_end = breath_start + trace_duration*FrameRate -1;

        raw_breath = breathRate_interp_smooth(breath_start:breath_end);

        % breath_zscored = (raw_breath - mean(raw_breath(1:pre_bout_duration*FrameRate))) / std(raw_breath(1:pre_bout_duration*FrameRate));

        breath_trials(II,:) = raw_breath;

    end

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

    MinPeakProminence=MinPeakProminences_ctrl(I);

    [pksVal, pksLocs]=findpeaks([pixDiff(1); pixDiff],FrameRate,'MaxPeakWidth',0.2,'MinPeakHeight',MinPeakProminence,'Threshold',1e-4,'WidthReference','halfheight');

    BB_intervals = diff(pksLocs);
    breathRate=1./BB_intervals;

    breath_outlier = filloutliers(breathRate,"nearest","percentiles",[0 100]);

    breathRate_interp = interp1(pksLocs(2:end)*FrameRate,breath_outlier,1:size(breath_data,1),'nearest','extrap');

    breathRate_interp_smooth = movmean(breathRate_interp,[smooth_window*FrameRate 0]);

    %%

    if I>2
        stim_num = length(idxFreqs);
    else
        stim_num = 2;
    end
    
    breath_trials = NaN(stim_num,trace_duration*FrameRate);

    for II = 1:stim_num

        breath_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        breath_end = breath_start + trace_duration*FrameRate -1;

        raw_breath = breathRate_interp_smooth(breath_start:breath_end);

        % breath_zscored = (raw_breath - mean(raw_breath(1:100*FrameRate))) / std(raw_breath(1:100*FrameRate));

        breath_trials(II,:) = raw_breath;

    end

    breath_ctrl(I,:) = mean(breath_trials);

end

%%

xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

fig = figure(1);
set(fig, 'Position', [2561 49 922 1315]);

subplot(4,2,1);
hold on
patch('XData',[0, 0, 60, 60],'YData',[1, 6.5, 6.5, 1],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(breath)+std(breath)/sqrt(trials_num) fliplr(mean(breath)-std(breath)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',1);
plot(xlims,mean(breath),'Color',[64 145 48]./255,'LineWidth',2)
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([2.5 5])
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title({'Respiratory Rate'},'FontSize',15,'FontWeight','bold')
ylabel('Hz','FontSize',12,'FontWeight','bold');
hold off

%%

peak_breath = mean(breath(:,100*FrameRate+1:160*FrameRate),2) - mean(breath(:,40*FrameRate+1:100*FrameRate),2);

peak_breath_ctrl = mean(breath_ctrl(:,100*FrameRate+1:160*FrameRate),2) - mean(breath_ctrl(:,40*FrameRate+1:100*FrameRate),2);

figure(1);

axes('Position', [0.5703 0.7673 0.15 0.1561]);
hold on

for k = 1:size(peak_breath,1)
    plot(1.3,peak_breath(k),'marker','o','markersize',4,...
        'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,...
        'linestyle','none');
end

for k = 1:size(peak_breath_ctrl,1)
    plot(1.7,peak_breath_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_breath,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_breath,1));

[S_chow,M_chow] = std(peak_breath_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_breath_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[64 145 48]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[0.5 0.5 0.5],'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[64 145 48]./255,'linewidth',3,'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[0.5 0.5 0.5],'linewidth',3,'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5);

line([1.3 1.7], [1.85, 1.85], 'Color', 'k', 'LineWidth', 2);
text(1.5, 1.9, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 1.1, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -0.29, 'Opsin', 'Color', [64 145 48]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -0.29, 'Control', 'Color', [0.5 0.5 0.5], 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-0.2 2])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Breath Rate (Hz)','FontSize',12,'FontWeight','bold');

hold off

[h_Breath, p_Breath, ~, stats_Breath] = ttest2(peak_breath,peak_breath_ctrl,'Tail','right')
[h_breath, p_breath, ~, stats_breath] = ttest(peak_breath,0,'Tail','right')
[h_breath_ctrl, p_breath_ctrl, ~, stats_breath_ctrl] = ttest(peak_breath_ctrl,0,'Tail','right')
