% clc
% clear
% close all

%% set the path for data.

Directory = 'P:\Yihui\data\';    % Main directory\

HR_folder = {'m1821\Jun_13_2024';
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

fpass1_trials=[4 5.5;
    4 6;
    4 6;
    3 6;
    3 6;
    3.5 5.5;
    1.5 6.5;
    4 6;
    2 6;
    4 6;
    1 5;
    2 7];

fpass2_trials=[4 5.5;
    4 6;
    4 6;
    3 6;
    3 6;
    3.5 5.5;
    1.5 6.5;
    4 9;
    3 10;
    4 9;
    3 8;
    4 9];

splitIndex=[1000000 1000000 1000000 1000000 1000000 1000000 1000000 1170000 1190000 750000 1180000 1600000];

fpass_ctrl=[3 4;
    3 4.5;
    1 4;
    2 4];

Sample_Rate = 1000;    % 1000 scans per second.

stimFreqs = 20;

pre_bout_duration = 100;

trace_duration = 500;   % 90 seconds.

smooth_window = 5;

trials_num = size(HR_folder,1);

ctrl_trials_num = size(ctrl_folder,1);

HR = [];
HR_ctrl = [];

HRV = [];
HRV_ctrl = [];

RMSSD_on = [];
RMSSD_on_ctrl = [];

for I=1:size(HR_folder,1)

    Data_Folder = [Directory HR_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs)-step_timepoint(1);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass1=fpass1_trials(I,:);
    fpass2=fpass2_trials(I,:);

    splitIdx = splitIndex(I);

    bufferLen = round(300 * Sample_Rate);

    segment1 = ECG_raw(1:splitIdx + bufferLen);
    segment2 = ECG_raw(splitIdx - bufferLen + 1:end);

    ECG_seg1 = bandpass(segment1, fpass1, Sample_Rate);
    ECG_seg2 = bandpass(segment2, fpass2, Sample_Rate);

    valid_seg1 = ECG_seg1(1:splitIdx);
    valid_seg2 = ECG_seg2(bufferLen+1:end);

    ECG_Bandpass = [valid_seg1, valid_seg2];

    % find peaks

    if I>4
        minPeakPromVal=0.07;
    else
        minPeakPromVal=0.007;
    end

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.3,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",[0 100]);

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm_outlier,1:size(datas,1),'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*Sample_Rate 0]);

    %% HRV

    SD = abs(diff(RR_intervals)).*1000;
    SD_interp = interp1(pksLocs(3:end)*Sample_Rate,SD,1:size(datas,1),'nearest','extrap');
    SD_interp_smooth = movmean(SD_interp,[smooth_window*Sample_Rate 0]);

    %% Optogenetic stimulation bouts

    if I>4
        stim_num = length(idxFreqs);
    else
        stim_num = 2;
    end

    HR_trials = [];
    HRV_trials = [];
    laser_on_RMSSD_trials = [];

    for II = 1:stim_num

        t_start = round(laser_onsets(II)*Sample_Rate) - pre_bout_duration*Sample_Rate +1;
        t_end = t_start + trace_duration*Sample_Rate -1;

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR_trials = [HR_trials; raw_HR];

        raw_HRV = SD_interp_smooth(t_start:t_end);
        HRV_trials = [HRV_trials; raw_HRV];

        baseline_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)-60) & pksLocs(2:end)<laser_onsets(II));
        baseline_RR_intervals = RR_intervals(baseline_heartbeat_index);
        baseline_RMSSD = sqrt(mean((diff(baseline_RR_intervals)).^2));

        laser_on_heartbeat_index = find(pksLocs(2:end)>=laser_onsets(II) & pksLocs(2:end)<(laser_onsets(II)+60));
        laser_on_RR_intervals = RR_intervals(laser_on_heartbeat_index);
        laser_on_RMSSD_trials = [laser_on_RMSSD_trials; sqrt(mean((diff(laser_on_RR_intervals)).^2))-baseline_RMSSD];

    end

    HR = [HR; mean(HR_trials,1)];

    HRV = [HRV; mean(HRV_trials,1)];

    RMSSD_on = [RMSSD_on; mean(laser_on_RMSSD_trials,1)];

end

for I=1:size(ctrl_folder,1)

    Data_Folder = [Directory ctrl_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs)-step_timepoint(1);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_ctrl(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    if I>2
        minPeakPromVal=0.07;
    else
        minPeakPromVal=0.007;
    end

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.3,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",[0 100]);

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm_outlier,1:size(datas,1),'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*Sample_Rate 0]);

    %% HRV

    SD = abs(diff(RR_intervals)).*1000;
    SD_interp = interp1(pksLocs(3:end)*Sample_Rate,SD,1:size(datas,1),'nearest','extrap');
    SD_interp_smooth = movmean(SD_interp,[smooth_window*Sample_Rate 0]);

    %% Optogenetic stimulation bouts

    if I>2
        stim_num = length(idxFreqs);
    else
        stim_num = 2;
    end

    HR_trials = [];
    HRV_trials = [];
    laser_on_RMSSD_trials = [];

    for II = 1:stim_num

        t_start = round(laser_onsets(II)*Sample_Rate) - pre_bout_duration*Sample_Rate +1;
        t_end = t_start + trace_duration*Sample_Rate -1;

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR_trials = [HR_trials; raw_HR];

        raw_HRV = SD_interp_smooth(t_start:t_end);
        HRV_trials = [HRV_trials; raw_HRV];

        baseline_heartbeat_index = find(pksLocs(2:end)>=(laser_onsets(II)-60) & pksLocs(2:end)<laser_onsets(II));
        baseline_RR_intervals = RR_intervals(baseline_heartbeat_index);
        baseline_RMSSD = sqrt(mean((diff(baseline_RR_intervals)).^2));

        laser_on_heartbeat_index = find(pksLocs(2:end)>=laser_onsets(II) & pksLocs(2:end)<(laser_onsets(II)+60));
        laser_on_RR_intervals = RR_intervals(laser_on_heartbeat_index);
        laser_on_RMSSD_trials = [laser_on_RMSSD_trials; sqrt(mean((diff(laser_on_RR_intervals)).^2))-baseline_RMSSD];

    end

    HR_ctrl = [HR_ctrl; mean(HR_trials,1)];

    HRV_ctrl = [HRV_ctrl; mean(HRV_trials,1)];

    RMSSD_on_ctrl = [RMSSD_on_ctrl; mean(laser_on_RMSSD_trials,1)];

end

%%

xlims = (-pre_bout_duration*Sample_Rate+1:(trace_duration-pre_bout_duration)*Sample_Rate)/Sample_Rate;

figure(1);

subplot(4,2,5)
hold on
patch('XData',[0, 0, 60, 60],'YData',[110, 340, 340, 110],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR)+std(HR)/sqrt(trials_num) fliplr(mean(HR)-std(HR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',1);
plot(xlims,mean(HR),'Color',[255 128 128]./255,'LineWidth',2)
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([240 340])
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Heart Rate','FontSize',15,'FontWeight','bold')
ylabel('bpm','FontSize',12,'FontWeight','bold');
hold off

subplot(4,2,7);
hold on
patch('XData',[0, 0, 60, 60],'YData',[0, 80, 80, 0],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HRV)+std(HRV)/sqrt(trials_num) fliplr(mean(HRV)-std(HRV)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',1);
plot(xlims,mean(HRV),'Color',[255 128 128]./255,'LineWidth',2)
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([10 70])
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Heart Rate Variability','FontSize',15,'FontWeight','bold')
ylabel('SD (ms)','FontSize',12,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold')
hold off

%%

peak_HRV = mean(HR(:,100*Sample_Rate+1:160*Sample_Rate),2) - mean(HR(:,40*Sample_Rate+1:100*Sample_Rate),2);

peak_HRV_ctrl = mean(HR_ctrl(:,100*Sample_Rate+1:160*Sample_Rate),2) - mean(HR_ctrl(:,40*Sample_Rate+1:100*Sample_Rate),2);

figure(1);

axes('Position', [0.5703 0.3291 0.15 0.1561]);
hold on

for k = 1:size(peak_HRV,1)
    plot(1.3,peak_HRV(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HRV_ctrl,1)
    plot(1.7,peak_HRV_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HRV,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HRV,1));

[S_chow,M_chow] = std(peak_HRV_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HRV_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[0.5 0.5 0.5],'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[0.5 0.5 0.5],'linewidth',3,'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5);

line([1.3 1.7], [44, 44], 'Color', 'k', 'LineWidth', 2);
text(1.5, 46, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 26, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -23, 'Opsin', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -23, 'Control', 'Color', [0.5 0.5 0.5], 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-20 50])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Heart Rate (bpm)','FontSize',12,'FontWeight','bold');

hold off

[h_HR, p_HR, ~, stats_HR] = ttest2(peak_HRV,peak_HRV_ctrl,'Tail','right')
[h_HR_trials, p_HR_trials, ~, stats_HR_trials] = ttest(peak_HRV,0,'Tail','right')
[h_HR_ctrl, p_HR_ctrl, ~, stats_HR_ctrl] = ttest(peak_HRV_ctrl,0,'Tail','right')

%%

peak_HRV = RMSSD_on.*1000;

peak_HRV_ctrl = RMSSD_on_ctrl.*1000;

figure(1);

axes('Position', [0.5703 0.1100 0.15 0.1561]);
hold on

for k = 1:size(peak_HRV,1)
    plot(1.3,peak_HRV(k),'marker','o','markersize',4,...
        'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,...
        'linestyle','none');
end

for k = 1:size(peak_HRV_ctrl,1)
    plot(1.7,peak_HRV_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_HRV,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_HRV,1));

[S_chow,M_chow] = std(peak_HRV_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_HRV_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[255 128 128]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[0.5 0.5 0.5],'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[255 128 128]./255,'linewidth',3,'markeredgecolor',[255 128 128]./255,'markerfacecolor',[255 128 128]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[0.5 0.5 0.5],'linewidth',3,'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5);

line([1.3 1.7], [23, 23], 'Color', 'k', 'LineWidth', 2);
text(1.5, 25, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, -30, '**', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -53, 'Opsin', 'Color', [255 128 128]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -53, 'Control', 'Color', [0.5 0.5 0.5], 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-50 30])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ RMSSD (ms)','FontSize',12,'FontWeight','bold');

hold off

[h_HRV, p_HRV, ~, stats_HRV] = ttest2(peak_HRV,peak_HRV_ctrl,'Tail','left')
[h_HRV_trials, p_HRV_trials, ~, stats_HRV_trials] = ttest(peak_HRV,0,'Tail','left')
[h_HRV_ctrl, p_HRV_ctrl, ~, stats_HRV_ctrl] = ttest(peak_HRV_ctrl,0,'Tail','right')

%%

saveas(gcf, 'Orexin_OptoStim_Anaesthetic_Run_Pupil_HR_HRV.svg')

saveas(gcf, 'Orexin_OptoStim_Anaesthetic_Run_Pupil_HR_HRV.png')
