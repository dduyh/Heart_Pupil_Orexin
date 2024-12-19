clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

% HR_folder = {'m1826\May_10_2024';
%     'm1821\May_10_2024';
%     'm1822\May_10_2024';
%     'm1825\May_14_2024'};
% 
% ctrl_folder = {'m1840\May_14_2024';
%     'm1841\May_14_2024'};

HR_folder = {'m1825\May_23_2024';
    'm1826\May_23_2024';
    'm1821\May_24_2024_2';
    'm1822\May_24_2024'};

ctrl_folder = {'m1840\May_23_2024';
    'm1821\May_24_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

stimFreqs = 5;

trace_duration = 500;   % 90 seconds.

trials_num = size(HR_folder,1);

ctrl_trials_num = size(ctrl_folder,1);

trials = NaN(trials_num,trace_duration*Sample_Rate);

trials_ctrl = NaN(ctrl_trials_num,trace_duration*Sample_Rate);

fpass_trials=[8 11;
    9 10;
    9 10.5;
    9 10];

fpass_ctrl=[8.5 9;
    9 11];

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

    fpass=fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.007;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    %%

    HR_trials = NaN(length(idxFreqs),trace_duration*Sample_Rate);

    for II = 1:length(laser_onsets)

        t_start = laser_onsets(II)*Sample_Rate - 100*Sample_Rate +1;
        t_end = laser_onsets(II)*Sample_Rate + 400*Sample_Rate;

        heartbeat_index = find(pksLocs>(laser_onsets(II)-100) & pksLocs<(laser_onsets(II)+400));
        raw_trial = interp1(pksLocs(heartbeat_index)*Sample_Rate,heartRate_bpm(heartbeat_index-1),t_start:t_end,'linear','extrap');
        raw_trial_smooth = movmean(raw_trial,5*Sample_Rate);

        HR_trial_zscored = (raw_trial_smooth - mean(raw_trial_smooth(1:100*Sample_Rate))) / std(raw_trial_smooth(1:100*Sample_Rate));

        HR_trials(II,:) = HR_trial_zscored;

    end

    trials(I,:) = mean(HR_trials);

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

    minPeakPromVal=0.007;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    %%

    HR_trials = NaN(length(idxFreqs),trace_duration*Sample_Rate);

    for II = 1:length(laser_onsets)

        t_start = laser_onsets(II)*Sample_Rate - 100*Sample_Rate +1;
        t_end = laser_onsets(II)*Sample_Rate + 400*Sample_Rate;

        heartbeat_index = find(pksLocs>(laser_onsets(II)-100) & pksLocs<(laser_onsets(II)+400));
        raw_trial = interp1(pksLocs(heartbeat_index)*Sample_Rate,heartRate_bpm(heartbeat_index-1),t_start:t_end,'linear','extrap');
        raw_trial_smooth = movmean(raw_trial,5*Sample_Rate);

        HR_trial_zscored = (raw_trial_smooth - mean(raw_trial_smooth(1:100*Sample_Rate))) / std(raw_trial_smooth(1:100*Sample_Rate));

        HR_trials(II,:) = HR_trial_zscored;

    end

    trials_ctrl(I,:) = mean(HR_trials);

end

%%

xlims = (-100*Sample_Rate+1:400*Sample_Rate)/Sample_Rate;

figure(1);

hold on
patch('XData',[0, 0, 60, 60],'YData',[-12, 5, 5, -12],'EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(trials),'Color',[255 128 128]./255,'LineWidth',2)

title('Heart Rate (Chrimson)','FontSize',20,'FontWeight','bold')
ylabel('HR z-score','FontSize',15,'FontWeight','bold');

xlim([-100 400])
ylim([-12 5])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


figure(2);

hold on
patch('XData',[0, 0, 60, 60],'YData',[-12, 5, 5, -12],'EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials_ctrl)+std(trials_ctrl)/sqrt(ctrl_trials_num) fliplr(mean(trials_ctrl)-std(trials_ctrl)/sqrt(ctrl_trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',0.7);
plot(xlims,mean(trials_ctrl),'Color','k','LineWidth',2);

title('Heart Rate (Control)','FontSize',20,'FontWeight','bold')
ylabel('HR z-score','FontSize',15,'FontWeight','bold');

xlim([-100 400])
ylim([-12 5])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off