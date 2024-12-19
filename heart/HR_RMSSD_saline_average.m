clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

ECG_folder = {'m1746\Mar_06_2024\session_1';
    'm1747\Mar_06_2024\session_1';
    'm1750\Mar_06_2024\session_1';
    'm1773\Mar_06_2024\session_1'};

pre_atropine_heartRate = [];
post_atropine_heartRate = [];

pre_atropine_RMSSD = [];
post_atropine_RMSSD = [];

%%
for I=1:size(ECG_folder,1)

    Data_Folder = [Directory ECG_folder{I} '\'];

    % Load data

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    % Analyse ECG signals

    Sample_Rate = 1000;    % 1000 scans per second.

    trace_duration = 360;   % 360 seconds.

    vid_start = ceil(step_timepoint(1))*Sample_Rate+1;
    pre_trace_end = vid_start+trace_duration*Sample_Rate-1;

    % Raw ECG signals

    pre_ECG_raw = datas(vid_start:pre_trace_end,2)';
    post_ECG_raw = datas((end-trace_duration*Sample_Rate+1):end,2)';

    % Remove baseline wandering

    pre_fpass=[7 14];
    pre_ECG_Bandpass = bandpass(pre_ECG_raw,pre_fpass,Sample_Rate);

    post_fpass=[7 14];
    post_ECG_Bandpass = bandpass(post_ECG_raw,post_fpass,Sample_Rate);

    % find peaks

    pre_minPeakPromVal=0.01;

    [pre_pksVal, pre_pksLocs]=findpeaks(pre_ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',pre_minPeakPromVal);

    post_minPeakPromVal=0.01;

    [post_pksVal, post_pksLocs]=findpeaks(post_ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',post_minPeakPromVal);

    % heart rate in time

    pre_RR_intervals = diff(pre_pksLocs);
    pre_heartRate=1./pre_RR_intervals;

    post_RR_intervals = diff(post_pksLocs);
    post_heartRate=1./post_RR_intervals;

    % Analyse running signals

    pre_running = datas(vid_start:pre_trace_end,1)';
    signedThreshold = 2^(32-1);
    pre_running(pre_running > signedThreshold) = pre_running(pre_running > signedThreshold) - 2^32;
    pre_speedDeg = diff(pre_running);
    pre_Abs_speedDeg = abs(pre_speedDeg);
    pre_speed = movmean(pre_Abs_speedDeg,500);

    speedThreshold = 0.002;
    pre_running_index = zeros(1,trace_duration*Sample_Rate);
    pre_running_index(pre_speed > speedThreshold) = 1;

    pre_heartbeat_index = round(pre_pksLocs * 1000);
    pre_movement_state = pre_running_index(pre_heartbeat_index);
    pre_movement_state(1) = [];

    pre_resting_RR_intervals = pre_RR_intervals(pre_movement_state==0);

    pre_atropine_heartRate = [pre_atropine_heartRate; mean(1./pre_resting_RR_intervals)];

    pre_atropine_RMSSD = [pre_atropine_RMSSD; sqrt(mean((diff(pre_resting_RR_intervals)).^2))];

    post_running = datas((end-trace_duration*Sample_Rate+1):end,1)';
    signedThreshold = 2^(32-1);
    post_running(post_running > signedThreshold) = post_running(post_running > signedThreshold) - 2^32;
    post_speedDeg = diff(post_running);
    post_Abs_speedDeg = abs(post_speedDeg);
    post_speed = movmean(post_Abs_speedDeg,500);

    speedThreshold = 0.002;
    post_running_index = zeros(1,trace_duration*Sample_Rate);
    post_running_index(post_speed > speedThreshold) = 1;

    post_heartbeat_index = round(post_pksLocs * 1000);
    post_movement_state = post_running_index(post_heartbeat_index);
    post_movement_state(1) = [];

    post_resting_RR_intervals = post_RR_intervals(post_movement_state==0);

    post_atropine_heartRate = [post_atropine_heartRate; mean(1./post_resting_RR_intervals)];

    post_atropine_RMSSD = [post_atropine_RMSSD; sqrt(mean((diff(post_resting_RR_intervals)).^2))];
end

%%

heart_rate_data = [pre_atropine_heartRate; post_atropine_heartRate];

state = [ones(4,1); repmat(2,4,1)];
state = categorical(state,[1 2],{'Pre','Post'});

figure;
hold on

b1 = boxchart(state,heart_rate_data,'GroupByColor',state,'BoxFaceAlpha',1,'BoxWidth',2,'LineWidth',5,'BoxMedianLineColor','black','BoxEdgeColor','none','MarkerStyle','none');
b1(1).BoxFaceColor = '#D5E4A8';
b1(2).BoxFaceColor = '#9BC985';

HR = [pre_atropine_heartRate post_atropine_heartRate];

for k = 1:size(HR,1)
    plot([0.8,2.2],HR(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

ylim([9.5 13])
ylabel('HR (Hz)','FontSize',20,'FontWeight','bold');
title({'Heart Rate Before/After Saline Injection','',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_HR, p_HR, ci_HR, stats_HR] = ttest(pre_atropine_heartRate, post_atropine_heartRate)

%%

heart_rate_variability = [pre_atropine_RMSSD; post_atropine_RMSSD].*1000;

state = [ones(4,1); repmat(2,4,1)];
state = categorical(state,[1 2],{'Pre','Post'});

figure;
hold on

b2 = boxchart(state,heart_rate_variability,'GroupByColor',state,'BoxFaceAlpha',1,'BoxWidth',2,'LineWidth',5,'BoxMedianLineColor','black','BoxEdgeColor','none','MarkerStyle','none');
b2(1).BoxFaceColor = '#D5E4A8';
b2(2).BoxFaceColor = '#9BC985';

HRV = [pre_atropine_RMSSD post_atropine_RMSSD].*1000;

for k = 1:size(HRV,1)
    plot([0.8,2.2],HRV(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

% ylim([0 50])
ylabel('RMSSD (ms)','FontSize',20,'FontWeight','bold');
title({'HRV Before/After Saline Injection','',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_HRV, p_HRV, ci_HRV, stats_HRV] = ttest(pre_atropine_RMSSD, post_atropine_RMSSD)

