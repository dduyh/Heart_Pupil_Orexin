clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

ECG_folder = {'m1746\Mar_05_2024\session_1';
    'm1746\Mar_06_2024\session_1';
    'm1746\Mar_07_2024\session_1';
    'm1747\Mar_05_2024\session_1';
    'm1747\Mar_06_2024\session_1';
    'm1747\Mar_07_2024\session_1';
    'm1748\Mar_05_2024\session_1';
    'm1748\Mar_07_2024\session_1';
    'm1749\Mar_05_2024\session_1';
    'm1750\Mar_05_2024\session_1';
    'm1750\Mar_06_2024\session_1';
    'm1750\Mar_07_2024\session_1';
    'm1772\Mar_05_2024\session_1';
    'm1772\Mar_07_2024\session_1';
    'm1773\Mar_05_2024\session_1';
    'm1773\Mar_06_2024\session_1'};

running_heartRate = [];
resting_heartRate = [];

running_RMSSD = [];
resting_RMSSD = [];

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

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

trace_end = vid_start+trace_duration*Sample_Rate-1;

% Raw ECG signals

ECG_raw = datas(vid_start:trace_end,2)';

% Remove baseline wandering

fpass=[7 14];

ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

% find peaks

minPeakPromVal=0.01;

[pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

% heart rate in time

RR_intervals = diff(pksLocs);
heartRate=1./RR_intervals;
heartRate_smooth = movmean(heartRate,100);

% Analyse running signals

running = datas(vid_start:trace_end,1)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
Abs_speedDeg = abs(speedDeg);
speed = movmean(Abs_speedDeg,500);

speedThreshold = 0.002;
running_index = zeros(1,trace_duration*Sample_Rate);
running_index(speed > speedThreshold) = 1;

heartbeat_index = round(pksLocs * 1000);
movement_state = running_index(heartbeat_index);
movement_state(1) = [];

running_RR_intervals = RR_intervals(movement_state==1);
resting_RR_intervals = RR_intervals(movement_state==0);

running_heartRate = [running_heartRate; mean(1./running_RR_intervals)];
resting_heartRate = [resting_heartRate; mean(1./resting_RR_intervals)];

running_RMSSD = [running_RMSSD; sqrt(mean((diff(running_RR_intervals)).^2))];
resting_RMSSD = [resting_RMSSD; sqrt(mean((diff(resting_RR_intervals)).^2))];

end

%%

heart_rate_data = [resting_heartRate; running_heartRate];

state = [ones(16,1); repmat(2,16,1)];
state = categorical(state,[1 2],{'Resting','Running'});

figure(1);
hold on

b1 = boxchart(state,heart_rate_data,'GroupByColor',state,'BoxWidth',1.5,'LineWidth',3.5,'MarkerStyle','none');

HR = [resting_heartRate running_heartRate];

for k = 1:size(HR,1)
    plot([0.8,2.2],HR(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

% ylim([9.5 15])
ylabel('HR (Hz)','FontSize',20,'FontWeight','bold');
title({'Heart Rate','',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2; 
ax.FontSize = 15; 
ax.FontWeight = 'bold'; 

hold off

[h_HR, p_HR, ci_HR, stats_HR] = ttest(resting_heartRate, running_heartRate)

%%

heart_rate_variability = [resting_RMSSD; running_RMSSD].*1000;

state = [ones(16,1); repmat(2,16,1)];
state = categorical(state,[1 2],{'Resting','Running'});

figure(2);
hold on

b2 = boxchart(state,heart_rate_variability,'GroupByColor',state,'BoxWidth',1.5,'LineWidth',3.5,'MarkerStyle','none');

HRV = [resting_RMSSD running_RMSSD].*1000;

for k = 1:size(HRV,1)
    plot([0.8,2.2],HRV(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

% ylim([0 50])
ylabel('RMSSD (ms)','FontSize',20,'FontWeight','bold');
title({'Heart Rate Variability','',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2; 
ax.FontSize = 15; 
ax.FontWeight = 'bold'; 

hold off

[h_HRV, p_HRV, ci_HRV, stats_HRV] = ttest(resting_RMSSD, running_RMSSD)


