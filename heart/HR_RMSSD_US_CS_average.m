clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

ECG_folder = {'m485_L_VTA_R_NAc_oxLight\Apr_05_2024';
    'm487_L_VTA_nLightR\Apr_05_2024';
    'm1772\Apr_05_2024';
    'm485_L_VTA_R_NAc_oxLight\Apr_08_2024';
    'm487_L_VTA_nLightR\Apr_08_2024';
    'm485_L_VTA_R_NAc_oxLight\Apr_10_2024';
    'm486_L_NAc_nLightR\Apr_10_2024';
    'm1772\Apr_10_2024'};


Sample_Rate = 200;    % 200 scans per second.

trace_duration = 60;   % 60 seconds.

US_HR = [];
CS_HR = [];

US_RMSSD = [];
CS_RMSSD = [];

%%
for I=1:size(ECG_folder,1)

    Data_Folder = [Directory ECG_folder{I} '\'];

    % Load data

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    % Analyse ECG signals


    US_start = ceil(step_timepoint(3)*Sample_Rate)+1;
    US_end = US_start+trace_duration*Sample_Rate-1;

    CS_start = ceil(step_timepoint(4)*Sample_Rate)+1;
    CS_end = CS_start+trace_duration*Sample_Rate-1;

    % Raw ECG signals

    US_ECG_raw = datas(US_start:US_end,2)';
    CS_ECG_raw = datas(CS_start:CS_end,2)';

    % Remove baseline wandering

    fpass=[10 13];

    US_ECG_Bandpass = bandpass(US_ECG_raw,fpass,Sample_Rate);

    CS_ECG_Bandpass = bandpass(CS_ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.007;

    [US_pksVal, US_pksLocs]=findpeaks(US_ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    [CS_pksVal, CS_pksLocs]=findpeaks(CS_ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    US_RR_intervals = diff(US_pksLocs);
    US_heartRate=1./US_RR_intervals;
    US_heartRate_bpm=US_heartRate*60;

    CS_RR_intervals = diff(CS_pksLocs);
    CS_heartRate=1./CS_RR_intervals;
    CS_heartRate_bpm=CS_heartRate*60;

    % Analyse running signals

    US_HR = [US_HR; mean(US_heartRate_bpm)];

    US_RMSSD = [US_RMSSD; sqrt(mean((diff(US_RR_intervals)).^2))];

    CS_HR = [CS_HR; mean(CS_heartRate_bpm)];

    CS_RMSSD = [CS_RMSSD; sqrt(mean((diff(CS_RR_intervals)).^2))];
end

%%

heart_rate_data = [US_HR; CS_HR];

state = [ones(8,1); repmat(2,8,1)];
state = categorical(state,[1 2],{'CS-','CS+'});

figure(1);
hold on

b1 = boxchart(state,heart_rate_data,'GroupByColor',state,'BoxWidth',1.5,'LineWidth',3.5,'MarkerStyle','none');

HR = [US_HR CS_HR];

for k = 1:size(HR,1)
    plot([0.8,2.2],HR(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

ylim([660 780])
ylabel('HR (bpm)','FontSize',20,'FontWeight','bold');
title({'Heart Rate',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_HR, p_HR, ci_HR, stats_HR] = ttest(US_HR, CS_HR, 'Tail','right')

%%

heart_rate_variability = [US_RMSSD; CS_RMSSD].*1000;

state = [ones(8,1); repmat(2,8,1)];
state = categorical(state,[1 2],{'CS-','CS+'});

figure;
hold on

b2 = boxchart(state,heart_rate_variability,'GroupByColor',state,'BoxWidth',1.5,'LineWidth',3.5,'MarkerStyle','none');

HRV = [US_RMSSD CS_RMSSD].*1000;

for k = 1:size(HRV,1)
    plot([0.8,2.2],HRV(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

ylim([3 11])
ylabel('RMSSD (ms)','FontSize',20,'FontWeight','bold');
title({'Heart Rate Variability','',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_HRV, p_HRV, ci_HRV, stats_HRV] = ttest(US_RMSSD, CS_RMSSD)




