% clc
% clear
% close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

ECG_folder = {'m39\Aug_06_2024';
    'm1026\Aug_06_2024';
    'm1028\Aug_06_2024';
    'm49\Aug_06_2024';
    'm2126\Jan_22_2026';
    'm2161\Jan_22_2026';
    'm2199\Jan_22_2026';
    'm2200\Jan_22_2026';
    'm2214\Jan_23_2026';
    'm2215\Jan_23_2026';
    'm2228\Jan_23_2026'};

Sample_Rate = 1000;    % 1000 scans per second.

trace_duration = 30;   % 60 seconds.

fpass_trials=[9 13;
    9 13;
    8.5 12;
    10 13;
    9 13;
    6 12;
    10 13;
    8.5 13.5;
    10 13;
    8.5 13;
    10 13];

HR = NaN(size(ECG_folder,1),4);

%%
for I=1:size(ECG_folder,1)

    Data_Folder = [Directory ECG_folder{I} '\'];

    % Load data

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    % Analyse ECG signals

    if I==2||I==3||I==4

        US_1_start = ceil(step_timepoint(4)*Sample_Rate)+1;
        US_1_end = US_1_start+trace_duration*Sample_Rate-1;
        US_1_baseline = US_1_start-trace_duration*Sample_Rate;


        CS_1_start = ceil(step_timepoint(3)*Sample_Rate)+1;
        CS_1_end = CS_1_start+trace_duration*Sample_Rate-1;
        CS_1_baseline = CS_1_start-trace_duration*Sample_Rate;

        US_2_start = ceil(step_timepoint(6)*Sample_Rate)+1;
        US_2_end = US_2_start+trace_duration*Sample_Rate-1;
        US_2_baseline = US_2_start-trace_duration*Sample_Rate;

        CS_2_start = ceil(step_timepoint(5)*Sample_Rate)+1;
        CS_2_end = CS_2_start+trace_duration*Sample_Rate-1;
        CS_2_baseline = CS_2_start-trace_duration*Sample_Rate;

    else

        US_1_start = ceil(step_timepoint(3)*Sample_Rate)+1;
        US_1_end = US_1_start+trace_duration*Sample_Rate-1;
        US_1_baseline = US_1_start-trace_duration*Sample_Rate;

        CS_1_start = ceil(step_timepoint(4)*Sample_Rate)+1;
        CS_1_end = CS_1_start+trace_duration*Sample_Rate-1;
        CS_1_baseline = CS_1_start-trace_duration*Sample_Rate;

        US_2_start = ceil(step_timepoint(5)*Sample_Rate)+1;
        US_2_end = US_2_start+trace_duration*Sample_Rate-1;
        US_2_baseline = US_2_start-trace_duration*Sample_Rate;

        CS_2_start = ceil(step_timepoint(6)*Sample_Rate)+1;
        CS_2_end = CS_2_start+trace_duration*Sample_Rate-1;
        CS_2_baseline = CS_2_start-trace_duration*Sample_Rate;

    end

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass=fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    if I>4
        minPeakPromVal=0.07;
    else
        minPeakPromVal=0.007;
    end

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm,1:size(datas,1),'nearest','extrap');
    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[5*Sample_Rate 0]);

    US_1_HR_baseline = heartRate_bpm_interp_smooth(US_1_baseline:US_1_start-1);
    US_1_HR = heartRate_bpm_interp_smooth(US_1_start:US_1_end);

    CS_1_HR_baseline = heartRate_bpm_interp_smooth(CS_1_baseline:CS_1_start-1);
    CS_1_HR = heartRate_bpm_interp_smooth(CS_1_start:CS_1_end);

    US_2_HR_baseline = heartRate_bpm_interp_smooth(US_2_baseline:US_2_start-1);
    US_2_HR = heartRate_bpm_interp_smooth(US_2_start:US_2_end);

    CS_2_HR_baseline = heartRate_bpm_interp_smooth(CS_2_baseline:CS_2_start-1);
    CS_2_HR = heartRate_bpm_interp_smooth(CS_2_start:CS_2_end);

    % Analyse running signals

    HR(I,:) = [mean(US_1_HR)-mean(US_1_HR_baseline) mean(CS_1_HR)-mean(CS_1_HR_baseline) mean(US_2_HR)-mean(US_2_HR_baseline) mean(CS_2_HR)-mean(CS_2_HR_baseline)];

end

%%

figure
hold on

for k = 1:size(HR,1)
    plot(1:4,HR(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(HR,'omitnan');
SEM = S/sqrt(size(HR,1));
errorbar(1:4, M, SEM, "Color","black",'LineWidth',3);
plot(1:4, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:4)
xticklabels({'CS -','CS +','CS -','CS +'})

ylim([-60 100])
xlim([0.5 4.5])

ylabel('ΔHR (bpm)','FontSize',20,'FontWeight','bold');
title({'Heart Rate',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

%%

figure(2)
hold on

line([0.5 4.5],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S,M] = std(HR,'omitnan');
SEM = S/sqrt(size(HR,1));
errorbar(1:4, M, SEM, "Color",[45 127 184]/255,'LineWidth',3);
plot4 = plot(1:4, M,'o-','color',[45 127 184]/255,'linewidth',4,'markeredgecolor',[45 127 184]/255,'markerfacecolor',[45 127 184]/255,'markersize',4);

xticks(1:4)
xticklabels({'CS -','CS +','CS -','CS +'})

ylim([-10 60])
xlim([0.5 4.5])

ylabel('ΔHR (bpm)','FontSize',20,'FontWeight','bold');
title('Heart Rate on Cue Testing Phase','FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

legend([plot1, plot2, plot3, plot4], {'DTR+ Fasted', 'DTR+ Fed', 'DTR- Fasted', 'DTR- Fed'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
legend('Location','best')

hold off




