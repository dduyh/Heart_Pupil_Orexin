%% ECG data processing

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

ECG_folder = {'m1746\Mar_05_2024\session_1';
    'm1747\Mar_07_2024\session_1';
    'm1748\Mar_05_2024\session_1';
    'm1750\Mar_05_2024\session_1';
    'm1750\Mar_07_2024\session_1';
    'm1772\Mar_07_2024\session_1';
    'm1773\Mar_05_2024\session_1'};

Injection_onsets = [363 442 483 481 357 454 547]-35;

Sample_Rate = 1000;    % 1000 scans per second.

trace_duration = 1200;   % 1260 seconds.

trials_num = size(ECG_folder,1);
trials = NaN(trials_num,trace_duration*Sample_Rate);    % 15 seconds per trial.

%%
for I=1:size(ECG_folder,1)

    Data_Folder = [Directory ECG_folder{I} '\'];
    % Load data

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    % Analyse ECG signals

    vid_start = ceil(step_timepoint(1)*Sample_Rate)+1;
    trace_start = vid_start+Injection_onsets(I)*Sample_Rate-325*Sample_Rate;
    trace_end = trace_start+trace_duration*Sample_Rate-1;

    % Raw ECG signals

    ECG_raw = datas(trace_start:trace_end,2)';

    % Remove baseline wandering

    fpass=[5 14];

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.015;
%     figure;
%     findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

%     figure
%     plot(pksLocs(2:end),heartRate_smooth)

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm,1:trace_duration*Sample_Rate,'spline','extrap');
    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,5*Sample_Rate);
    trials(I,:) = heartRate_bpm_interp_smooth;

end

%%

xlims = (-325*Sample_Rate+1:(trace_duration-325)*Sample_Rate)/Sample_Rate;

figure(1);

hold on
patch('XData',[0, 0, 60, 60],'YData',[400 780 780 400],'EdgeColor','none','FaceColor','#8030f5','FaceAlpha',0.2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[231 216 250]./255,'FaceAlpha',0.7);
plot(xlims,mean(trials),'Color','#9392BE','LineWidth',2)
% line([360,360],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);

xticks(-200:200:800);
yticks(400:50:780);

xlim([-325 trace_duration-325])
ylim([400 780])
% title({'Heart Rate Before/After Propranolol Injection',''},'FontSize',20,'FontWeight','bold','color','k')
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
ylabel('HR (bpm)','FontSize',15,'FontWeight','bold');

% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off


