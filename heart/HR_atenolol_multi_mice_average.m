%% ECG data processing

clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';    % Main directory\

ECG_folder = {'m1747\Jun_04_2024';
    'm1749\Jun_04_2024';
    'm1772\Jun_04_2024';
    'm1774\Jun_04_2024'};

Injection_onsets = [400 400 410 420];

Sample_Rate = 1000;    % 1000 scans per second.

trace_duration = 1600;   % 1260 seconds.

fpass_trials=[8 12;
    8 12.5;
    11 13.5;
    9 13];

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
    trace_start = vid_start+Injection_onsets(I)*Sample_Rate-360*Sample_Rate;
    trace_end = trace_start+trace_duration*Sample_Rate-1;

    % Raw ECG signals

    ECG_raw = datas(trace_start:trace_end,2)';

    % Remove baseline wandering

    fpass=fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.007;
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

xlims = (-360*Sample_Rate+1:(trace_duration-360)*Sample_Rate)/Sample_Rate;

figure(1);

hold on
patch('XData',[0, 0, 60, 60],'YData',[400 780 780 400],'EdgeColor','none','FaceColor','#6DBB86','FaceAlpha',0.5);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#D1E8E6','FaceAlpha',0.7);
plot(xlims,mean(trials),'Color','#7CBEAE','LineWidth',2)
% line([360,360],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);

xticks(-200:200:1200);
yticks(400:50:780);

xlim([-360 trace_duration-360])
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

