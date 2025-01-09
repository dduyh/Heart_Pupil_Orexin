%% ECG data processing

% clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

ECG_folder = {'m55\Jul_01_2024';
    'm39\Jun_20_2024';
    'm64\Jun_24_2024';
    'm58\Jul_01_2024';
    'm65\Jun_24_2024';
    'm1028\Jun_27_2024';
    'm49\Jun_27_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

trace_duration = 2400;   % 1260 seconds.

pre_duration = 1100;

fpass_trials=[7.5 12;
    9 11;
    6.5 10.5;
    7 12;
    8 11;
    9 11;
    8 12];

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

    Injection_onset = find(datas(:,3),1);
    trace_start = Injection_onset-pre_duration*Sample_Rate;
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

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm,1:trace_duration*Sample_Rate,'nearest','extrap');
    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[5*Sample_Rate 0]);
    HR_trial_zscored = (heartRate_bpm_interp_smooth - mean(heartRate_bpm_interp_smooth(1:pre_duration*Sample_Rate))) / std(heartRate_bpm_interp_smooth(1:pre_duration*Sample_Rate));
    trials(I,:) = HR_trial_zscored;

end

%%

xlims = (-pre_duration*Sample_Rate+1:(trace_duration-pre_duration)*Sample_Rate)/Sample_Rate;

figure;

hold on
patch('XData',[0, 0, 60, 60],'YData',[-5 2 2 -5],'EdgeColor','none','FaceColor','#8030f5','FaceAlpha',0.2);
% patch('XData',[0, 0, 60, 60],'YData',[500 660 660 500],'EdgeColor','none','FaceColor','#F6D7B4');

patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',0.5);
plot(xlims,mean(trials),'Color','#547DB1','LineWidth',2)
% line([360,360],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);

% xticks(-200:200:800);
% yticks(400:50:780);

xlim([-pre_duration trace_duration-pre_duration])
ylim([-4 2])
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


