clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

HR_folder = {'m1027\Jul_07_2024';
    'm37\Jul_08_2024';
    'm38\Jul_08_2024';
    'm40\Jul_08_2024';
    'm41\Jul_08_2024';
    'm53\Jul_08_2024';
    'm54\Jul_08_2024';
    'm63\Jul_09_2024';
    'm50\Jul_09_2024';
    'm59\Jul_09_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

trace_duration = 62;   % 90 seconds.

trials_num = size(HR_folder,1);

trials = NaN(trials_num,trace_duration*Sample_Rate);

fpass_trials=[11 13.5;
    10.5 13;
    11 13;
    11 12.5;
    11 13;
    9.5 13;
    11.5 14;
    10 13;
    11 13;
    11 13];

smooth_window = 5;


for I=1:size(HR_folder,1)

    Data_Folder = [Directory HR_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    tone_onsets = step_timepoint([3 4 5]);

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

    heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",[0 100]);

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm_outlier,1:size(datas,1),'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*Sample_Rate 0]);

    %%

    HR_trials = NaN(length(tone_onsets),trace_duration*Sample_Rate);


    for II = 1:length(tone_onsets)

        t_start = round(tone_onsets(II)*Sample_Rate) - 30*Sample_Rate +1;
        t_end = round(tone_onsets(II)*Sample_Rate) + 32*Sample_Rate;

        raw_trial = heartRate_bpm_interp_smooth(t_start:t_end);

        HR_trial_zscored = (raw_trial - mean(raw_trial(1:30*Sample_Rate))) / std(raw_trial(1:30*Sample_Rate));

        %         HR_trials(II,:) = HR_trial_zscored;

        HR_trials(II,:) = raw_trial;

    end

    trials(I,:) = mean(HR_trials);

end


%%

xlims = (-30*Sample_Rate+1:32*Sample_Rate)/Sample_Rate;

figure(1);

hold on

patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(trials),'Color',[255 128 128]./255,'LineWidth',2)
line([0,0],[690, 745],'Color','k','linestyle','--','LineWidth',2);

title({'Heart Rate',''},'FontSize',20,'FontWeight','bold')
ylabel('bpm','FontSize',15,'FontWeight','bold');

xlim([-30 32])
ylim([690 745])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

