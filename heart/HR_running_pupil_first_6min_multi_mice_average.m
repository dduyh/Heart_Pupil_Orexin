clc
clear
close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

folder = {'m1747\Mar_05_2024\session_1';
    'm1747\Mar_07_2024\session_1';
    'm1750\Mar_05_2024\session_1';
    'm1772\Mar_07_2024\session_1';
    'm1773\Mar_05_2024\session_1';
    'm1773\Mar_06_2024\session_1'};

Run_onsets = {[298.469]; % m1747\Mar_05_2024
    [37.312]; % m1747\Mar_07_2024
    [5.533 207.317 298.294]; % m1750\Mar_05_2024
    [203.897 301.08]; % m1772\Mar_07_2024
    [11.745 72.379 319.017]; % m1773\Mar_05_2024
    [51.123]};

trials_num = 0;
for j = 1:size(folder,1)
    trials_num = trials_num + numel(Run_onsets{j});
end

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

trace_duration = 360;   % 360 seconds.

run = NaN(trials_num,15*Sample_Rate);
trials = NaN(trials_num,15*Sample_Rate);    % 15 seconds per trial.
pupil = NaN(trials_num,15*FrameRate);

k = 1;

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    vid_start = ceil(step_timepoint(1)*Sample_Rate)+1;

    trace_end = vid_start+trace_duration*Sample_Rate-1;

    timepoint = times(vid_start:end,1)';
    time = timepoint(1,:)-timepoint(1,1);

    smooth_window = 0.5;

    %% Analyse running signals

    running = datas(vid_start:trace_end,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(vid_start:trace_end,2)';

    % Remove baseline wandering

    fpass=[9 13];

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.01;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_smooth = movmean(heartRate,100);


    %%
    Run_bout_onsets = Run_onsets{I};

    for II = 1:length(Run_bout_onsets)

        t_start = Run_bout_onsets(II)*Sample_Rate - 5*Sample_Rate +1;
        t_end = Run_bout_onsets(II)*Sample_Rate + 10*Sample_Rate;

        run(k,:) = speed(t_start:t_end);

        heartbeat_index = find(pksLocs>(Run_bout_onsets(II)-5) & pksLocs<(Run_bout_onsets(II)+10));
        raw_trial = interp1(pksLocs(heartbeat_index)*Sample_Rate,heartRate_bpm(heartbeat_index-1),t_start:t_end,'linear','extrap');

        % Standardize signals
        %         trial_zscored = (raw_trial - mean(raw_trial(1:5*Sample_Rate))) / std(raw_trial(1:5*Sample_Rate));
        %         trials(k,:) = trial_zscored;
        trials(k,:) = raw_trial;

        pupil_start = round(Run_bout_onsets(II)*FrameRate) - 5*FrameRate +1;
        pupil_end = round(Run_bout_onsets(II)*FrameRate) + 10*FrameRate;

        raw_pupil = areas(pupil_start:pupil_end);

        pupil_zscored = (raw_pupil - mean(raw_pupil(1:5*FrameRate))) / std(raw_pupil(1:5*FrameRate));
        pupil(k,:) = pupil_zscored;


        k = k+1;

    end

end

%%
figure;
xlims = (1:15*Sample_Rate)/Sample_Rate;
pupil_xlims = (1:15*FrameRate)/FrameRate;

subplot(3,1,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(run)+std(run)/sqrt(trials_num) fliplr(mean(run)-std(run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',0.7);
plot(xlims,mean(run),'Color','k','LineWidth',2)
line([5,5],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
axis off
title('Abs. Running Speed','FontSize',20,'FontWeight','bold')
hold off

subplot(3,1,2)
hold on
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2)
line([5,5],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
% axis off
title('Pupil Size','FontSize',20,'FontWeight','bold','color',[229 114 190]./255)
hold off

subplot(3,1,3)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(trials),'Color',[255 128 128]./255,'LineWidth',2)
line([5,5],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
% ylim([10.7 12.7])
title('Heart Rate','FontSize',20,'FontWeight','bold','color',[255 128 128]./255)
ylabel('HR (Hz)','FontSize',15,'FontWeight','bold');

% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


