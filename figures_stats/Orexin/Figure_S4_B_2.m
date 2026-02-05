%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

folder = {'m55\Jul_01_2024';
    'm39\Jun_20_2024';
    'm64\Jun_24_2024';
    'm58\Jul_01_2024';
    'm65\Jun_24_2024';
    'm1028\Jun_27_2024';
    'm49\Jun_27_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

fpass_trials=[7.5 12;
    9 11;
    6.5 10.5;
    7 12;
    8 11;
    9 11;
    8 12];

outlier_pupil=[0 99.95;
    0.1 98.8;
    0.3 98.6;
    0 99.5;
    2 92;
    0.2 97.4;
    0.1 99.8];

pre_bout_duration = 5;   % 5 seconds.

trace_duration = 30;   % 30 seconds.

speed_thresh = 0.02;

Run = [];
HR = [];
Pupil = [];

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    Injection_onset = find(datas(:,3),1);
    vid_start = ceil(step_timepoint(1)*Sample_Rate)+1;
    trace_end = Injection_onset;

    timepoint = times(vid_start:trace_end,1)';
    time = timepoint(1,:)-timepoint(1,1);

    smooth_window = 1;

    %% Analyse running signals

    running = datas(vid_start:trace_end,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed = [speed(1) speed];

    %% Analyse pupil size

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(vid_start:trace_end,2)';

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

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm,0:time(end)*Sample_Rate,'nearest','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*Sample_Rate 0]);

    %% Running bouts

    binarized_speed = speed > speed_thresh;

    running_onsets = find(diff(binarized_speed)==1);

    running_onsets(running_onsets<pre_bout_duration*Sample_Rate)=[];
    running_onsets(running_onsets>(length(speed)-trace_duration*Sample_Rate+pre_bout_duration*Sample_Rate))=[];

    good_running_onsets = [];

    for i = 1:length(running_onsets)

        ii = sum(binarized_speed(running_onsets(i)-pre_bout_duration*Sample_Rate+1:running_onsets(i)));

        if ii == 0
            good_running_onsets = [good_running_onsets,running_onsets(i)];
        end
    end

    %     run = NaN(length(good_running_onsets),trace_duration*Sample_Rate);


    for i = 1:length(good_running_onsets)

        t_start = good_running_onsets(i) - pre_bout_duration*Sample_Rate +1;
        t_end = t_start + trace_duration*Sample_Rate -1;

        Run = [Run; speed(t_start:t_end)];

        HR = [HR; heartRate_bpm_interp_smooth(t_start:t_end)];

        pupil_start = round(time(good_running_onsets(i))*FrameRate) - pre_bout_duration*FrameRate +1;
        pupil_end = pupil_start + trace_duration*FrameRate -1;

        raw_pupil = pupil(pupil_start:pupil_end)';

        pupil_zscored = (raw_pupil - mean(raw_pupil(1:pre_bout_duration*FrameRate))) / std(raw_pupil(1:pre_bout_duration*FrameRate));
        Pupil = [Pupil; raw_pupil];

    end

end

%% k-means clustering running bouts

% run_input_zscore = zscore(Run);
% 
% km = kmeans(run_input_zscore,2);

km = kmeans(Run(:,1:15000),2);

%%
xlims = (-5*Sample_Rate+1:25*Sample_Rate)/Sample_Rate;
pupil_xlims = (-5*FrameRate+1:25*FrameRate)/FrameRate;

figure;

run = Run(km==1,:);
pupil = Pupil(km==1,:);
trials = HR(km==1,:);

trials_num = size(run,1);

subplot(3,1,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(run)+std(run)/sqrt(trials_num) fliplr(mean(run)-std(run)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',0.5);
plot(xlims,mean(run),'Color','#547DB1','LineWidth',2)

subplot(3,1,2)
hold on
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',0.5);
plot(pupil_xlims,mean(pupil),'Color','#547DB1','LineWidth',2)

subplot(3,1,3)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',0.5);
plot(xlims,mean(trials),'Color','#547DB1','LineWidth',2)

%%

run = Run(km==2,:);
pupil = Pupil(km==2,:);
trials = HR(km==2,:);

trials_num = size(run,1);

subplot(3,1,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(run)+std(run)/sqrt(trials_num) fliplr(mean(run)-std(run)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#D1E8E6','FaceAlpha',0.5);
plot(xlims,mean(run),'Color',[124 190 174]./255,'LineWidth',2)
line([0,0],[0,0.1],'Color','k','linestyle','--','LineWidth',2);
ylim([0 0.1])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Abs. Running Speed','FontSize',20,'FontWeight','bold')
% hold off

subplot(3,1,2)
hold on
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#D1E8E6','FaceAlpha',0.5);
plot(pupil_xlims,mean(pupil),'Color',[124 190 174]./255,'LineWidth',2)
line([0,0],[120,520],'Color','k','linestyle','--','LineWidth',2);
ylim([120 520])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Pupil Size','FontSize',20,'FontWeight','bold')
ylabel('z-score','FontSize',15,'FontWeight','bold');
hold off

subplot(3,1,3)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#D1E8E6','FaceAlpha',0.5);
plot(xlims,mean(trials),'Color',[124 190 174]./255,'LineWidth',2)
line([0,0],[500,710],'Color','k','linestyle','--','LineWidth',2);
ylim([500 710])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Heart Rate','FontSize',20,'FontWeight','bold')
ylabel('HR (bpm)','FontSize',15,'FontWeight','bold');
hold off
