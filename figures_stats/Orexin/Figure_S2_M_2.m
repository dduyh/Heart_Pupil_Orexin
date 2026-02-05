clc
clear
close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

folder = {'m55\Jun_25_2024';
    'm1026\Jun_21_2024';
    'm39\Jun_28_2024';
    'm1028\Jun_21_2024';
    'm49\Jun_21_2024';
    'm58\Jun_25_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

fpass_trials=[8 11;
    7 10.5;
    7.5 11.5;
    5 8.5;
    6.5 10.5;
    8.5 11];

outlier_pupil=[0.2 98.8;
    0 98.6;
    0.3 99.6;
    0.1 99.1;
    0 99;
    0.3 98.1];

pre_bout_duration = 5;   % 5 seconds.

trace_duration = 30;   % 30 seconds.

speed_thresh = 0.02;

Run = [];
HR = [];
Pupil = [];

xcorrwind = 5*Sample_Rate; % in samples

r_RunPup = [];
r_RunHR = [];
r_HRPup = [];

smooth_window = 1;

%%

k = 1;

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

    %% Analyse running signals

    running = datas(vid_start:trace_end,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed = [Abs_speedDeg(1) Abs_speedDeg];
    speed_smooth = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];

    %% Analyse pupil size

    pupil_areas = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

    video_time = (1/FrameRate):(1/FrameRate):(size(pupil_areas,1)/FrameRate);

    pupil_interp = interp1(video_time*Sample_Rate,pupil_areas',1:length(running),'nearest','extrap');

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

    %% Running bouts

    binarized_speed = speed_smooth > speed_thresh;

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

        raw_running = speed(t_start:t_end);
        raw_speed = speed_smooth(t_start:t_end);
        Run = [Run; raw_speed];

        raw_pupil = pupil_interp(t_start:t_end);
        Pupil = [Pupil; raw_pupil];

        raw_HR = heartRate_bpm_interp(t_start:t_end);
        HR = [HR; raw_HR];

        [r_RunPup(k,:),lags] = xcorr(raw_running-mean(raw_running),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');

        [r_RunHR(k,:),lags] = xcorr(raw_running-mean(raw_running),raw_HR-mean(raw_HR),xcorrwind,'coeff');

        [r_HRPup(k,:),lags] = xcorr(raw_HR-mean(raw_HR),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');

        k = k+1;

    end

end

%% k-means clustering running bouts

% run_input_zscore = zscore(Run);
% 
% km = kmeans(run_input_zscore,2);

km = kmeans(Run(:,1:15000),2);

%%

xlims = (-5*Sample_Rate:5*Sample_Rate)/Sample_Rate;

figure;

r_runpup = r_RunPup(km==1,:);
r_runHR = r_RunHR(km==1,:);
r_HRpup = r_HRPup(km==1,:);

trials_num = size(r_runpup,1);

subplot(3,1,1);
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runpup)+std(r_runpup)/sqrt(trials_num) fliplr(mean(r_runpup)-std(r_runpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',0.5);
plot(xlims,mean(r_runpup),'Color','#547DB1','LineWidth',2)
[maxcorr_runpup,l] = max(mean(r_runpup));
plot([xlims(l) xlims(l)],[-0.05 maxcorr_runpup],'--','Color','#547DB1','LineWidth',2)

subplot(3,1,2)
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runHR)+std(r_runHR)/sqrt(trials_num) fliplr(mean(r_runHR)-std(r_runHR)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',0.5);
plot(xlims,mean(r_runHR),'Color','#547DB1','LineWidth',2)
[maxcorr_runHR,l] = max(mean(r_runHR));
plot([xlims(l) xlims(l)],[-0.02 maxcorr_runHR],'--','Color','#547DB1','LineWidth',2)

subplot(3,1,3)
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_HRpup)+std(r_HRpup)/sqrt(trials_num) fliplr(mean(r_HRpup)-std(r_HRpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#C6D6EA','FaceAlpha',0.5);
plot(xlims,mean(r_HRpup),'Color','#547DB1','LineWidth',2)
[maxcorr_HRpup,l] = max(mean(r_HRpup));
plot([xlims(l) xlims(l)],[-0.05 maxcorr_HRpup],'--','Color','#547DB1','LineWidth',2)

%%

r_runpup = r_RunPup(km==2,:);
r_runHR = r_RunHR(km==2,:);
r_HRpup = r_HRPup(km==2,:);

trials_num = size(r_runpup,1);

subplot(3,1,1);
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runpup)+std(r_runpup)/sqrt(trials_num) fliplr(mean(r_runpup)-std(r_runpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#D1E8E6','FaceAlpha',0.5);
plot(xlims,mean(r_runpup),'Color',[124 190 174]./255,'LineWidth',2)
[maxcorr_runpup,l] = max(mean(r_runpup));
plot([xlims(l) xlims(l)],[-0.05 maxcorr_runpup],'--','Color',[124 190 174]./255,'LineWidth',2)
ylim([-0.05 0.15])
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');
xlim([-5 5])
xticks(-5:1:5)
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

subplot(3,1,2)
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runHR)+std(r_runHR)/sqrt(trials_num) fliplr(mean(r_runHR)-std(r_runHR)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#D1E8E6','FaceAlpha',0.5);
plot(xlims,mean(r_runHR),'Color',[124 190 174]./255,'LineWidth',2)
[maxcorr_runHR,l] = max(mean(r_runHR));
plot([xlims(l) xlims(l)],[-0.02 maxcorr_runHR],'--','Color',[124 190 174]./255,'LineWidth',2)
ylim([-0.02 0.04])
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');
xlim([-5 5])
xticks(-5:1:5)
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off

subplot(3,1,3)
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_HRpup)+std(r_HRpup)/sqrt(trials_num) fliplr(mean(r_HRpup)-std(r_HRpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor','#D1E8E6','FaceAlpha',0.5);
plot(xlims,mean(r_HRpup),'Color',[124 190 174]./255,'LineWidth',2)
[maxcorr_HRpup,l] = max(mean(r_HRpup));
plot([xlims(l) xlims(l)],[-0.05 maxcorr_HRpup],'--','Color',[124 190 174]./255,'LineWidth',2)
ylim([-0.05 0.3])
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');
xlim([-5 5])
xticks(-5:1:5)
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off
