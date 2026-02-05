% clc
% clear
% close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

folder = {'m40\Jun_24_2024';
    'm54\Jul_01_2024';
    'm63\Jun_24_2024';
    'm1027\Jun_27_2024';
    'm53\Jul_01_2024';
    'm37\Jun_20_2024';
    'm59\Jul_01_2024';
    'm38\Jun_20_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

fpass_trials=[7 12;
    9 12.5;
    9.5 12;
    10 12;
    9 13;
    11 13;
    8 13;
    10 12];

outlier_pupil=[0 99.1;
    0 99.7;
    0.3 98.2;
    0.1 98.8;
    0 98.9;
    0.3 98.6;
    0 97.8;
    0 99.7];

pre_bout_duration = 5;   % 5 seconds.

trace_duration = 30;   % 30 seconds.

speed_thresh = 0.02;

Run = [];
HR = [];
Pupil = [];

xcorrwind = 5*FrameRate; % in samples

r_RunPup = [];
r_RunHR = [];
r_HRPup = [];

peak_r_RunPup_values_DTR = zeros(size(folder,1), 1);
peak_lags_RunPup_values_DTR = zeros(size(folder,1), 1);

peak_r_RunHR_values_DTR = zeros(size(folder,1), 1);
peak_lags_RunHR_values_DTR = zeros(size(folder,1), 1);

peak_r_HRPup_values_DTR = zeros(size(folder,1), 1);
peak_lags_HRPup_values_DTR = zeros(size(folder,1), 1);


Peak_r_RunPup_values_DTR = [];
Peak_lags_RunPup_values_DTR = [];

Peak_r_RunHR_values_DTR = [];
Peak_lags_RunHR_values_DTR = [];

Peak_r_HRPup_values_DTR = [];
Peak_lags_HRPup_values_DTR = [];

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
    speed_smooth = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);

    %% Analyse pupil size

    pupil_areas = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

    pupil_smooth = movmean(pupil_areas,[smooth_window*FrameRate 0]);

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

    heartRate_bpm_interp = interp1(pksLocs(2:end)*FrameRate,heartRate_bpm,0:time(end)*FrameRate,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %% Running bouts

    binarized_speed = speed_smooth_resampled > speed_thresh;

    running_onsets = find(diff(binarized_speed)==1);

    running_onsets(running_onsets<pre_bout_duration*FrameRate)=[];
    running_onsets(running_onsets>(length(speed_smooth_resampled)-trace_duration*FrameRate+pre_bout_duration*FrameRate))=[];

    good_running_onsets = [];

    for i = 1:length(running_onsets)

        ii = sum(binarized_speed(running_onsets(i)-pre_bout_duration*FrameRate+1:running_onsets(i)));

        if ii == 0
            good_running_onsets = [good_running_onsets,running_onsets(i)];
        end
    end

    %     run = NaN(length(good_running_onsets),trace_duration*Sample_Rate);

    peak_r_RunPup = zeros(length(good_running_onsets), 1);
    peak_lags_RunPup = zeros(length(good_running_onsets), 1); 

    peak_r_RunHR = zeros(length(good_running_onsets), 1);
    peak_lags_RunHR = zeros(length(good_running_onsets), 1); 

    peak_r_HRPup = zeros(length(good_running_onsets), 1);
    peak_lags_HRPup = zeros(length(good_running_onsets), 1);

    for i = 1:length(good_running_onsets)

        t_start = good_running_onsets(i) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        raw_speed = speed_smooth_resampled(t_start:t_end);
        Run = [Run; raw_speed];

        raw_pupil = pupil_smooth(t_start:t_end)';
        Pupil = [Pupil; raw_pupil];

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR = [HR; raw_HR];

        [r_RunPup(k,:),lags] = xcorr(raw_speed-mean(raw_speed),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');
        [peak_r_RunPup(i), max_idx] = max(r_RunPup(k,:));
        peak_lags_RunPup(i) = lags(max_idx)/FrameRate;
        [Peak_r_RunPup_values_DTR(k), max_idx] = max(r_RunPup(k,:));
        Peak_lags_RunPup_values_DTR(k) = lags(max_idx)/FrameRate;

        [r_RunHR(k,:),lags] = xcorr(raw_speed-mean(raw_speed),raw_HR-mean(raw_HR),xcorrwind,'coeff');
        [peak_r_RunHR(i), max_idx] = max(r_RunHR(k,:));
        peak_lags_RunHR(i) = lags(max_idx)/FrameRate;
        [Peak_r_RunHR_values_DTR(k), max_idx] = max(r_RunHR(k,:));
        Peak_lags_RunHR_values_DTR(k) = lags(max_idx)/FrameRate;

        [r_HRPup(k,:),lags] = xcorr(raw_HR-mean(raw_HR),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');
        [peak_r_HRPup(i), max_idx] = max(r_HRPup(k,:));
        peak_lags_HRPup(i) = lags(max_idx)/FrameRate;
        [Peak_r_HRPup_values_DTR(k), max_idx] = max(r_HRPup(k,:));
        Peak_lags_HRPup_values_DTR(k) = lags(max_idx)/FrameRate;

        k = k+1;

    end

    peak_r_RunPup_values_DTR(I) = mean(peak_r_RunPup);
    peak_lags_RunPup_values_DTR(I) = mean(peak_lags_RunPup);

    peak_r_RunHR_values_DTR(I) = mean(peak_r_RunHR);
    peak_lags_RunHR_values_DTR(I) = mean(peak_lags_RunHR);

    peak_r_HRPup_values_DTR(I) = mean(peak_r_HRPup);
    peak_lags_HRPup_values_DTR(I) = mean(peak_lags_HRPup);

end

%%

[h_RunPup_corr, p_RunPup_corr, ~, stats_RunPup_corr] = ttest(Peak_r_RunPup_values_DTR, 0)
[h_RunPup_lag, p_RunPup_lag, ~, stats_RunPup_lag] = ttest(Peak_lags_RunPup_values_DTR, 0)

[h_RunHR_corr, p_RunHR_corr, ~, stats_RunHR_corr] = ttest(Peak_r_RunHR_values_DTR, 0)
[h_RunHR_lag, p_RunHR_lag, ~, stats_RunHR_lag] = ttest(Peak_lags_RunHR_values_DTR, 0)

[h_HRPup_corr, p_HRPup_corr, ~, stats_HRPup_corr] = ttest(Peak_r_HRPup_values_DTR, 0)
[h_HRPup_lag, p_HRPup_lag, ~, stats_HRPup_lag] = ttest(Peak_lags_HRPup_values_DTR, 0)

%% k-means clustering running bouts

run_input_zscore = zscore(Run);

km = kmeans(run_input_zscore,2);

% km = kmeans(Run(:,1:15000),2);

%%

xlims = (-5*FrameRate+1:25*FrameRate)/FrameRate;

figure;

% run = Run(km==1,:);
% pupil = Pupil(km==1,:);
% trials = HR(km==1,:);
% 
% trials_num = size(run,1);
% 
% subplot(3,1,1);
% hold on
% patch('XData',[xlims fliplr(xlims)],'YData',[mean(run)+std(run)/sqrt(trials_num) fliplr(mean(run)-std(run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.5);
% plot(xlims,mean(run),'Color',[139 92 158]./255,'LineWidth',2)
% 
% subplot(3,1,2)
% hold on
% patch('XData',[xlims fliplr(xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.5);
% plot(xlims,mean(pupil),'Color',[139 92 158]./255,'LineWidth',2)
% 
% subplot(3,1,3)
% hold on
% patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.5);
% plot(xlims,mean(trials),'Color',[139 92 158]./255,'LineWidth',2)

%%

run = Run(km==2,:);
pupil = Pupil(km==2,:);
trials = HR(km==2,:);

trials_num = size(run,1);

subplot(3,1,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(run)+std(run)/sqrt(trials_num) fliplr(mean(run)-std(run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.5);
plot(xlims,mean(run),'Color',[229 114 190]./255,'LineWidth',2)
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
patch('XData',[xlims fliplr(xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.5);
plot(xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2)
line([0,0],[120,520],'Color','k','linestyle','--','LineWidth',2);
ylim([120 520])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Pupil Size','FontSize',20,'FontWeight','bold')
ylabel('pixels','FontSize',15,'FontWeight','bold');
hold off

subplot(3,1,3)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.5);
plot(xlims,mean(trials),'Color',[229 114 190]./255,'LineWidth',2)
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

%%

xlims = (-5*FrameRate:5*FrameRate)/FrameRate;

figure;

% r_runpup = r_RunPup(km==1,:);
% r_runHR = r_RunHR(km==1,:);
% r_HRpup = r_HRPup(km==1,:);
% 
% trials_num = size(r_runpup,1);
% 
% subplot(3,1,1);
% hold on
% xline(0,'--','LineWidth',2);
% patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runpup)+std(r_runpup)/sqrt(trials_num) fliplr(mean(r_runpup)-std(r_runpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.5);
% plot(xlims,mean(r_runpup),'Color',[139 92 158]./255,'LineWidth',2)
% [maxcorr_runpup,l] = max(mean(r_runpup));
% plot([xlims(l) xlims(l)],[-0.05 maxcorr_runpup],'--','Color',[139 92 158]./255,'LineWidth',2)
% 
% subplot(3,1,2)
% hold on
% xline(0,'--','LineWidth',2);
% patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runHR)+std(r_runHR)/sqrt(trials_num) fliplr(mean(r_runHR)-std(r_runHR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.5);
% plot(xlims,mean(r_runHR),'Color',[139 92 158]./255,'LineWidth',2)
% [maxcorr_runHR,l] = max(mean(r_runHR));
% plot([xlims(l) xlims(l)],[-0.02 maxcorr_runHR],'--','Color',[139 92 158]./255,'LineWidth',2)
% 
% subplot(3,1,3)
% hold on
% xline(0,'--','LineWidth',2);
% patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_HRpup)+std(r_HRpup)/sqrt(trials_num) fliplr(mean(r_HRpup)-std(r_HRpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.5);
% plot(xlims,mean(r_HRpup),'Color',[139 92 158]./255,'LineWidth',2)
% [maxcorr_HRpup,l] = max(mean(r_HRpup));
% plot([xlims(l) xlims(l)],[-0.05 maxcorr_HRpup],'--','Color',[139 92 158]./255,'LineWidth',2)

%%

r_runpup = r_RunPup(km==2,:);
r_runHR = r_RunHR(km==2,:);
r_HRpup = r_HRPup(km==2,:);

trials_num = size(r_runpup,1);

subplot(3,1,1);
hold on
xline(0,'--','LineWidth',2);
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runpup)+std(r_runpup)/sqrt(trials_num) fliplr(mean(r_runpup)-std(r_runpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.5);
plot(xlims,mean(r_runpup),'Color',[229 114 190]./255,'LineWidth',2)
[maxcorr_runpup,l] = max(mean(r_runpup));
plot([xlims(l) xlims(l)],[-0.2 maxcorr_runpup],'--','Color',[229 114 190]./255,'LineWidth',2)
ylim([-0.2 0.8])
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
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runHR)+std(r_runHR)/sqrt(trials_num) fliplr(mean(r_runHR)-std(r_runHR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.5);
plot(xlims,mean(r_runHR),'Color',[229 114 190]./255,'LineWidth',2)
[maxcorr_runHR,l] = max(mean(r_runHR));
plot([xlims(l) xlims(l)],[-0.1 maxcorr_runHR],'--','Color',[229 114 190]./255,'LineWidth',2)
ylim([-0.1 0.4])
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
patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_HRpup)+std(r_HRpup)/sqrt(trials_num) fliplr(mean(r_HRpup)-std(r_HRpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.5);
plot(xlims,mean(r_HRpup),'Color',[229 114 190]./255,'LineWidth',2)
[maxcorr_HRpup,l] = max(mean(r_HRpup));
plot([xlims(l) xlims(l)],[-0.1 maxcorr_HRpup],'--','Color',[229 114 190]./255,'LineWidth',2)
ylim([-0.1 0.5])
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');
xlim([-5 5])
xticks(-5:1:5)
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold off
