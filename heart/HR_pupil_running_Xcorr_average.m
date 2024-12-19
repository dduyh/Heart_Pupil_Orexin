clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

folder = {'m1821\May_30_2024';
    'm1822\May_30_2024';
    'm1825\May_31_2024';
    'm1826\May_31_2024'};

stim = 'Fear_Conditioning';

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

stimFreqs = 20;

trace_duration = 500;   % 90 seconds.

trials_num = 3*size(folder,1);

run = NaN(trials_num,trace_duration*Sample_Rate);
pupil = NaN(trials_num,trace_duration*Sample_Rate);
HR = NaN(trials_num,trace_duration*Sample_Rate);

xcorrwind = 5*Sample_Rate; % in samples

r_runpup = NaN(trials_num,2*xcorrwind+1);
r_runHR = NaN(trials_num,2*xcorrwind+1);
r_HRpup = NaN(trials_num,2*xcorrwind+1);

fpass_trials=[9 12;
    9 13;
    9.5 12;
    8 12];

outlier_running=[0 99.9;
    0 99.9;
    0 99.56;
    0 99.9];

outlier_pupil=[0 100;
    0 99.9;
    0 99.8;
    0 100];

outlier_HR=[0 100;
    0.1 100;
    0 100;
    0 100];



%%

k = 1;

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    data_start = round(step_timepoint(1)*FrameRate)+1;
    data_end = round(step_timepoint(1)*FrameRate + size(datas,1)/Sample_Rate*FrameRate);

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs)-step_timepoint(1);

    %% Analyse running signals

    running = datas(:,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    Abs_speedDeg_outlier = filloutliers(Abs_speedDeg,"nearest","percentiles",outlier_running(I,:));
    speed = [Abs_speedDeg_outlier(1) Abs_speedDeg_outlier];

    %% Analyse pupil size

    Pupil_up_x = pupil_data(data_start:data_end,1);
    Pupil_up_y = pupil_data(data_start:data_end,2);
    Pupil_left_x = pupil_data(data_start:data_end,4);
    Pupil_left_y = pupil_data(data_start:data_end,5);
    Pupil_down_x = pupil_data(data_start:data_end,7);
    Pupil_down_y = pupil_data(data_start:data_end,8);
    Pupil_right_x = pupil_data(data_start:data_end,10);
    Pupil_right_y = pupil_data(data_start:data_end,11);

    center_x = zeros(size(Pupil_up_x,1),1);
    center_y = zeros(size(Pupil_up_x,1),1);
    radii = zeros(size(Pupil_up_x,1),1);
    areas = zeros(size(Pupil_up_x,1),1);

    for i = 1:size(Pupil_up_x,1)
        X1(1) = Pupil_up_x(i,1);
        X1(2) = Pupil_up_y(i,1);

        X2(1) = Pupil_left_x(i,1);
        X2(2) = Pupil_left_y(i,1);

        Y1(1) = Pupil_down_x(i,1);
        Y1(2) = Pupil_down_y(i,1);

        Y2(1) = Pupil_right_x(i,1);
        Y2(2) = Pupil_right_y(i,1);

        MinorAxis = [X1(1), X1(2); Y1(1), Y1(2)];
        MinorAxisLength = pdist(MinorAxis, 'euclidean');

        MajorAxis = [X2(1), X2(2); Y2(1), Y2(2)];
        MajorAxisLength = pdist(MajorAxis, 'euclidean');

        diameters = mean([MajorAxisLength, MinorAxisLength]);
        radii(i,1) = diameters/2;
        areas(i,1) = pi*radii(i,1).^2;
    end

    pupil_areas = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

    video_time = (1/FrameRate):(1/FrameRate):(size(pupil_areas,1)/FrameRate);

    pupil_interp = interp1(video_time*Sample_Rate,pupil_areas,1:size(datas,1),'nearest','extrap');

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

    heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",outlier_HR(I,:));

    heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm_outlier,1:size(datas,1),'nearest','extrap');

    %%

    for II = 1:length(laser_onsets)

        t_start = round(laser_onsets(II)*Sample_Rate) - 100*Sample_Rate +1;
        t_end = round(laser_onsets(II)*Sample_Rate) + 400*Sample_Rate;

        raw_running = speed(t_start:t_end);
        run(k,:) = raw_running;

        raw_pupil = pupil_interp(t_start:t_end);
        pupil(k,:) = raw_pupil;

        raw_HR = heartRate_bpm_interp(t_start:t_end);
        HR(k,:) = raw_HR;

        [r_runpup(k,:),lags] = xcorr(raw_running-mean(raw_running),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');

        [r_runHR(k,:),lags] = xcorr(raw_running-mean(raw_running),raw_HR-mean(raw_HR),xcorrwind,'coeff');

        [r_HRpup(k,:),lags] = xcorr(raw_HR-mean(raw_HR),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');

        k = k+1;

    end

end

%%

xlims = (-5*Sample_Rate:5*Sample_Rate)/Sample_Rate;

figure;

hold on

xline(0,'--','LineWidth',2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runpup)+std(r_runpup)/sqrt(trials_num) fliplr(mean(r_runpup)-std(r_runpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',0.7);
plot(xlims,mean(r_runpup),'Color','k','LineWidth',2)

[maxcorr_runpup,l] = max(mean(r_runpup));
plot([xlims(l) xlims(l)],[0.04 maxcorr_runpup],'--','Color','k','LineWidth',2)

ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');

xlim([-5 5])
% ylim([540 740])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


figure;

hold on

xline(0,'--','LineWidth',2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_runHR)+std(r_runHR)/sqrt(trials_num) fliplr(mean(r_runHR)-std(r_runHR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(r_runHR),'Color',[255 128 128]./255,'LineWidth',2)

[maxcorr_runHR,l] = max(mean(r_runHR));
plot([xlims(l) xlims(l)],[0.02 maxcorr_runHR],'--','Color',[255 128 128]./255,'LineWidth',2)

ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');

xlim([-5 5])
% ylim([540 740])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


figure;

hold on

xline(0,'--','LineWidth',2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_HRpup)+std(r_HRpup)/sqrt(trials_num) fliplr(mean(r_HRpup)-std(r_HRpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(xlims,mean(r_HRpup),'Color',[229 114 190]./255,'LineWidth',2)

[maxcorr_HRpup,l] = max(mean(r_HRpup));
plot([xlims(l) xlims(l)],[0.15 maxcorr_HRpup],'--','Color',[229 114 190]./255,'LineWidth',2)

ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');

xlim([-5 5])
% ylim([540 740])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off
