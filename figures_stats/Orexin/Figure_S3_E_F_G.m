clc
clear
close all

%% set the path for data.

Directory = 'P:\Yihui\data\';    % Main directory\

folder = {'m1821\Jun_13_2024';
    'm1822\Jun_14_2024';
    'm1825\Jun_15_2024';
    'm1826\Jun_15_2024';
    'm2178\Jun_19_2025';
    'm2179\Jun_19_2025';
    'm2144\Jun_19_2025';
    'm2145\Jun_19_2025';
    'm2146\Jun_19_2025';
    'm2140\Jun_19_2025';
    'm2148\Jun_19_2025';
    'm2169\Jun_19_2025'};

stim = 'Fear_Conditioning';

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

stimFreqs = 20;

pre_bout_duration = 100;

trace_duration = 500;

HR = [];
pupil = [];
breath = [];

xcorrwind = 60*FrameRate; % in samples

r_HRpup = [];
r_pupBR = [];
r_HRBR = [];

MinPeakProminences=[0.3 0.05 0.2 0.06 0.3 0.25 0.3 0.3 0.2 0.2 0.2 0.08];

fpass1_trials=[4 5.5;
    4 6;
    4 6;
    3 6;
    3 6;
    3.5 5.5;
    1.5 6.5;
    4 6;
    2 6;
    4 6;
    1 5;
    2 7];

fpass2_trials=[4 5.5;
    4 6;
    4 6;
    3 6;
    3 6;
    3.5 5.5;
    1.5 6.5;
    4 9;
    3 10;
    4 9;
    3 8;
    4 9];

splitIndex=[1000000 1000000 1000000 1000000 1000000 1000000 1000000 1170000 1190000 750000 1180000 1600000];

outlier_pupil=[0.5 99.5;
    0.5 99.5;
    0.2 98.3;
    0 100;
    0 100;
    0 100;
    0 100;
    0.5 97.8;
    0 100;
    0.1 99.36;
    0 100;
    0 100];

smooth_window = 1;

%%

k = 1;

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    if I>4
        pupil_data = csvread([Data_Folder folder{I}(1:5) '_' folder{I}(7:end) 'DLC_resnet50_Pupil_trackingJun20shuffle1_1000000_filtered.csv'],3,1);
    else
        pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);
    end

    breath_data = csvread([Data_Folder 'Values.csv'],1,1);

    data_start = round(step_timepoint(1)*FrameRate)+1;
    data_end = round(step_timepoint(1)*FrameRate + size(datas,1)/Sample_Rate*FrameRate);

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs)-step_timepoint(1);

    %% Analyse pupil size

    Pupil_up_x = pupil_data(data_start:data_end,1);
    Pupil_up_y = pupil_data(data_start:data_end,2);
    Pupil_left_x = pupil_data(data_start:data_end,4);
    Pupil_left_y = pupil_data(data_start:data_end,5);
    Pupil_down_x = pupil_data(data_start:data_end,7);
    Pupil_down_y = pupil_data(data_start:data_end,8);
    Pupil_right_x = pupil_data(data_start:data_end,10);
    Pupil_right_y = pupil_data(data_start:data_end,11);

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

    pupil_smooth = movmean(pupil_areas,[smooth_window*FrameRate 0]);

    %% Analyse breathing rate

    breath_brightness = breath_data(data_start:data_end);

    pixDiff = diff(breath_brightness);

    MinPeakProminence=MinPeakProminences(I);

    [pksVal, pksLocs]=findpeaks([pixDiff(1); pixDiff],FrameRate,'MaxPeakWidth',0.2,'MinPeakProminence',MinPeakProminence,'Threshold',1e-4,'WidthReference','halfheight');

    BB_intervals = diff(pksLocs);
    breathRate=1./BB_intervals;

    breath_outlier = filloutliers(breathRate,"nearest","percentiles",[0 100]);

    breathRate_interp = interp1(pksLocs(2:end)*FrameRate,breath_outlier,0:(size(datas,1)-1)*FrameRate/Sample_Rate,'spline','extrap');

    breath_smooth = movmean(breathRate_interp,[smooth_window*FrameRate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(:,2)';

    % Remove baseline wandering

    fpass1=fpass1_trials(I,:);
    fpass2=fpass2_trials(I,:);

    splitIdx = splitIndex(I);

    bufferLen = round(300 * Sample_Rate);

    segment1 = ECG_raw(1:splitIdx + bufferLen);
    segment2 = ECG_raw(splitIdx - bufferLen + 1:end);

    ECG_seg1 = bandpass(segment1, fpass1, Sample_Rate);
    ECG_seg2 = bandpass(segment2, fpass2, Sample_Rate);

    valid_seg1 = ECG_seg1(1:splitIdx);
    valid_seg2 = ECG_seg2(bufferLen+1:end);

    ECG_Bandpass = [valid_seg1, valid_seg2];

    % find peaks

    if I>4
        minPeakPromVal=0.07;
    else
        minPeakPromVal=0.007;
    end

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.3,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",[0 100]);

    heartRate_bpm_interp = interp1(pksLocs(2:end)*FrameRate,heartRate_bpm_outlier,0:(size(datas,1)-1)*FrameRate/Sample_Rate,'spline','extrap');

    heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[smooth_window*FrameRate 0]);

    %%

    if I>4
        stim_num = length(idxFreqs);
    else
        stim_num = 2;
    end

    for II = 1:stim_num

        t_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        raw_breath = breath_smooth(t_start:t_end);
        breath(k,:) = raw_breath;

        raw_pupil = pupil_smooth(t_start:t_end);
        pupil(k,:) = raw_pupil;

        raw_HR = heartRate_bpm_interp_smooth(t_start:t_end);
        HR(k,:) = raw_HR;

        [r_HRpup(k,:),lags] = xcorr(raw_HR-mean(raw_HR),raw_pupil-mean(raw_pupil),xcorrwind,'coeff');

        [r_pupBR(k,:),lags] = xcorr(raw_pupil-mean(raw_pupil),raw_breath-mean(raw_breath),xcorrwind,'coeff');

        [r_HRBR(k,:),lags] = xcorr(raw_HR-mean(raw_HR),raw_breath-mean(raw_breath),xcorrwind,'coeff');

        k = k+1;

    end

end

%%

xlims = (-xcorrwind:xcorrwind)/FrameRate;

trials_num = size(r_HRpup,1);

figure;

hold on

xline(0,'--','LineWidth',2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_HRpup)+std(r_HRpup)/sqrt(trials_num) fliplr(mean(r_HRpup)-std(r_HRpup)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(xlims,mean(r_HRpup),'Color',[229 114 190]./255,'LineWidth',2)

[maxcorr_HRpup,l] = max(mean(r_HRpup));
plot([xlims(l) xlims(l)],[-0.4 maxcorr_HRpup],'--','Color',[229 114 190]./255,'LineWidth',2)

xlabel('Lag (s)','FontSize',15,'FontWeight','bold')
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');

% xlim([-5 5])
ylim([0 0.3])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


figure;

hold on

xline(0,'--','LineWidth',2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_pupBR)+std(r_pupBR)/sqrt(trials_num) fliplr(mean(r_pupBR)-std(r_pupBR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',0.7);
plot(xlims,mean(r_pupBR),'Color',[64 145 48]./255,'LineWidth',2)

[maxcorr_pupBR,l] = max(mean(r_pupBR));
plot([xlims(l) xlims(l)],[-0.4 maxcorr_pupBR],'--','Color',[64 145 48]./255,'LineWidth',2)

xlabel('Lag (s)','FontSize',15,'FontWeight','bold')
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');

% xlim([-5 5])
ylim([-0.2 0.4])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


figure;

hold on

xline(0,'--','LineWidth',2);

patch('XData',[xlims fliplr(xlims)],'YData',[mean(r_HRBR)+std(r_HRBR)/sqrt(trials_num) fliplr(mean(r_HRBR)-std(r_HRBR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(r_HRBR),'Color',[255 128 128]./255,'LineWidth',2)

[maxcorr_HRBR,l] = max(mean(r_HRBR));
plot([xlims(l) xlims(l)],[-0.4 maxcorr_HRBR],'--','Color',[255 128 128]./255,'LineWidth',2)

xlabel('Lag (s)','FontSize',15,'FontWeight','bold')
ylabel('Correlation coefficient','FontSize',15,'FontWeight','bold');

% xlim([-5 5])
ylim([-0.1 0.3])

ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


