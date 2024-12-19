%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

folder = {'m2070\Oct_16_2024';
    'm2070\Nov_09_2024';
    'm2070\Nov_15_2024';
    'm2071\Oct_16_2024';
    'm2071\Nov_09_2024';
    'm2071\Nov_15_2024';
    'm2072\Oct_16_2024';
    'm2072\Nov_09_2024';
    'm2072\Nov_15_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

fpass_trials=[9 13;
    9.5 13;
    10 13.5;
    11 13;
    9 13;
    10.5 13.5;
    13 14;
    9.5 14.5;
    10.5 13.5];

outlier_pupil=[0 99.95;
    0 99.95;
    0 100;
    0 99.91;
    0 99.99;
    0 100;
    0 100;
    0 99.99;
    0 98];

pre_bout_duration = 5;   % 5 seconds.

trace_duration = 30;   % 30 seconds.

speed_thresh = 0.02;

smooth_window = 1;

Run = [];
HR = [];
Pupil = [];
Cell = [];

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    % Load data
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);

    pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    Cell_Status_data = readtable([Data_Folder folder{I}(1:5) '_cell_traces.csv'],'Format','auto');
    Cell_Status = Cell_Status_data(1,2:end);
    Cell_Status_transposed = rows2vars(Cell_Status);
    Cell_Status_column = Cell_Status_transposed{:, 2};
    cell_status = strcmp(Cell_Status_column, 'accepted');

    cell_traces_data = readmatrix([Data_Folder folder{I}(1:5) '_cell_traces.csv'],'Range',[3 1]);
    cell_times = cell_traces_data(:,1);
    cell_traces = cell_traces_data(:,2:end);
    cell_traces_accepted = cell_traces(:,cell_status);

    nanRows = any(isnan(cell_traces_accepted), 2);
    cell_times = cell_times(~nanRows, :);
    cell_traces_accepted = cell_traces_accepted(~nanRows, :);

    data_onset = find(datas(:,3),1)-50;
    data_offset = find(datas(:,3),1,'last')-10;

    vid_onset = floor((data_onset/Sample_Rate-step_timepoint(1))*FrameRate)+1;
    vid_offset = ceil((data_offset/Sample_Rate-step_timepoint(1))*FrameRate);

    timepoint = times(data_onset:data_offset,1)';
    time = timepoint(1,:)-timepoint(1,1);

    frame_time = 0:(1/FrameRate):time(end);

    %% Analyse cell traces

    numCells = size(cell_traces_accepted, 2);
    cell_traces_interpolated = zeros(numCells,length(frame_time));

    for i = 1:numCells
        cell_traces_interpolated(i,:) = interp1(cell_times,cell_traces_accepted(:, i),frame_time,'nearest','extrap');
    end

    cell_traces_interpolated_smooth = movmean(cell_traces_interpolated,[smooth_window*FrameRate 0],2);

    %% Analyse running signals

    running = datas(data_onset:data_offset,1)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed_smooth = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);
    speed_smooth = [speed_smooth(1) speed_smooth];
    speed_smooth_resampled = resample(speed_smooth, FrameRate, Sample_Rate);
    speed = [Abs_speedDeg(1) Abs_speedDeg];
    speed_resampled = resample(speed, FrameRate, Sample_Rate);

    %% Analyse pupil size

    Pupil_up_x = pupil_data(vid_onset:min(end,vid_offset),1);
    Pupil_up_y = pupil_data(vid_onset:min(end,vid_offset),2);
    Pupil_left_x = pupil_data(vid_onset:min(end,vid_offset),4);
    Pupil_left_y = pupil_data(vid_onset:min(end,vid_offset),5);
    Pupil_down_x = pupil_data(vid_onset:min(end,vid_offset),7);
    Pupil_down_y = pupil_data(vid_onset:min(end,vid_offset),8);
    Pupil_right_x = pupil_data(vid_onset:min(end,vid_offset),10);
    Pupil_right_y = pupil_data(vid_onset:min(end,vid_offset),11);

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

        [center_x(i,1), center_y(i,1)]= node(X1,Y1,X2,Y2);

        MinorAxis = [X1(1), X1(2); Y1(1), Y1(2)];
        MinorAxisLength = pdist(MinorAxis, 'euclidean');

        MajorAxis = [X2(1), X2(2); Y2(1), Y2(2)];
        MajorAxisLength = pdist(MajorAxis, 'euclidean');

        diameters = mean([MajorAxisLength, MinorAxisLength]);
        radii(i,1) = diameters/2;
        areas(i,1) = pi*radii(i,1).^2;
    end

    pupil = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

    pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

    %% Analyse ECG signals

    % Raw ECG signals

    ECG_raw = datas(data_onset:data_offset,2)';

    % Remove baseline wandering

    fpass=fpass_trials(I,:);

    ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

    % find peaks

    minPeakPromVal=0.07;

    [pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

    % heart rate in time

    RR_intervals = diff(pksLocs);
    heartRate=1./RR_intervals;
    heartRate_bpm=heartRate*60;

    heartRate_bpm_interp = interp1(pksLocs(2:end),heartRate_bpm,frame_time,'spline','extrap');

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

    for i = 1:length(good_running_onsets)

        t_start = good_running_onsets(i) - pre_bout_duration*FrameRate +1;
        t_end = t_start + trace_duration*FrameRate -1;

        Run = [Run; speed_resampled(t_start:t_end)];

        Pupil = [Pupil; pupil(t_start:t_end)'];

        HR = [HR; heartRate_bpm_interp(t_start:t_end)];

        cell_traces_bouts = cell_traces_interpolated(:,t_start:t_end);
        cell_traces_bouts_zscored = (cell_traces_bouts - mean(cell_traces_bouts(:,1:5*FrameRate),2))./std(cell_traces_bouts(:,1:5*FrameRate),0,2);
        Cell = [Cell; cell_traces_bouts_zscored];

    end

end

%%

xlims = (-5*FrameRate+1:25*FrameRate)/FrameRate;

trials_num = size(Run,1);
traces_num = size(Cell,1);

figure;

subplot(4,1,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Run)+std(Run)/sqrt(trials_num) fliplr(mean(Run)-std(Run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',0.7);
plot(xlims,mean(Run),'Color','k','LineWidth',2)
line([0,0],[0,0.05],'Color','k','linestyle','--','LineWidth',2);
ylim([0 0.05])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Abs. Running Speed','FontSize',20,'FontWeight','bold')
hold off

subplot(4,1,2);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Cell)+std(Cell)/sqrt(traces_num) fliplr(mean(Cell)-std(Cell)/sqrt(traces_num))],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',0.7);
plot(xlims,mean(Cell),'Color',[64 145 48]./255,'LineWidth',2)
line([0,0],[-0.3,2],'Color','k','linestyle','--','LineWidth',2);
ylim([-0.3,2])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Orx-GCaMP6s','FontSize',20,'FontWeight','bold')
hold off

subplot(4,1,3)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(Pupil)+std(Pupil)/sqrt(trials_num) fliplr(mean(Pupil)-std(Pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(xlims,mean(Pupil),'Color',[229 114 190]./255,'LineWidth',2)
line([0,0],[100,180],'Color','k','linestyle','--','LineWidth',2);
ylim([100 180])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Pupil Size','FontSize',20,'FontWeight','bold')
ylabel('z-score','FontSize',15,'FontWeight','bold');
hold off

subplot(4,1,4)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(HR)+std(HR)/sqrt(trials_num) fliplr(mean(HR)-std(HR)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(HR),'Color',[255 128 128]./255,'LineWidth',2)
line([0,0],[700,780],'Color','k','linestyle','--','LineWidth',2);
ylim([700 780])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
title('Heart Rate','FontSize',20,'FontWeight','bold')
ylabel('HR (bpm)','FontSize',15,'FontWeight','bold');
hold off

%%

figure;

hold on

imagesc(Cell,[0 20]);
colormap('hot');

hbar = colorbar;
set(hbar, 'Position', [0.9230    0.7825    0.0377    0.1433])
hbar.FontSize = 15;
hbar.Label.String = 'z-score';
hbar.FontWeight = 'bold';

line([5*FrameRate,5*FrameRate],[1,traces_num],'Color','white','linestyle','--','LineWidth',2);
ylim([1,traces_num])

xlim([1,30*FrameRate])
xticks([1:5*FrameRate:30*FrameRate 30*FrameRate])
xticklabels({'-5','0','5','10','15','20','25'})

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

title('Orx-GCaMP6s Neural Activities','FontSize',20,'FontWeight','bold')
ylabel('cell traces','FontSize',15,'FontWeight','bold');
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')

hold off
