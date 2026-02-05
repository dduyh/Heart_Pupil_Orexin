% clc
% clear
% close all

%% set the path for data.

Directory = 'P:\Yihui\data\';    % Main directory\

pupil_folder = {'m1821\Jun_13_2024';
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

ctrl_folder = {'m1840\Jun_13_2024';
    'm1841\Jun_13_2024';
    'm2155\Jun_20_2025';
    'm2156\Jun_20_2025'};

stim = 'Fear_Conditioning';

FrameRate = 20;

stimFreqs = 20;

pre_bout_duration = 100;

trace_duration = 500;   % 90 seconds.

trials_num = size(pupil_folder,1);

ctrl_trials_num = size(ctrl_folder,1);

pupil = NaN(trials_num,trace_duration*FrameRate);

pupil_ctrl = NaN(ctrl_trials_num,trace_duration*FrameRate);

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

outlier_pupil_ctrl=[0 100;
    0 100;
    0 100;
    0 100];

smooth_window = 5;

for I=1:size(pupil_folder,1)

    Data_Folder = [Directory pupil_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    if I>4
        pupil_data = csvread([Data_Folder pupil_folder{I}(1:5) '_' pupil_folder{I}(7:end) 'DLC_resnet50_Pupil_trackingJun20shuffle1_1000000_filtered.csv'],3,1);
    else
        pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);
    end

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs);

    if I>4
        stim_num = length(idxFreqs);
    else
        stim_num = 2;
    end

    pupil_trials = NaN(stim_num,trace_duration*FrameRate);

    %% Analyse pupil size

    Pupil_up_x = pupil_data(:,1);
    Pupil_up_y = pupil_data(:,2);
    Pupil_left_x = pupil_data(:,4);
    Pupil_left_y = pupil_data(:,5);
    Pupil_down_x = pupil_data(:,7);
    Pupil_down_y = pupil_data(:,8);
    Pupil_right_x = pupil_data(:,10);
    Pupil_right_y = pupil_data(:,11);

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

    for II = 1:stim_num

        pupil_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        pupil_end = pupil_start + trace_duration*FrameRate -1;

        raw_pupil = pupil_areas(pupil_start:pupil_end);

        pupil_smooth = movmean(raw_pupil,[smooth_window*FrameRate 0]);

        pupil_zscored = (pupil_smooth - mean(pupil_smooth(1:pre_bout_duration*FrameRate))) / std(pupil_smooth(1:pre_bout_duration*FrameRate));

        pupil_trials(II,:) = pupil_smooth;

    end

    pupil(I,:) = mean(pupil_trials);

end

for I=1:size(ctrl_folder,1)

    Data_Folder = [Directory ctrl_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);

    if I>2
        pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingJun20shuffle1_1000000_filtered.csv'],3,1);
    else
        pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);
    end

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs);

    if I>2
        stim_num = length(idxFreqs);
    else
        stim_num = 2;
    end

    pupil_trials = NaN(stim_num,trace_duration*FrameRate);

    %% Analyse pupil size

    Pupil_up_x = pupil_data(:,1);
    Pupil_up_y = pupil_data(:,2);
    Pupil_left_x = pupil_data(:,4);
    Pupil_left_y = pupil_data(:,5);
    Pupil_down_x = pupil_data(:,7);
    Pupil_down_y = pupil_data(:,8);
    Pupil_right_x = pupil_data(:,10);
    Pupil_right_y = pupil_data(:,11);

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

    pupil_areas = filloutliers(areas,"nearest","percentiles",outlier_pupil_ctrl(I,:));

    for II = 1:stim_num

        pupil_start = round(laser_onsets(II)*FrameRate) - pre_bout_duration*FrameRate +1;
        pupil_end = pupil_start + trace_duration*FrameRate -1;

        raw_pupil = pupil_areas(pupil_start:pupil_end);

        pupil_smooth = movmean(raw_pupil,[smooth_window*FrameRate 0]);

        pupil_zscored = (pupil_smooth - mean(pupil_smooth(1:pre_bout_duration*FrameRate))) / std(pupil_smooth(1:pre_bout_duration*FrameRate));

        pupil_trials(II,:) = pupil_smooth;

    end

    pupil_ctrl(I,:) = mean(pupil_trials);

end
%%

pupil_xlims = (-pre_bout_duration*FrameRate+1:(trace_duration-pre_bout_duration)*FrameRate)/FrameRate;

figure(1);

subplot(4,2,3)
hold on
patch('XData',[0, 0, 60, 60],'YData',[600, 900, 900, 600],'EdgeColor','none','FaceColor',[255, 204, 204]./255,'FaceAlpha',1);
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',1);
plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2);
xlim([-pre_bout_duration (trace_duration-pre_bout_duration)])
ylim([600 850])
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('Pupil Size','FontSize',15,'FontWeight','bold','color','k')
ylabel('pixels','FontSize',12,'FontWeight','bold');
hold off

%%

peak_pupil = mean(pupil(:,100*FrameRate+1:160*FrameRate),2) - mean(pupil(:,40*FrameRate+1:100*FrameRate),2);

peak_pupil_ctrl = mean(pupil_ctrl(:,100*FrameRate+1:160*FrameRate),2) - mean(pupil_ctrl(:,40*FrameRate+1:100*FrameRate),2);

figure(1);

axes('Position', [0.5703 0.5482 0.15 0.1561]);
hold on

for k = 1:size(peak_pupil,1)
    plot(1.3,peak_pupil(k),'marker','o','markersize',4,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

for k = 1:size(peak_pupil_ctrl,1)
    plot(1.7,peak_pupil_ctrl(k),'marker','o','markersize',4,...
        'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],...
        'linestyle','none');
end

line([0.9 2.1],[0 0],'Color','k','linestyle','--','LineWidth',2);

[S_PB,M_PB] = std(peak_pupil,'omitnan');
SEM_PB = S_PB/sqrt(size(peak_pupil,1));

[S_chow,M_chow] = std(peak_pupil_ctrl,'omitnan');
SEM_chow = S_chow/sqrt(size(peak_pupil_ctrl,1));

errorbar(1.1, M_PB, SEM_PB, "Color",[229 114 190]./255,'LineWidth',2);
errorbar(1.9, M_chow, SEM_chow, "Color",[0.5 0.5 0.5],'LineWidth',2);
plot(1.1, M_PB,'marker','o','color',[229 114 190]./255,'linewidth',3,'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,'markersize',5);
plot(1.9, M_chow,'marker','o','color',[0.5 0.5 0.5],'linewidth',3,'markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5);

line([1.3 1.7], [130, 130], 'Color', 'k', 'LineWidth', 2);
text(1.5, 133, '*', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

text(1.1, 70, '***', 'FontSize', 20, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.XTickLabel = {};

xticks([1.2 1.8])

text(1.2, -27, 'Opsin', 'Color', [229 114 190]./255, 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');
text(1.8, -27, 'Control', 'Color', [0.5 0.5 0.5], 'FontSize', 12, 'FontWeight','bold', 'HorizontalAlignment', 'center');

ylim([-20 140])
xlim([0.9 2.1])

title({'Laser ON'},'FontSize',15,'FontWeight','bold')
ylabel('Δ Pupil (pixels)','FontSize',12,'FontWeight','bold');

hold off

[h_Pupil, p_Pupil, ~, stats_Pupil] = ttest2(peak_pupil,peak_pupil_ctrl,'Tail','right')
[h_pupil, p_pupil, ~, stats_pupil] = ttest(peak_pupil,0,'Tail','right')
[h_pupil_ctrl, p_pupil_ctrl, ~, stats_pupil_ctrl] = ttest(peak_pupil_ctrl,0,'Tail','right')
