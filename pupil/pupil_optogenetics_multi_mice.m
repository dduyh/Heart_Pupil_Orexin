clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

% pupil_folder = {'m1825\May_08_2024';
%     'm1826\May_10_2024';
%     'm1821\May_10_2024';
%     'm1825\May_14_2024'};
% 
% ctrl_folder = {'m1840\May_14_2024';
%     'm1841\May_14_2024'};

pupil_folder = {'m1825\May_23_2024';
    'm1826\May_23_2024';
    'm1821\May_24_2024_2';
    'm1822\May_24_2024'};

ctrl_folder = {'m1840\May_23_2024';
    'm1821\May_24_2024'};

FrameRate = 20;

stimFreqs = 20;

trace_duration = 120;   % 90 seconds.

trials_num = size(pupil_folder,1);

ctrl_trials_num = size(ctrl_folder,1);

pupil = NaN(trials_num,trace_duration*FrameRate);

pupil_ctrl = NaN(ctrl_trials_num,trace_duration*FrameRate);

pupil_xlims = (1:trace_duration*FrameRate)/FrameRate;

for I=1:size(pupil_folder,1)

    Data_Folder = [Directory pupil_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs);

    pupil_trials = NaN(length(idxFreqs),trace_duration*FrameRate);

    for II = 1:length(laser_onsets)

        pupil_start = round(laser_onsets(II)*FrameRate) - 30*FrameRate +1;
        pupil_end = round(laser_onsets(II)*FrameRate) + 90*FrameRate;

        raw_pupil = areas(pupil_start:pupil_end);

        pupil_zscored = (raw_pupil - mean(raw_pupil(1:30*FrameRate))) / std(raw_pupil(1:30*FrameRate));

        pupil_trials(II,:) = pupil_zscored;

    end

    figure(1);
    hold on
    plot(pupil_xlims,mean(pupil_trials),'LineWidth',1)

    pupil(I,:) = mean(pupil_trials);

end

for I=1:size(ctrl_folder,1)

    Data_Folder = [Directory ctrl_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs);

    pupil_trials = NaN(length(idxFreqs),trace_duration*FrameRate);

    for II = 1:length(laser_onsets)

        pupil_start = round(laser_onsets(II)*FrameRate) - 30*FrameRate +1;
        pupil_end = round(laser_onsets(II)*FrameRate) + 90*FrameRate;

        raw_pupil = areas(pupil_start:pupil_end);

        pupil_zscored = (raw_pupil - mean(raw_pupil(1:30*FrameRate))) / std(raw_pupil(1:30*FrameRate));

        pupil_trials(II,:) = pupil_zscored;

    end

    figure(2);
    hold on
    plot(pupil_xlims,mean(pupil_trials),'LineWidth',1)

    pupil_ctrl(I,:) = mean(pupil_trials);

end
%%

figure(1);

patch('XData',[30, 30, 90, 90],'YData',[-10, 270, 270, -10],'EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

plot(pupil_xlims,mean(pupil),'Color','k','LineWidth',2)

% axis off
ylabel('Pupil Area z-score','FontSize',15,'FontWeight','bold');
title({'Pupil Size (Chrimson)',''},'FontSize',20,'FontWeight','bold','color','k')

xlim([0 120])
ylim([-10 270])

% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


figure(2);

patch('XData',[30, 30, 90, 90],'YData',[-10, 270, 270, -10],'EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

plot(pupil_xlims,mean(pupil_ctrl),'Color','k','LineWidth',2)

% axis off
ylabel('Pupil Area z-score','FontSize',15,'FontWeight','bold');
title({'Pupil Size (Control)',''},'FontSize',20,'FontWeight','bold','color','k')

xlim([0 120])
ylim([-10 270])

% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off