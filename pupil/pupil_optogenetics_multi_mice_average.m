clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

pupil_folder = {'m1825\May_08_2024';
    'm1826\May_10_2024';
    'm1821\May_10_2024';
    'm1825\May_14_2024'};

% ctrl_folder = {'m1840\May_14_2024';
%     'm1841\May_14_2024'};

ctrl_folder = {'m1840\May_14_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

stimFreqs = 1;

trace_duration = 90;   % 90 seconds.

trials_num = 5*size(pupil_folder,1);

ctrl_trials_num = 5*size(ctrl_folder,1);

pupil = NaN(trials_num,trace_duration*FrameRate);

pupil_ctrl = NaN(ctrl_trials_num,trace_duration*FrameRate);

k = 1;

for I=1:size(pupil_folder,1)

    Data_Folder = [Directory pupil_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs);

    for II = 1:length(laser_onsets)

        pupil_start = round(laser_onsets(II)*FrameRate) - 30*FrameRate +1;
        pupil_end = round(laser_onsets(II)*FrameRate) + 60*FrameRate;

        raw_pupil = areas(pupil_start:pupil_end);

        pupil_zscored = (raw_pupil - mean(raw_pupil(1:30*FrameRate))) / std(raw_pupil(1:30*FrameRate));
        %         pupil_zscored = (raw_pupil - mean(raw_pupil)) / std(raw_pupil);
        pupil(k,:) = pupil_zscored;

        k = k+1;

    end

end

k = 1;

for I=1:size(ctrl_folder,1)

    Data_Folder = [Directory ctrl_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'delay.mat']);
    load([Data_Folder 'stimSchedule.mat'])
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    idxFreqs = find(stimSchedule == stimFreqs);
    laser_onsets = step_timepoint(idxFreqs)+delay(idxFreqs);

    for II = 1:length(laser_onsets)

        pupil_start = round(laser_onsets(II)*FrameRate) - 30*FrameRate +1;
        pupil_end = round(laser_onsets(II)*FrameRate) + 60*FrameRate;

        raw_pupil = areas(pupil_start:pupil_end);

        pupil_zscored = (raw_pupil - mean(raw_pupil(1:30*FrameRate))) / std(raw_pupil(1:30*FrameRate));
        %         pupil_zscored = (raw_pupil - mean(raw_pupil)) / std(raw_pupil);
        pupil_ctrl(k,:) = pupil_zscored;

        k = k+1;

    end

end

%%

pupil_xlims = (1:trace_duration*FrameRate)/FrameRate;

figure;

hold on
patch('XData',[30, 30, 60, 60],'YData',[-10, 110, 110, -10],'EdgeColor','none','FaceColor','green','FaceAlpha',0.2);

patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil_ctrl)+std(pupil_ctrl)/sqrt(ctrl_trials_num) fliplr(mean(pupil_ctrl)-std(pupil_ctrl)/sqrt(ctrl_trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',0.7);
p1 = plot(pupil_xlims,mean(pupil_ctrl),'Color','k','LineWidth',2);

patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
p2 = plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2);

legend([p2 p1],{'Chrimson','Control'},'Location','northwest')
legend('boxoff')

% axis off
ylabel('Pupil Area z-score','FontSize',15,'FontWeight','bold');
title({'Pupil Size',''},'FontSize',20,'FontWeight','bold','color','k')

xlim([0 90])
ylim([-10 90])

% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


