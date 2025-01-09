clc
clear
close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

pupil_folder = {'m1026\Jul_12_2024';
    'm1028\Jul_12_2024';
    'm39\Jul_12_2024';
    'm55\Jul_13_2024';
    'm64\Jul_14_2024';
    'm65\Jul_14_2024';
    'm49\Jul_13_2024';
    'm58\Jul_13_2024'};

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

% trace_duration = 62;   % 90 seconds.

% trials_num = 3*size(pupil_folder,1);
% 
% pupil = NaN(trials_num,trace_duration*FrameRate);
% 
% k = 1;

pupil = [];

for I=1:size(pupil_folder,1)

    Data_Folder = [Directory pupil_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    pupil = [pupil; mean(areas)];

%     tone_onsets = step_timepoint([3 4 5])-step_timepoint(1);
% 
%     for II = 1:length(tone_onsets)
% 
%         pupil_start = round(tone_onsets(II)*FrameRate) - 30*FrameRate +1;
%         pupil_end = round(tone_onsets(II)*FrameRate) + 32*FrameRate;
% 
%         raw_pupil = areas(pupil_start:pupil_end);
% 
%         pupil_zscored = (raw_pupil - mean(raw_pupil(1:30*FrameRate))) / std(raw_pupil(1:30*FrameRate));
%         %         pupil_zscored = (raw_pupil - mean(raw_pupil)) / std(raw_pupil);
%         pupil(k,:) = pupil_zscored;
% 
%         k = k+1;
% 
%     end

end

%%

% pupil_xlims = (-30*FrameRate+1:32*FrameRate)/FrameRate;
% 
% figure;
% 
% hold on
% patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
% plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2)
% line([0,0],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
% 
% % axis off
% ylabel('Pupil Area z-score','FontSize',15,'FontWeight','bold');
% title({'Pupil Size',''},'FontSize',20,'FontWeight','bold','color','k')
% 
% xlim([-30 32])
% ylim([-5 15])
% 
% % axis off
% ax = gca;
% % ax.XAxis.Visible = 'off'; % remove x-axis
% ax.LineWidth = 1.2;
% ax.FontSize = 15;
% ax.FontWeight = 'bold';
% 
% hold off


