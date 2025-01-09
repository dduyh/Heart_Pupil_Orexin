% clc
% clear
% close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

pupil_folder = {'m1027\Jul_07_2024';
    'm37\Jul_08_2024';
    'm38\Jul_08_2024';
    'm40\Jul_08_2024';
    'm41\Jul_08_2024';
    'm53\Jul_08_2024';
    'm54\Jul_08_2024';
    'm63\Jul_09_2024';
    'm50\Jul_09_2024';
    'm59\Jul_09_2024'};

outlier_pupil=[0 99.9;
    0 100;
    0 99.9;
    0 99.9;
    0 100;
    0 99.9;
    0 99.9;
    0 99.9;
    0 100;
    0 99.995];

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

trace_duration = 62;   % 90 seconds.

trials_num = size(pupil_folder,1);

pupil = NaN(trials_num,trace_duration*FrameRate);

k = 1;

for I=1:size(pupil_folder,1)

    Data_Folder = [Directory pupil_folder{I} '\'];

    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    tone_onsets = step_timepoint([5])-step_timepoint(1);

    for II = 1:length(tone_onsets)

        pupil_start = round(tone_onsets(II)*FrameRate) - 30*FrameRate +1;
        pupil_end = round(tone_onsets(II)*FrameRate) + 32*FrameRate;

        pupil_areas = filloutliers(areas,"nearest","percentiles",outlier_pupil(I,:));

        pupil_zscored = (pupil_areas - mean(pupil_areas)) / std(pupil_areas);

        pupil(k,:) = pupil_zscored(pupil_start:pupil_end);

        k = k+1;

    end

end

%%

pupil_xlims = (-30*FrameRate+1:32*FrameRate)/FrameRate;

figure(1);

hold on
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[215 195 223]./255,'FaceAlpha',0.5);
plot(pupil_xlims,mean(pupil),'Color',[139 92 158]./255,'LineWidth',2)
line([0,0],[-1 3],'Color','k','linestyle','--','LineWidth',2);

% axis off
% ylabel('Pupil Area z-score','FontSize',15,'FontWeight','bold');
% title({'Pupil Size',''},'FontSize',20,'FontWeight','bold','color','k')

xlim([-30 32])
ylim([-1 3])

% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off


