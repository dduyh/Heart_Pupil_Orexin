clc
clear
close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

pupil_folder = {'m1746\Mar_07_2024\session_1';
    'm1747\Mar_05_2024\session_1';
    'm1748\Mar_07_2024\session_1';
    'm1749\Mar_05_2024\session_1';
    'm1772\Mar_05_2024\session_1'};

stim = 'heart';

Injection_onsets = [430 360 400 500 500]+20;

FrameRate = 20;

trace_duration = 1200;   % 90 seconds.

pre_duration = 360;

trials_num = size(pupil_folder,1);

pupil = NaN(trials_num,trace_duration*FrameRate);

smooth_window = 1;

for I=1:size(pupil_folder,1)

    Data_Folder = [Directory pupil_folder{I} '\'];

    pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    %% Analyse pupil size

    Pupil_up_x = pupil_data(:,1);
    Pupil_up_y = pupil_data(:,2);
    Pupil_left_x = pupil_data(:,4);
    Pupil_left_y = pupil_data(:,5);
    Pupil_down_x = pupil_data(:,7);
    Pupil_down_y = pupil_data(:,8);
    Pupil_right_x = pupil_data(:,10);
    Pupil_right_y = pupil_data(:,11);

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

    pupil_areas = filloutliers(areas,"nearest","percentiles",[0 99.95]);

    pupil_smooth = movmean(pupil_areas,[smooth_window*FrameRate 0]);

    trace_start = Injection_onsets(I)*FrameRate-pre_duration*FrameRate;
    trace_end = trace_start+trace_duration*FrameRate-1;

    pupil(I,:) = pupil_smooth(trace_start:trace_end);

end

%%

pupil_xlims = (-pre_duration*FrameRate+1:(trace_duration-pre_duration)*FrameRate)/FrameRate;

figure;

hold on
patch('XData',[0, 0, 60, 60],'YData',[200, 1500, 1500, 200],'EdgeColor','none','FaceColor','#DF553F','FaceAlpha',0.5);

patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2)

% axis off
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
ylabel('Pupil Size (Pixels)','FontSize',15,'FontWeight','bold');
% title({'Pupil Size (Atropine)',''},'FontSize',20,'FontWeight','bold','color','k')

xlim([-pre_duration trace_duration-pre_duration])
ylim([200 1500])

% axis off
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

