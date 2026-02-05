% clc
% clear
% close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

pupil_folder = {'m37\Aug_06_2024';
    'm38\Aug_06_2024';
    'm40\Aug_06_2024';
    'm41\Aug_06_2024';
    'm50\Aug_06_2024';
    'm2147\Jan_22_2026';
    'm2148\Jan_22_2026';
    'm2160\Jan_22_2026';
    'm2198\Jan_22_2026';
    'm2216\Jan_23_2026';
    'm2220\Jan_23_2026';
    'm2225\Jan_23_2026'};

FrameRate = 20;

trials_num = size(pupil_folder,1);

pupil = NaN(trials_num,4);

outlier_pupil=[0 99.98;
    0 100;
    0 99.98;
    0 100;
    0 99.91;
    0 100;
    0 100;
    0 99.98;
    0 100;
    0 100;
    0 100;
    0 99.99];

%%

for I=1:size(pupil_folder,1)

    Data_Folder = [Directory pupil_folder{I} '\'];

    load([Data_Folder 'step_timepoint.mat']);

    if I>5
        pupil_data = csvread([Data_Folder pupil_folder{I}(1:5) '_' pupil_folder{I}(7:end) 'DLC_resnet50_Pupil_trackingJan29shuffle1_1000000_filtered.csv'],3,1);
    else
        pupil_data = csvread([Data_Folder 'Fear_ConditioningDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);
    end

    pupil_start = round((step_timepoint(2)-step_timepoint(1))*FrameRate)+1;

    %% Analyse pupil size

    Pupil_up_x = pupil_data(pupil_start:end,1);
    Pupil_up_y = pupil_data(pupil_start:end,2);
    Pupil_left_x = pupil_data(pupil_start:end,4);
    Pupil_left_y = pupil_data(pupil_start:end,5);
    Pupil_down_x = pupil_data(pupil_start:end,7);
    Pupil_down_y = pupil_data(pupil_start:end,8);
    Pupil_right_x = pupil_data(pupil_start:end,10);
    Pupil_right_y = pupil_data(pupil_start:end,11);

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

    pupil_zscored = zscore(pupil_areas);

    if I==1||I==2||I==3||I==4

        pupil(I,:) = [mean(pupil_zscored(120*FrameRate+1:150*FrameRate)) mean(pupil_zscored(360*FrameRate+1:390*FrameRate)) mean(pupil_zscored(600*FrameRate+1:630*FrameRate)) mean(pupil_zscored(840*FrameRate+1:870*FrameRate))];

    elseif I==5

        pupil(I,:) = [mean(pupil_zscored(360*FrameRate+1:390*FrameRate)) mean(pupil_zscored(120*FrameRate+1:150*FrameRate)) mean(pupil_zscored(840*FrameRate+1:870*FrameRate)) mean(pupil_zscored(600*FrameRate+1:630*FrameRate))];

    else

        pupil(I,:) = [mean(pupil_zscored(240*FrameRate+1:270*FrameRate)) mean(pupil_zscored(600*FrameRate+1:630*FrameRate)) mean(pupil_zscored(960*FrameRate+1:990*FrameRate)) mean(pupil_zscored(1320*FrameRate+1:1350*FrameRate))];

    end

end

%%

figure
hold on

for k = 1:size(pupil,1)
    plot(1:4,pupil(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(pupil,'omitnan');
SEM = S/sqrt(size(pupil,1));
errorbar(1:4, M, SEM, "Color","black",'LineWidth',3);
plot(1:4, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:4)
xticklabels({'CS -','CS +','CS -','CS +'})

ylim([-2 4])
xlim([0.5 4.5])

ylabel('z score (s.d.)','FontSize',20,'FontWeight','bold');
title('Pupil Size','FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

%%

figure(2)
hold on

[S,M] = std(pupil,'omitnan');
SEM = S/sqrt(size(pupil,1));
errorbar(1:4, M, SEM, "Color",[232 44 43]/255,'LineWidth',3);
plot2 = plot(1:4, M,'o-','color',[232 44 43]/255,'linewidth',4,'markeredgecolor',[232 44 43]/255,'markerfacecolor',[232 44 43]/255,'markersize',4);

xticks(1:4)
xticklabels({'CS -','CS +','CS -','CS +'})

ylim([-2 4])
xlim([0.5 4.5])

ylabel('z score (s.d.)','FontSize',20,'FontWeight','bold');
title('Pupil Size on Cue Testing Phase','FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off



