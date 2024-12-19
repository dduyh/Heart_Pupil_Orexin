clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

ECG_folder = {'m1746\Mar_07_2024\session_1';
    'm1747\Mar_05_2024\session_1';
    'm1748\Mar_07_2024\session_1';
    'm1749\Mar_05_2024\session_1';
    'm1772\Mar_05_2024\session_1'};

pre_atropine_pupil = [];
post_atropine_pupil = [];

%%
for I=1:size(ECG_folder,1)

    Data_Folder = [Directory ECG_folder{I} '\'];

    % Load data

    pupil_data = csvread([Data_Folder 'heartDLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

    % Analyse pupil size

    FrameRate = 20;

    trace_duration = 360;   % 360 seconds.

    Pupil_up_x = pupil_data(:,1);
    Pupil_up_y = pupil_data(:,2);
    Pupil_left_x = pupil_data(:,4);
    Pupil_left_y = pupil_data(:,5);
    Pupil_down_x = pupil_data(:,7);
    Pupil_down_y = pupil_data(:,8);
    Pupil_right_x = pupil_data(:,10);
    Pupil_right_y = pupil_data(:,11);

    center_x = zeros(size(pupil_data,1),1);
    center_y = zeros(size(pupil_data,1),1);
    radii = zeros(size(pupil_data,1),1);
    areas = zeros(size(pupil_data,1),1);

    for i = 1:size(pupil_data,1)
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

    pupil_zscored = (areas - mean(areas)) / std(areas);

    pre_pupil = pupil_zscored(1:trace_duration*FrameRate,1);
    pre_atropine_pupil = [pre_atropine_pupil; mean(pre_pupil)];

    post_pupil = pupil_zscored((900*FrameRate+1):1260*FrameRate,1);
    post_atropine_pupil = [post_atropine_pupil; mean(post_pupil)];
end

%%

pupil_data = [pre_atropine_pupil; post_atropine_pupil];

state = [ones(5,1); repmat(2,5,1)];
state = categorical(state,[1 2],{'Pre','Post'});

figure;
hold on

b1 = boxchart(state,pupil_data,'GroupByColor',state,'BoxFaceAlpha',1,'BoxWidth',2,'LineWidth',5,'BoxMedianLineColor','black','BoxEdgeColor','none','MarkerStyle','none');
b1(1).BoxFaceColor = '#D5E4A8';
b1(2).BoxFaceColor = '#F5D9E6';

HR = [pre_atropine_pupil post_atropine_pupil];

for k = 1:size(HR,1)
    plot([0.8,2.2],HR(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

% ylim([9.5 13])
ylabel('z score (s.d.)','FontSize',20,'FontWeight','bold');
title({'Pupil Size Before/After Atropine Injection','',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h, p, ci, stats] = ttest(pre_atropine_pupil, post_atropine_pupil)





