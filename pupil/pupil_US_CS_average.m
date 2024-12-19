clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

ECG_folder = {'m485_L_VTA_R_NAc_oxLight\Apr_05_2024';
    'm487_L_VTA_nLightR\Apr_05_2024';
    'm1772\Apr_05_2024';
    'm485_L_VTA_R_NAc_oxLight\Apr_08_2024';
    'm487_L_VTA_nLightR\Apr_08_2024';
    'm485_L_VTA_R_NAc_oxLight\Apr_10_2024';
    'm486_L_NAc_nLightR\Apr_10_2024';
    'm1772\Apr_10_2024'};


FrameRate = 20;

trace_duration = 600;   % 600 seconds.

pupil = zeros(size(ECG_folder,1),5);

%%
for I=1:size(ECG_folder,1)

    Data_Folder = [Directory ECG_folder{I} '\'];

    % Load data

    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');

    % pupil size

    pupil_start = round((step_timepoint(4)-step_timepoint(1))*FrameRate)-360*FrameRate+1;
    pupil_end = pupil_start+trace_duration*FrameRate-1;

    raw_pupil = areas(pupil_start:pupil_end);

    for i = 1 : trace_duration
        medianArea(i) = median(raw_pupil(((i-1)*FrameRate+1):i*FrameRate));
    end

    pupil_zscored = (medianArea - mean(medianArea)) / std(medianArea);

    pupil(I,:) = [mean(pupil_zscored(91:120)) mean(pupil_zscored(121:150)) mean(pupil_zscored(241:270)) mean(pupil_zscored(361:390)) mean(pupil_zscored(481:510))];

end

%%

figure
hold on

for k = 1:size(pupil,1)
    plot(1:5,pupil(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(pupil,'omitnan');
SEM = S/sqrt(8);
errorbar(1:5, M, SEM, "Color","black",'LineWidth',3);
plot(1:5, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:5)
xticklabels({'baseline','CS -','baseline','CS +','baseline',})

% ylim([0 100])
xlim([0.5 5.5])

ylabel('z score (s.d.)','FontSize',20,'FontWeight','bold');
title({'Pupil Size',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off



