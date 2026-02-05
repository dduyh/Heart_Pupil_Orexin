% %%
% clear; close all; clc;

%%
freezing =[32.3333 31.25 41.5417 41.1667;     % m53
    28 95.5 50 98.0833;                       % m54
    72.4583 89.25 82.5 69.7917;               % m59
    51.5833 74.375 62.125 77.0833;            % m63
    26.375 79.7083 58.4167 54.4167;           % m2227
    40.2917 97.2083 87.0833 97.125;           % m2230
    62.125 97.625 36 92.25;                   % m2212
    74.75 83.75 89.125 93.0833;               % m2201
    37.625 87.25 74.75 98.0417;               % m2204
    48.0833 93.9167 54.5417 94.125];          % m2205

figure
hold on

for k = 1:size(freezing,1)
    plot(1:4,freezing(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(freezing,'omitnan');
SEM = S/sqrt(size(freezing,1));
errorbar(1:4, M, SEM, "Color","black",'LineWidth',3);
plot(1:4, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:4)
xticklabels({'CS -','CS +','CS -','CS +'})

ylim([0 100])
xlim([0.5 4.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');
%
title({'Freezing on Extinction day 1 (DTR+ fasted)',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

% [h1, p1, ci1, stats1] = ttest(freezing(:,1), freezing(:,2))
%
% [h2, p2, ci2, stats2] = ttest(freezing(:,2), freezing(:,3))


%%

figure(2)
hold on

[S,M] = std(freezing,'omitnan');
SEM = S/sqrt(size(freezing,1));
errorbar(1:4, M, SEM, "Color",[233 102 94]/255,'LineWidth',3);
plot1 = plot(1:4, M,'o-','color',[233 102 94]/255,'linewidth',4,'markeredgecolor',[233 102 94]/255,'markerfacecolor',[233 102 94]/255,'markersize',4);

xticks(1:4)
xticklabels({'CS -','CS +','CS -','CS +'})

ylim([0 100])
xlim([0.5 4.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');
title('Freezing on Cue Testing Phase', 'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

