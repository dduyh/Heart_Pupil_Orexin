% %%
% clear; close all; clc;

%%
freezing =[64.875 98.375 11.375 97.25;       % m37
    22.7917 39.2083 38.7917 24.6667;         % m38
    55.1667 73.875 63.9167 82.625;           % m40
    33.9167 47.9167 56.8333 63.375;          % m41
    74.7083 95.7917 70.3333 93.9583;         % m50
    55.9583 96.4167 38.4583 96.5417;         % m2147
    25.875 95.875 79.375 92;                 % m2148
    54.5833 74.5833 76.9583 69.9167;         % m2160
    91.4583 95.5 40.7083 90.75;              % m2198
    80.3333 99.8333 82.5417 99.9167;         % m2216
    86.9167 99.7917 90.75 100;               % m2220
    82.3333 99.875 68.2917 99.5833];         % m2225

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
title({'Freezing on Extinction day 1 (DTR+ fed)',''},'FontSize',24,'FontWeight','bold')

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
errorbar(1:4, M, SEM, "Color",[232 44 43]/255,'LineWidth',3);
plot2 = plot(1:4, M,'o-','color',[232 44 43]/255,'linewidth',4,'markeredgecolor',[232 44 43]/255,'markerfacecolor',[232 44 43]/255,'markersize',4);

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
