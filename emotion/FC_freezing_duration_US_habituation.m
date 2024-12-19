%%
clear; close all; clc;

%% 
freezing =[72.1667 98.8333 89.5 38.1667 64 57.6667 49.8333 29.5 40.5;
45.5 NaN NaN 0 71.3333 66.5 0 85.3333 0;
50.8333 90.5 24.8333 0 22.8333 0 18.6667 0 63;
7.16667 0 30.3333 61.8333 44.8333 31 0 23.6667 20.6667;
0 0 0 0 0 0 0 0 0;
0 0 9.5 0 22.5 0 47 45.3333 73.5;
7.33333 0 0 0 0 0 0 11 0];

figure
hold on

for k = 1:size(freezing,1)
    plot(1:9,freezing(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(freezing,'omitnan');
SEM = S/sqrt(7);
errorbar(1:9, M, SEM, "Color","black",'LineWidth',3);
plot(1:9, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:9)
xticklabels({'1','2','3','1','2','3','1','2','3'})

% ylim([0 100])
xlim([0.5 9.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');

% title({'Freezing on Habituation days',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

% [h1, p1, ci1, stats1] = ttest(freezing(:,1), freezing(:,2))
% 
% [h2, p2, ci2, stats2] = ttest(freezing(:,2), freezing(:,3))
