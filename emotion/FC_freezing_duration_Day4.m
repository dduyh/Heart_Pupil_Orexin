%%
clear; close all; clc;

%% 
freezing =[4.38889 11.2222 6.83333 17.445 16.3327 10.6667 25.2231 32.7768;
           45.278 45.9445 21.6662 44.3894 54.6105 29.2222 67.9999 20.2222;
           48.0558 32.5001 43.5001 57.5 59.1112 51.1112 21.9436 51.1676;
           17.5555 27.1111 40.6666 10.7222 42.3333 29.6111 32.3333 18.5;
           80.5002 66.4441 26.1666 19.1111 19.5 23.7778 28.1111 57.2222;
           82.3331 49.667 89.6111 55.2223 62.2778 49.0556 63.2223 9.72131;
           79.1113 72.4996 67.2222 72.8339 74.6104 65.9452 67.1658 34.0565;
           87.1113 82.2778 96.4444 85.9994 34.7229 57.1112 26.2769 23.1676];

figure
hold on

for k = 1:size(freezing,1)
    plot(1:8,freezing(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(freezing);
SEM = S/sqrt(8);
errorbar(1:8, M, SEM, "Color","black",'LineWidth',3);
plot(1:8, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:8)
xticklabels({'1min','2min','3min','4min','5min','6min','7min','8min'})

% ylim([0 100])
xlim([0.5 8.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');

title({'Freezing on the Context day',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

% [h1, p1, ci1, stats1] = ttest(freezing(:,1), freezing(:,2))
% 
% [h2, p2, ci2, stats2] = ttest(freezing(:,2), freezing(:,3))
