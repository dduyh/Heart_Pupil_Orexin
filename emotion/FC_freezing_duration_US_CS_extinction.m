%%
clear; close all; clc;

%% 
freezing_day1 =[mean([34.4167 0]) mean([17.1667 37.4167]);
mean([3.33333 1.25]) mean([6.75 10.3333]);
mean([61.5 16.0833]) mean([45 62.5]);
mean([34.3333 28.25]) mean([53.3333 62.0833]);
mean([0 0]) mean([0 0]);
mean([32.5 16.3333]) mean([81.3333 76]);
mean([3.16667 0]) mean([8.25 0])];

freezing_day2 =[mean([23.4167 18.3333]) mean([35.4167 48.6667]);
mean([3.33333 0]) mean([50.3333 46.0833]);
mean([31.3333 77.9167]) mean([41.5 88.0833]);
mean([36.1667 21.4167]) mean([43.5833 61.9167]);
mean([0 0]) mean([0 0]);
mean([0 0]) mean([0 0]);
mean([0 0]) mean([0 0])];

freezing_day3 =[mean([55.9167 27.3333]) mean([43.5833 78.6667]);
mean([0 0]) mean([67.8333 53.4167]);
mean([76.0833 70.0833]) mean([67.6667 76.75]);
mean([45.3333 11.1667]) mean([73 66]);
mean([0 0]) mean([0 0]);
mean([0 0]) mean([0 0]);
mean([0 0]) mean([0 0])];

freezing_day4 =[22.9583 68.625;
2.04167 4.16667;
7.625 8.375;
39.6667 43.3333;
18.375 10.8333;
8.91667 28.7917;
0 0];

freezing_day5 =[41.8333 15.5833;
8.75 58.5417;
77.0833 75.2917;
19.5833 8.04167;
0 0;
8.91667 5.125;
0 0];

figure
hold on

[S_day1,M_day1] = std(freezing_day1,'omitnan');
SEM_day1 = S_day1/sqrt(7);
errorbar(1:2, M_day1, SEM_day1, "Color",'#CEE55D','LineWidth',3); %#D5E4A8
plot(1:2, M_day1,'o-','color','#CEE55D','linewidth',4,'markeredgecolor','#CEE55D','markerfacecolor','#CEE55D','markersize',4);

[S_day2,M_day2] = std(freezing_day2,'omitnan');
SEM_day2 = S_day2/sqrt(7);
errorbar(1:2, M_day2, SEM_day2, "Color",'#85B243','LineWidth',3);%#9BC985
plot(1:2, M_day2,'o-','color','#85B243','linewidth',4,'markeredgecolor','#85B243','markerfacecolor','#85B243','markersize',4);

[S_day3,M_day3] = std(freezing_day3,'omitnan');
SEM_day3 = S_day3/sqrt(7);
errorbar(1:2, M_day3, SEM_day3, "Color",'#33a02d','LineWidth',3);
plot(1:2, M_day3,'o-','color','#33a02d','linewidth',4,'markeredgecolor','#33a02d','markerfacecolor','#33a02d','markersize',4);

[S_day4,M_day4] = std(freezing_day4,'omitnan');
SEM_day4 = S_day4/sqrt(7);
errorbar(1:2, M_day4, SEM_day4, "Color",'#1B450E','LineWidth',3);
plot(1:2, M_day4,'o-','color','#1B450E','linewidth',4,'markeredgecolor','#1B450E','markerfacecolor','#1B450E','markersize',4);

[S_day5,M_day5] = std(freezing_day5,'omitnan');
SEM_day5 = S_day5/sqrt(7);
errorbar(1:2, M_day5, SEM_day5, "Color",'k','LineWidth',3);
plot(1:2, M_day5,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

legend('Day 1','Day 1','Day 2','Day 2','Day 3','Day 3','Day 4','Day 4','Day 5','Day 5')
 
xticks(1:2)
xticklabels({'CS -','CS +'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');
% 
title({'Freezing on Extinction Days',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

% [h1, p1, ci1, stats1] = ttest(freezing(:,1), freezing(:,2))
% 
% [h2, p2, ci2, stats2] = ttest(freezing(:,2), freezing(:,3))
