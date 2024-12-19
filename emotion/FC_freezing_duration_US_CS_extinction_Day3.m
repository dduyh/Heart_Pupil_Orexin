%%
clear; close all; clc;

%% 
freezing =[mean([55.9167 27.3333]) mean([43.5833 78.6667]);
mean([0 0]) mean([67.8333 53.4167]);
mean([76.0833 70.0833]) mean([67.6667 76.75]);
mean([45.3333 11.1667]) mean([73 66]);
mean([0 0]) mean([0 0]);
mean([0 0]) mean([0 0]);
mean([0 0]) mean([0 0])];

figure
hold on

for k = 1:size(freezing,1)
    plot(1:2,freezing(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(freezing,'omitnan');
SEM = S/sqrt(7);
errorbar(1:2, M, SEM, "Color","black",'LineWidth',3);
plot(1:2, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'CS -','CS +'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');

title({'Freezing on Extinction Day 3',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

% [h1, p1, ci1, stats1] = ttest(freezing(:,1), freezing(:,2))
% 
% [h2, p2, ci2, stats2] = ttest(freezing(:,2), freezing(:,3))
