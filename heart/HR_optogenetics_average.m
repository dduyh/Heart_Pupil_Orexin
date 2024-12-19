%%
clear; close all; clc;

%% 
HR =[2.22939 -2.21968;
2.44458 -15.8444;
-4.09989 -4.24933;
-3.60726 -20.9214];

figure
hold on

for k = 1:size(HR,1)
    plot(1:2,HR(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(HR,'omitnan');
SEM = S/sqrt(4);
errorbar(1:2, M, SEM, "Color","black",'LineWidth',3);
plot(1:2, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'5 Hz','20 Hz'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('HR z-score','FontSize',20,'FontWeight','bold');

title({'Heart Rate (Chrimson)',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h1, p1, ci1, stats1] = ttest(HR(:,1), HR(:,2),'Tail','right')