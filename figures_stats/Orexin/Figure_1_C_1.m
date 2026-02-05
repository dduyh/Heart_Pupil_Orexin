% %%
% clear; close all; clc;

%%
freezing =[60.2083 69.5 68.0417 71.4167;       % m39
    65.9167 98.9167 91.5833 99.125;            % m1026
    91.5833 99.0833 87.3333 97.9583;           % m1028
    43.625 52.6667 47.8333 67.125;             % m49
    45.2917 99.5417 86.5 99.25;                % m2126
    69.1667 92.75 98.5417 99.7083;             % m2161
    64.9167 88.0417 66.9583 74.375;            % m2199
    68.5 99.2083 79.375 98.0833;               % m2200
    92.9583 90.9167 80.7083 99.5417;           % m2214
    75.6667 99.875 69.7917 99.75;              % m2215
    74.25 100 80.375 99.7917];                 % m2228

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
title({'Freezing on Extinction day 1 (DTR- fed)',''},'FontSize',24,'FontWeight','bold')

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
errorbar(1:4, M, SEM, "Color",[45 127 184]/255,'LineWidth',3);
plot4 = plot(1:4, M,'o-','color',[45 127 184]/255,'linewidth',4,'markeredgecolor',[45 127 184]/255,'markerfacecolor',[45 127 184]/255,'markersize',4);

xticks(1:4)
xticklabels({'CS -','CS +','CS -','CS +'})

ylim([30 100])
xlim([0.5 4.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');
title('Freezing on Cue Testing Phase', 'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

legend([plot1, plot2, plot3, plot4], {'DTR+ Fasted', 'DTR+ Fed', 'DTR- Fasted', 'DTR- Fed'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
legend('Location','best')

hold off

