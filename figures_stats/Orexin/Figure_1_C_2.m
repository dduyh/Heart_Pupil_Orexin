% %%
% clear; close all; clc;

%%
freezing =[84.6667 91.0833 96.0417 96.125;     % m55
    98.0417 98.2917 98.7917 100;               % m58
    78.375 87.125 86.2083 88.9583;             % m64
    97.0833 99.6667 90.625 99.875;             % m65
    48.7083 94.2083 73.125 96.2917;            % m2208
    63.4167 99.8333 73.75 100;                 % m2229
    68.2083 73.75 96.3333 95.9583;             % m2231
    52.5 71.1667 72.4167 44.125;               % m2209
    82.625 83.7083 94.125 91.9167;             % m2211
    97.7083 99.4167 81 99.7917;                % m2202
    75.7917 99.2917 88.0417 99.7083];          % m2203

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
title({'Freezing on Extinction day 1 (DTR- fasted)',''},'FontSize',24,'FontWeight','bold')

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
errorbar(1:4, M, SEM, "Color",[128 204 187]/255,'LineWidth',3);
plot3 = plot(1:4, M,'o-','color',[128 204 187]/255,'linewidth',4,'markeredgecolor',[128 204 187]/255,'markerfacecolor',[128 204 187]/255,'markersize',4);

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

