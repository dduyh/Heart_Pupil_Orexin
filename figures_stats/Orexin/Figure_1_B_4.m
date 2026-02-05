% %%
% clear; close all; clc;

%%
freezing_CS =[92.13333 98.26667 98.66667;   % m53
    98.4 99.86667 99.6;                     % m54
    73.46667 100 100;                       % m59
    97.6 99.2 100;                          % m63
    35.6 79.06667 62.93333;                 % m2230
    21.73333 65.86667 71.33333;             % m2212
    45.6 51.73333 63.2;                     % m2201
    72.4 96.4 71.46667;                     % m2204
    77.33333 63.6 90.26667];                % m2205

freezing_US =[86.4 59.3333 36.1333;             % m53
    98.2667 99.3333 99.3333;                    % m54
    87.8667 66.8 68.2667;                       % m59
    36.5333 55.4667 71.7333;                    % m63
    48.9333 30.8 26.4;                          % m2230
    75.2 52.9333 96;                            % m2212
    41.4667 33.0667 34.5333;                    % m2201
    86.1333 78.8 45.4667;                       % m2204
    81.7333 60.2667 8.13333];                   % m2205

%%
figure
hold on

for k = 1:size(freezing_CS,1)
    plot(1:3,freezing_CS(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S_CS,M_CS] = std(freezing_CS);
SEM_CS = S_CS/sqrt(size(freezing_CS,1));
errorbar(1:3, M_CS, SEM_CS, "Color",[227 26 28]/255,'LineWidth',3);
plot(1:3, M_CS,'o-','color',[227 26 28]/255,'linewidth',4,'markeredgecolor',[227 26 28]/255,'markerfacecolor',[227 26 28]/255,'markersize',4);

xticks([1 2 3])
xticklabels({'Trial 1','Trial 2','Trial 3'})

ylim([0 100])
xlim([0.5 3.5])

ylabel('CS+ Freezing Duration (%)','FontSize',20,'FontWeight','bold','color','k');

title({'Freezing during CS+ Tone',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

% [h1, p1, ci1, stats1] = ttest(freezing(:,1), freezing(:,2))
%
% [h2, p2, ci2, stats2] = ttest(freezing(:,2), freezing(:,3))

%%
figure
hold on

for k = 1:size(freezing_US,1)
    plot(1:3,freezing_US(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S_US,M_US] = std(freezing_US);
SEM_US = S_US/sqrt(size(freezing_US,1));
errorbar(1:3, M_US, SEM_US, "Color",[33 113 181]/255,'LineWidth',3);
plot(1:3, M_US,'o-','color',[33 113 181]/255,'linewidth',4,'markeredgecolor',[33 113 181]/255,'markerfacecolor',[33 113 181]/255,'markersize',4);

xticks([1 2 3])
xticklabels({'Trial 1','Trial 2','Trial 3'})

ylim([0 100])
xlim([0.5 3.5])

ylabel('CS- Freezing Duration (%)','FontSize',20,'FontWeight','bold','color','k');

title({'Freezing during CS- Tone',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

%%
figure(3)
hold on

errorbar(1:3, M_CS, SEM_CS, "Color",[234 167 118]/255,'LineWidth',3);
plot1 = plot(1:3, M_CS,'o-','color',[234 167 118]/255,'linewidth',4,'markeredgecolor',[234 167 118]/255,'markerfacecolor',[234 167 118]/255,'markersize',4);

errorbar(1:3, M_US, SEM_US, "Color",[199 226 179]/255,'LineWidth',3);
plot2 = plot(1:3, M_US,'o-','color',[199 226 179]/255,'linewidth',4,'markeredgecolor',[199 226 179]/255,'markerfacecolor',[199 226 179]/255,'markersize',4);

xticks([1 2 3])
xticklabels({'Trial 1','Trial 2','Trial 3'})

ylim([20 100])
xlim([0.5 3.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold','color','k');

title({'Freezing during Tone of the Conditioning phase',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

