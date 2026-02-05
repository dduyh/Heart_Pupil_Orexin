% %%
% clear; close all; clc;

%%
freezing_CS =[81.86667 88.53333 86.26667;   % m55
    99.46667 79.86667 86.93333;             % m58
    90.8 100 50.13333;                      % m64
    92.53333 79.73333 96.93333;             % m65
    23.33333 43.2 82.13333;                 % m2229
    27.46667 58.8 80.4;                     % m2231
    36 26.13333 86.53333;                   % m2209
    20.66667 66.66667 97.33333;             % m2211
    77.33333 58.53333 100;                  % m2202
    74 88 82.8];                            % m2203
    
freezing_US =[38.9333 34.1333 46.8;             % m55
    87.0667 62.9333 62.4;                       % m58
    77.0667 50.2667 32.2667;                    % m64
    80.8 64.5333 46.5333;                       % m65
    60.6667 54.5333 16.5333;                    % m2229
    18.1333 7.33333 15.4667;                    % m2231
    14.2667 20.8 22.1333;                       % m2209
    55.7333 67.2 41.2;                          % m2211
    78.4 84 68.6667;                            % m2202
    88.1333 75.7333 79.2];                      % m2203

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

errorbar(1:3, M_CS, SEM_CS, "Color",[232 44 43]/255,'LineWidth',3);
plot5 = plot(1:3, M_CS,'o-','color',[232 44 43]/255,'linewidth',4,'markeredgecolor',[232 44 43]/255,'markerfacecolor',[232 44 43]/255,'markersize',4);

errorbar(1:3, M_US, SEM_US, "Color",[45 127 184]/255,'LineWidth',3);
plot6 = plot(1:3, M_US,'o-','color',[45 127 184]/255,'linewidth',4,'markeredgecolor',[45 127 184]/255,'markerfacecolor',[45 127 184]/255,'markersize',4);

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

