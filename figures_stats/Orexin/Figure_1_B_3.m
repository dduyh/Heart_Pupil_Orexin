% %%
% clear; close all; clc;

%%
freezing_CS =[79.86667 99.6 100;         % m37
    72.8 77.73333 94.26667;              % m38
    91.6 99.46667 100;                   % m40
    99.73333 94.8 100;                   % m41
    86.93333 54 77.86667;                % m1027
    68.66667 75.6 94.53333;              % m50
    26.8 71.06667 93.33333;              % m2147
    20.66667 65.46667 93.73333;          % m2148
    15.46667 54.53333 89.73333;          % m2160
    10.53333 69.73333 90.66667;          % m2198
    49.6 65.73333 98.13333;              % m2216
    50.13333 91.86667 76.26667;          % m2220
    56.8 69.33333 83.6;                  % m2225
    50.8 53.2 60.66667];                 % m2227 

freezing_US =[99.6 93.6 51.8667;             % m37
    99.0667 70.5333 49.8667;                 % m38
    59.2 49.0667 51.2;                       % m40
    100 99.8667 74.5333;                     % m41
    72.9333 57.7333 57.4667;                 % m1027
    43.8667 39.2 63.3333;                    % m50
    24.5333 41.6 26.5333;                    % m2148
    54.5333 17.6 84.2667;                    % m2147
    62.1333 28.4 33.0667;                    % m2160
    75.3333 23.6 25.2;                       % m2198
    60.5333 46.2667 71.4667;                 % m2216
    68.8 49.4667 33.6;                       % m2220
    42.2667 58.5333 22;                      % m2225
    72.5333 41.2 62.6667];                   % m2227

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

errorbar(1:3, M_CS, SEM_CS, "Color",[233 102 94]/255,'LineWidth',3);
plot3 = plot(1:3, M_CS,'o-','color',[233 102 94]/255,'linewidth',4,'markeredgecolor',[233 102 94]/255,'markerfacecolor',[233 102 94]/255,'markersize',4);

errorbar(1:3, M_US, SEM_US, "Color",[128 204 187]/255,'LineWidth',3);
plot4 = plot(1:3, M_US,'o-','color',[128 204 187]/255,'linewidth',4,'markeredgecolor',[128 204 187]/255,'markerfacecolor',[128 204 187]/255,'markersize',4);

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

