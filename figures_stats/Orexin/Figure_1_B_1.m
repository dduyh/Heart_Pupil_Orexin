% %%
% clear; close all; clc;

%%
freezing_CS =[57.33333 71.33333 96.53333;      % m39
    96.26667 95.06667 100;                     % m1026
    77.73333 87.46667 95.73333;                % m1028
    63.33333 64.93333 96.8;                    % m49
    20.66667 67.06667 90.13333;                % m2126
    28.26667 60.93333 87.86667;                % m2161
    20.26667 63.6 55.86667;                    % m2199
    20.4 58.66667 53.73333;                    % m2200
    26.66667 79.2 69.73333;                    % m2214
    69.33333 73.6 85.73333;                    % m2215
    10.26667 54.8 66.13333;                    % m2208
    20.4 59.6 46.4];                           % m2228
    

freezing_US =[63.7333 57.3333 50.6667;             % m39
    65.2 75.8667 40.1333;                          % m1026
    71.3333 42.8 38.9333;                          % m1028
    84.2667 55.3333 37.7333;                       % m49
    53.4667 35.7333 66.9333;                       % m2126
    72.5333 44.2667 69.0667;                       % m2161
    61.2 16.6667 30.6667;                          % m2199
    20.6667 27.6 28.4;                             % m2200
    20 38.6667 26.4;                               % m2214
    72.2667 46.8 38.5333;                          % m2215
    35.0667 32.4 23.3333;                          % m2208
    14.8 27.4667 18.2667];                         % m2228

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

fig = figure(3);
set(fig, 'Position', [2739 374 883 577]);

hold on

errorbar(1:3, M_CS, SEM_CS, "Color",[153 40 34]/255,'LineWidth',3);
plot7 = plot(1:3, M_CS,'o-','color',[153 40 34]/255,'linewidth',4,'markeredgecolor',[153 40 34]/255,'markerfacecolor',[153 40 34]/255,'markersize',4);

errorbar(1:3, M_US, SEM_US, "Color",[41 57 143]/255,'LineWidth',3);
plot8 = plot(1:3, M_US,'o-','color',[41 57 143]/255,'linewidth',4,'markeredgecolor',[41 57 143]/255,'markerfacecolor',[41 57 143]/255,'markersize',4);

xticks([1 2 3])
xticklabels({'Trial 1','Trial 2','Trial 3'})

ylim([30 100])
xlim([0.5 4.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold','color','k');

title('Freezing on Conditioning Phase','FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

legend([plot1, plot3, plot5, plot7, plot2, plot4, plot6, plot8], {'CS+ (DTR+ Fasted)', 'CS+ (DTR+ Fed)', 'CS+ (DTR- Fasted)', 'CS+ (DTR- Fed)', 'CS- (DTR+ Fasted)', 'CS- (DTR+ Fed)', 'CS- (DTR- Fasted)', 'CS- (DTR- Fed)'},'FontSize',12,'FontWeight','bold');
legend('boxoff')
legend('Location','best')

hold off

