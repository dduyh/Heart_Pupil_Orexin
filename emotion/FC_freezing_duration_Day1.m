%%
clear; close all; clc;

%%
freezing =[0 7.77777 30.5555;
    12.4438 71.4445 73.3334
    13.8889 34.6666 80.002
    56.4451 36.6653 48.9979
    0 4 52.6667
    31.8896 46.3334 61.2243
    16 53.2236 63.4465];

% freezing =[15.8349 27.4985 8.61112 5.83334 20 14.1651 12.2238 14.4461 24.4428 30.2794 37.2205 15.0017 28.8906 23.6129 48.3316;
%     33.902 35.8268 33.0555 27.2156 23.6178 42.2256 50.559 35.5488 55.8367 52.7778 43.3264 55.5625 60.8333 53.3333 47.4965;
% 55.5656 66.6717 73.0505 82.2324 69.9898 61.3837 63.071 89.1718 77.2119 66.6667 46.3785 47.5052 33.3439 81.6667 82.2117]';

figure
hold on

for k = 1:size(freezing,1)
    plot(1:3,freezing(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(freezing);
SEM = S/sqrt(7);
errorbar(1:3, M, SEM, "Color","black",'LineWidth',3);
plot(1:3, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks([1 2 3])
xticklabels({'Tone1','Tone2','Tone3'})

ylim([0 100])
xlim([0.5 3.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');
% 
% title({'Freezing during Tone on the Conditioning day',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h1, p1, ci1, stats1] = ttest(freezing(:,1), freezing(:,2))

[h2, p2, ci2, stats2] = ttest(freezing(:,2), freezing(:,3))
