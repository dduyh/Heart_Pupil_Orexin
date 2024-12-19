%%
clear; close all; clc;

%% 
% pupil =[1.23111 97.9125;
% 4.43752 105.779;
% -0.771938 52.4837;
% -0.884926 32.5404];

pupil =[12.8815 92.6003;
39.2629 105.779;
13.058 94.8296;
27.3134 258.609];

figure
hold on

for k = 1:size(pupil,1)
    plot(1:2,pupil(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

[S,M] = std(pupil,'omitnan');
SEM = S/sqrt(4);
errorbar(1:2, M, SEM, "Color","black",'LineWidth',3);
plot(1:2, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'5 Hz','20 Hz'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('Pupil z-score','FontSize',20,'FontWeight','bold');

title({'Pupil Size (Chrimson)',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h1, p1, ci1, stats1] = ttest(pupil(:,1), pupil(:,2),'Tail','left')
