clc
clear
close all

%%

DTR_fed = [0.962 0.537 0.397 0.843 0.461 -0.127 0.624 0.907 -0.955 0.999];

DTR_fasted = [1.253 1.131 1.076 -0.294 -0.755 -0.119 0.760 1.098 0.308 1.320];

control_fed = [1.11345 0.240 -0.453 1.160 -0.137 -0.149 0.137 0.360];

control_fasted = [1.020 1.442 0.563 0.197 0.935 1.082 -0.318 1.434];

%%

figure;
hold on

control = [control_fed' control_fasted'];
DTR = [DTR_fed' DTR_fasted'];

model_series = [mean(DTR_fed) mean(DTR_fasted); mean(control_fed) mean(control_fasted)];
model_error = [std(DTR_fed)/sqrt(10) std(DTR_fasted)/sqrt(10); std(control_fed)/sqrt(8) std(control_fasted)/sqrt(8)];

b = bar(model_series,'linewidth',3);
b(1).FaceColor = '#C6D6EA';
b(2).FaceColor = '#719DC9';

for k = 1:size(DTR,1)
    plot([0.8571, 1.1429],DTR(k,:),'marker','o','markersize',2,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

for k = 1:size(control,1)
    plot([1.8571, 2.1429],control(k,:),'marker','o','markersize',2,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

errorbar([0.8571 1.1429; 1.8571 2.1429],model_series,model_error,'k','linestyle','none','linewidth',3,'CapSize',15);

ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

ylabel('HR z-score','FontSize',15,'FontWeight','bold');
title({'Tone 1 response','',''},'FontSize',20,'FontWeight','bold','color','k')

xticks([1 2])
xticklabels({'DTR+','DTR-'})

%%

[h_HR_1, p_HR_1, ci_HR_1, stats_HR_1] = ttest2(DTR_fed, DTR_fasted)

[h_HR_2, p_HR_2, ci_HR_2, stats_HR_2] = ttest2(control_fed, control_fasted)

[h_HR_3, p_HR_3, ci_HR_3, stats_HR_3] = ttest2(DTR_fed,control_fed)

[h_HR_4, p_HR_4, ci_HR_4, stats_HR_4] = ttest2(DTR_fasted,control_fasted)





