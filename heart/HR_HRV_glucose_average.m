clc
clear
close all

%%

DTR_fed = [664.528484398650 661.434748408894 709.621265072213 649.079850588092 659.666189110814 712.924602238765 698.683922953942 666.126579602840];

% DTR_fed = [664.528484398650 661.434748408894 709.621265072213 667.939869445211 649.079850588092 659.666189110814 712.924602238765 698.683922953942 666.126579602840];

DTR_fasted = [621.235266430550 682.211564851960 608.670759041570 601.769082276858 583.614874595042 610.304653971274 591.440318188470 689.155922703436];

% DTR_fasted = [682.211564851960 621.235266430550 591.440318188470 608.670759041570 583.614874595042 601.769082276858 631.321676590102 610.304653971274 689.155922703436];

control_fed = [684.167078816214 613.416970918848 642.835089854995 588.793984765402 582.758295936355];

% control_fed = [684.167078816214 613.416970918848 537.967699383744 642.835089854995 577.990249971475 588.793984765402 582.758295936355];

control_fasted = [573.604179884029 600.218329673847 617.068354016641 395.424123305519 533.756903935411];

% control_fasted = [573.604179884029 544.522220664496 600.218329673847 395.424123305519 533.756903935411 617.068354016641];

%%

figure;
hold on

control = [control_fed' control_fasted'];
DTR = [DTR_fed' DTR_fasted'];

model_series = [mean(DTR_fed) mean(DTR_fasted); mean(control_fed) mean(control_fasted)];
model_error = [std(DTR_fed)/sqrt(9) std(DTR_fasted)/sqrt(9); std(control_fed)/sqrt(7) std(control_fasted)/sqrt(6)];

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

ylabel('bpm','FontSize',15,'FontWeight','bold');
title({'Heart Rate','',''},'FontSize',20,'FontWeight','bold','color','k')

xticks([1 2])
xticklabels({'DTR+','DTR-'})


%%

DTR_fed = [7.24310531330512 5.96297907979153 12.7718015304881 24.9999607879444 32.6925745718038 12.4724476901758 13.2191964647053 16.5075753901392];

% DTR_fed = [7.24310531330512 5.96297907979153 12.7718015304881 14.5447820686904 24.9999607879444 32.6925745718038 12.4724476901758 13.2191964647053 16.5075753901392];

DTR_fasted = [12.7407562915327 6.86105416746326 22.4181141286098 19.5095262240248 39.6940037116450 16.0023856246178 34.9891218529106 23.1875469647140];

% DTR_fasted = [6.86105416746326 12.7407562915327 34.9891218529106 22.4181141286098 39.6940037116450 19.5095262240248 9.96922551225401 16.0023856246178 23.1875469647140];

control_fed = [4.44642428234689 4.21302843695219 12.6565655984505 20.1307304123513 25.5085398266219];

% control_fed = [4.44642428234689 4.21302843695219 14.2626718555203 12.6565655984505 13.7138580682068 20.1307304123513 25.5085398266219];

control_fasted = [9.39040279062771 26.9390303431543 6.20146262016181 38.6412744532207 15.5103958542485];

% control_fasted = [9.39040279062771 12.5869784501544 26.9390303431543 38.6412744532207 15.5103958542485 6.20146262016181];

%%

figure;
hold on

control = [control_fed' control_fasted'];
DTR = [DTR_fed' DTR_fasted'];

model_series = [mean(DTR_fed) mean(DTR_fasted); mean(control_fed) mean(control_fasted)];
model_error = [std(DTR_fed)/sqrt(9) std(DTR_fasted)/sqrt(9); std(control_fed)/sqrt(7) std(control_fasted)/sqrt(6)];

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

ylabel('RMSSD (ms)','FontSize',15,'FontWeight','bold');
title({'Heart Rate Variability','',''},'FontSize',20,'FontWeight','bold','color','k')

xticks([1 2])
xticklabels({'DTR+','DTR-'})

%%

[h_HR_1, p_HR_1, ci_HR_1, stats_HR_1] = ttest(DTR_fed, DTR_fasted)

[h_HR_2, p_HR_2, ci_HR_2, stats_HR_2] = ttest(control_fed, control_fasted)

[h_HR_3, p_HR_3, ci_HR_3, stats_HR_3] = ttest2(DTR_fed,control_fed)

[h_HR_4, p_HR_4, ci_HR_4, stats_HR_4] = ttest2(DTR_fasted,control_fasted)


