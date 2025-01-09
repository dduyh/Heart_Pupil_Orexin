clc
clear
close all

%% set the path for data.

% DTR- fed
pre_glucose_heartRate = [650; 634; 590; 732; 541; 666; 635];
post_glucose_heartRate = [478; 537; 431; 456; 517; 569; 525];

% % DTR+ fed
% pre_glucose_heartRate = [686; 680; 708; 671; 652; 712; 714; 680; 700];
% post_glucose_heartRate = [455; 460; 568; 562; 598; 600; 695; 612; 671];

% % DTR- fasted
% pre_glucose_heartRate = [554; 594; 633; 332; 483; 618];
% post_glucose_heartRate = [513; 479; 514; 319; 500; 580];

% % DTR+ fasted
% pre_glucose_heartRate = [647; 661; 557; 577; 598; 594; 631; 621; 716];
% post_glucose_heartRate = [519; 531; 535; 519; 579; 583; 641; 599; 675];

%%

heart_rate_data = [pre_glucose_heartRate; post_glucose_heartRate];

state = [ones(7,1); repmat(2,7,1)];
state = categorical(state,[1 2],{'Pre','Post'});

figure(1);
hold on

b1 = boxchart(state,heart_rate_data,'GroupByColor',state,'BoxFaceAlpha',1,'BoxWidth',2,'LineWidth',5,'BoxMedianLineColor','black','BoxEdgeColor','black','MarkerStyle','none');
b1(1).BoxFaceColor = '#DECAE0';
b1(2).BoxFaceColor = '#8B5C9E';

HR = [pre_glucose_heartRate post_glucose_heartRate];

for k = 1:size(HR,1)
    plot([0.8,2.2],HR(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

ylim([300 750])
ylabel('HR (bpm)','FontSize',20,'FontWeight','bold');
% title({'Heart Rate Before/After Propranolol Injection','',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.YAxis.Visible = 'off'; % remove y-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_HR, p_HR, ci_HR, stats_HR] = ttest(pre_glucose_heartRate, post_glucose_heartRate,"Tail","right")





