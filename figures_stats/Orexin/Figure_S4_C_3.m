figure
hold on

for k = 1:size(peak_r_RunPup_values,1)
    plot(1,peak_r_RunPup_values(k),'marker','o','markersize',5,...
        'markeredgecolor',[124 190 174]./255,'markerfacecolor',[124 190 174]./255,...
        'linestyle','none');
end

for k = 1:size(peak_r_RunPup_values_DTR,1)
    plot(2,peak_r_RunPup_values_DTR(k),'marker','o','markersize',5,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

[S,M] = std(peak_r_RunPup_values,'omitnan');
SEM = S/sqrt(size(peak_r_RunPup_values,1));

[S_DTR,M_DTR] = std(peak_r_RunPup_values_DTR,'omitnan');
SEM_DTR = S_DTR/sqrt(size(peak_r_RunPup_values_DTR,1));

errorbar(1:2, [M M_DTR], [SEM SEM_DTR], "Color","black",'LineWidth',3);
plot(1:2, [M M_DTR],'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'DTR-','DTR+'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('Running-pupil c.c.','FontSize',20,'FontWeight','bold');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_RunPup_r, p_RunPup_r, ~, stats_RunPup_r] = ttest2(peak_r_RunPup_values,peak_r_RunPup_values_DTR)

%%
figure
hold on

for k = 1:size(peak_r_RunHR_values,1)
    plot(1,peak_r_RunHR_values(k),'marker','o','markersize',5,...
        'markeredgecolor',[124 190 174]./255,'markerfacecolor',[124 190 174]./255,...
        'linestyle','none');
end

for k = 1:size(peak_r_RunHR_values_DTR,1)
    plot(2,peak_r_RunHR_values_DTR(k),'marker','o','markersize',5,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

[S,M] = std(peak_r_RunHR_values,'omitnan');
SEM = S/sqrt(size(peak_r_RunHR_values,1));

[S_DTR,M_DTR] = std(peak_r_RunHR_values_DTR,'omitnan');
SEM_DTR = S_DTR/sqrt(size(peak_r_RunHR_values_DTR,1));

errorbar(1:2, [M M_DTR], [SEM SEM_DTR], "Color","black",'LineWidth',3);
plot(1:2, [M M_DTR],'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'DTR-','DTR+'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('Running-HR c.c.','FontSize',20,'FontWeight','bold');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_RunHR_r, p_RunHR_r, ~, stats_RunHR_r] = ttest2(peak_r_RunHR_values,peak_r_RunHR_values_DTR)

%%
figure
hold on

for k = 1:size(peak_r_HRPup_values,1)
    plot(1,peak_r_HRPup_values(k),'marker','o','markersize',5,...
        'markeredgecolor',[124 190 174]./255,'markerfacecolor',[124 190 174]./255,...
        'linestyle','none');
end

for k = 1:size(peak_r_HRPup_values_DTR,1)
    plot(2,peak_r_HRPup_values_DTR(k),'marker','o','markersize',5,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

[S,M] = std(peak_r_HRPup_values,'omitnan');
SEM = S/sqrt(size(peak_r_HRPup_values,1));

[S_DTR,M_DTR] = std(peak_r_HRPup_values_DTR,'omitnan');
SEM_DTR = S_DTR/sqrt(size(peak_r_HRPup_values_DTR,1));

errorbar(1:2, [M M_DTR], [SEM SEM_DTR], "Color","black",'LineWidth',3);
plot(1:2, [M M_DTR],'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'DTR-','DTR+'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('HR-Pupil c.c.','FontSize',20,'FontWeight','bold');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_HRPup_r, p_HRPup_r, ~, stats_HRPup_r] = ttest2(peak_r_HRPup_values,peak_r_HRPup_values_DTR)

%%
figure
hold on

for k = 1:size(peak_lags_RunPup_values,1)
    plot(1,peak_lags_RunPup_values(k),'marker','o','markersize',5,...
        'markeredgecolor',[124 190 174]./255,'markerfacecolor',[124 190 174]./255,...
        'linestyle','none');
end

for k = 1:size(peak_lags_RunPup_values_DTR,1)
    plot(2,peak_lags_RunPup_values_DTR(k),'marker','o','markersize',5,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

[S,M] = std(peak_lags_RunPup_values,'omitnan');
SEM = S/sqrt(size(peak_lags_RunPup_values,1));

[S_DTR,M_DTR] = std(peak_lags_RunPup_values_DTR,'omitnan');
SEM_DTR = S_DTR/sqrt(size(peak_lags_RunPup_values_DTR,1));

errorbar(1:2, [M M_DTR], [SEM SEM_DTR], "Color","black",'LineWidth',3);
plot(1:2, [M M_DTR],'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'DTR-','DTR+'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('Running-pupil lag (s)','FontSize',20,'FontWeight','bold');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_RunPup_lag, p_RunPup_lag, ~, stats_RunPup_lag] = ttest2(peak_lags_RunPup_values,peak_lags_RunPup_values_DTR)

%%
figure
hold on

for k = 1:size(peak_lags_RunHR_values,1)
    plot(1,peak_lags_RunHR_values(k),'marker','o','markersize',5,...
        'markeredgecolor',[124 190 174]./255,'markerfacecolor',[124 190 174]./255,...
        'linestyle','none');
end

for k = 1:size(peak_lags_RunHR_values_DTR,1)
    plot(2,peak_lags_RunHR_values_DTR(k),'marker','o','markersize',5,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

[S,M] = std(peak_lags_RunHR_values,'omitnan');
SEM = S/sqrt(size(peak_lags_RunHR_values,1));

[S_DTR,M_DTR] = std(peak_lags_RunHR_values_DTR,'omitnan');
SEM_DTR = S_DTR/sqrt(size(peak_lags_RunHR_values_DTR,1));

errorbar(1:2, [M M_DTR], [SEM SEM_DTR], "Color","black",'LineWidth',3);
plot(1:2, [M M_DTR],'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'DTR-','DTR+'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('Running-HR lag (s)','FontSize',20,'FontWeight','bold');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_RunHR_lag, p_RunHR_lag, ~, stats_RunHR_lag] = ttest2(peak_lags_RunHR_values,peak_lags_RunHR_values_DTR)

%%
figure
hold on

for k = 1:size(peak_lags_HRPup_values,1)
    plot(1,peak_lags_HRPup_values(k),'marker','o','markersize',5,...
        'markeredgecolor',[124 190 174]./255,'markerfacecolor',[124 190 174]./255,...
        'linestyle','none');
end

for k = 1:size(peak_lags_HRPup_values_DTR,1)
    plot(2,peak_lags_HRPup_values_DTR(k),'marker','o','markersize',5,...
        'markeredgecolor',[229 114 190]./255,'markerfacecolor',[229 114 190]./255,...
        'linestyle','none');
end

[S,M] = std(peak_lags_HRPup_values,'omitnan');
SEM = S/sqrt(size(peak_lags_HRPup_values,1));

[S_DTR,M_DTR] = std(peak_lags_HRPup_values_DTR,'omitnan');
SEM_DTR = S_DTR/sqrt(size(peak_lags_HRPup_values_DTR,1));

errorbar(1:2, [M M_DTR], [SEM SEM_DTR], "Color","black",'LineWidth',3);
plot(1:2, [M M_DTR],'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:2)
xticklabels({'DTR-','DTR+'})

% ylim([0 100])
xlim([0.5 2.5])

ylabel('HR-Pupil lag (s)','FontSize',20,'FontWeight','bold');

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

[h_HRPup_lag, p_HRPup_lag, ~, stats_HRPup_lag] = ttest2(peak_lags_HRPup_values,peak_lags_HRPup_values_DTR)

