%%

clear; close all; clc;

%% Set the path for output data

Directory = 'E:\data\';    % Main directory\

folder = {'m2070\Oct_16_2024';
    'm2070\Nov_09_2024';
    'm2070\Nov_15_2024';
    'm2071\Oct_16_2024';
    'm2071\Nov_09_2024';
    'm2071\Nov_15_2024';
    'm2072\Nov_09_2024';
    'm2072\Nov_15_2024'};

r2_full_all = [];
relRun = [];
relPup = [];
relHR = [];

for I=1:size(folder,1)

    Data_Folder = [Directory folder{I} '\'];

    load([Data_Folder 'regressions.mat']);

    r2_full_all = [r2_full_all r2_full];

    relRun = [relRun percContrib(:,1)'];

    relPup = [relPup percContrib(:,2)'];

    relHR = [relHR percContrib(:,3)'];



end

%% Regression plotings

figure

t = tiledlayout(4, 9, 'TileSpacing', 'compact', 'Padding', 'compact');

% sort by r squared
[~,idx]  = sort(r2_full_all);

% plot the relative contributions 
nexttile(t, 7, [1, 3]);
bar(r2_full_all(idx))
title('full model r^2','FontSize',15,'FontWeight','bold')

nexttile(t, 16, [1, 3]);
bar(relRun)
title('run contribution','FontSize',15,'FontWeight','bold')
ylim([0 100])

nexttile(t, 25, [1, 3]);
bar(relPup)
title('pupil contribution','FontSize',15,'FontWeight','bold')
ylim([0 100])

nexttile(t, 34, [1, 3]);
bar(relHR)
title('HR contribution','FontSize',15,'FontWeight','bold')
ylim([0 100])

% relative pupil vs relative heartrate 
nexttile(t, 1, [2, 2]);
scatter(relRun,relPup)
xlabel('relative Run contribution','FontSize',15,'FontWeight','bold')
ylabel('relative Pupil contribution','FontSize',15,'FontWeight','bold')
xlim([-10 110])
ylim([-10 110])
xlim([-10 110])
ylim([-10 110])
axis square

nexttile(t, 3, [2, 2]);
scatter(relPup,relHR)
xlabel('relative pupil contribution','FontSize',15,'FontWeight','bold')
ylabel('relative HR contribution','FontSize',15,'FontWeight','bold')
xlim([-10 110])
ylim([-10 110])
xlim([-10 110])
ylim([-10 110])
axis square

nexttile(t, 5, [2, 2]);
scatter(relHR,relRun)
xlabel('relative HR contribution','FontSize',15,'FontWeight','bold')
ylabel('relative Run contribution','FontSize',15,'FontWeight','bold')
xlim([-10 110])
ylim([-10 110])
xlim([-10 110])
ylim([-10 110])
axis square

% histograms
edges = [0 10 20 30 40 50 60 70 80 90 100];
nexttile(t, 19, [2, 2]);
histogram(relRun(idx),edges)
xlabel('% contribution','FontSize',15,'FontWeight','bold')
ylabel('# cells','FontSize',15,'FontWeight','bold')
title('Run Contributions','FontSize',15,'FontWeight','bold')
axis square

nexttile(t, 21, [2, 2]);
histogram(relPup(idx),edges)
xlabel('% contribution','FontSize',15,'FontWeight','bold')
ylabel('# cells','FontSize',15,'FontWeight','bold')
title('Pupil Contributions','FontSize',15,'FontWeight','bold')
axis square

nexttile(t, 23, [2, 2]);
histogram(relHR(idx),edges)
xlabel('% contribution','FontSize',15,'FontWeight','bold')
ylabel('# cells','FontSize',15,'FontWeight','bold')
title('HR Contributions','FontSize',15,'FontWeight','bold')
axis square


