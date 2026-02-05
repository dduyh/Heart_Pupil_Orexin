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
    'm2072\Nov_15_2024';
    'm2151\Dec_30_2024';
    'm2151\Jan_02_2025';
    'm2152\Dec_30_2024';
    'm2152\Jan_02_2025';
    'm2154\Dec_30_2024';
    'm2154\Jan_02_2025'};

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

%%

threshold = 20;

noneCells = [];   
pupilCells = [];   
heartCells = [];    
bothCells = [];  

noneCells = find(relRun < threshold & relPup < threshold & relHR < threshold); 
pupilCells = find(relPup > threshold & relRun < threshold & relHR < threshold); % Pupil only
heartCells = find(relHR > threshold & relRun < threshold & relPup < threshold); % HR only
runCells = find(relRun > threshold & relPup < threshold & relHR < threshold); % Running only
bothCells = find(relRun > threshold | relPup > threshold | relHR > threshold); % Mixed
bothCells = setdiff(bothCells, [noneCells pupilCells heartCells runCells]); % Exclude others

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
bar(relRun(idx))
title('run contribution','FontSize',15,'FontWeight','bold')
ylim([0 100])

nexttile(t, 25, [1, 3]);
bar(relPup(idx))
title('pupil contribution','FontSize',15,'FontWeight','bold')
ylim([0 100])

nexttile(t, 34, [1, 3]);
bar(relHR(idx))
title('HR contribution','FontSize',15,'FontWeight','bold')
ylim([0 100])

nexttile(t, 1, [2, 2]);
hold on
scatter(relRun(noneCells), relPup(noneCells), [], [0.5, 0.5, 0.5],'filled');
scatter(relRun(pupilCells), relPup(pupilCells), [], [229 114 190]./255,'filled');
scatter(relRun(heartCells), relPup(heartCells), [], [255 128 128]./255,'filled');
scatter(relRun(runCells), relPup(runCells), [], [0.0, 0.0, 1.0], 'filled'); 
scatter(relRun(bothCells), relPup(bothCells), [], 'black','filled');
xline(threshold, '--k');
yline(threshold, '--k');
xlabel('relative Run contribution','FontSize',15,'FontWeight','bold')
ylabel('relative Pupil contribution','FontSize',15,'FontWeight','bold')
xlim([-10 110])
ylim([-10 110])
axis square
hold off

% relative pupil vs relative heartrate 
nexttile(t, 3, [2, 2]);
hold on
scatter(relHR(noneCells), relPup(noneCells), [], [0.5, 0.5, 0.5],'filled');
scatter(relHR(pupilCells), relPup(pupilCells), [], [229 114 190]./255,'filled');
scatter(relHR(heartCells), relPup(heartCells), [], [255 128 128]./255,'filled');
scatter(relHR(runCells), relPup(runCells), [], [0.0, 0.0, 1.0], 'filled'); 
scatter(relHR(bothCells), relPup(bothCells), [], 'black','filled');
xline(threshold, '--k');
yline(threshold, '--k');
xlabel('relative HR contribution','FontSize',15,'FontWeight','bold')
ylabel('relative pupil contribution','FontSize',15,'FontWeight','bold')
xlim([-10 110])
ylim([-10 110])
axis square
hold off

nexttile(t, 5, [2, 2]);
cellCounts = [length(noneCells), length(pupilCells), length(heartCells), length(runCells), length(bothCells)];
totalCells = sum(cellCounts);
percentages = round((cellCounts / totalCells) * 100);
pieColors = [0.5, 0.5, 0.5; 
    229./255 114./255 190./255; 
    255./255 128./255 128./255; 
    0.0, 0.0, 1.0;
    0.0, 0.0, 0.0];
cellLabels = {'none', 'Pupil', 'HR', 'Run', 'Both'};
p = pie(cellCounts);
for i = 1:2:length(p) 
    pieIndex = ceil(i / 2); 
    set(p(i), 'FaceColor', pieColors(pieIndex, :)); 
end
for i = 2:2:length(p) 
    pieIndex = ceil(i / 2); 

    centroid = mean(p(i-1).Vertices, 1); 

    text(centroid(1), centroid(2), sprintf('%s: %d%%', cellLabels{pieIndex}, percentages(pieIndex)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 15, 'Color', 'white','FontWeight','bold'); 
end
for i = 2:2:length(p)
    delete(p(i)); 
end
title(sprintf('n = %d cells', totalCells),'FontSize',15,'FontWeight','bold');

% histograms
edges = [0 10 20 30 40 50 60 70 80 90 100];
nexttile(t, 19, [2, 2]);
histogram(relRun(idx),edges)
ylim([0 300])
xlabel('% contribution','FontSize',15,'FontWeight','bold')
ylabel('# cells','FontSize',15,'FontWeight','bold')
title('Run Contributions','FontSize',15,'FontWeight','bold')
axis square

nexttile(t, 21, [2, 2]);
histogram(relPup(idx),edges)
ylim([0 300])
xlabel('% contribution','FontSize',15,'FontWeight','bold')
ylabel('# cells','FontSize',15,'FontWeight','bold')
title('Pupil Contributions','FontSize',15,'FontWeight','bold')
axis square

nexttile(t, 23, [2, 2]);
histogram(relHR(idx),edges)
ylim([0 300])
xlabel('% contribution','FontSize',15,'FontWeight','bold')
ylabel('# cells','FontSize',15,'FontWeight','bold')
title('HR Contributions','FontSize',15,'FontWeight','bold')
axis square

