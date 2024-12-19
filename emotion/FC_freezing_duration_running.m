%%
clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm484_L_NAc_R_VTA_oxLight';            % Mouse name\
date = 'Mar_21_2024';                             % Date\

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

%%

Sample_Rate = 200;    % 200 scans per second.

vid_start = ceil(step_timepoint(1)*Sample_Rate)+1;

Tone_start = ceil(step_timepoint(2)*Sample_Rate)+1;
Tone_end = ceil(step_timepoint(2)*Sample_Rate)+480*Sample_Rate;

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

%% Analyse running signals

running = datas(vid_start:end,1)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
Abs_speedDeg = abs(speedDeg);
speed = movmean(Abs_speedDeg,100);

freezing_timepoint = find(speed<0.004);

c1 = 1;
arrset = cell(0,0);
while(c1<numel(freezing_timepoint))
    c2 = 0;
    while (c1+c2+1<=numel(freezing_timepoint)&&freezing_timepoint(c1)+c2+1==freezing_timepoint(c1+c2+1))
        c2 = c2+1;
    end
    if(c2>=1)
        arrset= [arrset;(freezing_timepoint(c1:1:c1+c2))];
    end
    c1 = c1 + c2 +1;
end
fprintf('有%d组连续数：\n',numel(arrset))

duration = zeros(1,length(arrset));

for i = 1:length(arrset)
    duration(i) = size(arrset{i},2);
end

duration_thershold = 0.75 * Sample_Rate;
arrset(duration<duration_thershold)=[];

freezing_timepoint = [];
for i = 1:length(arrset)
    freezing_timepoint = [freezing_timepoint arrset{i}];
end

duration_new = zeros(1,length(arrset));

for i = 1:length(arrset)
    duration_new(i) = size(arrset{i},2);
end

freezing = zeros(1,length(timepoint));
freezing(freezing_timepoint) = 1;

freezing_extinction_trace = freezing(Tone_start:Tone_end);

freezing_duration = zeros(1,8);
for i = 1:8
    freezing_trace = freezing_extinction_trace(((i-1)*60*Sample_Rate+1):i*60*Sample_Rate);
    freezing_duration(i) = length(find(freezing_trace == 1))/(60*Sample_Rate);
end
%% Plot running, freezing states

figure(1)
subplot(2,1,1);
plot(time,[speed speed(end)],'k')
xlim([0 time(end)])
ylim([0 0.5])
title('Running Speed','FontSize',15,'FontWeight','bold')
ylabel('Speed','FontSize',15,'FontWeight','bold')
axis off

h2 = subplot(2,1,2);
mymap = [235 145 132
    128 172 249];
mymap = mymap./255;
imagesc(freezing)
colormap(h2,mymap);
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
title('Freezing States','FontSize',15,'FontWeight','bold','color',[128 172 249]./255)
axis off

%%
figure
hold on

for k = 1:size(freezing_duration,1)
    plot(1:8,freezing_duration(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

% [S,M] = std(freezing);
% SEM = S/sqrt(8);
% errorbar(1:8, M, SEM, "Color","black",'LineWidth',3);
% plot(1:8, M,'o-','color','k','linewidth',4,'markeredgecolor','k','markerfacecolor','k','markersize',4);

xticks(1:8)
xticklabels({'1min','2min','3min','4min','5min','6min','7min','8min'})

% ylim([0 100])
% xlim([0.5 3.5])

ylabel('Freezing Duration (%)','FontSize',20,'FontWeight','bold');

title({'Freezing on the Extinction day',''},'FontSize',24,'FontWeight','bold')

ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

hold off

