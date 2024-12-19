clc
clear
close all

%%
load('P:\Yihui\m1271\19092023m127 file1_superRigMainMat.mat')

dataLength=size(superRigMainMat,1);

licking = superRigMainMat(:,1);
lick = zeros(dataLength,1);
% lick(licking>0.5) = 1;
lick(licking>4) = 1;

sucrose = superRigMainMat(:,5);
sweet = zeros(dataLength,1);
sweet(sucrose>3) = 1;

water = superRigMainMat(:,6);
neutral = zeros(dataLength,1);
neutral(water>3) = 1;

qunine = superRigMainMat(:,7);
bitter = zeros(dataLength,1);
bitter(qunine>0.3) = 1;

%%

sucrose_duration = (1+[9702-8512, 33601-31886, 40944-40165, 64718-62731, 76959-71824, 83255-81184, 90820-88481, 95457-93044, 125045-122699, 130121-128133])./400; %seconds

water_duration = (1+[0, 48679-46481, 105461-103922, 110579-108752, 138796-136895])./400; %seconds

qunine_duration = (1+[26582-25321, 56111-53456, 115712-112794, 147051-145265])./400; %seconds

duration_mean = [mean(sucrose_duration); mean(water_duration); mean(qunine_duration)];
duration_error = [std(sucrose_duration)./sqrt(length(sucrose_duration)); std(water_duration)./sqrt(length(water_duration)); std(qunine_duration)./sqrt(length(qunine_duration))];

figure
b = bar(duration_mean,'grouped','FaceColor','none','linewidth',2);
hold on

x = 1:3;
errorbar(x,duration_mean,duration_error,'k','linestyle','none','linewidth',2);

xticklabels({'Sucrose','Water','Qunine'})
hold off

hy = ylabel('Licking Bout Duration (s)');
set(hy,'fontname','Times New Roman','fontsize',14)

ht = title('Licking Bout Duration of different Liquids');
set(ht,'fontname','Times New Roman','fontsize',18)

%%

sucrose_onset = [];
for k = 1:length(sweet)-2
        if sweet(k)==0 && sweet(k+1)==0 && sweet(k+2)==1
            sucrose_onset = [sucrose_onset k+2];
        end
end
    
sucrose_lick_offset = [];
% latency

sucrose_lick_timepoints = cell(size(sucrose_onset,2),1);
% sucrose_lick_rate = cell(size(sucrose_onset,2),1);
sucrose_lick_times = zeros(1,size(sucrose_onset,2));

for i = 1:size(sucrose_onset,2)
    
    sucrose_lick_series = lick(sucrose_onset(i):sucrose_lick_offset(i));
    
    lick_timepoints = [];
    for j = 1:length(sucrose_lick_series)
        if sucrose_lick_series(j)==0 && sucrose_lick_series(j+1)==1
            lick_timepoints = [lick_timepoints j+1];
        end
    end
    
    sucrose_lick_timepoints{i} = lick_timepoints;
    sucrose_lick_times(i) = length(lick_timepoints);
    
%     sucrose_lick_rate{i} = 
end

water_onset = [];
for k = 1:length(neutral)-2
        if neutral(k)==0 && neutral(k+1)==0 && neutral(k+2)==1
            water_onset = [water_onset k+2];
        end
end

water_lick_offset = [0 48679 110579 138796];

water_lick_timepoints = cell(size(water_onset,2),1);
% water_lick_rate = cell(size(water_onset,2),1);
water_lick_times = zeros(1,size(water_onset,2));

for i = 1:size(water_onset,2)
    
    water_lick_series = lick(water_onset(i):water_lick_offset(i));
    
    lick_timepoints = [];
    for j = 1:length(water_lick_series)
        if water_lick_series(j)==0 && water_lick_series(j+1)==1
            lick_timepoints = [lick_timepoints j+1];
        end
    end
    
    water_lick_timepoints{i} = lick_timepoints;
    water_lick_times(i) = length(lick_timepoints);
    
%     water_lick_rate{i} = 
end

qunine_onset = [];
for k = 1:length(bitter)-2
        if bitter(k)==0 && bitter(k+1)==0 && bitter(k+2)==1
            qunine_onset = [qunine_onset k+2];
        end
end

qunine_lick_offset = [26582 56111 115712 147051];

qunine_lick_timepoints = cell(size(qunine_onset,2),1);
% qunine_lick_rate = cell(size(qunine_onset,2),1);
qunine_lick_times = zeros(1,size(qunine_onset,2));

for i = 1:size(qunine_onset,2)
    
    qunine_lick_series = lick(qunine_onset(i):qunine_lick_offset(i));
    
    lick_timepoints = [];
    for j = 1:length(qunine_lick_series)
        if qunine_lick_series(j)==0 && qunine_lick_series(j+1)==1
            lick_timepoints = [lick_timepoints j+1];
        end
    end
    
    qunine_lick_timepoints{i} = lick_timepoints;
    qunine_lick_times(i) = length(lick_timepoints);
    
%     qunine_lick_rate{i} = 
end

times_mean = [mean(sucrose_lick_times); mean(water_lick_times); mean(qunine_lick_times)];
times_error = [std(sucrose_lick_times)./sqrt(length(sucrose_lick_times)); std(water_lick_times)./sqrt(length(water_lick_times)); std(qunine_lick_times)./sqrt(length(qunine_lick_times))];

figure
b = bar(times_mean,'grouped','FaceColor','none','linewidth',2);
hold on

x = 1:3;
errorbar(x,times_mean,times_error,'k','linestyle','none','linewidth',2);

xticklabels({'Sucrose','Water','Qunine'})
hold off

hy = ylabel('Bout size (licks/bout)');
set(hy,'fontname','Times New Roman','fontsize',14)

ht = title('Licking Times of different Liquids');
set(ht,'fontname','Times New Roman','fontsize',18)

%%
%ILI
%Lick duration


