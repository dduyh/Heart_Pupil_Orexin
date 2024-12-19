clc
clear
close all

%% set the path for data.

Directory = 'P:\Yihui\data\';    % Main directory\

sucrose_folder = {'m483_L_VTA_R_NAc_oxLight\Dec_09_2023\session_1';
    'm483_L_VTA_R_NAc_oxLight\Dec_12_2023\session_2';
    'm484_L_NAc_R_VTA_oxLight\Dec_09_2023\session_1';
    'm484_L_NAc_R_VTA_oxLight\Dec_12_2023\session_2';
    'm485_L_VTA_R_NAc_oxLight\Dec_12_2023\session_2';
    'm486_L_NAc_nLightR\Dec_10_2023\session_2';
    'm486_R_VTA_nLightR\Dec_09_2023\session_1';
    'm487_L_VTA_nLightR\Dec_10_2023\session_2';
    'm487_R_NAc_nLightR\Dec_09_2023\session_1';
    'm751_L_Orx_GCaMP6s_R_Orx_Chrimson\Dec_08_2023\session_1';
    'm751_L_Orx_GCaMP6s_R_Orx_Chrimson\Dec_10_2023\session_2';
    'm758_L_MCH_R_Orx_GCaMP\Dec_08_2023\session_1';
    'm758_L_MCH_R_Orx_GCaMP\Dec_10_2023\session_2'};

quinine_folder = {'m483_L_VTA_R_NAc_oxLight\Dec_09_2023\session_2';
    'm483_L_VTA_R_NAc_oxLight\Dec_12_2023\session_1';
    'm484_L_NAc_R_VTA_oxLight\Dec_09_2023\session_2';
    'm484_L_NAc_R_VTA_oxLight\Dec_12_2023\session_1';
    'm485_L_VTA_R_NAc_oxLight\Dec_12_2023\session_1';
    'm486_L_NAc_nLightR\Dec_10_2023\session_1';
    'm486_R_VTA_nLightR\Dec_09_2023\session_2';
    'm487_L_VTA_nLightR\Dec_10_2023\session_1';
    'm487_R_NAc_nLightR\Dec_09_2023\session_2';
    'm751_L_Orx_GCaMP6s_R_Orx_Chrimson\Dec_08_2023\session_2';
    'm751_L_Orx_GCaMP6s_R_Orx_Chrimson\Dec_10_2023\session_1';
    'm758_L_MCH_R_Orx_GCaMP\Dec_08_2023\session_2';
    'm758_L_MCH_R_Orx_GCaMP\Dec_10_2023\session_1'};

Sucrose = [];
Quinine = [];

for I=1:size(sucrose_folder,1)
    
    Sucrose_Folder = [Directory sucrose_folder{I} '\'];
    
    Sucrose_Data = load([Sucrose_Folder 'datas.mat']).datas;
    Sucrose_times = load([Sucrose_Folder 'times.mat']).times;
    Sucrose_step_timepoint = load([Sucrose_Folder 'step_timepoint.mat']).step_timepoint;
    
    Sample_Rate = 200; % 200 scans per second.
    Sucrose_vid_start = ceil(Sucrose_step_timepoint(1))*Sample_Rate+1;
    
    Sucrose_timepoint = Sucrose_times(Sucrose_vid_start:end,1)';
    Sucrose_time = Sucrose_timepoint(1,:)-Sucrose_timepoint(1,1);
    Sucrose_total_time = floor(Sucrose_time(end));
    
    Sucrose_running = Sucrose_Data(Sucrose_vid_start:end,2)';
    signedThreshold = 2^(32-1);
    Sucrose_running(Sucrose_running > signedThreshold) = Sucrose_running(Sucrose_running > signedThreshold) - 2^32;
    Sucrose_speedDeg = diff(Sucrose_running);
    Sucrose_Abs_Deg = abs(Sucrose_speedDeg);
    Sucrose_Abs_speed = sum(Sucrose_Abs_Deg)/Sucrose_total_time;
    Sucrose = [Sucrose; Sucrose_Abs_speed];
    
    Quinine_Folder = [Directory quinine_folder{I} '\'];
    
    Quinine_Data = load([Quinine_Folder 'datas.mat']).datas;
    Quinine_times = load([Quinine_Folder 'times.mat']).times;
    Quinine_step_timepoint = load([Quinine_Folder 'step_timepoint.mat']).step_timepoint;
    
    Sample_Rate = 200; % 200 scans per second.
    Quinine_vid_start = ceil(Quinine_step_timepoint(1))*Sample_Rate+1;
    
    Quinine_timepoint = Quinine_times(Quinine_vid_start:end,1)';
    Quinine_time = Quinine_timepoint(1,:)-Quinine_timepoint(1,1);
    Quinine_total_time = floor(Quinine_time(end));
    
    
    Quinine_running = Quinine_Data(Quinine_vid_start:end,2)';
    signedThreshold = 2^(32-1);
    Quinine_running(Quinine_running > signedThreshold) = Quinine_running(Quinine_running > signedThreshold) - 2^32;
    Quinine_speedDeg = diff(Quinine_running);
    Quinine_Abs_Deg = abs(Quinine_speedDeg);
    Quinine_Abs_speed = sum(Quinine_Abs_Deg)/Quinine_total_time;
    Quinine = [Quinine; Quinine_Abs_speed];
end

%%

figure;
hold on

speed = [Sucrose Quinine];

model_series = [mean(Sucrose); mean(Quinine)];
model_error = [std(Sucrose); std(Quinine)];

for k = 1:size(speed,1)
    plot(1:2,speed(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],...
        'linestyle','-','color',[0.75 0.75 0.75 0.9],'linewidth',3);
end

bar(1,model_series(1), 'grouped','FaceColor','none','EdgeColor',[227 26 28]/255,'linewidth',4);
bar(2,model_series(2), 'grouped','FaceColor','none','EdgeColor',[33 113 181]/255,'linewidth',4);

errorbar(1,model_series(1),model_error(1),'k','linestyle','none','linewidth',2,'color',[227 26 28]/255,'CapSize',15);
errorbar(2,model_series(2),model_error(2),'k','linestyle','none','linewidth',2,'color',[33 113 181]/255,'CapSize',15);

% xticklabels({'Sucrose','Quinine'})
hold off
% axis off
hy = ylabel('Speed (\circ/s)','FontSize',15,'FontWeight','bold');
set(hy,'fontname','Times New Roman','fontsize',18)
xlim([0.3 2.8])
