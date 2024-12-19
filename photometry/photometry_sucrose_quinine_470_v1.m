clc
clear
close all

%% set the path for data.

Directory = 'P:\Yihui\data\';    % Main directory\

sucrose_folder = {'m483_L_VTA_R_NAc_oxLight\Dec_09_2023\session_1';
    'm483_L_VTA_R_NAc_oxLight\Dec_12_2023\session_2';
    'm484_L_NAc_R_VTA_oxLight\Dec_09_2023\session_1';
    'm484_L_NAc_R_VTA_oxLight\Dec_12_2023\session_2';
    'm485_L_VTA_R_NAc_oxLight\Dec_12_2023\session_2'};

quinine_folder = {'m483_L_VTA_R_NAc_oxLight\Dec_09_2023\session_2';
    'm483_L_VTA_R_NAc_oxLight\Dec_12_2023\session_1';
    'm484_L_NAc_R_VTA_oxLight\Dec_09_2023\session_2';
    'm484_L_NAc_R_VTA_oxLight\Dec_12_2023\session_1';
    'm485_L_VTA_R_NAc_oxLight\Dec_12_2023\session_1'};

order = [0 0 1 1 0];

Sample_Rate = 200; % 200 scans per second.

Sucrose_max_zscore_ch1 = [];
Sucrose_max_zscore_ch2 = [];

Quinine_max_zscore_ch1 = [];
Quinine_max_zscore_ch2 = [];

for I=1:size(sucrose_folder,1)
    %% Sucrose
    
    Sucrose_Folder = [Directory sucrose_folder{I} '\'];
    
    Sucrose_Data = load([Sucrose_Folder 'datas.mat']).datas;
    Sucrose_times = load([Sucrose_Folder 'times.mat']).times;
    Sucrose_step_timepoint = load([Sucrose_Folder 'step_timepoint.mat']).step_timepoint;
    
    Sucrose_vid_start = ceil(Sucrose_step_timepoint(1))*Sample_Rate+1;
    
    Sucrose_timepoint = Sucrose_times(Sucrose_vid_start:end,1)';
    Sucrose_time = Sucrose_timepoint(1,:)-Sucrose_timepoint(1,1);
    Sucrose_total_time = floor(Sucrose_time(end));
    
    %% Analyse running signals
    
    Sucrose_running = Sucrose_Data(Sucrose_vid_start:end,2)';
    signedThreshold = 2^(32-1);
    Sucrose_running(Sucrose_running > signedThreshold) = Sucrose_running(Sucrose_running > signedThreshold) - 2^32;
    Sucrose_speedDeg = diff(Sucrose_running);
    Sucrose_Abs_Deg = abs(Sucrose_speedDeg);
    Sucrose_Abs_speed = movmean(Sucrose_Abs_Deg,100);
    Sucrose_index = find(~Sucrose_Abs_speed);
    
    %% Analyse photometry signals
    
    if order(I)
        raw_ch1 = Sucrose_Data(Sucrose_vid_start:end,3)';
        raw_ch2 = Sucrose_Data(Sucrose_vid_start:end,4)';
    else
        raw_ch1 = Sucrose_Data(Sucrose_vid_start:end,4)';
        raw_ch2 = Sucrose_Data(Sucrose_vid_start:end,3)';
    end
    
    % Smooth data
    smooth_win = 10;
    smooth_ch1 = movmean(raw_ch1,smooth_win);
    smooth_ch2 = movmean(raw_ch2,smooth_win);
    
    % Remove bleaching slope
    F = @(x,t)x(1)*exp(-t/x(2)) + x(3)*exp(-t/(x(4)*x(2))) + x(5);
    
    x0_ch1 = [max(smooth_ch1)/4, 3600, max(smooth_ch1)/4, 0.1, max(smooth_ch1)/2] ;
    lb_ch1 = [0, 600, 0, 0, 0];
    ub_ch1 = [max(smooth_ch1), 36000, max(smooth_ch1), 1, max(smooth_ch1)];
    x_ch1 = lsqcurvefit(F, x0_ch1, Sucrose_time, smooth_ch1, lb_ch1, ub_ch1);
    base_ch1 = F(x_ch1,Sucrose_time);
    ch1 = smooth_ch1 - base_ch1;
    
    x0_ch2 = [max(smooth_ch2)/4, 3600, max(smooth_ch2)/4, 0.1, max(smooth_ch2)/2] ;
    lb_ch2 = [0, 600, 0, 0, 0];
    ub_ch2 = [max(smooth_ch2), 36000, max(smooth_ch2), 1, max(smooth_ch2)];
    x_ch2 = lsqcurvefit(F, x0_ch2, Sucrose_time, smooth_ch2, lb_ch2, ub_ch2);
    base_ch2 = F(x_ch2,Sucrose_time);
    ch2 = smooth_ch2 - base_ch2;
    
    % Exclude running timepoints
    Sucrose_ch1_static = ch1(Sucrose_index);
    Sucrose_ch2_static = ch2(Sucrose_index);
    
    % Standardize signals
    Sucrose_ch1_zscored = (Sucrose_ch1_static - mean(Sucrose_ch1_static)) / std(Sucrose_ch1_static);
    Sucrose_ch2_zscored = (Sucrose_ch2_static - mean(Sucrose_ch2_static)) / std(Sucrose_ch2_static);
    
    Sucrose_max_zscore_ch1 = [Sucrose_max_zscore_ch1; max(Sucrose_ch1_zscored)];
    Sucrose_max_zscore_ch2 = [Sucrose_max_zscore_ch2; max(Sucrose_ch2_zscored)];
    
    
    %% Quinine
    
    Quinine_Folder = [Directory quinine_folder{I} '\'];
    
    Quinine_Data = load([Quinine_Folder 'datas.mat']).datas;
    Quinine_times = load([Quinine_Folder 'times.mat']).times;
    Quinine_step_timepoint = load([Quinine_Folder 'step_timepoint.mat']).step_timepoint;
    
    Quinine_vid_start = ceil(Quinine_step_timepoint(1))*Sample_Rate+1;
    
    Quinine_timepoint = Quinine_times(Quinine_vid_start:end,1)';
    Quinine_time = Quinine_timepoint(1,:)-Quinine_timepoint(1,1);
    Quinine_total_time = floor(Quinine_time(end));
    
    %% Analyse running signals
    
    Quinine_running = Quinine_Data(Quinine_vid_start:end,2)';
    signedThreshold = 2^(32-1);
    Quinine_running(Quinine_running > signedThreshold) = Quinine_running(Quinine_running > signedThreshold) - 2^32;
    Quinine_speedDeg = diff(Quinine_running);
    Quinine_Abs_Deg = abs(Quinine_speedDeg);
    Quinine_Abs_speed = movmean(Quinine_Abs_Deg,100);
    Quinine_index = find(~Quinine_Abs_speed);
    
    %% Analyse photometry signals
    
    if order(I)
        raw_ch1 = Quinine_Data(Quinine_vid_start:end,3)';
        raw_ch2 = Quinine_Data(Quinine_vid_start:end,4)';
    else
        raw_ch1 = Quinine_Data(Quinine_vid_start:end,4)';
        raw_ch2 = Quinine_Data(Quinine_vid_start:end,3)';
    end
    
    % Smooth data
    smooth_win = 10;
    smooth_ch1 = movmean(raw_ch1,smooth_win);
    smooth_ch2 = movmean(raw_ch2,smooth_win);
    
    % Remove bleaching slope
    F = @(x,t)x(1)*exp(-t/x(2)) + x(3)*exp(-t/(x(4)*x(2))) + x(5);
    
    x0_ch1 = [max(smooth_ch1)/4, 3600, max(smooth_ch1)/4, 0.1, max(smooth_ch1)/2] ;
    lb_ch1 = [0, 600, 0, 0, 0];
    ub_ch1 = [max(smooth_ch1), 36000, max(smooth_ch1), 1, max(smooth_ch1)];
    x_ch1 = lsqcurvefit(F, x0_ch1, Quinine_time, smooth_ch1, lb_ch1, ub_ch1);
    base_ch1 = F(x_ch1,Quinine_time);
    ch1 = smooth_ch1 - base_ch1;
    
    x0_ch2 = [max(smooth_ch2)/4, 3600, max(smooth_ch2)/4, 0.1, max(smooth_ch2)/2] ;
    lb_ch2 = [0, 600, 0, 0, 0];
    ub_ch2 = [max(smooth_ch2), 36000, max(smooth_ch2), 1, max(smooth_ch2)];
    x_ch2 = lsqcurvefit(F, x0_ch2, Quinine_time, smooth_ch2, lb_ch2, ub_ch2);
    base_ch2 = F(x_ch2,Quinine_time);
    ch2 = smooth_ch2 - base_ch2;
    
    % Exclude running timepoints
    Quinine_ch1_static = ch1(Quinine_index);
    Quinine_ch2_static = ch2(Quinine_index);
    
    % Standardize signals
    Quinine_ch1_zscored = (Quinine_ch1_static - mean(Quinine_ch1_static)) / std(Quinine_ch1_static);
    Quinine_ch2_zscored = (Quinine_ch2_static - mean(Quinine_ch2_static)) / std(Quinine_ch2_static);
    
    Quinine_max_zscore_ch1 = [Quinine_max_zscore_ch1; max(Quinine_ch1_zscored)];
    Quinine_max_zscore_ch2 = [Quinine_max_zscore_ch2; max(Quinine_ch2_zscored)];
    
end

%%

figure;

subplot(2,1,1);
hold on

speed_ch1 = [Sucrose_max_zscore_ch1 Quinine_max_zscore_ch1];

model_series = [mean(Sucrose_max_zscore_ch1); mean(Quinine_max_zscore_ch1)];
model_error = [std(Sucrose_max_zscore_ch1)/sqrt(size(speed_ch1,1)); std(Quinine_max_zscore_ch1)/sqrt(size(speed_ch1,1))];

for k = 1:size(speed_ch1,1)
    plot(1:2,speed_ch1(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[121 191 112]./255,'markerfacecolor',[121 191 112]./255,...
        'linestyle','-','color',[121./255 191./255 112./255 0.9],'linewidth',3);
end

bar(1,model_series(1), 'grouped','FaceColor','none','EdgeColor',[227 26 28]/255,'linewidth',4);
bar(2,model_series(2), 'grouped','FaceColor','none','EdgeColor',[33 113 181]/255,'linewidth',4);

errorbar(1,model_series(1),model_error(1),'linestyle','none','linewidth',2,'color',[227 26 28]/255,'CapSize',15);
errorbar(2,model_series(2),model_error(2),'linestyle','none','linewidth',2,'color',[33 113 181]/255,'CapSize',15);

xticklabels({})
xlim([0.3 2.8])
ylabel('Max z score (s.d.)','FontSize',18,'FontWeight','bold');
title('oxLight @ NAc','FontSize',20,'FontWeight','bold','color',[121 191 112]./255)
hold off

subplot(2,1,2);
hold on

speed_ch2 = [Sucrose_max_zscore_ch2 Quinine_max_zscore_ch2];

model_series = [mean(Sucrose_max_zscore_ch2); mean(Quinine_max_zscore_ch2)];
model_error = [std(Sucrose_max_zscore_ch2)/sqrt(size(speed_ch2,1)); std(Quinine_max_zscore_ch2)/sqrt(size(speed_ch2,1))];

for k = 1:size(speed_ch2,1)
    plot(1:2,speed_ch2(k,:),'marker','o','markersize',4,...
        'markeredgecolor',[64 145 48]./255,'markerfacecolor',[64 145 48]./255,...
        'linestyle','-','color',[64./255 145./255 48./255 0.9],'linewidth',3);
end

bar(1,model_series(1), 'grouped','FaceColor','none','EdgeColor',[227 26 28]/255,'linewidth',4);
bar(2,model_series(2), 'grouped','FaceColor','none','EdgeColor',[33 113 181]/255,'linewidth',4);

errorbar(1,model_series(1),model_error(1),'linestyle','none','linewidth',2,'color',[227 26 28]/255,'CapSize',15);
errorbar(2,model_series(2),model_error(2),'linestyle','none','linewidth',2,'color',[33 113 181]/255,'CapSize',15);

xticklabels({})
xlim([0.3 2.8])
ylabel('Max z score (s.d.)','FontSize',18,'FontWeight','bold');
title('oxLight @ VTA','FontSize',20,'FontWeight','bold','color',[64 145 48]./255)
hold off







