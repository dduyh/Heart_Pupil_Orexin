clc
clear
close all

%% set the path for data.

Directory = 'P:\Yihui\data\';    % Main directory\

sucrose_folder = {'m751_L_Orx_GCaMP6s_R_Orx_Chrimson\Dec_08_2023\session_1';
    'm751_L_Orx_GCaMP6s_R_Orx_Chrimson\Dec_10_2023\session_2';
    'm758_L_MCH_R_Orx_GCaMP\Dec_08_2023\session_1';
    'm758_L_MCH_R_Orx_GCaMP\Dec_10_2023\session_2'};

quinine_folder = {'m751_L_Orx_GCaMP6s_R_Orx_Chrimson\Dec_08_2023\session_2';
    'm751_L_Orx_GCaMP6s_R_Orx_Chrimson\Dec_10_2023\session_1';
    'm758_L_MCH_R_Orx_GCaMP\Dec_08_2023\session_2';
    'm758_L_MCH_R_Orx_GCaMP\Dec_10_2023\session_1'};

order = [0 0 1 0];
start_index = [2 2 2 2; 2 2 1 2];

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
    
    Sucrose_vid_start = ceil(Sucrose_step_timepoint(1))*Sample_Rate + start_index(1,I);
    
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
    Sucrose_index = find(~Sucrose_Abs_speed(1:end-200));
    
    %% Analyse photometry signals
    
    if order(I)
        raw_ch1 = Sucrose_Data(Sucrose_vid_start:end,3)';
        raw_ch2 = Sucrose_Data(Sucrose_vid_start:end,4)';
    else
        raw_ch1 = Sucrose_Data(Sucrose_vid_start:end,4)';
        raw_ch2 = Sucrose_Data(Sucrose_vid_start:end,3)';
    end
    
    % Seperate two channels
    index_470 = find(~repmat([ones(1,10) zeros(1,10)],1,10*Sucrose_total_time));
    index_405 = find(repmat([ones(1,10) zeros(1,10)],1,10*Sucrose_total_time));
    
    raw_470_ch1 = raw_ch1(index_470);
    raw_405_ch1 = raw_ch1(index_405);
    raw_470_ch2 = raw_ch2(index_470);
    raw_405_ch2 = raw_ch2(index_405);
    
    time_470 = Sucrose_time(index_470);
    time_405 = Sucrose_time(index_405);
    
    % Smooth data
    smooth_win = 10;
    smooth_470_ch1 = movmean(raw_470_ch1,smooth_win);
    smooth_405_ch1 = movmean(raw_405_ch1,smooth_win);
    smooth_470_ch2 = movmean(raw_470_ch2,smooth_win);
    smooth_405_ch2 = movmean(raw_405_ch2,smooth_win);
    
    % Remove bleaching slope
    F = @(x,t)x(1)*exp(-t/x(2)) + x(3)*exp(-t/(x(4)*x(2))) + x(5);
    
    x0_470_ch1 = [max(smooth_470_ch1)/4, 3600, max(smooth_470_ch1)/4, 0.1, max(smooth_470_ch1)/2] ;
    lb_470_ch1 = [0, 600, 0, 0, 0];
    ub_470_ch1 = [max(smooth_470_ch1), 36000, max(smooth_470_ch1), 1, max(smooth_470_ch1)];
    x_470_ch1 = lsqcurvefit(F, x0_470_ch1, time_470, smooth_470_ch1, lb_470_ch1, ub_470_ch1);
    base_470_ch1 = F(x_470_ch1,time_470);
    ch1_470 = smooth_470_ch1 - base_470_ch1;
    
    x0_405_ch1 = [max(smooth_405_ch1)/4, 3600, max(smooth_405_ch1)/4, 0.1, max(smooth_405_ch1)/2] ;
    lb_405_ch1 = [0, 600, 0, 0, 0];
    ub_405_ch1 = [max(smooth_405_ch1), 36000, max(smooth_405_ch1), 1, max(smooth_405_ch1)];
    x_405_ch1 = lsqcurvefit(F, x0_405_ch1, time_405, smooth_405_ch1, lb_405_ch1, ub_405_ch1);
    base_405_ch1 = F(x_405_ch1,time_405);
    ch1_405 = smooth_405_ch1 - base_405_ch1;
    
    x0_470_ch2 = [max(smooth_470_ch2)/4, 3600, max(smooth_470_ch2)/4, 0.1, max(smooth_470_ch2)/2] ;
    lb_470_ch2 = [0, 600, 0, 0, 0];
    ub_470_ch2 = [max(smooth_470_ch2), 36000, max(smooth_470_ch2), 1, max(smooth_470_ch2)];
    x_470_ch2 = lsqcurvefit(F, x0_470_ch2, time_470, smooth_470_ch2, lb_470_ch2, ub_470_ch2);
    base_470_ch2 = F(x_470_ch2,time_470);
    ch2_470 = smooth_470_ch2 - base_470_ch2;
    
    x0_405_ch2 = [max(smooth_405_ch2)/4, 3600, max(smooth_405_ch2)/4, 0.1, max(smooth_405_ch2)/2] ;
    lb_405_ch2 = [0, 600, 0, 0, 0];
    ub_405_ch2 = [max(smooth_405_ch2), 36000, max(smooth_405_ch2), 1, max(smooth_405_ch2)];
    x_405_ch2 = lsqcurvefit(F, x0_405_ch2, time_405, smooth_405_ch2, lb_405_ch2, ub_405_ch2);
    base_405_ch2 = F(x_405_ch2,time_405);
    ch2_405 = smooth_405_ch2 - base_405_ch2;
    
    % Calculate difference between 470nm and 405nm signals
    fitdata_ch1 = fit(ch1_405',ch1_470',fittype('poly1'),'Robust','on');
    fitdata_ch2 = fit(ch2_405',ch2_470',fittype('poly1'),'Robust','on');
    
    if fitdata_ch1.p1 > 0
        ch1_405 = fitdata_ch1(ch1_405)';
        ch1 = ch1_470 - ch1_405;
    else
        ch1 = ch1_470;
    end
    
    if fitdata_ch2.p1 > 0
        ch2_405 = fitdata_ch2(ch2_405)';
        ch2 = ch2_470 - ch2_405;
    else
        ch2 = ch2_470;
    end
    
    % Exclude running timepoints
    Sucrose_ch1_static = ch1(round(Sucrose_index/2));
    Sucrose_ch2_static = ch2(round(Sucrose_index/2));
    
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
    
    Quinine_vid_start = ceil(Quinine_step_timepoint(1))*Sample_Rate + start_index(2,I);
    
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
    Quinine_index = find(~Quinine_Abs_speed(1:end-200));
    
    %% Analyse photometry signals
    
    if order(I)
        raw_ch1 = Quinine_Data(Quinine_vid_start:end,3)';
        raw_ch2 = Quinine_Data(Quinine_vid_start:end,4)';
    else
        raw_ch1 = Quinine_Data(Quinine_vid_start:end,4)';
        raw_ch2 = Quinine_Data(Quinine_vid_start:end,3)';
    end
    
        % Seperate two channels
    index_470 = find(~repmat([ones(1,10) zeros(1,10)],1,10*Quinine_total_time));
    index_405 = find(repmat([ones(1,10) zeros(1,10)],1,10*Quinine_total_time));
    
    raw_470_ch1 = raw_ch1(index_470);
    raw_405_ch1 = raw_ch1(index_405);
    raw_470_ch2 = raw_ch2(index_470);
    raw_405_ch2 = raw_ch2(index_405);
    
    time_470 = Quinine_time(index_470);
    time_405 = Quinine_time(index_405);
    
    % Smooth data
    smooth_win = 10;
    smooth_470_ch1 = movmean(raw_470_ch1,smooth_win);
    smooth_405_ch1 = movmean(raw_405_ch1,smooth_win);
    smooth_470_ch2 = movmean(raw_470_ch2,smooth_win);
    smooth_405_ch2 = movmean(raw_405_ch2,smooth_win);
    
    % Remove bleaching slope
    F = @(x,t)x(1)*exp(-t/x(2)) + x(3)*exp(-t/(x(4)*x(2))) + x(5);
    
    x0_470_ch1 = [max(smooth_470_ch1)/4, 3600, max(smooth_470_ch1)/4, 0.1, max(smooth_470_ch1)/2] ;
    lb_470_ch1 = [0, 600, 0, 0, 0];
    ub_470_ch1 = [max(smooth_470_ch1), 36000, max(smooth_470_ch1), 1, max(smooth_470_ch1)];
    x_470_ch1 = lsqcurvefit(F, x0_470_ch1, time_470, smooth_470_ch1, lb_470_ch1, ub_470_ch1);
    base_470_ch1 = F(x_470_ch1,time_470);
    ch1_470 = smooth_470_ch1 - base_470_ch1;
    
    x0_405_ch1 = [max(smooth_405_ch1)/4, 3600, max(smooth_405_ch1)/4, 0.1, max(smooth_405_ch1)/2] ;
    lb_405_ch1 = [0, 600, 0, 0, 0];
    ub_405_ch1 = [max(smooth_405_ch1), 36000, max(smooth_405_ch1), 1, max(smooth_405_ch1)];
    x_405_ch1 = lsqcurvefit(F, x0_405_ch1, time_405, smooth_405_ch1, lb_405_ch1, ub_405_ch1);
    base_405_ch1 = F(x_405_ch1,time_405);
    ch1_405 = smooth_405_ch1 - base_405_ch1;
    
    x0_470_ch2 = [max(smooth_470_ch2)/4, 3600, max(smooth_470_ch2)/4, 0.1, max(smooth_470_ch2)/2] ;
    lb_470_ch2 = [0, 600, 0, 0, 0];
    ub_470_ch2 = [max(smooth_470_ch2), 36000, max(smooth_470_ch2), 1, max(smooth_470_ch2)];
    x_470_ch2 = lsqcurvefit(F, x0_470_ch2, time_470, smooth_470_ch2, lb_470_ch2, ub_470_ch2);
    base_470_ch2 = F(x_470_ch2,time_470);
    ch2_470 = smooth_470_ch2 - base_470_ch2;
    
    x0_405_ch2 = [max(smooth_405_ch2)/4, 3600, max(smooth_405_ch2)/4, 0.1, max(smooth_405_ch2)/2] ;
    lb_405_ch2 = [0, 600, 0, 0, 0];
    ub_405_ch2 = [max(smooth_405_ch2), 36000, max(smooth_405_ch2), 1, max(smooth_405_ch2)];
    x_405_ch2 = lsqcurvefit(F, x0_405_ch2, time_405, smooth_405_ch2, lb_405_ch2, ub_405_ch2);
    base_405_ch2 = F(x_405_ch2,time_405);
    ch2_405 = smooth_405_ch2 - base_405_ch2;
    
    % Calculate difference between 470nm and 405nm signals
    fitdata_ch1 = fit(ch1_405',ch1_470',fittype('poly1'),'Robust','on');
    fitdata_ch2 = fit(ch2_405',ch2_470',fittype('poly1'),'Robust','on');
    
    if fitdata_ch1.p1 > 0
        ch1_405 = fitdata_ch1(ch1_405)';
        ch1 = ch1_470 - ch1_405;
    else
        ch1 = ch1_470;
    end
    
    if fitdata_ch2.p1 > 0
        ch2_405 = fitdata_ch2(ch2_405)';
        ch2 = ch2_470 - ch2_405;
    else
        ch2 = ch2_470;
    end
    
    % Exclude running timepoints
    Quinine_ch1_static = ch1(round(Quinine_index/2));
    Quinine_ch2_static = ch2(round(Quinine_index/2));
    
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

speed_ch1 = [Sucrose_max_zscore_ch1(3:4) Quinine_max_zscore_ch1(3:4)];

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
title('MCH-GCaMP6s','FontSize',20,'FontWeight','bold','color',[121 191 112]./255)
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
title('Orx-GCaMP6s','FontSize',20,'FontWeight','bold','color',[64 145 48]./255)
hold off







