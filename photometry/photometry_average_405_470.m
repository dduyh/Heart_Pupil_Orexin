clc
clear
close all

%% set the path for data.

Directory = 'D:\data\';    % Main directory\

folder = {'m758_L_MCH_R_Orx_GCaMP\Dec_10_2023\session_2';
    'm758_L_MCH_R_Orx_GCaMP\Dec_10_2023\session_3';
    'm758_L_MCH_R_Orx_GCaMP\Dec_08_2023\session_2'};

Run_onsets = {[12270 46740 101600];
    [24820 54560];
    [3327 23390 45360]};

order = [0 0 1];
start_index = [2 2 1];

trials_num = 0;
for j = 1:size(folder,1)
    trials_num = trials_num + numel(Run_onsets{j});
end

Sample_Rate = 200;    % 200 scans per second.
FrameRate = 20;

run = NaN(trials_num,15*Sample_Rate);
trials_ch1 = NaN(trials_num,15*Sample_Rate/2);    % 15 seconds per trial.
trials_ch2 = NaN(trials_num,15*Sample_Rate/2);    % 15 seconds per trial.
pupil = NaN(trials_num,15*FrameRate);

k = 1;

for I=1:size(folder,1)
    
    Data_Folder = [Directory folder{I} '\'];
    
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');
    
    vid_start = ceil(step_timepoint(1))*Sample_Rate + start_index(I);
    
    timepoint = times(vid_start:end,1)';
    time = timepoint(1,:)-timepoint(1,1);
    total_time = floor(time(end));
    
    %% Analyse running signals
    
    running = datas(vid_start:end,2)';
    signedThreshold = 2^(32-1);
    running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
    speedDeg = diff(running);
    Abs_speedDeg = abs(speedDeg);
    speed = movmean(Abs_speedDeg,100);
    
    %% Analyse photometry signals
    
    if order(I)
        raw_ch1 = datas(vid_start:end,3)';
        raw_ch2 = datas(vid_start:end,4)';
    else
        raw_ch1 = datas(vid_start:end,4)';
        raw_ch2 = datas(vid_start:end,3)';
    end
    
    % Seperate two channels
    index_470 = find(~repmat([ones(1,10) zeros(1,10)],1,10*total_time));
    index_405 = find(repmat([ones(1,10) zeros(1,10)],1,10*total_time));
    
    raw_470_ch1 = raw_ch1(index_470);
    raw_405_ch1 = raw_ch1(index_405);
    raw_470_ch2 = raw_ch2(index_470);
    raw_405_ch2 = raw_ch2(index_405);
    
    time_470 = time(index_470);
    time_405 = time(index_405);
    
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
    
    %%
    Run_bout_onsets = Run_onsets{I};
    
    for II = 1:length(Run_bout_onsets)
        
        t_start = Run_bout_onsets(II) - 5*Sample_Rate +1;
        t_end = Run_bout_onsets(II) + 10*Sample_Rate;
        
        pupil_start = round(Run_bout_onsets(II)*FrameRate/Sample_Rate) - 5*FrameRate +1;
        pupil_end = round(Run_bout_onsets(II)*FrameRate/Sample_Rate) + 10*FrameRate;
        
        run(k,:) = speed(t_start:t_end);
        
        raw_trial_ch1 = ch1(ceil(t_start/2):(t_end/2));
        raw_trial_ch2 = ch2(ceil(t_start/2):(t_end/2));
        
        raw_pupil = areas(pupil_start:pupil_end);
        
        % Standardize signals
        trial_ch1_zscored = (raw_trial_ch1 - mean(raw_trial_ch1(1:2.5*Sample_Rate))) / std(raw_trial_ch1(1:2.5*Sample_Rate));
        trials_ch1(k,:) = trial_ch1_zscored;
        
        trial_ch2_zscored = (raw_trial_ch2 - mean(raw_trial_ch2(1:2.5*Sample_Rate))) / std(raw_trial_ch2(1:2.5*Sample_Rate));
        trials_ch2(k,:) = trial_ch2_zscored;
        
        pupil_zscored = (raw_pupil - mean(raw_pupil(1:5*FrameRate))) / std(raw_pupil(1:5*FrameRate));
        pupil(k,:) = pupil_zscored;
        
        k = k+1;
        
%                 figure(1)
%                 hold on
%                 plot(trial_ch1_zscored,'Color',[64 145 48]./255,'LineWidth',2)
    end
    
end

%%
figure;
xlims = (1:15*Sample_Rate)/Sample_Rate;
xlims_half = (1:15*Sample_Rate/2)/(Sample_Rate/2);
pupil_xlims = (1:15*FrameRate)/FrameRate;

subplot(4,1,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(run)+std(run)/sqrt(trials_num) fliplr(mean(run)-std(run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',0.7);
plot(xlims,mean(run),'Color','k','LineWidth',2)
line([5,5],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
axis off
title('Abs. Running Speed','FontSize',20,'FontWeight','bold')
hold off

subplot(4,1,2)
hold on
patch('XData',[xlims_half fliplr(xlims_half)],'YData',[mean(trials_ch1)+std(trials_ch1)/sqrt(trials_num) fliplr(mean(trials_ch1)-std(trials_ch1)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[211 229 209]./255,'FaceAlpha',1);
plot(xlims_half,mean(trials_ch1),'Color',[121 191 112]./255,'LineWidth',2)
line([5,5],[min(mean(trials_ch1)-std(trials_ch1)/sqrt(trials_num)),max(mean(trials_ch1)+std(trials_ch1)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
% axis off
title('MCH-GCaMP6s','FontSize',20,'FontWeight','bold','color',[121 191 112]./255)
hold off

subplot(4,1,3)
hold on
patch('XData',[xlims_half fliplr(xlims_half)],'YData',[mean(trials_ch2)+std(trials_ch2)/sqrt(trials_num) fliplr(mean(trials_ch2)-std(trials_ch2)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',0.7);
plot(xlims_half,mean(trials_ch2),'Color',[64 145 48]./255,'LineWidth',2)
line([5,5],[min(mean(trials_ch2)-std(trials_ch2)/sqrt(trials_num)),max(mean(trials_ch2)+std(trials_ch2)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
% axis off
title('Orx-GCaMP6s','FontSize',20,'FontWeight','bold','color',[64 145 48]./255)
hold off

subplot(4,1,4)
hold on
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2)
line([5,5],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
% axis off
title('Pupil Size','FontSize',20,'FontWeight','bold','color',[229 114 190]./255)
hold off


