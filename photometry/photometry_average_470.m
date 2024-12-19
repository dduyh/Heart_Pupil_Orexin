clc
clear
close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

folder = {'m483_L_VTA_R_NAc_oxLight\Dec_12_2023\session_1';
    'm483_L_VTA_R_NAc_oxLight\Dec_12_2023\session_2';
    'm483_L_VTA_R_NAc_oxLight\Dec_09_2023\session_1';
    'm483_L_VTA_R_NAc_oxLight\Dec_09_2023\session_2';
    'm484_L_NAc_R_VTA_oxLight\Dec_09_2023\session_1';
    'm484_L_NAc_R_VTA_oxLight\Dec_09_2023\session_2'};

Run_onsets = {[14770 38460 88500];
    [22830 59950];
    [2407 13930 38130];
    [48870 89240];
    [24720];
    [20570 47860];};

order = [0 0 0 0 1 1];

trials_num = 0;
for j = 1:size(folder,1)
    trials_num = trials_num + numel(Run_onsets{j});
end

Sample_Rate = 200;    % 200 scans per second.
FrameRate = 20;

run = NaN(trials_num,15*Sample_Rate);
trials_ch1 = NaN(trials_num,15*Sample_Rate);    % 15 seconds per trial.
trials_ch2 = NaN(trials_num,15*Sample_Rate);    % 15 seconds per trial.
pupil = NaN(trials_num,15*FrameRate);

k = 1;

for I=1:size(folder,1)
    
    Data_Folder = [Directory folder{I} '\'];
    
    load([Data_Folder 'datas.mat']);
    load([Data_Folder 'times.mat']);
    load([Data_Folder 'step_timepoint.mat']);
    load([Data_Folder 'Pupil.mat'], 'areas');
    
    vid_start = ceil(step_timepoint(1))*Sample_Rate+1;
    
    timepoint = times(vid_start:end,1)';
    time = timepoint(1,:)-timepoint(1,1);
    
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
    
    % Smooth data
    smooth_win = 10;
    smooth_ch1 = movmean(raw_ch1,smooth_win);
    smooth_ch2 = movmean(raw_ch2,smooth_win);
    
    % Remove bleaching slope
    F = @(x,t)x(1)*exp(-t/x(2)) + x(3)*exp(-t/(x(4)*x(2))) + x(5);
    
    x0_ch1 = [max(smooth_ch1)/4, 3600, max(smooth_ch1)/4, 0.1, max(smooth_ch1)/2] ;
    lb_ch1 = [0, 600, 0, 0, 0];
    ub_ch1 = [max(smooth_ch1), 36000, max(smooth_ch1), 1, max(smooth_ch1)];
    x_ch1 = lsqcurvefit(F, x0_ch1, time, smooth_ch1, lb_ch1, ub_ch1);
    base_ch1 = F(x_ch1,time);
    ch1 = smooth_ch1 - base_ch1;
    
    x0_ch2 = [max(smooth_ch2)/4, 3600, max(smooth_ch2)/4, 0.1, max(smooth_ch2)/2] ;
    lb_ch2 = [0, 600, 0, 0, 0];
    ub_ch2 = [max(smooth_ch2), 36000, max(smooth_ch2), 1, max(smooth_ch2)];
    x_ch2 = lsqcurvefit(F, x0_ch2, time, smooth_ch2, lb_ch2, ub_ch2);
    base_ch2 = F(x_ch2,time);
    ch2 = smooth_ch2 - base_ch2;
    
    %%
    Run_bout_onsets = Run_onsets{I};
    
    for II = 1:length(Run_bout_onsets)
        
        t_start = Run_bout_onsets(II) - 5*Sample_Rate +1;
        t_end = Run_bout_onsets(II) + 10*Sample_Rate;
        
        pupil_start = round(Run_bout_onsets(II)*FrameRate/Sample_Rate) - 5*FrameRate +1;
        pupil_end = round(Run_bout_onsets(II)*FrameRate/Sample_Rate) + 10*FrameRate;
        
        run(k,:) = speed(t_start:t_end);
        
        raw_trial_ch1 = ch1(t_start:t_end);
        raw_trial_ch2 = ch2(t_start:t_end);
        
        raw_pupil = areas(pupil_start:pupil_end);
        
        % Standardize signals
        trial_ch1_zscored = (raw_trial_ch1 - mean(raw_trial_ch1(1:5*Sample_Rate))) / std(raw_trial_ch1(1:5*Sample_Rate));
        trials_ch1(k,:) = trial_ch1_zscored;
        
        trial_ch2_zscored = (raw_trial_ch2 - mean(raw_trial_ch2(1:5*Sample_Rate))) / std(raw_trial_ch2(1:5*Sample_Rate));
        trials_ch2(k,:) = trial_ch2_zscored;
        
        pupil_zscored = (raw_pupil - mean(raw_pupil(1:5*FrameRate))) / std(raw_pupil(1:5*FrameRate));
        pupil(k,:) = pupil_zscored;
        
        k = k+1;
        
        %         figure(1)
        %         hold on
        %         plot(trial_ch1_zscored,'Color',[64 145 48]./255,'LineWidth',2)
    end
    
end

%%
figure;
xlims = (1:15*Sample_Rate)/Sample_Rate;
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
patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials_ch1)+std(trials_ch1)/sqrt(trials_num) fliplr(mean(trials_ch1)-std(trials_ch1)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[211 229 209]./255,'FaceAlpha',1);
plot(xlims,mean(trials_ch1),'Color',[121 191 112]./255,'LineWidth',2)
line([5,5],[min(mean(trials_ch1)-std(trials_ch1)/sqrt(trials_num)),max(mean(trials_ch1)+std(trials_ch1)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
axis off
title('oxLight @ NAc','FontSize',20,'FontWeight','bold','color',[121 191 112]./255)
hold off

subplot(4,1,3)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials_ch2)+std(trials_ch2)/sqrt(trials_num) fliplr(mean(trials_ch2)-std(trials_ch2)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[193 219 187]./255,'FaceAlpha',0.7);
plot(xlims,mean(trials_ch2),'Color',[64 145 48]./255,'LineWidth',2)
line([5,5],[min(mean(trials_ch2)-std(trials_ch2)/sqrt(trials_num)),max(mean(trials_ch2)+std(trials_ch2)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
axis off
title('oxLight @ VTA','FontSize',20,'FontWeight','bold','color',[64 145 48]./255)
hold off

subplot(4,1,4)
hold on
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2)
line([5,5],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
axis off
title('Pupil Size','FontSize',20,'FontWeight','bold','color',[229 114 190]./255)
hold off


