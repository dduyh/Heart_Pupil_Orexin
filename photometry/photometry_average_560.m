clc
clear
close all

%% set the path for data.

Directory = 'E:\data\';    % Main directory\

folder = {'m487_L_VTA_nLightR\Dec_10_2023\session_1'};

Run_onsets = {[10370 35040 46830 58920 71120]};

trials_num = 0;
for j = 1:size(folder,1)
    trials_num = trials_num + numel(Run_onsets{j});
end

Sample_Rate = 200;    % 200 scans per second.
FrameRate = 20;

run = NaN(trials_num,15*Sample_Rate);    
trials = NaN(trials_num,15*Sample_Rate);    % 15 seconds per trial.
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
    
    raw = datas(vid_start:end,4)';
    
    % Smooth data
    smooth_win = 10;
    smooth = movmean(raw,smooth_win);
    
    % Remove bleaching slope
    F = @(x,t)x(1)*exp(-t/x(2)) + x(3)*exp(-t/(x(4)*x(2))) + x(5);
    x0 = [max(smooth)/4, 3600, max(smooth)/4, 0.1, max(smooth)/2] ;
    lb = [0, 600, 0, 0, 0];
    ub = [max(smooth), 36000, max(smooth), 1, max(smooth)];
    x = lsqcurvefit(F, x0, time, smooth, lb, ub);
    base = F(x,time);
    ch1 = smooth - base;
    
    %%
    Run_bout_onsets = Run_onsets{I};
    
    for II = 1:length(Run_bout_onsets)
        
        t_start = Run_bout_onsets(II) - 5*Sample_Rate +1;
        t_end = Run_bout_onsets(II) + 10*Sample_Rate;
        
        pupil_start = round(Run_bout_onsets(II)*FrameRate/Sample_Rate) - 5*FrameRate +1;
        pupil_end = round(Run_bout_onsets(II)*FrameRate/Sample_Rate) + 10*FrameRate;
        
        run(k,:) = speed(t_start:t_end);
        
        raw_trial = ch1(t_start:t_end);
        
        raw_pupil = areas(pupil_start:pupil_end);

        % Standardize signals
        trial_zscored = (raw_trial - mean(raw_trial(1:5*Sample_Rate))) / std(raw_trial(1:5*Sample_Rate));
        trials(k,:) = trial_zscored;
        
        pupil_zscored = (raw_pupil - mean(raw_pupil(1:5*FrameRate))) / std(raw_pupil(1:5*FrameRate));
        pupil(k,:) = pupil_zscored;
        
        k = k+1;
        
    end
    
end

%%
figure;
xlims = (1:15*Sample_Rate)/Sample_Rate;
pupil_xlims = (1:15*FrameRate)/FrameRate;

subplot(3,1,1);
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(run)+std(run)/sqrt(trials_num) fliplr(mean(run)-std(run)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[0.75 0.75 0.75],'FaceAlpha',0.7);
plot(xlims,mean(run),'Color','k','LineWidth',2)
line([5,5],[min(mean(run)-std(run)/sqrt(trials_num)),max(mean(run)+std(run)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
axis off
title('Abs. Running Speed','FontSize',20,'FontWeight','bold')
hold off

subplot(3,1,2)
hold on
patch('XData',[xlims fliplr(xlims)],'YData',[mean(trials)+std(trials)/sqrt(trials_num) fliplr(mean(trials)-std(trials)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[255 217 217]./255,'FaceAlpha',0.7);
plot(xlims,mean(trials),'Color',[255 128 128]./255,'LineWidth',2)
line([5,5],[min(mean(trials)-std(trials)/sqrt(trials_num)),max(mean(trials)+std(trials)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
axis off
title('nLightR @ VTA','FontSize',20,'FontWeight','bold','color',[255 128 128]./255)
hold off

subplot(3,1,3)
hold on
patch('XData',[pupil_xlims fliplr(pupil_xlims)],'YData',[mean(pupil)+std(pupil)/sqrt(trials_num) fliplr(mean(pupil)-std(pupil)/sqrt(trials_num))],'EdgeColor','none','FaceColor',[246 211 235]./255,'FaceAlpha',0.7);
plot(pupil_xlims,mean(pupil),'Color',[229 114 190]./255,'LineWidth',2)
line([5,5],[min(mean(pupil)-std(pupil)/sqrt(trials_num)),max(mean(pupil)+std(pupil)/sqrt(trials_num))],'Color','k','linestyle','--','LineWidth',2);
axis off
title('Pupil Size','FontSize',20,'FontWeight','bold','color',[229 114 190]./255)
hold off


