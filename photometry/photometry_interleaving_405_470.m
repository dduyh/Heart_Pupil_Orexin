%% Fiber photometry data processing

clear; close all; clc;

%% set the path for output data

Directory = 'P:\Yihui\data\';                     % Main directory\
mouse_name = 'm751_L_Orx_GCaMP6s_R_Orx_Chrimson';            % Mouse name\
date = 'Dec_08_2023';                             % Date\
session = 'session_2';

Data_Folder = [Directory mouse_name '\' date '\' session '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

Sample_Rate = 200; % 200 scans per second.
vid_start = ceil(step_timepoint(1))*Sample_Rate+2;

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

raw_ch1 = datas(vid_start:end,3)';
raw_ch2 = datas(vid_start:end,4)';

% Plot raw data
figure(1);
subplot(2,1,1)
plot(time,raw_ch1,'b')
ylabel('Channel 1 (V)','FontSize',10,'FontWeight','bold')
title('Raw signals (m751, Dec 08 2023, session 2, Quinine)','FontSize',10,'FontWeight','bold')

subplot(2,1,2)
plot(time,raw_ch2,'b')
xlabel('Time (seconds)','FontSize',10,'FontWeight','bold')
ylabel('Channel 2 (V)','FontSize',10,'FontWeight','bold')

%% Seperate two channels

index_470 = find(~repmat([ones(1,10) zeros(1,10)],1,10*total_time));
index_405 = find(repmat([ones(1,10) zeros(1,10)],1,10*total_time));

raw_470_ch1 = raw_ch1(index_470);
raw_405_ch1 = raw_ch1(index_405);
raw_470_ch2 = raw_ch2(index_470);
raw_405_ch2 = raw_ch2(index_405);

time_470 = time(index_470);
time_405 = time(index_405);

% Plot raw data
figure(2)

subplot(4,1,1)
plot(time_470,raw_470_ch1,'color',[85 196 212]./255)
ylabel('Ch 1 (V)','FontSize',10,'FontWeight','bold')
legend('470 nm','FontSize',10,'TextColor',[85 196 212]./255)
legend('boxoff')
title('Raw signals','FontSize',10,'FontWeight','bold')

subplot(4,1,2)
plot(time_405,raw_405_ch1,'color',[122 88 164]./255)
ylabel('Ch 1 (V)','FontSize',10,'FontWeight','bold')
legend('405 nm','FontSize',10,'TextColor',[122 88 164]./255)
legend('boxoff')

subplot(4,1,3)
plot(time_470,raw_470_ch2,'color',[85 196 212]./255)
ylabel('Ch 2 (V)','FontSize',10,'FontWeight','bold')
legend('470 nm','FontSize',10,'TextColor',[85 196 212]./255)
legend('boxoff')

subplot(4,1,4)
plot(time_405,raw_405_ch2,'color',[122 88 164]./255)
xlabel('Time (seconds)','FontSize',10,'FontWeight','bold')
ylabel('Ch 2 (V)','FontSize',10,'FontWeight','bold')
legend('405 nm','FontSize',10,'TextColor',[122 88 164]./255)
legend('boxoff')

%% Smooth data

smooth_win = 10;
smooth_470_ch1 = movmean(raw_470_ch1,smooth_win);
smooth_405_ch1 = movmean(raw_405_ch1,smooth_win);
smooth_470_ch2 = movmean(raw_470_ch2,smooth_win);
smooth_405_ch2 = movmean(raw_405_ch2,smooth_win);

%% Remove bleaching slope

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

% Plot smoothed signals with bleaching slopes
figure(3)
subplot(4,1,1)
hold on
title('Smoothed signals with double exponential bleaching fits','FontSize',10,'FontWeight','bold')
plot(time_470, smooth_470_ch1, 'color', [85 196 212]./255)
plot(time_470, base_470_ch1, 'k-')
ylabel('Ch 1 (V)','FontSize',10,'FontWeight','bold')
hold off

subplot(4,1,2)
hold on
plot(time_405, smooth_405_ch1, 'color', [122 88 164]./255)
plot(time_405, base_405_ch1,'k-')
ylabel('Ch 1 (V)','FontSize',10,'FontWeight','bold')
hold off

subplot(4,1,3)
hold on
plot(time_470, smooth_470_ch2, 'color', [85 196 212]./255)
plot(time_470, base_470_ch2,'k-')
ylabel('Ch 2 (V)','FontSize',10,'FontWeight','bold')
hold off

subplot(4,1,4)
hold on
plot(time_405, smooth_405_ch2, 'color', [122 88 164]./255)
plot(time_405, base_405_ch2,'k-')
xlabel('Time (seconds)','FontSize',10,'FontWeight','bold')
ylabel('Ch 2 (V)','FontSize',10,'FontWeight','bold')
hold off

% Plot bleaching baseline substracted signals
figure(4)
subplot(4,1,1)
plot(time_470,ch1_470,'color',[85 196 212]./255)
ylabel('Ch 1 (V)','FontSize',10,'FontWeight','bold')
legend('470 nm','FontSize',10,'TextColor',[85 196 212]./255)
legend('boxoff')
title('Bleaching Correction by Double Exponential Fit','FontSize',10,'FontWeight','bold')

subplot(4,1,2)
plot(time_405,ch1_405,'color',[122 88 164]./255)
ylabel('Ch 1 (V)','FontSize',10,'FontWeight','bold')
legend('405 nm','FontSize',10,'TextColor',[122 88 164]./255)
legend('boxoff')

subplot(4,1,3)
plot(time_470,ch2_470,'color',[85 196 212]./255)
ylabel('Ch 2 (V)','FontSize',10,'FontWeight','bold')
legend('470 nm','FontSize',10,'TextColor',[85 196 212]./255)
legend('boxoff')

subplot(4,1,4)
plot(time_405,ch2_405,'color',[122 88 164]./255)
xlabel('Time (seconds)','FontSize',10,'FontWeight','bold')
ylabel('Ch 2 (V)','FontSize',10,'FontWeight','bold')
legend('405 nm','FontSize',10,'TextColor',[122 88 164]./255)
legend('boxoff')

%% Plot robust linear regression fit

fitdata_ch1 = fit(ch1_405',ch1_470',fittype('poly1'),'Robust','on');
fitdata_ch2 = fit(ch2_405',ch2_470',fittype('poly1'),'Robust','on');

figure(5)
subplot(2,1,1)
hold on
plot(ch1_405,ch1_470,'k.')
plot(fitdata_ch1,'b')
xlabel('ch1 405','FontSize',10,'FontWeight','bold')
ylabel('ch1 470','FontSize',10,'FontWeight','bold')
hold off

subplot(2,1,2)
hold on
plot(ch2_405,ch2_470,'k.')
plot(fitdata_ch2,'b')
xlabel('ch2 405','FontSize',10,'FontWeight','bold')
ylabel('ch2 470','FontSize',10,'FontWeight','bold')
hold off

%% Calculate difference between 470nm and 405nm signals

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

% Standardize signals
ch1_zscored = (ch1 - mean(ch1)) / std(ch1);
ch2_zscored = (ch2 - mean(ch2)) / std(ch2);

% Plot normalized signals
figure(6)
subplot(2,1,1)
plot(time_470,ch1_zscored,'color',[56 181 72]./255)
title('Orx-GCaMP6s L','FontSize',10,'FontWeight','bold')
ylabel('z-score','FontSize',10,'FontWeight','bold')

subplot(2,1,2)
plot(time_470,ch2_zscored,'color',[56 181 72]./255)
title('Orx-Chrimson R','FontSize',10,'FontWeight','bold')
ylabel('z-score','FontSize',10,'FontWeight','bold')
xlabel('Time (seconds)','FontSize',10,'FontWeight','bold')

%% Plot running and licking signals

licking = datas(vid_start:end,1)';

running = datas(vid_start:end,2)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);

% Plot normalized signals
figure(7)
subplot(4,1,1);
plot(time,licking,'b')
ylim([-0.5 1.5])
title('Licking','FontSize',15,'FontWeight','bold')

subplot(4,1,2);
plot(time,[speedDeg speedDeg(end)],'k')
ylim([-2 2])
title('Running','FontSize',15,'FontWeight','bold')

subplot(4,1,3)
plot(time_470,ch1_zscored,'color',[56 181 72]./255)
title('Orx-GCaMP6s L','FontSize',15,'FontWeight','bold')
ylabel('z-score','FontSize',15,'FontWeight','bold')

subplot(4,1,4)
plot(time_470,ch2_zscored,'color',[56 181 72]./255)
title('Orx-Chrimson R','FontSize',15,'FontWeight','bold')
ylabel('z-score','FontSize',15,'FontWeight','bold')
xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')

