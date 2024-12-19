%% ECG data processing

clear; close all; clc;

%% Set the path for output data

Directory = 'D:\data\';                     % Main directory\
mouse_name = 'm37';            % Mouse name\
date = 'Jul_18_2024';                             % Date\
stim = 'Fear_Conditioning'; % Sucrose Quinine

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% Load data

load([Data_Folder 'datas.mat']);
load([Data_Folder 'times.mat']);
load([Data_Folder 'step_timepoint.mat']);

pupil_data = csvread([Data_Folder stim 'DLC_resnet50_Pupil_trackingDec25shuffle1_1000000_filtered.csv'],3,1);

%% Analyse ECG signals

Sample_Rate = 1000;    % 1000 scans per second.

FrameRate = 20;

vid_start = ceil(step_timepoint(1)*Sample_Rate)+1;

timepoint = times(vid_start:end,1)';
time = timepoint(1,:)-timepoint(1,1);
total_time = floor(time(end));

smooth_window = 1;

%% Analyse running signals

running = datas(vid_start:end,1)';
signedThreshold = 2^(32-1);
running(running > signedThreshold) = running(running > signedThreshold) - 2^32;
speedDeg = diff(running);
Abs_speedDeg = abs(speedDeg);
speed = movmean(Abs_speedDeg,[smooth_window*Sample_Rate 0]);

%% Analyse pupil size

Pupil_up_x = pupil_data(:,1);
Pupil_up_y = pupil_data(:,2);
Pupil_left_x = pupil_data(:,4);
Pupil_left_y = pupil_data(:,5);
Pupil_down_x = pupil_data(:,7);
Pupil_down_y = pupil_data(:,8);
Pupil_right_x = pupil_data(:,10);
Pupil_right_y = pupil_data(:,11);

center_x = zeros(size(Pupil_up_x,1),1);
center_y = zeros(size(Pupil_up_x,1),1);
radii = zeros(size(Pupil_up_x,1),1);
areas = zeros(size(Pupil_up_x,1),1);

for i = 1:size(pupil_data,1)
    X1(1) = Pupil_up_x(i,1);
    X1(2) = Pupil_up_y(i,1);

    X2(1) = Pupil_left_x(i,1);
    X2(2) = Pupil_left_y(i,1);

    Y1(1) = Pupil_down_x(i,1);
    Y1(2) = Pupil_down_y(i,1);

    Y2(1) = Pupil_right_x(i,1);
    Y2(2) = Pupil_right_y(i,1);

    MinorAxis = [X1(1), X1(2); Y1(1), Y1(2)];
    MinorAxisLength = pdist(MinorAxis, 'euclidean');

    MajorAxis = [X2(1), X2(2); Y2(1), Y2(2)];
    MajorAxisLength = pdist(MajorAxis, 'euclidean');

    diameters = mean([MajorAxisLength, MinorAxisLength]);
    radii(i,1) = diameters/2;
    areas(i,1) = pi*radii(i,1).^2;
end

save([Data_Folder 'Pupil.mat'], 'radii', 'areas');

pupil = filloutliers(areas,"nearest","percentiles",[0 100]);

pupil_smooth = movmean(pupil,[smooth_window*FrameRate 0]);

%% Raw ECG signals

ECG_raw = datas(vid_start:end,2)';

fpass=[9 14];

%% Plot spectrogram

window_length = 2 * Sample_Rate;  % Window length in samples (2 seconds)

% Calculate spectrogram
[S, F, T] = spectrogram(ECG_raw, hamming(window_length), round(window_length/2), [], Sample_Rate, 'yaxis');

% Limit the frequency range to 30 Hz
max_freq = 30;  % Maximum frequency to display (Hz)
freq_indices = F <= max_freq;  % Indices of frequencies <= 30 Hz
S_limited = S(freq_indices, :);  % Limit spectrogram to these frequencies
F_limited = F(freq_indices);  % Limit frequency vector

% Plot spectrogram
figure('Renderer', 'painters', 'Position', [600 50 1500 800]);

subplot(2,1,1)

imagesc(T/60, F_limited, 10*log10(abs(S_limited)));  % Convert time to minutes
axis xy;  % Flip y-axis to have low frequencies at the bottom
colormap('jet');
colorbar;

% Adjust color axis for better visibility (you may need to tweak these values)
caxis([0 30]);

% Add colorbar label
c = colorbar;
ylabel(c, 'Power/Frequency (dB/Hz)');

hold on;
h1 = yline(fpass(1), '--r', 'LineWidth', 2, 'Label', sprintf('x = %d', fpass(1)));
h2 = yline(fpass(2), '--b', 'LineWidth', 2, 'Label', sprintf('x = %d', fpass(2)));
hold off;

% Set labels and title
title('ECG Spectrogram');
xlabel('Time (minutes)');
ylabel('Frequency (Hz)');
grid on;

subplot(2,1,2)

plot(ECG_raw,'-k')

%% User prompt loop
while true
    % Create a dialog box for user input
    prompt = {'Enter new value for y-line 1:', 'Enter new value for y-line 2:'};
    dlgtitle = 'Adjust X-Lines';
    dims = [1 35];
    definput = {num2str(fpass(1)), num2str(fpass(2))};
    %answer = inputdlg(prompt, dlgtitle, dims, definput);
    answer = inputdlg(prompt,dlgtitle,1,{num2str(fpass(1)),num2str(fpass(2))},struct('WindowStyle','normal'));
    % Check if the user pressed cancel
    if isempty(answer)
        break;
    end
    
    % Parse user input
    newValues = str2double(answer);
    
    if length(newValues) == 2 && all(~isnan(newValues))
        % Update y-line positions
        fpass(1) = newValues(1);
        fpass(2) = newValues(2);
        
        % Update plot
        delete(h1);
        delete(h2);
        hold on;

        subplot(2,1,1)
        h1 = yline(fpass(1), '--r', 'LineWidth', 2, 'Label', sprintf('x = %d', fpass(1)));
        h2 = yline(fpass(2), '--b', 'LineWidth', 2, 'Label', sprintf('x = %d', fpass(2)));
        hold off;
    else
        % Display an error dialog if input is invalid
        errordlg('Invalid input. Please enter two numbers.', 'Input Error');
    end
    
    % Ask the user if they are satisfied with the y-line values
    
prompt = {'Are new fpass values correct? (0 = No, 1 = Yes)'};
    dlgtitle = 'Adjust X-Lines';
    dims = [1 35];
    definput = {'1'};
    answer = inputdlg(prompt,dlgtitle,1,{'0'},struct('WindowStyle','normal'));
    
    % Handle response
    if contains(answer,'1')
        break;
    end
end

%% Remove baseline wandering

ECG_Bandpass = bandpass(ECG_raw,fpass,Sample_Rate);

%% find peaks

minPeakPromVal=0.007;

% figure;
% findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

[pksVal, pksLocs]=findpeaks(ECG_Bandpass,Sample_Rate,'MaxPeakWidth',0.15,'MinPeakProminence',minPeakPromVal);

%% heart rate in time

RR_intervals = diff(pksLocs);
heartRate=1./RR_intervals;
heartRate_bpm=heartRate*60;

heartRate_bpm_outlier = filloutliers(heartRate_bpm,"nearest","percentiles",[0 100]);

heartRate_bpm_interp = interp1(pksLocs(2:end)*Sample_Rate,heartRate_bpm_outlier,0:time(end)*Sample_Rate,'nearest','extrap');

heartRate_bpm_interp_smooth = movmean(heartRate_bpm_interp,[5*Sample_Rate 0]);

%% HRV

for i = 1 : 5 : total_time
    heartbeat_index = find(pksLocs(2:end)>=(i-1) & pksLocs(2:end)<i+4);
    RR_intervals_trace = RR_intervals(heartbeat_index);
    RR_intervals_std(round(i/5)+1) = std(RR_intervals_trace);
end

SSD = abs(diff(RR_intervals));
SSD_smooth = movmean(SSD,50);

LF = bandpass(heartRate_bpm_interp,[0.4 0.8],Sample_Rate);

%% Plot running, heartRate, HRV and pupil size

video_time = (1/FrameRate):(1/FrameRate):(size(pupil,1)/FrameRate);

figure
subplot(6,1,1);
plot(time,[speed speed(end)],'k')
xlim([0 time(end)])
% ylim([0 0.5])
title('Running Speed','FontSize',15,'FontWeight','bold')
% ylabel('Speed','FontSize',15,'FontWeight','bold')
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(6,1,2)
plot(video_time, pupil_smooth,'color',[229 114 190]./255,'LineWidth',2)
xlim([0 time(end)])
% ylim([100 600])
% ylabel('pupil area (pixels)','FontSize',15,'FontWeight','bold')
title('Pupil Size','FontSize',15,'FontWeight','bold','color',[229 114 190]./255)
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(6,1,3)
plot(time,heartRate_bpm_interp_smooth,'LineWidth',2)
% plot(1:total_time, heartRate_median_bpm_smooth, 'LineWidth',2)
xlim([0 time(end)])
% ylim([600 780])
title('Heart Rate','FontSize',15,'FontWeight','bold','color',[34 75 160]./255)
% axis off
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(6,1,4)
% plot(2:total_time,SSD_smooth,'color','k','LineWidth',2)
plot(pksLocs(3:end),SSD_smooth,'color','k','LineWidth',2)
title('Root Mean Square of Successive Differences','FontSize',15,'FontWeight','bold','color','k')
xlim([0 time(end)])
% ylim([0 0.02])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(6,1,5)
plot(1:5:total_time,RR_intervals_std,'color','k','LineWidth',2)
% plot(pksLocs(3:end),SSD,'LineWidth',2)
title('Standard Deviation','FontSize',15,'FontWeight','bold','color','k')
xlim([0 time(end)])
% ylim([0 0.03])
ax = gca;
ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

subplot(6,1,6)
% plot(1:total_time,LF,'color','k','LineWidth',2)
plot(time,LF,'color','k','LineWidth',2)
title('LF Signal (0.4-0.8Hz)','FontSize',15,'FontWeight','bold','color','k')
xlim([0 time(end)])
% ylim([-100 100])
ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.FontWeight = 'bold';

xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')




