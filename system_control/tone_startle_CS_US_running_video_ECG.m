%%

clear; close all; clc;

%% set the path for output data

Directory = 'C:\Users\yihudu\Desktop\Yihui\data\';   % Main directory\
mouse_name = 'm833';            % Mouse name\
date = 'Oct_29_2025';           % Date\

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name
if ~exist(Data_Folder,'dir')                           % Create data folder if it doesn't exist
    mkdir(Data_Folder)
end

%%

PinkNoise = generateWhiteNoise([6000 8000],0.015,1); % volume = 90 dB

WhiteNoise = generateWhiteNoise([1 10000],1,1); % volume = 120 dB

%% Open video input (face)

imaqreset                                              % Disconnect and delete all image acquisition objects
imaqmex('feature','-limitPhysicalMemoryUsage',false);  % Set unlimited physical memory usage

vid = videoinput('gentl', 1, 'Mono8');
vid.FramesPerTrigger = Inf;
vid.ReturnedColorspace = 'grayscale';
% vid.ROIPosition = [600 338 400 400];

vid.LoggingMode = 'disk';
diskLogger = VideoWriter([Data_Folder mouse_name '_' date '.mp4'], 'MPEG-4');
vid.DiskLogger = diskLogger;
diskLogger.FrameRate = 20;

src = getselectedsource(vid);
src.AcquisitionFrameRate = 20;

% Turn off auto exposure and gain, and set exposure time to 45000
src.ExposureAuto = 'off';                             % Disable auto exposure
src.ExposureTime = 45000;                             % Set fixed exposure time to 40000

src.GainAuto = 'off';                                 % Turn off automatic gain control
src.Gain = 0;                                        % Disable gain by setting it to 0

%% Create a figure window to show traces

figure(1);

ax1 = subplot(3,1,1);
title('Running')
ylabel('Speed (deg.)');
hold(ax1,'off')

ax2 = subplot(3,1,2);
title('ECG')
hold(ax2,'off')

ax3 = subplot(3,1,3);
title('Camera')
hold(ax3,'off')

%% Save the data

global times datas step_timepoint

times = [];
datas = [];
step_timepoint = [];               % Record the exact time when each operation step was performed

%% Initialize the NI-DAQ USB-6343

daqreset

s = daq.createSession('ni'); % Create a session.

s.Rate = 1000; % Sessions run for one second at 1000 scans per second by default.

s.IsNotifyWhenDataAvailableExceedsAuto = true; % Listener is called ten times per second by default.

% s.NotifyWhenDataAvailableExceeds = 2000;

s.IsContinuous = true; % Configure the session to acquire continuously.

%% Signals Acquisition

% Control the running wheel rotory encoder
% add a counter input channel for |Position| measurement type.
% (A, B, and Z signal outputs to PFI8, PFI10, and PFI9).
ch1 = addCounterInputChannel(s, 'Dev1', 0, 'Position');

% Configure quadrature cycle encoding type (X1, X2, or X4).
ch1.EncoderType = 'X1';

ch2 = addAnalogInputChannel(s,'Dev1',3,'Voltage');

ch3 = addAnalogInputChannel(s,'Dev1',0,'Voltage');

% Configure a new listener to process the incoming data.
lh = addlistener(s,'DataAvailable', @(src, event)RunningPlot(src, event, ax1, ax2, ax3));

%% Start video preview

preview(vid);

%% Start signals generation and acquisition

s.startBackground()
tStart = tic; % Record the exact time for starting signals generation and acquisition

%% Start video recording

fprintf('Starting video...\n')
start(vid);
step_timepoint = [step_timepoint toc(tStart)]; % Record the exact time when video collection was started

fprintf('Video started. \n')

%% Control the speaker

t = timer;
t.StopFcn = @(~,~) disp('Session ended. \n'); 
t.StartDelay = 5; % 2400 secs = 40 mins
t.TimerFcn = {@speakerPlayer, PinkNoise, WhiteNoise, tStart};
% t.Period = 480;
t.TasksToExecute = 1;
t.ExecutionMode = 'fixedRate';

start(t)
step_timepoint = [step_timepoint toc(tStart)];
disp('Session started. \n'); 

%%

delete(t)  % Delete the timer

s.stop();  % Stop the Continuous Background Acquisition
fprintf('Acquisition finished successfully.\n')

delete(lh) % Delete the listener.

stoppreview(vid);                                   % Stop face video preview
stop(vid);                                          % Stop face video collection
tEnd = toc(tStart);                                 % Get the exact ending time
disp(['Elapsed time: ',num2str(tEnd),' seconds.']); % Display the exact ending time
step_timepoint = [step_timepoint tEnd];             % Record the exact ending time

close(diskLogger);                                  % Close output video object

% save data files
save([Data_Folder 'times.mat'],'times');
save([Data_Folder 'datas.mat'],'datas');
save([Data_Folder 'step_timepoint.mat'],'step_timepoint'); % Save the exact time for each operation step

%%

function RunningPlot(~, event, ax1, ax2, ax3)

global times datas

times = [times; event.TimeStamps];
datas = [datas; event.Data];

xlims = [max([0 ceil(event.TimeStamps(end))-10]) ceil(event.TimeStamps(end))];

speedDeg = diff(datas(max([0 end-10*1000])+1:end,1));
plot(ax1, times(max([0 end-10*1000])+1:end,1), [speedDeg; speedDeg(end)], 'b');
ax1.XLim = xlims;
ax1.YLim = [-2 2];

plot(ax2, times(max([0 end-10*1000])+1:end,1), datas(max([0 end-10*1000])+1:end,2), 'b')
ax2.XLim = xlims;

plot(ax3, times(max([0 end-10*1000])+1:end,1), datas(max([0 end-10*1000])+1:end,3), 'b')
ax3.XLim = xlims;

end

%%

function speakerPlayer(~, ~, PinkNoise, WhiteNoise, tStart)

global step_timepoint

step_timepoint = [step_timepoint toc(tStart)];
player = audioplayer(PinkNoise, 48000, 24);
play(player); % play tone through sound card (1 secs)
pause(5) % wait 300 secs

step_timepoint = [step_timepoint toc(tStart)];
player = audioplayer(WhiteNoise, 48000, 24);
play(player); % play tone through sound card (1 secs)
pause(5) % wait 300 secs

step_timepoint = [step_timepoint toc(tStart)];
player = audioplayer(WhiteNoise, 48000, 24);
play(player); % play tone through sound card (1 secs)
pause(5) % wait 300 secs

step_timepoint = [step_timepoint toc(tStart)];
player = audioplayer(WhiteNoise, 48000, 24);
play(player); % play tone through sound card (1 secs)
pause(5) % wait 300 secs

end

function [stereoNoise] = generateWhiteNoise(band,volume,duration)
%% band = [lowerfreq higherfreq]
% --- Parameters ---
fs = 48000;        % Sample rate (Hz)
rampDur = 0.05;    % Ramp duration (s)

fprintf('Generating %.2f s noise (%.0f-%.0f Hz)...\n', duration, band(1), band(2));

% --- Generate band-limited white noise ---
nSamples = round(duration * fs);
noise = randn(1, nSamples);  % White noise
[b, a] = butter(3, band / (fs/2), 'bandpass');
noise = filtfilt(b, a, noise);

% --- Normalize amplitude ---
noise = noise ./ max(abs(noise));

% --- Apply cosine onset/offset ramps ---
rampSamples = round(rampDur * fs);
ramp = 0.5 - 0.5 * cos(pi * (0:rampSamples-1) / rampSamples);  % half-cosine ramp
env = [ramp, ones(1, nSamples - 2*rampSamples), fliplr(ramp)];
noise = noise .* env;   % <-- Apply ramp envelope

% --- RMS limiter (protects amp/speakers) ---
targetRMS = 0.8;
currentRMS = rms(noise);
noise = noise * (targetRMS / currentRMS);

% --- Apply overall volume scaling ---
noise = noise * volume;

% --- Stereo (duplicate mono signal to both channels) ---
stereoNoise = [noise; noise];   % 2 x nSamples

end
