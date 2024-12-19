%%

clear; close all; clc;

%% set the path for output data

Directory = 'C:\Users\yihudu\Desktop\Yihui\data\';   % Main directory\
mouse_name = 'm1750';            % Mouse name\
date = 'May_25_2024';           % Date\

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name
if ~exist(Data_Folder,'dir')                           % Create data folder if it doesn't exist
    mkdir(Data_Folder)
end

%%

A = 2; % Amplitude
f_0 = 2900; % Frequency of sound (2.9k Hz)
fs = 29000;   % Sampling frequency (40k Hz)
N = 29000*120;    % Signal sampling points number, Playback duration (120 secs)
CS = A*sin(2*pi*f_0*(0:N-1)/fs); % Single frequency sine signal

Fs = 20000;  % sampling frequency in Hz
Ts = 1/Fs; % sampling time in seconds
te = 120; % signal duration in seconds
US = randn(te*Fs,1); % generate White noise of duration te

%% Open video input (face)

imaqreset                                              % Disconnect and delete all image acquisition objects
imaqmex('feature','-limitPhysicalMemoryUsage',false);  % Set unlimited physical memory usage

vid = videoinput('gentl', 1, 'Mono8');
vid.FramesPerTrigger = Inf;
vid.ReturnedColorspace = 'grayscale';
vid.ROIPosition = [151 291 1143 726];

vid.LoggingMode = 'disk';
diskLogger = VideoWriter([Data_Folder 'Fear_Conditioning.avi'], 'Grayscale AVI');
vid.DiskLogger = diskLogger;
diskLogger.FrameRate = 20;

src = getselectedsource(vid);
src.AcquisitionFrameRate = 20;
src.AutoExposureExposureTimeLowerLimit = 1000;%1000
src.AutoExposureExposureTimeUpperLimit = 1000;
src.AutoExposureGainLowerLimit = 0;
src.AutoExposureGainUpperLimit = 0.003;
src.AutoExposureGreyValueLowerLimit = 30;
src.AutoExposureGreyValueUpperLimit = 30.01;

%% Create a figure window to show traces

figure(1);

ax1 = subplot(3,1,1);
title('Running')
ylabel('Running Speed (deg.)');
hold(ax1,'off')

ax2 = subplot(3,1,2);
title('ECG')
hold(ax2,'off')

ax3 = subplot(3,1,3);
title('Injection')
xlabel('Time (s)');
hold(ax3,'off')

%% Save the data

global times datas step_timepoint

times = [];
datas = [];
step_timepoint = [];               % Record the exact time when each operation step was performed

%% Initialize the NI-DAQ USB-6343

daqreset

s = daq.createSession('ni'); % Create a session.

s.Rate = 200; % Sessions run for one second at 1000 scans per second by default.

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

ch3 = addDigitalChannel(s,'Dev1','port0/line0','InputOnly');

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
t.StopFcn = {@speakerStop, tStart}; 
t.StartDelay = 240;
t.TimerFcn = {@speakerPlayer, CS, US, fs, Fs, tStart};
% t.Period = 210;
t.TasksToExecute = 1;
t.ExecutionMode = 'fixedRate';

start(t)
step_timepoint = [step_timepoint toc(tStart)];
disp('Session started \n'); 

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

function speakerPlayer(~, ~, CS, US, fs, Fs, tStart)

global step_timepoint

step_timepoint = [step_timepoint toc(tStart)];

sound(US,Fs);  % play tone through sound card (30 secs)

pause(240)

step_timepoint = [step_timepoint toc(tStart)];

sound(CS,fs);  % play tone through sound card (30 secs)

end

function speakerStop(~, ~, tStart)

global step_timepoint

pause(780)    % wait 780 secs

step_timepoint = [step_timepoint toc(tStart)];

disp('Session ended \n'); 

end
