%% Facial expressions and locomotion of emotion states and their neuronal correlates in mice

% Author:  Yihui Du
% Date:    December 04th, 2023

%%

clear; close all; clc;

%% set the path for output data

% Directory = 'P:\Yihui\data\';   % Main directory\
Directory = 'C:\Users\yihudu\Desktop\Yihui\data\';   % Main directory\
% mouse_name = 'mouse_LLR';     % Mouse name\
mouse_name = 'm485_L_VTA_R_NAc_oxLight';            % Mouse name\
date = 'Dec_12_2023';           % Date\
session = 'session_2';  

TTL_trigger = 0;

Data_Folder = [Directory mouse_name '\' date '\' session '\'];     % Set data folder name
if ~exist(Data_Folder,'dir')                           % Create data folder if it doesn't exist
    mkdir(Data_Folder)
end

%% Initialize the Arduino

delete(instrfindall);                % Delete all existing serial port objects
a = serial('COM4');                  % Connect to serial port (COM4)
set(a,'BaudRate',9600);              % Set baud rate (communication speed in bits per second)
set(a,'Timeout',30);                 % Set time out (allowed time in seconds to complete read and write)
set(a,'InputBufferSize',8388608);    % Set input buffer size

fopen(a);                            % Open connection to Arduino board
if (exist('board1','var'))           % Stop the running program in Arduino
    board1.stop;pause(0);
end

%% Open video input (face)

imaqreset                                              % Disconnect and delete all image acquisition objects
imaqmex('feature','-limitPhysicalMemoryUsage',false);  % Set unlimited physical memory usage

vid = videoinput('gentl', 1, 'Mono8');
vid.FramesPerTrigger = Inf;
vid.ReturnedColorspace = 'grayscale';
vid.ROIPosition = [210 226 1017 669];

vid.LoggingMode = 'disk';
diskLogger = VideoWriter([Data_Folder 'Sucrose_Quinine.avi'], 'Grayscale AVI');
vid.DiskLogger = diskLogger;
diskLogger.FrameRate = 20;

src = getselectedsource(vid);
src.AcquisitionFrameRate = 20;
src.AutoExposureExposureTimeLowerLimit = 1000;
src.AutoExposureExposureTimeUpperLimit = 1000;
src.AutoExposureGainLowerLimit = 0;
src.AutoExposureGainUpperLimit = 0.003;
src.AutoExposureGreyValueLowerLimit = 30;
src.AutoExposureGreyValueUpperLimit = 30.01;

%% Open video input (body)



%% Create a figure window to show traces

figure(1);

ax1 = subplot(4,1,1);
title('Licking')
ylabel('Capacity');
hold(ax1,'off')

ax2 = subplot(4,1,2);
title('Running')
ylabel('Speed (deg.)');
hold(ax2,'off')

ax3 = subplot(4,1,3);
title('Photometry 1')
hold(ax3,'off')

ax4 = subplot(4,1,4);
title('Photometry 2')
xlabel('Time (s)');
hold(ax4,'off')

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

% Control the licking detectors
% ch1 = addAnalogInputChannel(s,'Dev1',0,'Voltage');
ch1 = addDigitalChannel(s,'Dev1','port0/line0','InputOnly');

% Control the running wheel rotory encoder
% add a counter input channel for |Position| measurement type.
% (A, B, and Z signal outputs to PFI8, PFI10, and PFI9).
ch2 = addCounterInputChannel(s, 'Dev1', 0, 'Position');

% Configure quadrature cycle encoding type (X1, X2, or X4).
ch2.EncoderType = 'X1';

% Control the Doric photometry setup
if TTL_trigger
    
    % Add counter output channel 1, 2 with |PulseGeneration| measurement type.
    addCounterOutputChannel(s,'Dev1', 1, 'PulseGeneration');
    addCounterOutputChannel(s,'Dev1', 2, 'PulseGeneration');
    
    % Clocked Counter Output
    ch3 = s.Channels(3);
    ch3.Frequency = 10;
    ch3.InitialDelay = 0;
    ch3.DutyCycle = 0.5;
    
    ch4 = s.Channels(4);
    ch4.Frequency = 10;
    ch4.InitialDelay = 0.05;
    ch4.DutyCycle = 0.5;
    
end

ch5 = addAnalogInputChannel(s,'Dev1',3,'Voltage');

ch6 = addAnalogInputChannel(s,'Dev1',19,'Voltage');
% ch6.TerminalConfig = 'SingleEnded'

% Configure a new listener to process the incoming data.
lh = addlistener(s,'DataAvailable', @(src, event)continuousPlot(src, event, ax1, ax2, ax3, ax4));

% Control the background LEDs

%% Start video preview

preview(vid);

%% Start signals generation and acquisition

s.startBackground()
tStart = tic; % Record the exact time for starting signals generation and acquisition

%% Start video recording

fprintf('Starting video...\n')
start(vid);
step_timepoint = [step_timepoint toc(tStart)]; % Record the exact time when video collection was started

fprintf('Video started. ')


%% Control the valve drivers

t = timer;
t.StopFcn = @(~,~)disp('valve timer stopped');
t.StartDelay = 50;
t.TimerFcn = {@valveTrigger, a, tStart, ax1};
t.Period = 50;
t.TasksToExecute = 6;
t.ExecutionMode = 'fixedSpacing';

start(t)
step_timepoint = [step_timepoint toc(tStart)];

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

