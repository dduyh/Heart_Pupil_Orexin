%%

clear; close all; clc;

%% set the path for output data

Directory = 'C:\Users\labadmin\Desktop\Yihui\data\';   % Main directory\
mouse_name = 'm841';            % Mouse name\
date = 'Dec_05_2025';           % Date\

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name
if ~exist(Data_Folder,'dir')                           % Create data folder if it doesn't exist
    mkdir(Data_Folder)
end

%% Open video input (face)

imaqreset                                              % Disconnect and delete all image acquisition objects
imaqmex('feature','-limitPhysicalMemoryUsage',false);  % Set unlimited physical memory usage

vid = videoinput('gentl', 1, 'Mono8');
vid.FramesPerTrigger = Inf;
vid.ReturnedColorspace = 'grayscale';
% vid.ROIPosition = [438 307 450 450];

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

ax1 = subplot(5,1,1);
title('Running')
ylabel('Running Speed (deg.)');
hold(ax1,'off')

ax2 = subplot(5,1,2);
title('ECG')
hold(ax2,'off')

ax3 = subplot(5,1,3);
title('Photometry')
hold(ax3,'off')

ax4 = subplot(5,1,4);
title('Pulse')
hold(ax4,'off')

ax5 = subplot(5,1,5);
title('Injection')
xlabel('Time (s)');
hold(ax5,'off')

%% Save the data

global datas

datas = [];

step_timepoint = [];               % Record the exact time when each operation step was performed

%% Initialize the NI-DAQ USB-6343

daqreset

s1 = daq.createSession('ni'); % Create a session.

s1.Rate = 1000; % Sessions run for one second at 1000 scans per second by default.

s1.IsNotifyWhenDataAvailableExceedsAuto = true; % Listener is called ten times per second by default.

% s.NotifyWhenDataAvailableExceeds = 2000;

s1.IsContinuous = true; % Configure the session to acquire continuously.

s2 = daq.createSession('ni'); % Create a session.

s2.Rate = 1000; % Sessions run for one second at 1000 scans per second by default.

s2.IsNotifyWhenDataAvailableExceedsAuto = true; % Listener is called ten times per second by default.

% s.NotifyWhenDataAvailableExceeds = 2000;

s2.IsContinuous = true; % Configure the session to acquire continuously.

%% Signals Acquisition

% Control the running wheel rotory encoder
% add a counter input channel for |Position| measurement type.
% (A, B, and Z signal outputs to PFI8, PFI10, and PFI9).
ch1 = addCounterInputChannel(s1, 'Dev2', 0, 'Position');

% Configure quadrature cycle encoding type (X1, X2, or X4).
ch1.EncoderType = 'X1';

ch2 = addAnalogInputChannel(s1,'Dev2',3,'Voltage');

% Add counter output channel 1, 2 with |PulseGeneration| measurement type.
ch3 = addCounterOutputChannel(s2,'Dev2','ctr1','PulseGeneration'); % make sure to set this to proper counter output and device

% Photometry
ch4 = addAnalogInputChannel(s1,'Dev2',19,'Voltage');

% Pulses detector
ch5 = addDigitalChannel(s1,'Dev2','port0/line2','InputOnly');

ch6 = addDigitalChannel(s1,'Dev2','port0/line0','InputOnly');

% Configure a new listener to process the incoming data.
lh = addlistener(s1,'DataAvailable', @(src, event)RunningECGPulsePlot(src, event, ax1, ax2, ax3, ax4, ax5));

%% Start video preview

preview(vid);

%% Start video recording

fprintf('Starting video...\n')
start(vid);
tStart = tic; % Record the exact time when video collection was started
fprintf('Video started. \n')

%% Control the laser

% setup stim schedule
stimSchedule = [20 5 5 20];
ntrials = length(stimSchedule);

% settings for trial timeline
stimLength = 60;%30
interTrialInterval = [900 420 420 420];%120

delay = zeros(1,ntrials);

% Start signals generation and acquisition

s1.startBackground
% step_timepoint = [step_timepoint toc(tStart)]; % Record the exact time for starting signals generation and acquisition

for trialInd = 1:ntrials

% select frequency from schedule and calculate duty cycle
ch3.Frequency = stimSchedule(trialInd);
ch3.DutyCycle = 5/(1000/ch3.Frequency);
% ch3.DutyCycle = 0.99;

initial_delay = interTrialInterval(trialInd); % constant
ch3.InitialDelay = initial_delay;
delay(trialInd) = initial_delay;

disp(['Trial number ' num2str(trialInd) ': ' num2str(stimSchedule(trialInd)) ' Hz ' ])

%%%%%%%%%% STIM %%%%%%%%%%%%%
s2.startBackground
step_timepoint = [step_timepoint toc(tStart)]; % Record the exact time for starting signals generation and acquisition

pause(stimLength + initial_delay)

s2.stop
%%%%%%%%%% END STIM %%%%%%%%%%

end

ch3.InitialDelay = 450; 
s2.startBackground
pause(420) % wait for 180
s2.stop
s1.stop

fprintf('Acquisition finished successfully.\n')

%%

delete(lh) % Delete the listener.

stoppreview(vid);                                   % Stop face video preview
stop(vid);                                          % Stop face video collection
tEnd = toc(tStart);                                 % Get the exact ending time
disp(['Elapsed time: ',num2str(tEnd),' seconds.']); % Display the exact ending time
step_timepoint = [step_timepoint tEnd];             % Record the exact ending time

close(diskLogger);                                  % Close output video object

% save data files
save([Data_Folder 'datas.mat'],'datas');
save([Data_Folder 'step_timepoint.mat'],'step_timepoint'); % Save the exact time for each operation step
save([Data_Folder 'delay.mat'],'delay');
save([Data_Folder 'stimSchedule.mat'],'stimSchedule');

%%

function RunningECGPulsePlot(~, event, ax1, ax2, ax3, ax4, ax5)

global datas

datas = [datas; event.Data];

data_length = size(datas,1);

xlims = [max([0 (ceil(data_length/1000)-10)*1000])+1 ceil(data_length/1000)*1000];

speedDeg = diff(datas(max([0 end-10*1000])+1:end,1));
plot(ax1, (max([0 data_length-10*1000])+1):data_length, [speedDeg; speedDeg(end)], 'b');
ax1.XLim = xlims;
ax1.YLim = [-2 2];

plot(ax2, (max([0 data_length-10*1000])+1):data_length, datas(max([0 end-10*1000])+1:end,2), 'b')
ax2.XLim = xlims;

plot(ax3, (max([0 data_length-10*1000])+1):data_length, datas(max([0 end-10*1000])+1:end,3), 'b')
ax3.XLim = xlims;

plot(ax4, (max([0 data_length-10*1000])+1):data_length, datas(max([0 end-10*1000])+1:end,4), 'b')
ax4.XLim = xlims;

plot(ax5, (max([0 data_length-10*1000])+1):data_length, datas(max([0 end-10*1000])+1:end,5), 'b')
ax5.XLim = xlims;

end