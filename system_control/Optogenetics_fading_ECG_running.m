%%

clear; close all; clc;

%% set the path for output data

Directory = 'C:\Users\yihudu\Desktop\Yihui\data\';   % Main directory\
mouse_name = 'm2179';            % Mouse name\
date = 'Mar_26_2025';           % Date\

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
vid.ROIPosition = [192 338 998 520];

vid.LoggingMode = 'disk';
diskLogger = VideoWriter([Data_Folder 'Fear_Conditioning.avi'], 'Grayscale AVI');
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
title('Pulse')
xlabel('Time (s)');
hold(ax3,'off')

%% Save the data

global datas

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

% Add counter output channel 1, 2 with |PulseGeneration| measurement type.
ch3 = addDigitalChannel(s,'Dev1','port0/line3', 'OutputOnly');

% Pulses detector
ch4 = addDigitalChannel(s,'Dev1','port0/line2','InputOnly');

% Configure a new listener to process the incoming data.
lh = addlistener(s,'DataAvailable', @(src, event)RunningECGPulsePlot(src, event, ax1, ax2, ax3));

%% Start video preview

preview(vid);

%% Start video recording

fprintf('Starting video...\n')
start(vid);
tStart = tic; % Record the exact time when video collection was started
fprintf('Video started. \n')

%% Generate TTL Stim

stimLength = 20;

Fs = 1000;

high_duration = 5e-3;

% settings for trial timeline
interTrialInterval = 420; % 500

% setup stim schedule
ntrials = 5;

ttl_signal = [];

for sec = 1:stimLength

    f = 21 - sec;
    T = 1 / f;
    low_duration = T - high_duration;

    ttl_1s = [];
    t_elapsed = 0;
    while t_elapsed + high_duration < 1
        ttl_1s = [ttl_1s; ones(round(high_duration * Fs), 1); zeros(round(low_duration * Fs), 1)];
        t_elapsed = t_elapsed + high_duration + low_duration;
    end

    if length(ttl_1s) < Fs
        ttl_1s = [ttl_1s; zeros(Fs - length(ttl_1s), 1)];
    elseif length(ttl_1s) > Fs
        ttl_1s = ttl_1s(1:Fs);
    end

    ttl_signal = [ttl_signal; ttl_1s];
end

zero_signal = zeros(Fs * interTrialInterval, 1);

final_signal = [zero_signal; repmat([ttl_signal; zero_signal], ntrials, 1)];

% Control the laser

fprintf('nTrials %d , Stim Duration %d secs + Rest %d secs, Total %d secs\n', ntrials, stimLength, interTrialInterval, length(final_signal)/Fs);

queueOutputData(s, final_signal);

s.startBackground();

step_timepoint = [step_timepoint toc(tStart)]; % Record the exact time for starting signals generation and acquisition

pause(interTrialInterval+ntrials*(stimLength+interTrialInterval))

s.stop();

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

%%

function RunningECGPulsePlot(~, event, ax1, ax2, ax3)

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

end