%%

clear; close all; clc;

%% Create a figure window to show traces

figure(1);

ax1 = subplot(2,1,1);
title('Running')
ylabel('Angular position (deg.)');
hold(ax1,'off')

ax2 = subplot(2,1,2);
xlabel('Time (s)');
ylabel('Speed (deg.)');
hold(ax2,'off')

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
ch1.EncoderType = 'X4';

ch2 = addAnalogInputChannel(s,'Dev1',3,'Voltage');



% Configure a new listener to process the incoming data.
lh = addlistener(s,'DataAvailable', @(src, event)continuousPlot(src, event, ax1, ax2));

%% Start signals generation and acquisition

s.startBackground()

%%
s.stop();  % Stop the Continuous Background Acquisition
fprintf('Acquisition finished successfully.\n')

delete(lh) % Delete the listener.

%%
function continuousPlot(src, event, ax1, ax2)

positionDataDeg = event.Data(:,1);

counterNBits = 32;
signedThreshold = 2^(counterNBits-1);
positionDataDeg(positionDataDeg > signedThreshold) = positionDataDeg(positionDataDeg > signedThreshold) - 2^counterNBits;

plot(ax1, event.TimeStamps, positionDataDeg, 'b');
title('Running')
ylabel('Angular position (deg.)');

plot(ax1, event.TimeStamps, positionDataDeg, 'b');
xlabel('Time (s)');
ylabel('Speed (deg.)');

speedDeg = diff(positionDataDeg);
plot(ax2, event.TimeStamps, [speedDeg; speedDeg(end)], 'b');
ax2.YLim = [-4 4];



end