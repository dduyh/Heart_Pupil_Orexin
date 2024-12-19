%% Acquire Continuous and Background Data Using NI Devices
%
% This example shows how to acquire analog input data using non-blocking 
% commands. This allows you to continue working in the MATLAB command window 
% during the acquisition. This is called *background acquisition*.
% Use *foreground acquisition* to cause MATLAB to wait for the
% entire acquisition to complete before you can execute your next command.
% Learn more about foreground acquisition in 
% <docid:daq_examples.example-ex39516307 Acquire Data Using NI devices>.
% 
%   Copyright 2010-2014 The MathWorks, Inc.

%% Create and Configure the Session Object
% * Use the |daq.createSession| command to create a session. 
% * Add channel 0 from our NI 9205 with the ID of cDAQ1Mod1.
% * Acquire 2000 scans per second for 5 seconds worth of data.

s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1', 0, 'Voltage');
s.Rate = 2000;
% s.DurationInSeconds = 5;
s

%% Plot Live Data as it is Acquired
% When you acquire data in the background, you can provide the acquisition
% session directions to handle the incoming data, using listeners and
% events. A DataAvailable event occurs when a specific amount of data is
% available to the session. A listener can respond to that event and
% initiate specified function.
%
% Use |addlistener| to add an anonymous function to the session. This
% function is called every time the DataAvailable event occurs, and plots
% the acquired data against time. By default this listener is called 10
% times per second.
%
% When you add a listener, a handle to the listener is returned. Save the
% handle in the variable |lh| and delete it later.

lh = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps, event.Data));

%%
% The default rate that the listener is called is 10 times a second.
% However, you can change this my modifying the
% |NotifyWhenDataAvailableExceeds| property.  The listener will be called
% when the number of points accumulated exceeds this value.  In this case,
% we set it to 2000, which results in a frequency of once per second.
s.NotifyWhenDataAvailableExceeds = 2000;

%%
% After you add the listener, start the background acquisition.  When you
% run this example in MATLAB a plot updates 10 times a second. 

s.IsContinuous = true;
s.startBackground();

%%
% Use |pause| in a loop to monitor the number of scans acquired for the
% duration of the acquisition. Note that the output is mixed in with output
% from the event listener.
while s.IsRunning
    pause(0.5)
    fprintf('While loop: Scans acquired = %d\n', s.ScansAcquired)
end

fprintf('Acquisition has terminated with %d scans acquired\n', s.ScansAcquired);

%%
% Delete the listener.
delete(lh)
