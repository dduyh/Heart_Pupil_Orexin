function [stereoNoise] = generateWhiteNoise(band,volume,duration)
%% band = [lowerfreq higherfreq]
% --- Parameters ---
fs = 48000;        % Sample rate (Hz)
rampDur = 0.05;    % Ramp duration (s)

fprintf('Generating %.2f s noise (%.0f–%.0f Hz)...\n', duration, band(1), band(2));

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



