


% Parameters
fs = 1000;  % Sampling frequency (Hz)
duration = 30 * 60;  % Duration in seconds (30 minutes)
window_length = 2 * fs;  % Window length in samples (2 seconds)

% Generate or load your ECG data here
% For this example, we'll generate random data
t = 0:1/fs:duration-1/fs;

% Calculate spectrogram
[S, F, T] = spectrogram(ECG_raw, hamming(window_length), round(window_length/2), [], fs, 'yaxis');

% Limit the frequency range to 30 Hz
max_freq = 30;  % Maximum frequency to display (Hz)
freq_indices = F <= max_freq;  % Indices of frequencies <= 30 Hz
S_limited = S(freq_indices, :);  % Limit spectrogram to these frequencies
F_limited = F(freq_indices);  % Limit frequency vector

%% Plot spectrogram
figure;

subplot(2,1,1)

imagesc(T/60, F_limited, 10*log10(abs(S_limited)));  % Convert time to minutes
axis xy;  % Flip y-axis to have low frequencies at the bottom
colormap('jet');
colorbar;

% Set labels and title
xlabel('Time (minutes)');
ylabel('Frequency (Hz)');
title('ECG Spectrogram (up to 30 Hz)');

% Adjust color axis for better visibility (you may need to tweak these values)
caxis([0 30]);

% Add colorbar label
c = colorbar;
ylabel(c, 'Power/Frequency (dB/Hz)');

yline(fpass(1),'-r')
yline(fpass(2),'-r')

subplot(2,1,2)

plot(ECG_raw,'-k')