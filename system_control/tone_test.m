%%
load gong.mat;
sound(y);

%%
A=2;%振幅
f_0=4000;%声音频率
fs=40000;   %采样频率
N=40000;    % 信�?�样点数，播放时长
y=A*sin(2*pi*f_0*(0:N-1)/fs); %�?�频信�?�
sound(y,fs);  %通过声�?�放音 

%%
Fs = 20000;  % sampling frequency in Hz
Ts = 1/Fs; % sampling time in seconds
te = 10; % signal duration in seconds
y = randn(te*Fs,1); % generate White noise of duration te
sound(y,Fs); % listen to the generated noise if you feel adventurous!