clc
clear
close all

%%
load('P:\Yihui\m1271\19092023m127 file1_superRigMainMat.mat')

dataLength=size(superRigMainMat,1);
numSubplots=size(superRigMainMat,2)+1;

%%

fileName = 'P:\Yihui\m1271\2023_09_07_m1271_Quinine_s3-09072023190146-0000.avi';
obj = VideoReader(fileName);
numFrames = obj.NumFrames;
FrameRate = obj.FrameRate;

IRsynch = zeros(numFrames,1);

figure(1);                                            
frame = read(obj,1);
IR=roipoly(frame);   

for i = 1:numFrames
    
    fprintf('running recording = %i\n',i)
    frame = read(obj,i);
    IR_region = frame.*uint8(IR);   
    IRsynch(i) = sum(IR_region,'all');
    
end

%%

figure(2)
titles = {'Licking', 'Running', 'StimType', 'Blank', 'Sucrose Valve', 'H2O Valve', 'Qunine Valve', 'Blank', 'IR synch'};

%%

for subplotLoop = 1:numSubplots-1
    subplot(numSubplots,1,subplotLoop)
    plot(1:dataLength,superRigMainMat(:,subplotLoop))
    xlim([0 dataLength])
    title(titles{subplotLoop})
end

%%
emotions = zeros(dataLength,1);
emotions(7311:20503) = 1;
emotions(28812:34452) = -1;
emotions(34452:57914) = 1;
emotions(65217:73542) = -1;
emotions(73542:79982) = 1;
emotions(96153:103344) = -1;
emotions(103344:119141) = 1;
emotions(127185:135799) = -1;
emotions(135799:150524) = 1;
emotions(159169:dataLength) = -1;

subplot(numSubplots,1,numSubplots)
map = [0 0 1
    0.8 0.8 0.8
    1 0 0];
imagesc(emotions')
colormap(map)
title('Emotion States')

%%

figure(3)
subplot(2,1,1)
% plot(1:dataLength,superRigMainMat(:,9))
plot(6254:157872,superRigMainMat(6254:157872,9))
xlim([6254 157872])

subplot(2,1,2)
% plot(1:numFrames,IRsynch)
plot(764:16496,IRsynch(764:16496))
xlim([764 16496])

%%

figure(4)
numSubplots = 7;

xlimits = [6254 157872];

licking = superRigMainMat(:,1);
h1 = subplot(numSubplots,1,1);
plot(xlimits(1):xlimits(2),licking(xlimits(1):xlimits(2)))
xlim(xlimits)
title('Licking')

% lick = zeros(dataLength,1);
% lick(licking>0.5) = 1;
% h2 = subplot(numSubplots,1,2);
% map = [0.8 0.8 0.8
%     1 1 0];
% imagesc(lick')
% colormap(h2,map)

sucrose = superRigMainMat(:,5);
h2 = subplot(numSubplots,1,2);
plot(xlimits(1):xlimits(2),sucrose(xlimits(1):xlimits(2)))
xlim(xlimits)
title('Sucrose')

water = superRigMainMat(:,6);
h3 = subplot(numSubplots,1,3);
plot(xlimits(1):xlimits(2),water(xlimits(1):xlimits(2)))
xlim(xlimits)
title('H2O')

qunine = superRigMainMat(:,7);
h4 = subplot(numSubplots,1,4);
plot(xlimits(1):xlimits(2),qunine(xlimits(1):xlimits(2)))
xlim(xlimits)
title('Qunine')

h5 = subplot(numSubplots,1,5);
map = [0 0 1
    0.8 0.8 0.8
    1 0 0];
imagesc(emotions(xlimits(1):xlimits(2))')
colormap(h5,map)
title('Emotion States')

running = superRigMainMat(:,2);
h6 = subplot(numSubplots,1,6);
plot(xlimits(1):xlimits(2),running(xlimits(1):xlimits(2)))
xlim(xlimits)
title('Running')

data = csvread('P:\Yihui\m1271\Pupil Labeling\Filtered_Sync.csv',1,0);
pupil_size = data(:,7);
h7 = subplot(numSubplots,1,7);
plot(1:size(pupil_size,1),pupil_size)
xlim([1 size(pupil_size,1)])
title('Pupil Size')

%%

emotion_offset = [7311 20503 28812 34452 57914 65217 73542 79982 96153 103344 119141 127185 135799 150524 157872]-6254;
emotion_onset = [0 emotion_offset(1:end-1)];
emotion_duration = (emotion_offset - emotion_onset)./151619.*5245;

% emotion = emotions(xlimits(1):xlimits(2));
% find(emotion==0);
% for 1 = 1:length(emotion)
% 
% end





