clc
clear
close all

%%
fileName = 'P:\Yihui\m1271\2023_09_07_m1271_Quinine_s3-09072023190146-0000.avi';
obj = VideoReader(fileName);
numFrames = obj.NumFrames;
FrameRate = obj.FrameRate;

%%
for i = 764:3:16496

    fprintf('running recording = %i\n',i)
    
    image_name = strcat('C:\Users\yihudu\Desktop\Yihui\data\m1271\frames_downsample_3\lick (', num2str(i), ').jpg');
    frame = read(obj,i);
    imwrite(frame,image_name,'jpg');
   
end