clc
clear
close all

%%
fileName = 'P:\Yihui\data\m484_L_NAc_R_VTA_oxLight\Dec_12_2023\session_1\Quinine.avi';
obj = VideoReader(fileName);
numFrames = obj.NumFrames;
FrameRate = obj.FrameRate;

%%
writerObj=VideoWriter('P:\Yihui\data\m484_L_NAc_R_VTA_oxLight\Dec_12_2023\session_1\Quinine_trim.avi');
writerObj.FrameRate = FrameRate;
open(writerObj);

%%
for j = 6011 : 6610
    fprintf('running frame = %i\n',j)
    frame = read(obj,j);
    writeVideo(writerObj,frame);
end

close(writerObj);