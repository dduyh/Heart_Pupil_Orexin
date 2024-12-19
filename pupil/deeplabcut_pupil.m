clc
clear
close all

%%
folderPath = 'D:\duyh\video\pupil video\';

data = csvread([folderPath 'WIN_20200814_16_17_29_ProDLC_mobnet_100_pupil_cropMar31shuffle1_300000.csv'],3,1);

Pupil_up_x = data(:,1);
Pupil_up_y = data(:,2);
Pupil_left_x = data(:,4);
Pupil_left_y = data(:,5);
Pupil_down_x = data(:,7);
Pupil_down_y = data(:,8);
Pupil_right_x = data(:,10);
Pupil_right_y = data(:,11);

%%
center_x = zeros(size(data,1),1);
center_y = zeros(size(data,1),1);
radii = zeros(size(data,1),1); 
areas = zeros(size(data,1),1); 

for i = 1:size(data,1)
    X1(1) = Pupil_up_x(i,1);
    X1(2) = Pupil_up_y(i,1);
    
    X2(1) = Pupil_left_x(i,1);
    X2(2) = Pupil_left_y(i,1);
    
    Y1(1) = Pupil_down_x(i,1);
    Y1(2) = Pupil_down_y(i,1);
    
    Y2(1) = Pupil_right_x(i,1);
    Y2(2) = Pupil_right_y(i,1);
    
    [center_x(i,1), center_y(i,1)]= node(X1,Y1,X2,Y2);
    
    MinorAxis = [X1(1), X1(2); Y1(1), Y1(2)];
    MinorAxisLength = pdist(MinorAxis, 'euclidean'); 
    
    MajorAxis = [X2(1), X2(2); Y2(1), Y2(2)];
    MajorAxisLength = pdist(MajorAxis, 'euclidean');
    
    diameters = mean([MajorAxisLength, MinorAxisLength]);
    radii(i,1) = diameters/2;
    areas(i,1) = pi*radii(i,1).^2;
end

save([folderPath 'Pupil_Tracking_1.mat'], 'radii', 'areas', 'center_x', 'center_y');

%% Read original video

fileName = [folderPath 'WIN_20200814_16_17_29_ProDLC_mobnet_100_pupil_cropMar31shuffle1_300000_labeled.mp4'];
obj = VideoReader(fileName);
numFrames = obj.NumFrames;
FrameRate = obj.FrameRate;

figure(1);

%% Create analysed video

writerObj=VideoWriter([folderPath 'Pupil_Tracking_1.mp4'], 'MPEG-4');
writerObj.FrameRate = FrameRate;
open(writerObj);

for i = 1 : numFrames
    
    frame = read(obj,i);
    disp(['processing...',num2str(i),'/',num2str(numFrames)]);
    
    
    figure(1);
    ROI_Im = frame(201:345,726:1000,:);
    imshow(ROI_Im);
    hold on
    
    set(gca,'YDir','reverse')
    plot(center_x(i,1), center_y(i,1), 'r.')
    viscircles([center_x(i,1), center_y(i,1)],radii(i,1),'LineWidth',0.05);
    hold off
    
    fig = getframe(gcf);
    writeVideo(writerObj,fig.cdata);
end

close(writerObj);

%% plot

figure(2)
time = 1 : numFrames/FrameRate;
medianArea = [];
stdArea = [];
for i = 1 : numFrames/FrameRate
medianArea(i) = median(areas(((i-1)*FrameRate+1):i*FrameRate));
stdArea(i) = std(areas(((i-1)*FrameRate+1):i*FrameRate));
end
fill([time fliplr(time)],[medianArea+stdArea fliplr(medianArea-stdArea)], [0.8 0.8 0.8])
hold on
plot(medianArea, 'LineWidth',2)
%errorbar(time,medianArea,stdArea);
hx = xlabel('recording time(s)');
set(hx,'fontname','Times New Roman','fontsize',18)
hy = ylabel('pupil area (pixels)');
set(hy,'fontname','Times New Roman','fontsize',18)
savefig([folderPath 'WIN_20200814_16_17_29_Pro.fig']);
saveas(gcf,[folderPath 'WIN_20200814_16_17_29_Pro.png']) ;

ht = title('Extract Mouse Pupil Size by Image Threshold');
set(ht,'fontname','Times New Roman','fontsize',18)

ylim([0 5000])