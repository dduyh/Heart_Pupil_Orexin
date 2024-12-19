%%
clear; close all; clc;

%% set the path for output data

Directory = 'P:\Yihui\data\';                     % Main directory\
mouse_name = 'm483_L_VTA_R_NAc_oxLight';            % Mouse name\
date = 'Mar_18_2024';                             % Date\

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% Load the data.

fileName = [Data_Folder 'Trial_22.mp4'];
obj = VideoReader(fileName);
numFrames = obj.NumFrames;
FrameRate= obj.FrameRate;

%%
figure(1);

frame1 = read(obj,1);

bw_in_right = roipoly(frame1);                        % Interactively create a region for inner boundary of right cage
[r_in_right,c_in_right]=find(bw_in_right==1);             % Get the row and column coordinates of all points in inner region

thresh = 0.5;                  % Threshold for converting image to binary image

% Save the tracking parameters.
centroids = zeros(1,2);                    % Centroids to temporarily save coordinates of mouse in each frame
centers = zeros(numFrames,2); % Centers to save coordinates of mouse for all frames

%% Tone 1.

writerObj=VideoWriter([Data_Folder 'Trial_22_Tone_1.mp4'],'MPEG-4');
writerObj.FrameRate = FrameRate;
open(writerObj);

for i = 180*FrameRate+1 : 210*FrameRate

    frame1 = read(obj,i);
    disp(['processing...',num2str(i-180*FrameRate),'/',num2str(numFrames)]);

    f=frame1(:,:,1);

    ROI_Im = f(min(r_in_right):max(r_in_right),min(c_in_right):max(c_in_right));   % Extract the left inner cage image

    figure(1);                                  % Open figure window 1
    imshow(ROI_Im);                          % Show the current frame
    hold on                                     % Hold on the plotted figure

    %% Find the mouse center

    I = im2bw(ROI_Im,thresh);           % Convert the subtracted image to binary image, based on threshold
    I = imcomplement(I);
    s = regionprops('table',I,'Area','Eccentricity');
    idx = find([s.Area] > 1000 & [s.Eccentricity] < 0.9);
    cc = bwconncomp(I);                 % Find and count all connected components in the binary image
    g = ismember(labelmatrix(cc), idx); % Get the binary image of only the largest region
    m = regionprops('table',g,'Area','Centroid','MajorAxisLength','MinorAxisLength'); % Measure the size properties of the largest region
    if(size(m,1)==1)                    % If only one region survives
        centroids = cat(1, m.Centroid); % Get the coordinates of mouse
        centers(i,:) = centroids;       % Save the coordinates of mouse for this frame in centers
        plot(centroids(:,1),centroids(:,2), 'r.')                  % Plot a red dot in the mouse center
        diameters = mean([m.MajorAxisLength m.MinorAxisLength],2); % Calculate the diameter of the mouse region
        radii = diameters/2;                                       % Calculate the radius of the mouse region
        viscircles(centroids,radii,'LineWidth',0.1);               % Plot a red circle as the mouse position
    end
    hold off                            % Hold off the plotted frame

    frame2 = getframe;                  % Get the plotted frame
    writeVideo(writerObj,frame2);       % Write the plotted frame into the processed video

end

close(writerObj);                       % Close output video object

%%
Tone_1_centers = centers((180*FrameRate+1 : 210*FrameRate),:);
% Tone_1_centers(all(Tone_1_centers==0,2),:) = [];
Tone_1_locomotion = zeros(1,size(Tone_1_centers,1)-1);

for i = 1 : size(Tone_1_centers,1)-1
    mouse_dist = [Tone_1_centers(i,1), Tone_1_centers(i,2); Tone_1_centers(i+1,1), Tone_1_centers(i+1,2)];
    Tone_1_locomotion(i) = pdist(mouse_dist, 'euclidean');
end

figure
plot(Tone_1_locomotion)
ylim([0 27])

%% Tone 2.

writerObj=VideoWriter([Data_Folder 'Trial_22_Tone_2.mp4'],'MPEG-4');
writerObj.FrameRate = FrameRate;
open(writerObj);

for i = 392*FrameRate+1 : 422*FrameRate

    frame1 = read(obj,i);
    disp(['processing...',num2str(i-392*FrameRate),'/',num2str(numFrames)]);

    f=frame1(:,:,1);

    ROI_Im = f(min(r_in_right):max(r_in_right),min(c_in_right):max(c_in_right));   % Extract the left inner cage image

    figure(1);                                  % Open figure window 1
    imshow(ROI_Im);                          % Show the current frame
    hold on                                     % Hold on the plotted figure

    %% Find the mouse center

    I = im2bw(ROI_Im,0.45);           % Convert the subtracted image to binary image, based on threshold
    I = imcomplement(I);
    s = regionprops('table',I,'Area','Eccentricity');
    idx = find([s.Area] > 2000 & [s.Eccentricity] < 0.9);
    cc = bwconncomp(I);                 % Find and count all connected components in the binary image
    g = ismember(labelmatrix(cc), idx); % Get the binary image of only the largest region
    m = regionprops('table',g,'Area','Centroid','MajorAxisLength','MinorAxisLength'); % Measure the size properties of the largest region
    if(size(m,1)==1)                    % If only one region survives
        centroids = cat(1, m.Centroid); % Get the coordinates of mouse
        centers(i,:) = centroids;       % Save the coordinates of mouse for this frame in centers
        plot(centroids(:,1),centroids(:,2), 'r.')                  % Plot a red dot in the mouse center
        diameters = mean([m.MajorAxisLength m.MinorAxisLength],2); % Calculate the diameter of the mouse region
        radii = diameters/2;                                       % Calculate the radius of the mouse region
        viscircles(centroids,radii,'LineWidth',0.1);               % Plot a red circle as the mouse position
    end
    hold off                            % Hold off the plotted frame

    frame2 = getframe;                  % Get the plotted frame
    writeVideo(writerObj,frame2);       % Write the plotted frame into the processed video

end

close(writerObj);                       % Close output video object

%%
Tone_2_centers = centers((392*FrameRate+1 : 422*FrameRate),:);
% Tone_2_centers(all(Tone_2_centers==0,2),:) = [];
Tone_2_locomotion = zeros(1,size(Tone_2_centers,1)-1);

for i = 1 : size(Tone_2_centers,1)-1
    mouse_dist = [Tone_2_centers(i,1), Tone_2_centers(i,2); Tone_2_centers(i+1,1), Tone_2_centers(i+1,2)];
    Tone_2_locomotion(i) = pdist(mouse_dist, 'euclidean');
end

figure
plot(Tone_2_locomotion)
ylim([0 27])

%% Tone 3

writerObj=VideoWriter([Data_Folder 'Trial_22_Tone_3.mp4'],'MPEG-4');
writerObj.FrameRate = FrameRate;
open(writerObj);

for i = 604*FrameRate+1 : 634*FrameRate

    frame1 = read(obj,i);
    disp(['processing...',num2str(i-604*FrameRate),'/',num2str(numFrames)]);

    f=frame1(:,:,1);

    ROI_Im = f(min(r_in_right):max(r_in_right),min(c_in_right):max(c_in_right));   % Extract the left inner cage image

    figure(1);                                  % Open figure window 1
    imshow(ROI_Im);                          % Show the current frame
    hold on                                     % Hold on the plotted figure

    %% Find the mouse center

    I = im2bw(ROI_Im,0.5);           % Convert the subtracted image to binary image, based on threshold
    I = imcomplement(I);
    s = regionprops('table',I,'Area','Eccentricity');
    idx = find([s.Area] > 1000 & [s.Eccentricity] < 0.9);
    cc = bwconncomp(I);                 % Find and count all connected components in the binary image
    g = ismember(labelmatrix(cc), idx); % Get the binary image of only the largest region
    m = regionprops('table',g,'Area','Centroid','MajorAxisLength','MinorAxisLength'); % Measure the size properties of the largest region
    if(size(m,1)==1)                    % If only one region survives
        centroids = cat(1, m.Centroid); % Get the coordinates of mouse
        centers(i,:) = centroids;       % Save the coordinates of mouse for this frame in centers
        plot(centroids(:,1),centroids(:,2), 'r.')                  % Plot a red dot in the mouse center
        diameters = mean([m.MajorAxisLength m.MinorAxisLength],2); % Calculate the diameter of the mouse region
        radii = diameters/2;                                       % Calculate the radius of the mouse region
        viscircles(centroids,radii,'LineWidth',0.1);               % Plot a red circle as the mouse position
    end
    hold off                            % Hold off the plotted frame

    frame2 = getframe;                  % Get the plotted frame
    writeVideo(writerObj,frame2);       % Write the plotted frame into the processed video

end

close(writerObj);                       % Close output video object

%%
Tone_3_centers = centers((604*FrameRate+1 : 634*FrameRate),:);
% Tone_3_centers(all(Tone_3_centers==0,2),:) = [];
Tone_3_locomotion = zeros(1,size(Tone_3_centers,1)-1);

for i = 1 : size(Tone_3_centers,1)-1
    mouse_dist = [Tone_3_centers(i,1), Tone_3_centers(i,2); Tone_3_centers(i+1,1), Tone_3_centers(i+1,2)];
    Tone_3_locomotion(i) = pdist(mouse_dist, 'euclidean');
end

figure
plot(Tone_3_locomotion)
ylim([0 27])

%%
figure
freezing =[0 68/900 123/900]*100;
plot(1:3,freezing,'o-','linewidth',3,'markerfacecolor','b')
ylim([0 40])
xlim([0 4])
ylabel('Freezing (%)','FontSize',20,'FontWeight','bold');
xticklabels({'Tone1','Tone2','Tone3'})
