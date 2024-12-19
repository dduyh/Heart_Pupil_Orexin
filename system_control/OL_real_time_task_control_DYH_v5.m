%Spatial and Reward coding during Observational learning in mice

%Real-time image acquisition and mice tracking workflow:

%Left cage: Demo Mouse
%Right cage: Observer Mouse

%Open the matlab data acquisition toolbox to preview,then close preview.
%Change the saving file folders’ names;
%Check the serial port number of Arduino board;
%Start the video recording, take the first frame as background;
%Place the mice into cages;

%Run the full script, which automatically runs the whole trial:
%Step 1: To start the trial: open the cage door to release the mouse;
%Step 2: Close the cage door when the mouse is detected leaving the cage;
%(Step 3: Provide food for the observer mice, when the demo mouse find the reward)
%Step 4: After the mouse finding the reward , open the cage door to let the mouse
%        back to cage, and provide food for the demo mouse;
%        Or after maximum 3 mins, open the cage door to let the mouse back to cage;
%Step 5: Close the cage door when the mouse is detected entering the cage;

%Author:  Yihui Du
%Date: 	November 22th, 2022
%%
clear; close all; clc;

%% set the path for output data
Directory = 'D:\yihui\CBM_data\';    % Main directory\
date = 'day7_Nov_23_2022';           % Date\
mouse_name = 'mouse_LLR';            % Mouse name\
type = 'demo';                       % Mouse type (demo or obs)\

%% Initialize the door and food dispenser controlled by Arduino
delete(instrfindall);                % Delete all existing serial port objects
s = serial('COM4');                  % Connect to serial port (COM4)
set(s,'BaudRate',9600);              % Set baud rate (communication speed in bits per second)
set(s,'Timeout',30);                 % Set time out (allowed time in seconds to complete read and write)
set(s,'InputBufferSize',8388608);    % Set input buffer size

fopen(s);                            % Open connection to Arduino board
if (exist('board1','var'))           % Stop the running program in Arduino
    board1.stop;pause(0);
end

%% OPEN VIDEO INPUT (maze)

imaqreset                            % Disconnect and delete all image acquisition objects
imaqmex('feature','-limitPhysicalMemoryUsage',false);  % Set unlimited physical memory usage

vid_maze = videoinput('winvideo', 1, 'RGB24_744x480'); % Create video input for the camera above maze
set(vid_maze,'Timeout',35);                            % Set time out (allowed time to wait for image)
vid_maze.FramesPerTrigger = Inf;                       % Set frames per trigger to infinity 
vid_maze.ReturnedColorspace = 'grayscale';             % Set color space to gray scale 
src_maze = getselectedsource(vid_maze);                % Get current video source object
src_maze.FrameRate = '15.0000';                        % Set frame rate to 15 fps 
src_maze.ExposureMode = 'manual';                      % Set exposure mode to manual 
src_maze.Exposure = -7;                                % Set exposure value to -7 
src_maze.GainMode = 'manual';                          % Set gain mode to manual 
src_maze.Gain = 16;                                    % Set gain value to 16 

%% Get the background image

start(vid_maze);                                       % Start maze video collection
pause(3)                                               % Wait 3 secs for video configuration
background = getsnapshot(vid_maze);                    % Get one background image of the maze
stop(vid_maze);                                        % Stop maze video collection

Data_Folder = [Directory date '\' mouse_name '_' type '\']; % Set data folder name
if ~exist(Data_Folder,'dir')                           % Create data folder if it doesn't exist
    mkdir(Data_Folder)
end
save([Data_Folder 'background.mat'],'background');     % Save the background image of the maze

%%
% set the area coordinates of inner boundary for left cage
figure(1);                                             % Create a figure window to show background image
disp('Please draw the inner boundary for left cage, and double click to save the coordinates    ');
bw_in_left=roipoly(background);                        % Interactively create a region for inner boundary of left cage  
[r_in_left,c_in_left]=find(bw_in_left==1);             % Get the row and column coordinates of all points in inner region

left_in_x1 = min(c_in_left);               % Get the horizontal coordinate of the first vertex of left region
left_in_y1 = min(r_in_left(c_in_left == min(c_in_left))); % Get the vertical coordinate of the first vertex of left region

left_in_x2 = min(c_in_left(r_in_left == min(r_in_left))); % Get the horizontal coordinate of the second vertex of left region
left_in_y2 = min(r_in_left);               % Get the vertical coordinate of the second vertex of left region

left_in_x3 = max(c_in_left);               % Get the horizontal coordinate of the third vertex of left region
left_in_y3 = max(r_in_left(c_in_left == max(c_in_left))); % Get the vertical coordinate of the third vertex of left region

left_in_x4 = max(c_in_left(r_in_left == max(r_in_left))); % Get the horizontal coordinate of the fourth vertex of left region
left_in_y4 = max(r_in_left);               % Get the vertical coordinate of the fourth vertex of left region

coordinates_in_left = [left_in_x1 left_in_x2 left_in_x3 left_in_x4; ...
    left_in_y1 left_in_y2 left_in_y3 left_in_y4]';  % Get the coordinates of all four vertices of left inner region

%%
% set the area coordinates of outer boundary for left cage
figure(1);                                             % open the figure window showing background image
disp('Please draw the outer boundary for left cage, and double click to save the coordinates    ');
bw_out_left=roipoly(background);                       % Interactively create a region for outer boundary of left cage  
[r_out_left,c_out_left]=find(bw_out_left==1);          % Get the row and column coordinates of all points in outer region

left_out_x1 = min(c_out_left);               % Get the horizontal coordinate of the first vertex of left region
left_out_y1 = min(r_out_left(c_out_left == min(c_out_left))); % Get the vertical coordinate of the first vertex of left region

left_out_x2 = min(c_out_left(r_out_left == min(r_out_left))); % Get the horizontal coordinate of the second vertex of left region
left_out_y2 = min(r_out_left);               % Get the vertical coordinate of the second vertex of left region

left_out_x3 = max(c_out_left);               % Get the horizontal coordinate of the third vertex of left region
left_out_y3 = max(r_out_left(c_out_left == max(c_out_left))); % Get the vertical coordinate of the third vertex of left region

left_out_x4 = max(c_out_left(r_out_left == max(r_out_left))); % Get the horizontal coordinate of the fourth vertex of left region
left_out_y4 = max(r_out_left);               % Get the vertical coordinate of the fourth vertex of left region

coordinates_out_left = [left_out_x1 left_out_x2 left_out_x3 left_out_x4; ...
    left_out_y1 left_out_y2 left_out_y3 left_out_y4]';  % Get the coordinates of all four vertices of left outer region

%% OPEN VIDEO INPUT (cages)

vid_cage = videoinput('tisimaq_r2013_64', 4, 'RGB24 (640x480) [Binning 2x]'); % Create video input for the camera above cage
set(vid_cage,'Timeout',35);                            % Set time out (allowed time to wait for image)
vid_cage.FramesPerTrigger = Inf;                       % Set frames per trigger to infinity 
vid_cage.ReturnedColorspace = 'grayscale';             % Set color space to gray scale 
src_cage = getselectedsource(vid_cage);                % Get current video source object
src_cage.ExposureAuto = 'Off';                         % Set exposure mode to manual 
src_cage.Exposure = 0.0666;                            % Set exposure value to 0.0666 
src_cage.GainAuto = 'Off';                             % Set gain mode to manual 
src_cage.Gain = 0;                                     % Set gain value to 0 

FrameRate = 15;                                        % Set frame rate parameter to 15 fps 

%% Set the center coordinates for rewards
figure(1);                                                % Open the figure window showing background image
imshow(background);
reward_center = ginput(1);                                % Interactively select a point for reward zone center  

reward_on_maze = 1;                                       % Whether food are placed on maze: true(1) / false(0)
reward_radius = 45;                                       % Radius of reward zone: 45 pixels (10 cm)
duration_range = 3;                                       % Search mouse position in previous 3 secs from current time point
reward_duration_range = duration_range*FrameRate;         % Number of frames in previous 3 secs
duration_threshold = 2;                                   % Mouse needs to stay in reward zone for 2 secs
reward_duration_threshold = duration_threshold*FrameRate; % Number of frames for the 2 secs threshold

%% Set the tracking Parameters

trial_length = 600;            % Set longest duration for one trial to 600 secs
exploration_threshold = 180;   % End the trial after 180 secs of exploration without reaching reward zone
initialization_threshold = 60; % End the trial after 60 secs without mouse's leaving the cage

thresh = 0.4;                  % Threshold for converting image to binary image
cage_in_thresh = 20000;        % Threshold for testing whether mouse enters the inner boundary of cage   
cage_out_thresh = 36000;       % Threshold for testing whether mouse leaves the outer boundary of cage

% start(vid_maze);
% pause(3)
% ROI_frame = getsnapshot(vid_maze);  % acquisite the maze image
% stop(vid_maze);
% 
% ROI_Im = imabsdiff(ROI_frame,background);
% cage_in_left = ROI_Im.*uint8(bw_in_left);
% cage_out_left = ROI_Im.*uint8(bw_out_left);
% 
% cage_in_thresh = sum(cage_in_left,'all')+1000;
% cage_out_thresh = sum(cage_out_left,'all')+1000;

%% Initiate the trial

trial_num = '1';               % Change the number for every trial

trial = ['trial' trial_num];   % Convert the trial number to a string

Data_Folder = [Directory date '\' mouse_name '_' type '\' trial '\']; % Set data folder name for every trial

if exist(Data_Folder,'dir')    % Ask user if they want to overwrite the data folder if it exists
    promptMessage = sprintf('This directory already exists:\n%s\nDo you want to overwrite it?', Data_Folder);
    titleBarCaption = 'Overwrite?';
    buttonText = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
    if strcmpi(buttonText, 'No')
        error('trial number already exists') % Display message if user does not want to overwrite
    else
        close;                               % Close the question dialog box if user want to overwrite
    end
else
    mkdir(Data_Folder)                       % Create data folder if it doesn't exist
end

% Open video output

writerObj=VideoWriter([Data_Folder 'processed_behavior.avi']); % Create object to write processed video
writerObj.FrameRate = FrameRate;                               % Set frame rate for output video to 15 fps 
open(writerObj);                                               % Open output video object

% Save the tracking parameters.

centroids = zeros(1,2);                    % Centroids to temporarily save coordinates of mouse in each frame 
centers = zeros(trial_length*FrameRate,2); % Centers to save coordinates of mouse for all frames 

i = 1;                                     % Parameter i records the currently analyzed frame number
T = zeros(1,trial_length*FrameRate);       % T record the exact time when each frame was captured

step_timepoint = zeros(1,6);               % Record the exact time when each operation step was performed

Step_2 = 1; % Set the initial state of step 2 to True 
Step_4 = 1; % Set the initial state of step 4 to True 
Step_5 = 1; % Set the initial state of step 5 to True 

reward_duration = 0;                                % Record the duration in number of frames for mouse staying in the reward zone
reward_timepoint = zeros(1,trial_length*FrameRate); % Record the frame number when mouse stays in the reward zone

% start the trial

start(vid_maze);        % Start maze video collection
t_vid_maze_Start = tic; % Record the exact time when maze video collection was started

start(vid_cage);        % Start cage video collection
t_vid_cage_Start = tic; % Record the exact time when cage video collection was started

tStart = tic; % Start timer to measure elapsed time

fprintf(s,'left_door_open/'); % Step 1: Open the left cage door to release the Demo mouse
step_timepoint(1)= toc(tStart); % Record the exact time for step 1
disp('Step 1: Open the left cage door to release the Demo mouse    ');

% Start capturing video frames in real time
while toc(tStart) < trial_length                % Capture video frames before exceeding longest trial duration
    
    ROI_frame = getsnapshot(vid_maze);          % Acquisite one frame from maze video
    T(i)= toc(tStart);                          % Record the exact time when every frame was captured
    i = i + 1;                                  % Parameter i records the currently analyzed frame number
    
    ROI_Im = imabsdiff(ROI_frame,background);   % Subtract the background image from the current frame to extract mouse shape
    
    cage_in_left = ROI_Im.*uint8(bw_in_left);   % Extract the left inner cage image
    cage_out_left = ROI_Im.*uint8(bw_out_left); % Extract the left outer cage image
    
    %     figure(2);
    %     imshow(ROI_Im);
    
    figure(1);                                  % Open figure window 1
    imshow(ROI_frame);                          % Show the current frame
    hold on                                     % Hold on the plotted figure
    
    %% Plot rectangular regions for two cages and circular area for reward
    
    patch('XData',[left_in_x1 left_in_x2 left_in_x3 left_in_x4],'YData',[left_in_y1 left_in_y2 left_in_y3 left_in_y4],...
        'EdgeColor','yellow','FaceColor','none','LineWidth',1); % Plot a yellow rectangular region for the left inner cage
    
    patch('XData',[left_out_x1 left_out_x2 left_out_x3 left_out_x4],'YData',[left_out_y1 left_out_y2 left_out_y3 left_out_y4],...
        'EdgeColor','green','FaceColor','none','LineWidth',1); % Plot a green rectangular region for the left outer cage
    
    plot(reward_center(:,1),reward_center(:,2), 'b.'); % Plot a blue dot in the reward center
    viscircles(reward_center,reward_radius,'Color','blue','LineWidth',0.4,'LineStyle','--'); % Plot a blue circle as the reward zone
    
    %% Detection of mouse in cage 
    if sum(cage_in_left,'all') > cage_in_thresh   % mouse is detected in the cage
        patch('XData',[left_in_x1 left_in_x2 left_in_x3 left_in_x4],'YData',[left_in_y1 left_in_y2 left_in_y3 left_in_y4],...
            'EdgeColor','yellow','FaceColor','yellow','FaceAlpha',0.5,'LineWidth',1); 
                                                  % Fill the inner rectangle with semitransparent yellow
    end
    
    if sum(cage_out_left,'all') < cage_out_thresh % mouse is detected out the cage
        patch('XData',[left_out_x1 left_out_x2 left_out_x3 left_out_x4],'YData',[left_out_y1 left_out_y2 left_out_y3 left_out_y4],...
            'EdgeColor','green','FaceColor','green','FaceAlpha',0.3,'LineWidth',1); 
                                                  % Fill the outer rectangle with semitransparent green
    end
    
    %% Find the mouse center
    
    I = im2bw(ROI_Im,thresh);           % Convert the subtracted image to binary image, based on threshold
    k = regionprops('table',I,'Area');  % Measure the area(number of pixels) of binary image regions
    idx = find(max([k.Area]));          % Get the index of the largest region
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
      
    %% Control the door and food dispenser

    if Step_2                                      % The step 2 operation has not been performed
        if sum(cage_out_left,'all') < cage_out_thresh % Mouse is detected out the cage
            fprintf(s,'left_door_close/');         % Step 2: Close the left cage door when Demo mouse has left
            step_timepoint(2)= toc(tStart);        % Record the exact time for step 2
            disp('Step 2: Close the left cage door, Demo mouse is leaving the left cage    ');
            Step_2 = 0;                            % Set the state of step 2 to False 
            
        elseif toc(tStart) > initialization_threshold % After waiting for maximum 1 mins
            fprintf(s,'left_door_close/');         % Close the cage door to end the trial
            step_timepoint(2)= toc(tStart);        % Record the exact time for step 2
            disp('Step 2: After 1 min, close the cage door to end the trial    ');
            Step_2 = 0;                            % Set the state of step 2 to False 
            break                                  % Break from the loop of capturing video frames
        end
        
    % (Step 3: Provide food for the observer mice, when the demo mouse found the reward)
          
    elseif Step_4                                  % The step 4 operation has not been performed
        
        if pdist([reward_center ; centroids],'euclidean') < reward_radius % The mouse is in the reward zone
            reward_timepoint(i) = 1;               % Mark the frame number when mouse stays in the reward zone
            reward_duration = length(find(reward_timepoint(max([i-reward_duration_range,1]):i)));
                                                   % Calculate the duration in number of frames for mouse staying in the reward zone
            if reward_duration > reward_duration_threshold % Mouse has stayed in reward zone for 2 secs
                fprintf(s,'left_door_open/');      % Step 4: Open the cage door to let the mouse back to cage
                step_timepoint(4)= toc(tStart);    % Record the exact time for step 4
                disp('Step 4: Open the left cage door, Demo mouse is in the reward zone for 2 secs    ');
                Step_4 = 0;                        % Set the state of step 4 to False 
                
                pause(2)                           % Wait for 2 secs
                fprintf(s,'left_food/');           % Provide food for the Demo mouse;
                disp('Provide food for the demo mouse    ');
            end
           
        elseif toc(tStart) > exploration_threshold % After exploring for maximum 3 mins
            fprintf(s,'left_door_open/');          % Step 4: Open the cage door to let the mouse back to cage
            step_timepoint(4)= toc(tStart);        % Record the exact time for step 4
            disp('Step 4: After 3 mins, open the cage door to let the mouse back to cage    ');
            Step_4 = 0;                            % Set the state of step 4 to False 
        end

    elseif Step_5                                  % The step 5 operation has not been performed
        if sum(cage_in_left,'all') > cage_in_thresh   % Mouse is detected in the cage
            fprintf(s,'left_door_close/');         % Step 5: Close the left cage door when Demo mouse is in cage
            step_timepoint(5)= toc(tStart);        % Record the exact time for step 5
            disp('Step 5: Close the left cage door, Demo mouse is in the left cage    ');
            Step_5 = 0;                            % Set the state of step 5 to False 
            pause(3)                               % Wait for 3 secs
            fprintf(s,'left_door_close/');         % Close the left cage door again in case the operation went wrong;
            pause(7)                               % Wait for 7 secs
            break                                  % Break from the loop of capturing video frames
        end 
    end
end
tEnd = toc(tStart);                                % Get the exact ending time
step_timepoint(6)= tEnd;                           % Record the exact ending time
disp(['Elapsed time: ',num2str(tEnd),' seconds.']);% Display the exact ending time

close(writerObj);                       % Close output video object

stop(vid_maze);                         % Stop maze video collection
t_vid_maze_End = toc(t_vid_maze_Start); % Record the exact time when maze video was ended

stop(vid_cage);                         % Stop cage video collection
t_vid_cage_End = toc(t_vid_cage_Start); % Record the exact time when cage video was ended

close(figure(1));                       % Close the figure window
% close(figure(2));

% save maze video
data = getdata(vid_maze, vid_maze.FramesAvailable);% Get all collected frames of maze video
numFrames = size(data, 4);              % Get number of frames
 
diskLogger = VideoWriter([Data_Folder 'shaping_to_maze.avi'], 'Uncompressed AVI'); % Create object to write the uncompressed maze video
diskLogger.FrameRate = numFrames./t_vid_maze_End; % Set the frame rate
open(diskLogger);                                 % Open output video object
for ii = 1:numFrames
    writeVideo(diskLogger, data(:,:,:,ii));       % Write each frame into the maze video
end
close(diskLogger);                                % Close output video object

% save cage video
data = getdata(vid_cage, vid_cage.FramesAvailable);% Get all collected frames of cage video
numFrames = size(data, 4);              % Get number of frames

diskLogger = VideoWriter([Data_Folder 'shaping_to_cage.avi'], 'Uncompressed AVI'); % Create object to write the uncompressed cage video
diskLogger.FrameRate = numFrames./t_vid_cage_End; % Set the frame rate
open(diskLogger);                                 % Open output video object
for ii = 1:numFrames
    writeVideo(diskLogger, data(:,:,:,ii));       % Write each frame into the cage video
end
close(diskLogger);                                % Close output video object

% save data files
save([Data_Folder 'centers.mat'],'centers'); % Save coordinates of mouse for all frames 
save([Data_Folder 'step_timepoint.mat'],'step_timepoint'); % Save the exact time for each operation step
save([Data_Folder 'parameters.mat'],'reward_center','reward_radius','duration_threshold','reward_on_maze', ...
    'reward_timepoint','exploration_threshold','coordinates_in_left','coordinates_out_left','T'); % Save other parameters

% imaqreset;
