
clear; close all; clc;

%% set the path for output data

Directory = 'P:\Yihui\data\';                     % Main directory\
mouse_name = 'm483_L_VTA_R_NAc_oxLight';            % Mouse name\
date = 'Dec_12_2023';                             % Date\

Sucrose = 0;

%% Load data

session_1_Folder = [Directory mouse_name '\' date '\session_1\'];     % Set data folder name
session_1_timepoint = load([session_1_Folder 'step_timepoint.mat']).step_timepoint;

session_2_Folder = [Directory mouse_name '\' date '\session_2\'];     % Set data folder name
session_2_timepoint = load([session_2_Folder 'step_timepoint.mat']).step_timepoint;

%%

Frame_Rate = 20;

session_1 = session_1_timepoint(1,:) - session_1_timepoint(1,1);
session_1_frames = ceil(session_1*Frame_Rate);
session_1_emotions = zeros(1,session_1_frames(end));

session_2 = session_2_timepoint(1,:) - session_2_timepoint(1,1);
session_2_frames = ceil(session_2*Frame_Rate);
session_2_emotions = zeros(1,session_2_frames(end));

if Sucrose
    
    session_1_emotions(session_1_frames(3):session_1_frames(4)) = 4;
    session_1_emotions(session_1_frames(5):session_1_frames(6)) = 4;
    session_1_emotions(session_1_frames(7):session_1_frames(8)) = 4;
    session_1_emotions(session_1_frames(9):session_1_frames(10)) = 4;
    session_1_emotions(session_1_frames(11):session_1_frames(12)) = 4;
    session_1_emotions(session_1_frames(13):session_1_frames(14)) = 4;
    
    session_1_emotions(1:session_1_frames(3)) = 3;
    session_1_emotions(session_1_frames(4):session_1_frames(5)) = 3;
    session_1_emotions(session_1_frames(6):session_1_frames(7)) = 3;
    session_1_emotions(session_1_frames(8):session_1_frames(9)) = 3;
    session_1_emotions(session_1_frames(10):session_1_frames(11)) = 3;
    session_1_emotions(session_1_frames(12):session_1_frames(13)) = 3;
    session_1_emotions(session_1_frames(14):session_1_frames(15)) = 3;
    
    session_2_emotions(session_2_frames(3):session_2_frames(4)) = 2;
    session_2_emotions(session_2_frames(5):session_2_frames(6)) = 2;
    session_2_emotions(session_2_frames(7):session_2_frames(8)) = 2;
    session_2_emotions(session_2_frames(9):session_2_frames(10)) = 2;
    session_2_emotions(session_2_frames(11):session_2_frames(12)) = 2;
    session_2_emotions(session_2_frames(13):session_2_frames(14)) = 2;
    
    session_2_emotions(1:session_2_frames(3)) = 1;
    session_2_emotions(session_2_frames(4):session_2_frames(5)) = 1;
    session_2_emotions(session_2_frames(6):session_2_frames(7)) = 1;
    session_2_emotions(session_2_frames(8):session_2_frames(9)) = 1;
    session_2_emotions(session_2_frames(10):session_2_frames(11)) = 1;
    session_2_emotions(session_2_frames(12):session_2_frames(13)) = 1;
    session_2_emotions(session_2_frames(14):session_2_frames(15)) = 1;
    
else
    
    session_1_emotions(session_1_frames(3):session_1_frames(4)) = 2;
    session_1_emotions(session_1_frames(5):session_1_frames(6)) = 2;
    session_1_emotions(session_1_frames(7):session_1_frames(8)) = 2;
    session_1_emotions(session_1_frames(9):session_1_frames(10)) = 2;
    session_1_emotions(session_1_frames(11):session_1_frames(12)) = 2;
    session_1_emotions(session_1_frames(13):session_1_frames(14)) = 2;
    
    session_1_emotions(1:session_1_frames(3)) = 1;
    session_1_emotions(session_1_frames(4):session_1_frames(5)) = 1;
    session_1_emotions(session_1_frames(6):session_1_frames(7)) = 1;
    session_1_emotions(session_1_frames(8):session_1_frames(9)) = 1;
    session_1_emotions(session_1_frames(10):session_1_frames(11)) = 1;
    session_1_emotions(session_1_frames(12):session_1_frames(13)) = 1;
    session_1_emotions(session_1_frames(14):session_1_frames(15)) = 1;
    
    session_2_emotions(session_2_frames(3):session_2_frames(4)) = 4;
    session_2_emotions(session_2_frames(5):session_2_frames(6)) = 4;
    session_2_emotions(session_2_frames(7):session_2_frames(8)) = 4;
    session_2_emotions(session_2_frames(9):session_2_frames(10)) = 4;
    session_2_emotions(session_2_frames(11):session_2_frames(12)) = 4;
    session_2_emotions(session_2_frames(13):session_2_frames(14)) = 4;
    
    session_2_emotions(1:session_2_frames(3)) = 3;
    session_2_emotions(session_2_frames(4):session_2_frames(5)) = 3;
    session_2_emotions(session_2_frames(6):session_2_frames(7)) = 3;
    session_2_emotions(session_2_frames(8):session_2_frames(9)) = 3;
    session_2_emotions(session_2_frames(10):session_2_frames(11)) = 3;
    session_2_emotions(session_2_frames(12):session_2_frames(13)) = 3;
    session_2_emotions(session_2_frames(14):session_2_frames(15)) = 3;
    
end

emotions = [session_1_emotions session_2_emotions];

%%

uiopen([Directory mouse_name '\' date '\cluster_index_8.fig'],1)

h=findobj(gcf,'type','axes');

h2 = h(2);
set(h2, 'Position', [0.298    0.83    0.44    0.0312])

mymap = [198 219 239
    33 113 181
    233 157 147
    227 26 28];
mymap = mymap./255;

imagesc(emotions)
colormap(h2,mymap);
axis off


