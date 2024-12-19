%% Analyse ECG signals

clear; close all; clc;

%% set the path for output data

Directory = 'P:\Yihui\data\';                     % Main directory\
mouse_name = 'm758_L_MCH_R_Orx_GCaMP';            % Mouse name\
date = 'Feb_20_2024';                             % Date\

Data_Folder = [Directory mouse_name '\' date '\'];     % Set data folder name

%% load data

datas_rest = csvread([Data_Folder 'ECG_shirt_rest.csv'],1,0);
datas_run = csvread([Data_Folder 'ECG_shirt_run.csv'],1,0);

%% Raw ECG signals

channelTime=1;
time_rest = datas_rest(:,channelTime);
time_run = datas_run(:,channelTime);

channelHeartRate=2;
ECG_rest_raw = datas_rest(:,channelHeartRate);
ECG_run_raw = datas_run(:,channelHeartRate);

figure(1);
plot(time_rest(1:1000),ECG_rest_raw(1:1000),'g')

%% Remove baseline wandering

Sample_Rate = 100;    % 100 scans per second.

fpass=10;
ECG_rest_Highpass = highpass(ECG_rest_raw,fpass,Sample_Rate);
ECG_run_Highpass = highpass(ECG_run_raw,fpass,Sample_Rate);

figure(2);
plot(time_rest(1:1000),ECG_rest_Highpass(1:1000),'g')
hold on
plot(time_run(1:1000)-time_run(1),ECG_run_Highpass(1:1000),'r')

%% find peaks

minPeakPromVal=200;
figure(3);
findpeaks(ECG_rest_Highpass,'MaxPeakWidth',50,'MinPeakProminence',minPeakPromVal);
findpeaks(ECG_run_Highpass,'MaxPeakWidth',50,'MinPeakProminence',minPeakPromVal);
[pksVal_rest, pksLocs_rest]=findpeaks(ECG_rest_Highpass,'MaxPeakWidth',50,'MinPeakProminence',minPeakPromVal);
[pksVal_run, pksLocs_run]=findpeaks(ECG_run_Highpass,'MaxPeakWidth',50,'MinPeakProminence',minPeakPromVal);

%% heart rate in time 

heartRate_rest=1./diff(pksLocs_rest);
heartRate_run=1./diff(pksLocs_run);
figure(4);
plot([heartRate_rest; heartRate_run])

