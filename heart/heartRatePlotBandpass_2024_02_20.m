


%%

piezoChannel=4;

injLoc=find(superRigMatMaster(:,piezoChannel)<-2,1);


selectDuration=2000;

selectStart1=injLoc-selectDuration*22;  % select regions to zoom in on
selectStop1=injLoc-selectDuration*21;

selectStart2=injLoc+selectDuration*121;
selectStop2=injLoc+selectDuration*122;

movMeanFact=10;

fpass=[6 15];
fs=400;

channelHeartRate=6;
channelRun=2;

minPeakPromVal=0.001;

%%

heartSensorBandpass = bandpass(superRigMatMaster(:,channelHeartRate),fpass,fs);


heartData=superRigMatMaster(:,channelHeartRate);

%% put axes in seconds

xLength=size(superRigMatMaster,1);

realTimeMsMax=xLength*2.5;  %1000ms/400Hz

xAxisMs=1:2.5:realTimeMsMax;
xAxisMs=xAxisMs/1000;

%% calculate running correctly

runningData=abs(diff(superRigMatMaster(:,channelRun)));


%%

figure;

% subplot(3,1,1)

plot(xAxisMs,superRigMatMaster(:,piezoChannel))
% 
xline((injLoc*2.5)/1000,'-k')

xline((selectStart1*2.5)/1000,'-c')
xline((selectStop1*2.5)/1000,'-c')
xline((selectStart2*2.5)/1000,'-r')
xline((selectStop2*2.5)/1000,'-r')



xlabel('time (s)');
ylabel ('piezo a.u.');

sgtitle('check injection time found')
% 
% subplot(3,1,2)
% 
% plot(heartSensorBandpass)
% 
% 
% subplot(3,1,3)
% 
% plot(movmean(heartData,movMeanFact))
% hold on;
% plot(heartSensorBandpass+medianHR,'-k')


%% plot pre and post propranolol

medianHR=median(heartData);
movMeanFactRun=100;



yRun(1)=min(min(movmean(runningData(selectStart1:selectStop1),movMeanFactRun)),min(movmean(runningData(selectStart2:selectStop2),movMeanFactRun)));
yRun(2)=max(max(movmean(runningData(selectStart1:selectStop1),movMeanFactRun)),max(movmean(runningData(selectStart2:selectStop2),movMeanFactRun)));


%% plot pre propranolol 

figure('Renderer', 'painters', 'Position', [50 100 1200 800]);


subplot(3,2,1)
plot(movmean(runningData(selectStart1:selectStop1),movMeanFactRun),'-k')

axis([1 size(selectStart1:selectStop1,2)+500 yRun(1) yRun(2)])

ylabel('Run speed')

subplot(3,2,3)

plot(heartData(selectStart1:selectStop1),'-c')
hold on;
plot(heartSensorBandpass(selectStart1:selectStop1)+medianHR,'-k')

legend('raw','bandpass')

% xline(100,'-r')
% xline(500,'-r')

ylabel('ecg raw')

subplot(3,2,5)
plot(movmean(heartData(selectStart1:selectStop1),movMeanFact),'-c')
hold on;
plot(heartSensorBandpass(selectStart1:selectStop1)+medianHR,'-k')

legend('movMean','bandpass')

xlabel('PRe time(s)')

ylabel('ecg movmean')

% plot post propranolol 


subplot(3,2,2)

plot(movmean(runningData(selectStart2:selectStop2),movMeanFactRun),'-k')

ylabel('Run speed')

xline(100,'-r')
xline(500,'-r')

axis([1 size(selectStart1:selectStop1,2)+500 yRun(1) yRun(2)])

subplot(3,2,4)

plot(heartData(selectStart2:selectStop2),'-r')
hold on;
plot(heartSensorBandpass(selectStart2:selectStop2)+medianHR,'-k')

legend('raw','bandpass')

ylabel('ecg raw')


ylabel('HRSensor + bandpass')

subplot(3,2,6)
plot(movmean(heartData(selectStart2:selectStop2),movMeanFact),'-r')
hold on;
plot(heartSensorBandpass(selectStart2:selectStop2)+medianHR,'-k')

legend('movMean','bandpass')

ylabel('ecg movmean')
xlabel('POST time(s)')

sgtitleTitle=strcat(fileName(12:17),'pre post propranolol');

sgtitle(sgtitleTitle);

saveas(gcf,(strcat(datestr(now, 'ddmmyyyy'),sgtitleTitle)),'svg');

%%


[pksVal pksLocs]=findpeaks(heartSensorBandpass,'MaxPeakWidth',80,'MinPeakProminence',minPeakPromVal);

%%

pksLocsRealTime=(pksLocs*2.5)/1000;


figure;

plot(xAxisMs,heartData-medianHR,'-m')
hold on;
plot(xAxisMs,heartSensorBandpass,'-k')
%plot(pksLocs,pksVal,'.r')
plot(pksLocsRealTime,pksVal,'or')

xline((selectStart1*2.5)/1000,'-c')
xline((selectStop1*2.5)/1000,'-c')
xline((selectStart2*2.5)/1000,'-r')
xline((selectStop2*2.5)/1000,'-r')
xline((injLoc*2.5)/1000,'-k','linewidth',3)


axis([0 max(xAxisMs) -0.05 0.05])

xlabel('time (s)')

sgtitle('check detected peaks')

%% get heart rate from peakLocs

heartRateVals=1./diff(pksLocsRealTime);

%% plot heart rate over time



figure('Renderer', 'painters', 'Position', [400 400 800 400]);

subplot(2,1,1)

plot(xAxisMs(2:end),movmean(abs(diff(superRigMatMaster(:,channelRun))),movMeanFactRun),'-b')
xline((injLoc*2.5)/1000,'-k','linewidth',3)
ylabel('Running Speed (dist/s)')


subplot(2,1,2)
plot(pksLocsRealTime(2:end),heartRateVals,'-m')
hold on
plot(pksLocsRealTime(2:end),movmean(heartRateVals,movMeanFact),'-k')


xline((selectStart1*2.5)/1000,'-c')
xline((selectStop1*2.5)/1000,'-c')
xline((selectStart2*2.5)/1000,'-r')
xline((selectStop2*2.5)/1000,'-r')
xline((injLoc*2.5)/1000,'-k','linewidth',3)

axis([0 max(xAxisMs) 0 20])

xlabel('Time (s)')
ylabel('Heart Rate (Hz)')

sgtitleTitle=strcat(fileName(12:17),'heart Rate');

saveas(gcf,(strcat(datestr(now, 'ddmmyyyy'),sgtitleTitle)),'svg');

%%

%powerSpec=pspectrum(superRigMatMaster(:,channelHeartRate),fs,spectogram);

% poweSpecMax=max(powerSpec);
% 
% window=1:10:900:

% 
% heartSpec=spectrogram(superRigMatMaster(:,hrChannel),fs);
% 
% figure;
% spectogram()

%% power spectrum

M=100; %number of samples
L=16; % overlap percent
lk=0.7;

OverlapPercent=(L/M)*100
TimeResolution=M/fs

figure

pspectrum(superRigMatMaster(:,channelHeartRate),fs,"spectrogram", ...
    TimeResolution=M/fs,OverlapPercent=(L/M)*100, ...
    Leakage=lk);

title("pspectrum")
cc = clim;
xl = xlim; 
yl = 20;


%%
% 
% rotThresh=((max(superRigMatMaster(:,2))-min(superRigMatMaster(:,2)))/2);
% 
% 
% for binarisedLoop=1:numel(superRigMatMaster(:,2))
% 
% 
%     if superRigMatMaster(binarisedLoop,2)>rotThresh
% 
%         binRunMat(binarisedLoop,1)=1;
% 
%     else
% 
%         binRunMat(binarisedLoop,1)=0;
% 
% 
%     end
% 
% end
% 
% %%
% 
% binRunMatBinned=[];
% 
% 
% runMatDiff=[(abs(diff(binRunMat))); NaN];
% 
% binSize=400;
% 
% for binLoop=1:binSize:numel(runMatDiff)-binSize %:size(masterTrialMat,1)
%     
% binRunMatBinned(binLoop:binLoop + (binSize-1),:)=nanmean(runMatDiff(binLoop:binLoop + (binSize-1)));
% 
% end


%% removes nan values

% idxToRemove = all(all(isnan(matHRPeaks),3),2);
% matHRPeaks(idxToRemove,:,:) = [];

      
%%
% 
% HRRealTime=pksLocs*2.5; %converted detected peaks to real time
% HrRtIntervl=diff(HRRealTime); % get inter-beat interval

% for HRconvertLoop=1:size(HrRtIntervl,1)
% 
%     HrRtHz(HRconvertLoop)=1000/HrRtIntervl(HRconvertLoop); %calculate num pulses per second
% 
% end

%% discretize and bin data
% 
% 
% binSizeDiscret=1000; %in ms (running is 400, becasue that's 1000 ms)
% 
% edges=1:binSizeDiscret:size(superRigMatMaster(:,channelHeartRate)*2.5,1);
% 
% HeartRateDiscret = discretize(HRRealTime,edges);
% 
% for findLoop=1:max(HeartRateDiscret)
% 
%    HrPerBin(1,findLoop)=numel(find(HeartRateDiscret==findLoop));
% 
% end
% 
% 
% %%
% 
% 
% binHRMat=[];
% 
% HRBin=1000; %1000 = 1 second 
% 
% binHRedges=1:HRBin:max(HRRealTime);
% 
% for binHRLoop=1:max(HRRealTime)
% 
% 
% binHRMat(binHRLoop,1)=numel(find((HRRealTime>binHRedges(binHRLoop)) & (HRRealTime<binHRedges(binHRLoop+1))));
% 
% 
% end
% 
% 
% binHRMat=binHRMat
% 
% figure;plot(binHRMat,'-ok')
% 
% %%
% 
% moveMeanFact=400;
% 
% figure;
% subplot(3,1,1)
% plot(binRunMatBinned)
% 
% subplot(3,1,2)
% plot(HrPerBin,'-ok')
% 
% subplot(3,1,3)
% plot(heartSensorBandpass,'-k')
% hold on;
% plot(pksLocs,pksVal,'.r')
% plot(pksLocs,pksVal,'or')
% 
% 
% 



