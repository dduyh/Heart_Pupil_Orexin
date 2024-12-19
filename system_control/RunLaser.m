% initialise NI box
s = daq.createSession('ni');
ch = addCounterOutputChannel(s,'Dev3','ctr0','PulseGeneration'); % make sure to set this to proper counter output and device

s.IsContinuous = true;

% setup stim schedule
nRepsPerFreq = 5;
stimFreqs = [1 5 10 20];

stimSchedule = repmat(stimFreqs,1,nRepsPerFreq);
stimSchedule = stimSchedule(randperm(length(stimSchedule)));

ntrials = length(stimSchedule);


% settings for trial timeline
stimLength = 30;
interTrialInterval = 120;
randITI = 60;

% wait for 180
pause(180)
for trialInd = 1:ntrials

disp(['Trial number ' num2str(trialInd)])

% select frequency from schedule and calculate duty cycle
ch.Frequency = stimSchedule(trialInd);
ch.DutyCycle = 5/(1000/ch.Frequency);

%%%%%%%%%% STIM %%%%%%%%%%%%%
s.startBackground

pause(stimLength)

s.stop
%%%%%%%%%% END STIM %%%%%%%%%%


% constant + random ITI 
pause(interTrialInterval + rand*randITI)


end

% save schedule and parameters 
save(['optoExp_' date])
