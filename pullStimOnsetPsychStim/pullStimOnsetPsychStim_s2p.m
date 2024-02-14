function pullStimOnsetPsychStim_s2p()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% code written by RC 2024%%%%%%%%%%%%%
%%%% based on code  written %%%%%%%%%%%%%
%%%% by Adam Ranson 2014 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extracts calcium data from protocols run with the vStimGui
% ch format = [1 2] or [1] or [2]
%% Set parameters (debug)
%expIDs = {{'2017-10-26_01_CFAP086'}};

% manually set red or green channel. Green = 1, Red = 2.
% ch = 2;

%% Get experimentdir manually.
% 
% expdir=uigetdir(); 
% expID= fullfile(expdir.name); 


%%get animalID
[aniID]=GetAnimalID2();

%localReposPath= "G:\\FfionData_toFindAnalysis"; 
%localReposPath= "E:\\FfionData_toFindAnalysis"; 
%%%%change path here
localReposPath="E:\\MyData";


%%if local copy of folder exists, cd to there, otherwise stop and print error
if isfolder(fullfile(localReposPath, aniID))
	cd(fullfile(localReposPath,aniID)); 
	
else 
	disp("Error: Local copy of Animal data not found")
	clearvars;
	return;
end 

%%get Exps
[Exps]=WhichExps4(localReposPath, aniID);

for jj=1:length(Exps)
    exp=string(Exps{jj}); 
    expID=fullfile(localReposPath, aniID, exp);

ch = [1 2];

%% Extract Data

%%%change preEventTime and postEventTime as appropriate
% parameters of trials to extract:
% eventName  = 'Stim onset';
% 
% 
% preEventTime = 2; %% =preeventtime from eventMatrix.events
% postEventTime = 5; %%= preeventtime+ drifttime+posteventime from eventMatrix.events
% baselineWindow = [-0.5 -0];
% responseWindow = [.2 1];
%         
%% Load signal traces and get into correct format:
%load event matrix
%if ~exist(fullfile(expID,'Events','psychStim.mat'))
 %   % if event file not generated then make it
%   eventGen.psychStimBasic(expID);
%end
        
load(fullfile(expID,'Events','psychStim.mat'));

load(fullfile(expID,'s2pData.mat'),'s2pData');
%neuralResampleFreq = 5;
%%find the sample frequency by finding the number of samples with
%%timestamp under 1s
%neuralResampleFrequency=ceil(find(eventMatrix.frameTimesTS<=1,1, 'last')/3);
%neuralResampleFreq=1;
%neuralResampleFreq=ceil(find(eventMatrix.frameTimesTS<=1, 1, 'last'));
%%% below works, but the data are not downsampled!!
%neuralResampleFreq=floor(find(s2pData.t<=1,1,'last'));
%%%find the first index where time is greater than 1 second and then divide
%%%length by 1 to get frequency. 
%%%%finished up to here.

neuralResampleFreq=5;

% eventName  = 'Stim onset';
% 
% 
% preEventTime = 2; %% =preeventtime from eventMatrix.events
% postEventTime = 5; %%= preeventtime+ drifttime+posteventime from eventMatrix.events
% baselineWindow = [-0.5 -0];
% responseWindow = [.2 1];
        
eventName=string(eventMatrix.eventTypes{1,1}); 

%%this refers to the drift time,... This needs to be included in the events
%%psychStim file at position 14

%%the iti/the pre event time
preEventTime=eventMatrix.events(1,29);
%postEventTime=eventMatrix.events(1,14)+ eventMatrix.events(1,15) + eventMatrix.events(1,16); 

%%this refers to the iti& poststationary time, which were at positions 15
%%and 16 in the event file, need to extract the iti from the psychstim file
% and get the post-stationary from the psych stim file too.
postEventTime= eventMatrix.events(1,25);
%%Rosie changed 25/07/2023
baselineWindow=[-1, 0]; 
responseWindow=[0.2, postEventTime];



% load timeline
%load(fullfile(expID,strcat(exp,'_Timeline.mat')));
% load expData
load(fullfile(expID,strcat(exp,'_psychstim.mat')));
stimParams = expData.stim.stims;
% load Ca2+ signals

        
        
for iCh = 1:length(s2pData.alldF)
            
outputFilename = ['Stim onset trials neural all ',num2str(ch(iCh)),'.mat'];
neuralData = s2pData.alldF{iCh};
neuralDataTime = s2pData.t;
            
% pull out signal snippets
cutSignalsNeural = extractAllTrials(neuralData,neuralDataTime,preEventTime,postEventTime,neuralResampleFreq,eventName,eventMatrix);
%cutSignalsWheel  = extractAllTrials(s2pData.wheelVelocity,neuralDataTime,preEventTime,postEventTime,neuralResampleFreq,eventName,eventMatrix);
%cutSignalsNeural.signalsWheel = cutSignalsWheel.signals;
% calculate response amplitude values for each neuron/trial
baselineSamples = find(cutSignalsNeural.timeVector>=baselineWindow(1),1);
baselineSamples = [baselineSamples find(cutSignalsNeural.timeVector>=baselineWindow(2),1)];
responseSamples = find(cutSignalsNeural.timeVector>=responseWindow(1),1);
responseSamples = [responseSamples find(cutSignalsNeural.timeVector>=responseWindow(2),1)];
            
			%%avg activity for each baseline and for each response period
			%%of a trial
for iCell = 1:length(cutSignalsNeural.signals)
    cutSignalsNeural.baselines(iCell,:) = mean(cutSignalsNeural.signals{iCell}(:,baselineSamples(1):baselineSamples(2)),2);
    cutSignalsNeural.responses(iCell,:) = mean(cutSignalsNeural.signals{iCell}(:,responseSamples(1):responseSamples(2)),2);
end
            
cutSignalsNeural.stimConfig.stimParams = stimParams;
cutSignalsNeural.stimConfig.stimParamNames = expData.stim.params;
            
% save data
if ~exist(fullfile(expID,'Processed'),'dir')
mkdir(fullfile(expID,'Processed'));
end
            
save(fullfile(expID,'Processed',outputFilename),'cutSignalsNeural');
            
end
end
end



