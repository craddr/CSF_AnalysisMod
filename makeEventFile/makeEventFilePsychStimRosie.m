function makeEventFilePsychStimRosie()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% written by RC 2024 %%%%%%%%%%%%%%%
%%%% methods for extraction of stim %%%
%%%% data based on codes of Adam %%%%%%
%%%% Ranson, 2014 %%%%%%%%%%%%%%%%%%%%%
%%%% all other codes RC original %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% select Animal whose experiments you want to create event files for
[aniID] =GetAnimalID2(); 
%%%%change path here
%localReposPath= "E:\\MyData"; 
localReposPath= "E:\\MyData"; 
%localReposPath="E:\\FfionData_toFindAnalysis";
%%if local copy of folder exists, cd to there, otherwise stop and print error
if isfolder(fullfile(localReposPath, aniID))
	cd(fullfile(localReposPath,aniID)); 
	
else 
	disp("Error: Local copy of Animal data not found")
	clearvars;
	return;
end 



%% select Experiments from local repository to create event files for
[Exps]= WhichExps2_New();

for j= 1:length(Exps) 
	%%get exp name
	expID=string(Exps{j}); 
	%%load the associated timeline file (holding timestamp data)
	%load(fullfile(localReposPath, aniID, expID, strcat(expID,'_Timeline.mat'))); 
    
         d= struct2cell(dir(fullfile(localReposPath, aniID, expID))); 
     d_timeline=d(1,4);
     %%change directory to the file containing the h5 file and then
     %%manually load this file
%      cd(fullfile(localReposPath, aniID, expID, d_timeline{1,1})); 
%      [FrameCounter, time]=LoadSyncEpisode;

%%try to do this using base codes instead of ThorSync Codes 

h5_info= h5info(fullfile(localReposPath, aniID, expID, d_timeline{1,1}, "Episode001.h5"));
clockRate = 20000000;
pathname=strcat(fullfile(localReposPath, aniID, expID, d_timeline{1,1}), '\');
sampleRate = LoadSyncXML(pathname);

props = {'start','length','interval'};
data = {[1,1],[1 Inf],[1 1]};

% if(~isempty(varargin))
%     assert(rem(length(varargin),2)==0 && iscellstr(varargin(1:2:end)), 'Inputs failed to conform to expected string-value pair format, eg> ''start'',1');
%     %foundProps = intersect(varargin(1:2:end), props); 
%     IdxCell = cellfun(@(x) strcmpi(x,props),varargin(1:2:end),'UniformOutput',false);
%     val = double(cell2mat(varargin(2:2:end)))*sampleRate;
%     for i=1:length(val)
%         data{cell2mat(IdxCell(i))>0} = [1 val(i)];
%     end
% end
time1= h5read(fullfile(localReposPath, aniID, expID, d_timeline{1,1}, "Episode001.h5"), '/Global/GCtr', data{1}, data{2}, data{3})';
time=double(time1)./clockRate;

%%get the frameCounterData
FrameCounter=h5read(fullfile(localReposPath, aniID, expID, d_timeline{1,1}, "Episode001.h5"), '/CI/FrameCounter', data{1}, data{2}, data{3})';

     %% get the timedata in a nice, usable format
     
     %% get the timestamp for the first sample of each new frame, by indices then times
         [C,ia,ic]=unique(FrameCounter, 'rows');
         %%indices are held in ia
         %% get rid of index for 0 frame count 
       %  frameIndices= ia(2:length(ia));
         %%get corresponding timestamps (wrt the ThorSync Timestamp)
        % frameTimesTS= time(frameIndices);
     	%% load the associated psychstim file (holding parameters of stimuli shown during experiment)
	load(fullfile(localReposPath, aniID, expID, strcat(expID,'_psychstim.mat'))); 
     %%get the frame indices for frames collected by the VS PC 
     
     %%now, get the TS data correct to the VS PC NI counter, so that the time=0, means the time of the first 2p frame recorded: 
  %   frameTimesTS=frameTimesTS-frameTimesTS(1);
     %% change this in Ffion's code also...
     %%% instead, we should change this so we've got all the timestamps,
     %%% but we 0 at frame=1, then we add expData.neuralFrames Data final time to to this time,
     %%% then we get rid of all the samples which are below t=0, then we
     %%% can take start2p away from timestamp data
     %%% should have our timestamp data zeroed to the neural Frames data

%      frameCounter=ic(ia(1):ia(end));
%      CounterTimer=time(ia(1):ia(end));
     %%% check the below is correct next week
     frameCounter=ic-1; 
     CounterTimer=double(time); 
     %% find the index of the first frame andn the counter timer for this frame.
     FirstFrameIndex=find(frameCounter>=1,1); 
     TimeAtFrame1=CounterTimer(FirstFrameIndex);
     %TimeAtLastFrameVS=seconds(max(expData.neuralFramesData.Time));
     %% find the index of the last frame given in VS, and the time of this last frame according to the VS
     LastFrameVS=max(expData.neuralFramesData.Dev1_ctr0);
     IndexLastFrameVS=find(expData.neuralFramesData.Dev1_ctr0==LastFrameVS,1);
     TimeLastFrameVS=seconds(expData.neuralFramesData.Time(IndexLastFrameVS));
     %%find the index for the frame which was the last frame collected by
     %%THORLABS thorsync.
     %%find the index of the last VS frame on the frame counter from 2P PC
     %%and the time for this frame
     LastVSFrameIndexOnTL=find(frameCounter>=LastFrameVS, 1);
     TimeLastFrameTL=CounterTimer(LastVSFrameIndexOnTL);
     %% change the counter timer so that the timestamps are around the same as those used for VS PC counter
     CounterTimer=CounterTimer-TimeLastFrameTL+TimeLastFrameVS;
     
	 %% get the timer data and then put it wrt the first frame collected (first frame=0), 
	 %%then put all times up by whatever the VS 
%      CounterTimer=CounterTimer-TimerLastFrame+TimeLastFrameTL;
%      
     
     %%testing
     CounterTimerLastFrameNew=CounterTimer(LastVSFrameIndexOnTL);

     
%      FirstSupZero= find(CounterTimer>=0,1); 
%      frameCounter=ic(FirstSupZero:ia(end)); 
%      CounterTimer=CounterTimer(FirstSupZero:ia(end));

              %frameIndices= ia(1:length(ia));
              
              frameIndices=diff(frameCounter);
              frameIndices=find(frameIndices==1);
         %%get corresponding timestamps (wrt the ThorSync Timestamp)
		 
		 %frameIndices=frameIndices-FirstSupZero;
		 %%get the frame indices correct to the new CounterTimer
         frameTimesTS= CounterTimer(frameIndices+1);
          %frameTimesTS=frameTimesTS-frameTimesTS(1);
     %%%% now we have the timers correct to the psych stim timers, and now
     %%%% we should, empirically, know what 0 in CounterTimer means, as it
     %%%% should EMPIRICALLY be at  the same point in time as when start2p
     %%%% command is given, give or take 0.2s
     
     
     
	% expData.startNeuralFrames=expData.startNeuralFrames+8;
	%%%try 4,3,2 for CFR030 first experiment
%  expData.startNeuralFrames=expData.startNeuralFrames+4;
	%write an empty matrix (later used to store a list of all stimulus events and properties for stimuli shown during experiment) 
	eventMatrix.events=[];
    %%offset all the timings by the checkFrames1 timestamp minus the Time
    %%at the last frame VS, to get everything correct to the VS counter
    %%times.
    timeOffset=expData.timeCheckFrames2-TimeLastFrameVS;
	%%for each stimulus presentation
	for iStim = 1:length(expData.trialData)
    expData.trialData{1,iStim}.timing.StartTrial= expData.trialData{1,iStim}.timing.StartTrial-timeOffset; 
    expData.trialData{1,iStim}.timing.StartIti= expData.trialData{1,iStim}.timing.StartIti-timeOffset;
    expData.trialData{1,iStim}.timing.EndIti= expData.trialData{1,iStim}.timing.EndIti-timeOffset;
    expData.trialData{1,iStim}.timing.StimulusStart= expData.trialData{1,iStim}.timing.StimulusStart-timeOffset;
    expData.trialData{1,iStim}.timing.StimulusEnd= expData.trialData{1,iStim}.timing.StimulusEnd-timeOffset;
	%%find the time the stim GO command was given (the timepoint where the stimulus came on, in s)
	%stimCommandTime=expData.trialData{iStim}.timing.StimulusStart(1);
	%stimCommandTime=expData.trialData{iStim}.timing.StartTrial(1);
	stimCommandTime=expData.trialData{iStim}.timing.StimulusStart(1);
    stimGoTime=stimCommandTime;
	%%set the stimulus stop time as =(stimGoTime+1)s
	stimStopTime= stimGoTime+1;
	%%get the index of the stimulus from the list of stimuli 
	stimulusID=expData.trialData{iStim}.stimID; 
	%%use index to find stimulus properties
	%(xpos, ypos, xsize, ysize, ori, sf, tf, duty, shape, prestationary, drifttime, poststationary, iti)
	stimulusParams= expData.stim.stims(stimulusID,:); 
	%%write vector containing the following
	newRow= [stimGoTime, stimStopTime, 1, stimulusParams, stimulusID];
	
	%%write this vector into the eventMatrix
	eventMatrix.events= [eventMatrix.events;newRow]; 
	
	end
	%%find the number of events (the number of trials run during the experiment)
	eventMatrix.eventCount=size(eventMatrix.events,1); 
	%%set all stimulus events as being type= onset
	eventMatrix.eventTypes={'Stim onset'}; 
	%%write list of parameter names into the matrix
	eventMatrix.paramNames= {'starttime','endtime','eventtype',expData.stim.params{:},'stimID', 'iti'}'; 
	%%write again as paramDefs (can be removed if no one uses this) 
	eventMatrix.paramDefs={'starttime','endtime','eventtype',expData.stim.params{:},'stimID', 'iti'}';
	
	eventMatrix.frameTimesTS=frameTimesTS; 
    eventMatrix.frameCounterSampleFreq=frameCounter; 
    eventMatrix.counterTimerSampleFreq=CounterTimer;
    eventMatrix.events(:,29)=3;
	eventDirectory= fullfile(localReposPath, aniID, expID, 'Events');
	%%if the events directory doesn't exist, make it 
	if ~isfolder(eventDirectory)
		mkdir(eventDirectory); 
	end
	%%filename to save eventMatrix to
	eventFile= fullfile(localReposPath, aniID, expID, 'Events', 'psychStim.mat'); 
	
	save(eventFile, 'eventMatrix'); 
    
    
    
    fprintf('Event File Created for %s\n', expID)
	
	end 
clearvars; 
end
