function collateDataBatchRosie()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% code written by Rosie Craddock 2024 %%
%%%%% extraction of and smoothing of dF/F %%
%%%%% uses codes based on those written %%%%
%%%%% by Adam Ranson 2014 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function opens the calcium trace, timestamps with timeline time, and
% resamples to time vector input if requested

%%basically, we're rewriting the fluorescence data output by suite2p to
%%list:
%all fluorescence traces per potential cell
%recording depth (usually just the 1) ie how deep into the cortex were we
%looking?
%timestamps which correspond to when each fluorescence value was obtained
%size and location of each potential cell in terms of pixels
%wheel velocity at each timestamp


%%use GUI to select expDir (1 by 1 only)
% expdir=uigetdir(); %%
% expID= fullfile(expdir);


%%get animalID
[aniID]=GetAnimalID2();

%localReposPath= "G:\\FfionData_toFindAnalysis";
%localReposPath="E:\\FfionData_toFindAnalysis";
%localReposPath="E:\\MyData";
%%%change path here
localReposPath="E:\\MyData";
%%if local copy of folder exists, cd to there, otherwise stop and print error
if isfolder(fullfile(localReposPath, aniID))
	cd(fullfile(localReposPath,aniID)); 
	
else 
	disp("Error: Local copy of Animal data not found")
	clearvars;
	return;
end 

%%get Exps for this batch
[Exps]=WhichExps3(localReposPath, aniID);


     resampleFreq = 5;

%% split the Fall files.

%SplitFall(localReposPath, aniID,Exps); 

    
for j=1:length(Exps)
    exp=string(Exps{j});
    expID=fullfile(localReposPath, aniID, exp);
    

%     resampleFreq = 5;

    load(fullfile(expID, 'Events', 'psychStim.mat')); 
%     

%     outputTimes = Timeline.rawDAQTimestamps(1):1/resampleFreq:Timeline.rawDAQTimestamps(end);
% 
%     if exist(fullfile(expID,'ch2'))
%         % then there are 2 functional channels
%         dataPath{1} = fullfile(expID,'suite2p');
%         dataPath{2} = fullfile(expID,'ch2','suite2p');
%     else
%         dataPath{1} = fullfile(expID, 'suite2p');
%     end
% 
%     depthCount = length(dir(fullfile(dataPath{1},'*plane*')));
% 
%     if depthCount==1
%         % then we might be doing frame averaging
%         load(fullfile(expID,'tifHeader.mat'));
%         acqNumAveragedFrames = header.acqNumAveragedFrames;
%     else
%         % then we assume no averaging
%         acqNumAveragedFrames = 1;
%     end

    % determine which channel has frame timing pulses      
%     neuralFramesIdx  = find(ismember({Timeline.hw.inputs().name},'neuralFrames'));
        
    % divide the frame counter by the number of depths & averaging factor.
    %Timeline.rawDAQData(:,neuralFramesIdx)=ceil(Timeline.rawDAQData(:,neuralFramesIdx)/depthCount/acqNumAveragedFrames);
%     frameCount=ceil(Timeline.rawDAQData(:,neuralFramesIdx)/acqNumAveragedFrames);
%% change here as needed (check the xml file if you forgot what settings you used)
nFramesAveraged= 6; 
channels=2; 
%%channel used for green/fluorescence indicator
channelused=1;
totalNumberFrames=length(eventMatrix.frameTimesTS);
%%frame count is for the green channel, not for the red or for both
    frameCount=ceil(totalNumberFrames/nFramesAveraged);
    
%     frameCounter=[1:1:frameCount];
%     frameCounter=frameCounter';
    frameTimes=    eventMatrix.frameTimesTS;
    %% downsample added 24/07/2023
    frameTimes=frameTimes(1:nFramesAveraged:length(frameTimes));
    
    
    %% frame times where we're taking every 6th Frame
%    frameTimes= eventMatrix.frameTimesTS(1:nFramesAveraged:totalNumberFrames);
%outputTimes= frameTimes;
%lastSample=length(eventMatrix.frameTimesTS);
outputTimes=eventMatrix.counterTimerSampleFreq(1):1/resampleFreq:eventMatrix.counterTimerSampleFreq(end);
    % determine tlTime of each frame
    
    %frameTimes = Timeline.rawDAQTimestamps(diff(frameCount)==1);
 
 
 
 dataPath{1}= fullfile(expID, 'suite2p'); 

        
% for each channel combine all valid rois
%% for ease of converting old scripts to new ones: 

iCh=1;
iDepth=0;
 %depthCount = length(dir(fullfile(dataPath{1},'*plane*')));

        alldF{iCh} = [];
        allDepths{iCh} = [];
        cellCount = 0;
            % load s2p data
            load(fullfile(dataPath{iCh},['plane',num2str(iDepth)],'Fall.mat'));
            % load numpy file containing cell classification
            cellValid = readNPY(fullfile(dataPath{iCh},['plane',num2str(iDepth)],'iscell.npy'));
            % neuropil subtraction
            F = F - (Fneu*ops.neucoeff);
            % dF/F calculation
            smoothingWindowSize = 100;
			%%%% RC: 2D convolution
            smoothed = conv2(F,ones(1,smoothingWindowSize)/smoothingWindowSize,'same');
        % remove edge effects
                
            smoothed(:,1:smoothingWindowSize)=repmat(smoothed(:,smoothingWindowSize+1),[1,smoothingWindowSize]);
            smoothed(:,end-smoothingWindowSize+1:end)=repmat(smoothed(:,end-smoothingWindowSize-1),[1,smoothingWindowSize]);
            % replace nans with large values (so they don't get picked up as mins)
            smoothed(isnan(smoothed))=max(smoothed(:))*2;
            baseline = imerode(smoothed,strel('rectangle',[1 smoothingWindowSize]));
			%%% i think this means that the dF value is actual fluroesnce-the smoothed fluorescence for the area over a window size of 100 frames
			%%%TODO: find a better way of dF/F probably with a larger sampled time
            % calc dF/F
            dF = (F-baseline)./baseline;
            % get times of each frame
%             depthFrameTimes = frameTimes(iDepth+1:depthCount:length(frameTimes));
             depthFrameTimes = frameTimes(1:size(dF,2));


% depthFrameTimes=outputTimes;
%             % resample to get desired sampling rate
             dF = interp1(depthFrameTimes,dF',outputTimes)';
            if size(dF,2)==1
                dF = dF';
            end
            
            %dF=dF';
            % pick out valid cells
            alldF{iCh} = [alldF{iCh};dF(cellValid(:,1)==1,:)];
            allDepths{iCh} = [allDepths{iCh};repmat(iDepth,[sum(cellValid(:,1)),1])];
                
            % store masks of cells (for longitidinal tracking etc)
            cellMask = zeros(size(ops.meanImg));
            validCellList = find(cellValid(:,1)==1);
            for iCell = 1:length(validCellList)
                cellID = validCellList(iCell);
                % roiPix = sub2ind(size(cellMask),int64(stat{cellID}.ypix)+int64(ops.yrange(1))-1,int64(stat{cellID}.xpix)+int64(ops.xrange(1))-1);
                roiPix = sub2ind(size(cellMask),int64(stat{cellID}.ypix),int64(stat{cellID}.xpix));
                cellMask(roiPix) = iCell+cellCount;
            end
            % set offset so that cell numbers match to the number of
            % rows in the matrix with all cell dF/f0 extracted
            cellCount = iCell+cellCount;
            allMeanFrames{iCh,iDepth+1} = ops.meanImg;
            allROIMasks{iCh,iDepth+1}   = cellMask;
                
        
            %outputTimes=outputTimes-outputTimes(1);
            
    
        % save all data
    s2pData.alldF = alldF;
    s2pData.allDepths = allDepths;
    s2pData.t = outputTimes;
    s2pData.allMeanFrames = allMeanFrames;
    s2pData.allROIMasks = allROIMasks;
   % s2pData.wheelVelocity = pullWheelVelocity(Timeline,outputTimes);
    save(fullfile(expID, 's2pData.mat'),'s2pData');
%%legacy code for when wheelData was also stored

% function [wheelVelocity] = pullWheelVelocity(Timeline,outputTimes)
% % returns wheel velocity from raw timeline data
% rotaryEncoderIdx  = find(ismember({Timeline.hw.inputs().name},'rotaryEncoder'));
% if ~isempty(rotaryEncoderIdx)
%     rotaryEncoderData = Timeline.rawDAQData(:,rotaryEncoderIdx);
%     % make rotary encoder value sensible by removing wraparound
%     midBound = (2^32)/2;
%     rotaryEncoderData(rotaryEncoderData>midBound)=rotaryEncoderData(rotaryEncoderData>midBound) - 2^32;
%     rotaryEncoderData = smooth(rotaryEncoderData,200);
%     rotaryEncoderSpeedData = [diff(rotaryEncoderData);0];
%     rotaryEncoderTime = Timeline.rawDAQTimestamps;
% else
%     rotaryEncoderTime = Timeline.rawDAQTimestamps;
%     rotaryEncoderSpeedData=zeros(size(rotaryEncoderTime));
%     msgbox('Warning: No rotary encoder data found!');
%     disp('Warning: No rotary encoder data found!')
% end
% 
% wheelVelocity = interp1(rotaryEncoderTime,rotaryEncoderSpeedData,outputTimes);
% end


end
end