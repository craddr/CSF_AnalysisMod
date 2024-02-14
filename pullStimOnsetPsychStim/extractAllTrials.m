
function cutSignals = extractAllTrials(inputData,inputDataTimes,preEventTime,postEventTime,inputFreq,eventType,eventMatrix)

%%%%%% written by Adam Ranson 2014 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to extract snippets of ca signal or other signals at 1000Hz
% sampling rate using event start and stop times in the event matrix.

% Arguments:
% inputData         - Ca2 signal etc, in cell x time matrix
% preEventTime      - time before event that should be extracted
% postEventTime     - time after event that should be extracted
% zeroTrace         - whether trace is set to zero at time of event start

% determine the length of the trace to extract in samples
samplesPerEvent = ceil((preEventTime+postEventTime)*inputFreq)+1;
preEventSamples = ceil(preEventTime*inputFreq);
timeVector = linspace(-preEventTime,postEventTime,samplesPerEvent);

% reduce the event matrix down to the type of events which has been
% specified:
eventID = searchCellArray(eventMatrix.eventTypes,eventType);
validEvents = find(eventMatrix.events(:,3)==eventID);

cutSignals.trialProperties = eventMatrix.events(validEvents,:);
cutSignals.paramNames = eventMatrix.paramNames;
cutSignals.eventType = eventType;

%preallocate
for iCell = 1:size(inputData,1)
    cutSignals.signals{iCell} = zeros(length(validEvents),samplesPerEvent);
end

for iEvent = 1:length(validEvents)
    currentEvent = validEvents(iEvent);
    eventStartTime = eventMatrix.events(currentEvent,1)-preEventTime;
    eventEndTime   = eventMatrix.events(currentEvent,1)+postEventTime;
    signalStart  = find(inputDataTimes > eventStartTime ,1);
    signalEnd    = signalStart + samplesPerEvent-1;
	

    for iCell = 1:size(inputData,1)
        signalTrace = inputData(iCell,signalStart:signalEnd);
        cutSignals.signals{iCell}(iEvent,:)=signalTrace;
    end
end
  
cutSignals.timeVector = timeVector;

end






