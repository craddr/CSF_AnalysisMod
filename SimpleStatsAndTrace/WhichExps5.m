function [Exps]=WhichExps5(localReposPath, aniID)
%%%% written by Rosie Craddock 2024, based on code written by Adam Ranson, 2014


AllExps=dir(); 
AllExps=struct2cell(AllExps);


for kk=3:length(AllExps(1,:))-1
    expName=AllExps{1,kk};  
    if isfile(fullfile(localReposPath, aniID, expName, 'Processed','Stim onset trials neural all 1.mat'))==1
        hasProcessed{kk-2}=1;
    else 
       hasProcessed{kk-2}=0;
        
    end 
end 

cd("log"); 
        

AllExpIDs=dir(); 

AllExpIDs=struct2cell(AllExpIDs);

%ExpLogs= zeros(length(AllExpIDs), 3); 

%%get log data including experiment comments

for i = 3:length(AllExpIDs(1,:))
    ThisExpID=AllExpIDs{1,i};

load(fullfile(ThisExpID)); 

%ExpLogs[i, 1]= expLog.expID; 
%ExpLogs[i, 2]= expLog.comments[1]; 
%ExpLogs[i,3]=expLog.expSummary[1];

ExpID{i-2}= ThisExpID(1:end-4); 
%ExperimentComments{i-2}= expLog.comments;
ExperimentSummary{i-2}= string(expLog.expSummary(1,1));

end

ExpID=ExpID'; 
%ExperimentComments=ExperimentComments';
ExperimentSummary=ExperimentSummary';
hasProcessed=hasProcessed';
%%print exp log info as a table
ExpTable= table(ExpID, hasProcessed,  ExperimentSummary); 

ExpTable(:,:) 

listOfExperiments=ExpID;

[indx, tf]= listdlg('Name', 'Select ExpIDs you want to analyse and plot', 'PromptString', {'Comments and summaries written about each', 'Experiment ID have been printed in the MATLAB', 'command window','Only works for experiments where basic grating stim', '(PsychStim) was used', 'AND where data has been fully processed'}, 'ListString', listOfExperiments, 'ListSize', [300, 300]);


% if user selected experiment files, Exps= those selected, otherwise, Exps=0; 
if tf==1 
for jj=1:length(indx)
    Exps{jj}=ExpID(indx(jj));
end
else 
Exps= 0; 

disp("No experiments selected")
clearvars; 
return;

end 



end