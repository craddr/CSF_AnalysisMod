function [Exps]=WhichExps3(localReposPath, aniID)

%%%% written by Rosie Craddock 2024, based on code written by Adam Ranson, 2014

AllExps=dir(); 

AllExps=struct2cell(AllExps);

for kk=3:length(AllExps(1,:))-1
    expName=AllExps{1,kk};  
    if isfolder(fullfile(localReposPath, aniID, expName, 'suite2p'))==1
        hasSuite2p{kk-2}=1;
    else 
        hasSuite2p{kk-2}=0;
        
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
ExperimentSummary{i-2}= string(expLog.expSummary(1,1));

end

ExpID=ExpID'; 
ExperimentSummary=ExperimentSummary';
hasSuite2p=hasSuite2p';

%%print exp log info as a table
ExpTable= table(ExpID, hasSuite2p, ExperimentSummary); 

ExpTable(:,:) 

listOfExperiments=ExpID;

[indx, tf]= listdlg('Name', 'Select ExpIDs you want to rewrite suite2p data for (all exps in current batch)', 'PromptString', {'Comments and summaries written about each', 'Experiment ID have been printed in the MATLAB', 'command window','Only works for experiments where basic grating stim', '(PsychStim) was used', 'AND for experiments which have been processed by suite2p'}, 'ListString', listOfExperiments, 'ListSize', [300, 300]);


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