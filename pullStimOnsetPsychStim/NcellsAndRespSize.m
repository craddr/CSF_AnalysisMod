%%%number of visually responsive cells 

%%%%% code written by RC 2024
function NcellsAndRespSize() 

[aniID]=GetAnimalID2();

%localReposPath= "E:\\FfionData_toFindAnalysis"; 
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

%%get Exps
[Exps]=WhichExps4(localReposPath, aniID);

for jj=1:length(Exps)
    exp=string(Exps{jj}); 
    expID=fullfile(localReposPath, aniID, exp);

%ch = [1 2];

load(fullfile(expID, 'Processed', 'Stim onset trials neural all 1.mat')); 

baselines= cutSignalsNeural.baselines; 
responses= cutSignalsNeural.responses; 

    for k= 1:max(cutSignalsNeural.trialProperties(:,28))
		visstim(k,:)=find(cutSignalsNeural.trialProperties(:,28) ==k); 
        contrast(k)= cutSignalsNeural.trialProperties(visstim(k,1),15);
        sf(k)=cutSignalsNeural.trialProperties(visstim(k,1),17);
    end
%%% find the trials where the contrast was ==1
    %contrast1=find(contrast==1|| contrast==0.5);
    contrast1=find(contrast==1);

visstimContrast1=visstim(contrast1, :);
sfList=sf(contrast1);
    
%%% each cell is 1 row, each trial is 1 column
for kk=1:height(responses)
for jk= 1:height(visstimContrast1)
    visstimContrast1ForStimx=visstimContrast1(jk,:);
    WilcoxonContrast1Resp{kk, jk}=ranksum(baselines(kk,visstimContrast1ForStimx), responses(kk,visstimContrast1ForStimx), 'tail', 'left');
    Increase{kk,jk}=responses(kk,visstimContrast1ForStimx)-baselines(kk,visstimContrast1ForStimx);
    IncreaseMean{kk,jk}=mean(responses(kk,visstimContrast1ForStimx)-baselines(kk,visstimContrast1ForStimx));
end
if any(cell2mat(WilcoxonContrast1Resp(kk,:))<=0.01)
    VisualResp{kk}=1; 
else 
    VisualResp{kk}=0; 
end 

if any(cell2mat(VisualResp(kk))==1)&any(cell2mat(Increase(kk, :))>=0.5)
    VisualResp2{kk}=1; 
else 
    VisualResp2{kk}=0;
end
NvisRespCells= sum(cell2mat(VisualResp(1,:)));
NvisRespCells2=sum(cell2mat(VisualResp2(1,:)));



%%%find where stimType ==x


end 

%%% so far we have found cells for which the dF/F is significantly
%%% different in baseline than for response period, and the dF/F is at
%%% least 0.5 higher in response than in baseline. 

%%%% for the visually responsive cells, find what is the highest SF they
%%%% respond to: 


for lm=1:height(WilcoxonContrast1Resp)
    %%% finding all visually responsive cells by significant difference
    %%% between baseline and response to 0.05
for ll=1:width(Increase)
    if cell2mat(WilcoxonContrast1Resp(lm,ll))<=0.01&any(cell2mat(Increase(lm, ll))>=0.5)



        allVisResp{lm, ll}=sfList(ll);
    else 
        allVisResp{lm, ll}=0;
    end 

highestSFIndex{lm}=find(cell2mat(allVisResp(lm, :))>0, 1, 'last');
highestSF{lm}=allVisResp(lm,cell2mat(highestSFIndex(lm)));

%%for highestSFIndex{lm} find all the trials for all contrasts tested, and
%%then find the lowest contrast for which a cell is visually responsive.
if (cell2mat(highestSF{1,lm}))>=0
    %% these are the stimulusNumbers which correspond to the highest SF responded to.
sfIndices{lm}=find(sf==(cell2mat(highestSF{1,lm})));

%%% use the stimulus numbers to find the contrasts for these stimulus IDs,
%%% then find all the trials for each of the contrasts, then test for each
%%% of these contrasts if the cell is responsive as done before.
end 
end


end 

end

