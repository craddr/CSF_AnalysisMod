%%%number of visually responsive cells 

function NcellsAndRespSize4() 

[aniID]=GetAnimalID2();

%localReposPath= "E:\\FfionData_toFindAnalysis"; 
localReposPath="G:\\MyData";


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
        %visstim: each row is for one trial type, each column is a repeat
		visstim(k,:)=find(cutSignalsNeural.trialProperties(:,28) ==k); 
        %%just a list of the contrasts
        contrast(k)= cutSignalsNeural.trialProperties(visstim(k,1),15);
        sf(k)=cutSignalsNeural.trialProperties(visstim(k,1),17);     
    end
    
            %%change where sf was input at 0.126 instead of 0.128
        IndicesForWrongSF=find(sf==0.126);
        sf(IndicesForWrongSF)=0.128;
%%% find the trials where the contrast was ==1
    %contrast1=find(contrast==1|| contrast==0.5);
    contrast1=find(contrast==1);
    
    contrastHalf=find(contrast==0.5);

visstimContrast1=visstim(contrast1, :);
visstimContrastHalf=visstim(contrastHalf,:);
sfListC1=sf(contrast1);
sfListChalf=sf(contrastHalf);
    
%%% each cell is 1 row, each trial is 1 column

%%kk=cell number
%%jk = number of spatial frequencies
for kk=1:height(responses)
for jk= 1:height(visstimContrast1)
    visstimContrast1ForStimx=visstimContrast1(jk,:);
    visstimContrastHalfForStimx=visstimContrastHalf(jk,:);
    WilcoxonContrast1Resp{kk, jk}=ranksum(baselines(kk,visstimContrast1ForStimx), responses(kk,visstimContrast1ForStimx), 'tail', 'left');
    WilcoxonContrastHalfResp{kk, jk}=ranksum(baselines(kk,visstimContrastHalfForStimx), responses(kk,visstimContrastHalfForStimx), 'tail', 'left');
    IncreaseC1{kk,jk}=responses(kk,visstimContrast1ForStimx)-baselines(kk,visstimContrast1ForStimx);
    IncreaseChalf{kk,jk}=responses(kk,visstimContrastHalfForStimx)-baselines(kk,visstimContrastHalfForStimx);
    IncreaseMeanC1{kk,jk}=mean(responses(kk,visstimContrast1ForStimx)-baselines(kk,visstimContrast1ForStimx));
    IncreaseMeanChalf{kk,jk}=mean(responses(kk,visstimContrastHalfForStimx)-baselines(kk,visstimContrastHalfForStimx));
    
    %% tag cells which respond to any full contrast or half contrast stimuli
    %% criteria for a responsive cell is that the stat test gives a P value for difference between baselines and responses of below 0.01, and the mean dF/F for the response period must be 0.3 greater than for baseline
    if cell2mat(WilcoxonContrast1Resp(kk,jk))<=0.01&cell2mat(IncreaseMeanC1(kk,jk))>=0.3
        VisualResp{kk, jk}=1;
        VisualResp2{kk,jk}=1;
    else
        VisualResp{kk,jk}=0;
        if cell2mat(WilcoxonContrastHalfResp(kk,jk))<=0.01&cell2mat(IncreaseMeanChalf(kk,jk))>=0.3
            VisualResp2{kk,jk}=1;
        else 
            VisualResp2{kk,jk}=0;
        end
    end
    

%% tag cells which respond to any full contrast stimuli


%%at this point, VisualResp2 contains an array of numberSFs testedxncells,
%%with 1s and 0s in each slot, where 1 means that the cell responded to
%%that SF, and 0 means it did not

if any(cell2mat(VisualResp2(kk,:))==1)
    VisRespList{kk}=1;
else 
    VisRespList{kk}=0;
end



end 
end
% NvisRespCells= sum(cell2mat(VisualResp(1,:)));
% NvisRespCells2=sum(cell2mat(VisualResp2(1,:)));
NvisRespCells=sum(cell2mat(VisRespList(1,:)));


    %%% only take into account where the contrast ==1, so the first 7,
    %%% rows, then compare the 5 values in the row against the 5 values in
    %%% the baseline row? see what we did for the HVA paper.

%%%find where stimType ==x



%%% so far we have found cells for which the dF/F is significantly
%%% different in baseline than for response period, and the dF/F is at
%%% least 0.5 higher in response than in baseline. 

%%%% for the visually responsive cells, find what is the highest SF they
%%%% respond to: 
%% a vector of indices listing cells which respond to at least one stimulus
VisResp2Indices=find(cell2mat(VisRespList)==1);

%% lm= number of cells
for lm=1:height(WilcoxonContrast1Resp)
    %%% finding all visually responsive cells by significant difference
    %%% between baseline and response to 0.3
    
    %% ll = number of SFs
for ll=1:width(IncreaseMeanC1)
    if cell2mat(VisualResp2(lm,ll))==1
    %if cell2mat(WilcoxonContrast1Resp(lm,ll))<=0.01&cell2mat(IncreaseMeanC1(lm, ll))>=0.3



        allVisResp{lm, ll}=sfListC1(ll);
        allVisRespReplacezeros{lm,ll}=sfListC1(ll);
    else 
        allVisResp{lm, ll}=0;
        allVisRespReplacezeros{lm,ll}=[];
    end 
end
    highestSF{lm}=max(cell2mat(allVisResp(lm,:)));
%%for highestSFIndex{lm} find all the trials for all contrasts tested, and
%%then find the lowest contrast for which a cell is visually responsive.
if (cell2mat(highestSF(1,lm)))>=0
    %% these are the stimulusNumbers which correspond to the highest SF responded to.
highestsfIndices{lm}=find(sf==(cell2mat(highestSF(1,lm))));

%% all SFs the cell responds to
SFsresponsiveTo{lm}=unique(cell2mat(allVisRespReplacezeros(lm,:)));

%%% use the stimulus numbers to find the contrasts for these stimulus IDs,
%%% then find all the trials for each of the contrasts, then test for each
%%% of these contrasts if the cell is responsive as done before.

end
end

%%for all visually responsive cells, get the indices of trials for stimuli
%%being of the highest SF that the cell responds to, get these trial
%%indices in a vector form
 nContrasts=length(unique(contrast));
 nRepeats=width(baselines)/width(contrast);
% 
 nStimToTest=nContrasts*nRepeats;

for hg=1:length(VisResp2Indices)
    %% highestsfIndicesForThisCell is actually only the indices for the highest
        highestsfIndicesForThisCell=highestsfIndices(:,VisResp2Indices(hg));
        visStimForMaxSF=visstim(cell2mat(highestsfIndicesForThisCell),:);

%         visStimForMaxSFVecAll(hg, :)=[visStimForMaxSFVector(1,:)];
visStimForMaxSFVArrAll{hg}=visStimForMaxSF;

end 

%%for each cell
respCellcounter=1;
for ln=1:height(baselines)
    %%if the cell is visually responsive
    if any(VisResp2Indices==ln)
        %%%then for each contrast, 
                    %%find the indices for the trials of the highest SF that this
            %%cell responds to
            highestSFforTHIScell=visStimForMaxSFVArrAll(respCellcounter);
            for jk=1:nContrasts
                %%trials for this contrast:
                trialsForThisContrast=highestSFforTHIScell{1,1}(jk,:);

                WilcoxonAllContrastMaxSF{ln, jk}=ranksum(baselines(ln,trialsForThisContrast), responses(ln,trialsForThisContrast), 'tail', 'left');
                IncreaseMeanAllContrastsMaxSF{ln,jk}=mean(responses(ln,trialsForThisContrast)-baselines(ln, trialsForThisContrast));
                if cell2mat(WilcoxonAllContrastMaxSF(ln,jk))<=0.01&cell2mat(IncreaseMeanAllContrastsMaxSF(ln,jk))>=0.3
                    VisualRespMaxSFByContrast{ln, jk}=1;
                else
                    VisualRespMaxSFByContrast{ln,jk}=0;
                end
            end
            respCellcounter=respCellcounter+1;
    
    else 
         for jk=1:nContrasts
            WilcoxonAllContrastMaxSF{ln,jk}=[];
            IncreaseMeanAllContrastsMaxSF{ln,jk}=[];
            VisualRespMaxSFByContrast{ln,jk}=0;
         end
     end 

end 

%%%find the absoluteHighestResponder: 
highestSFofAllCells=max(cell2mat(highestSF));
indexCellHighestSFresponder=find(cell2mat(highestSF)==highestSFofAllCells);

if length(indexCellHighestSFresponder)==1
    %%%find the lowest contrast responded to for that responder
    highestSF_respectiveContrast= find(cell2mat(VisualRespMaxSFByContrast(indexCellHighestSFresponder,:))==1,1,'last');
else 
    nCellsHighestSFresponder=length(indexCellHighestSFresponder);
    highestSF_respectiveContrast{nCellsHighestSFresponder}=find(cell2mat(VisualRespMaxSFByContrast(indexCellHighestSFresponder,:))==1,1,'last');
end


% highestSFResults= zeros(length(highestSF_respectiveContrast),2);
% highestSFResults(:,1)=highestSF_respectiveContrast; 
% highestSFResults(:,2)=highestSFofAllCells;
if ~exist(fullfile(expID,'ExtractedResults'),'dir')
    mkdir(fullfile(expID,'ExtractedResults'));
end

% finalResultsHighestSF= horzcat(highestSFResults, NvisRespCells);
colnames={ 'highestSFresponded to','Contrast for highest SF responded to', 'N visRespCells'};
resultsTableHighestSF=table(highestSFofAllCells, highestSF_respectiveContrast, NvisRespCells,'VariableNames', colnames);
writetable(resultsTableHighestSF, fullfile(expID, 'ExtractedResults', 'resultsHighestSFandContrastWithNVisRespCells0_05.csv'));
sfsTested=unique(sf);

listOfContrasts=unique(contrast, 'stable');

for gg=1:length(sfsTested)
    %%find the indicies for each  sf for trial types, put this in a nice format
    sfIndices{gg}=find(sf==sfsTested(gg));
    vectorListOfContrasts{1,gg}=listOfContrasts;
    for hj=1:NvisRespCells
        cellNumber=VisResp2Indices(1,hj);
        if any(cell2mat(SFsresponsiveTo(1,cellNumber))==sfsTested(gg))
        
            for ff=1:nContrasts
                
                

                %%visStimIndices
                %%get the trial numbers for all repeats for each contrast specified by (ff) of
                %%stimuli of SF specified by gg
%                 trialNumber=sfIndices{1,gg}(1,ff);
%                 visStimTrialIndicesEachSF{gg,ff}= visstim(trialNumber,:);
%                 visStimTrialIndicesThisTrialType=cell2mat(visStimTrialIndicesEachSF(gg,ff));
%                 WilcoxonAllContrastsTestedPerSF{hj,1}(gg,ff)=ranksum(baselines(cellNumber,visStimTrialIndicesThisTrialType), responses(cellNumber, visStimTrialIndicesThisTrialType), 'tail', 'left');
%                 IncreaseMeanAllContrastsTestedPerSF{hj,1}(gg,ff)=mean(responses(cellNumber,visStimTrialIndicesThisTrialType)-baselines(cellNumber, visStimTrialIndicesThisTrialType));
%                 
                trialNumber=sfIndices{1,gg}(1,ff);
                visStimTrialIndicesEachSF{gg,ff}= visstim(trialNumber,:);
                visStimTrialIndicesThisTrialType=cell2mat(visStimTrialIndicesEachSF(gg,ff));
                WilcoxonAllContrastsTestedPerSF{hj,1}(ff,gg)=ranksum(baselines(cellNumber,visStimTrialIndicesThisTrialType), responses(cellNumber, visStimTrialIndicesThisTrialType), 'tail', 'left');
                IncreaseMeanAllContrastsTestedPerSF{hj,1}(ff,gg)=mean(responses(cellNumber,visStimTrialIndicesThisTrialType)-baselines(cellNumber, visStimTrialIndicesThisTrialType));
                
                if WilcoxonAllContrastsTestedPerSF{hj,1}(ff,gg)<=0.01&IncreaseMeanAllContrastsTestedPerSF{hj,1}(ff,gg)>=0.3
                    VisRespCellsResponseArray{hj,1}(ff,gg)=1; 
                    VisRespCellsArrayFormatted{ff,gg}(1,hj)=1;
                    VisRespCellsArrayFormattedMagnitude{ff,gg}(1,hj)=IncreaseMeanAllContrastsTestedPerSF{hj,1}(ff,gg);
                else 
                    VisRespCellsResponseArray{hj,1}(ff,gg)=0;
                    VisRespCellsArrayFormatted{ff,gg}(1,hj)=0;
                    VisRespCellsArrayFormattedMagnitude{ff,gg}(1,hj)=NaN;
                end
            end
        else
            for ff=1:nContrasts
                
                WilcoxonAllContrastsTestedPerSF{hj,1}(ff,gg)=1;
                IncreaseMeanAllContrastsTestedPerSF{hj,1}(ff,gg)=0;
                VisRespCellsResponseArray{hj,1}(ff,gg)=0;
                VisRespCellsArrayFormatted{ff,gg}(1,hj)=0;
                VisRespCellsArrayFormattedMagnitude{ff,gg}(1,hj)=NaN;
            end 
        end
    end
    

end

%%find CSF taking into account the responses of all cells:

%%%find the minimum
%%%contrast responded to for each SF for ANY cell:

for gg=1:length(sfsTested)
    for ff=1:nContrasts
        CSF_arrayMax{ff,gg}=max(VisRespCellsArrayFormatted{ff,gg});
        CSF_arrayMeanMag{ff,gg}=mean(VisRespCellsArrayFormattedMagnitude{ff,gg}, "omitnan");
        
    end
    CSF_IndicesForContrasts{1,gg}=find(cell2mat(CSF_arrayMax(:,gg))==1,1,'last');
    CSF_arrayContrasts{1,gg}=vectorListOfContrasts{1,gg}(cell2mat(CSF_IndicesForContrasts(1,gg)));
    for hg=1:NvisRespCells
        CSF_IndicesForContrastsAll{hg,gg}=find(VisRespCellsResponseArray{hg,1}(:,gg)==1,1,'last');
        CSF_arrayContrastsAllCells{hg,gg}=vectorListOfContrasts{1,gg}(cell2mat(CSF_IndicesForContrastsAll(hg,gg)));
    end
end 

varNames={strcat('SF:',num2str(sf(1,1))),...
    strcat('SF:',num2str(sf(1,2))),...
    strcat('SF:',num2str(sf(1,3))),...
    strcat('SF:',num2str(sf(1,4))),...
    strcat('SF:',num2str(sf(1,5))),...
    strcat('SF:',num2str(sf(1,6))),...
    strcat('SF:',num2str(sf(1,7)))};
rowNames={strcat('Contrast:',num2str(contrast(1,1))),...
    strcat('Contrast:',num2str(contrast(1,8))),...
    strcat('Contrast:',num2str(contrast(1,15))),...
    strcat('Contrast:',num2str(contrast(1,22))),...
    strcat('Contrast:',num2str(contrast(1,29))),...
    strcat('Contrast:',num2str(contrast(1,36)))};
%%cet the minContrast required in a form suitable to use in R for plotting
for hj=1:NvisRespCells
    SFindices=[1:1:length(sfsTested)];
    IndicesForThisCell=[(hj*SFindices(end)-length(SFindices)+1):1:(hj*SFindices(end))];
    for gg=1:length(sfsTested)
        FormattedForR_minContrastArray{IndicesForThisCell(gg),1}=aniID;
        FormattedForR_minContrastArray{IndicesForThisCell(gg),2}=hj;
        FormattedForR_minContrastArray{IndicesForThisCell(gg),3}=sfsTested(gg);
        FormattedForR_minContrastArray{IndicesForThisCell(gg),4}=cell2mat(CSF_arrayContrastsAllCells(hj,gg));
    end
end

    

CSFtable_max=cell2table(CSF_arrayMax, 'VariableNames', varNames, 'RowNames', rowNames);
CSFtable_meanMag=cell2table(CSF_arrayMeanMag, 'VariableNames', varNames, 'RowNames', rowNames);
CSFtable_ContrastsVector=cell2table(CSF_arrayContrasts, 'VariableNames', varNames);
CSFtable_ContrastsAll=cell2table(CSF_arrayContrastsAllCells, 'VariableNames', varNames);
CSFtable_ContrastsAllCellsFormattedForR=cell2table(FormattedForR_minContrastArray, 'VariableNames', {'animalID','cellNumber','SF','MinContrast'});


writetable(CSFtable_max, fullfile(expID, 'ExtractedResults','ConstrastSensitivityFunction_maxResponses0_05.csv'), 'WriteRowNames', true);
writetable(CSFtable_meanMag,fullfile(expID, 'ExtractedResults','ConstrastSensitivityFunction_meanMagnitude0_05.csv'), 'WriteRowNames',true);
writetable(CSFtable_ContrastsVector,fullfile(expID, 'ExtractedResults','ConstrastSensitivityFunction_lowestContrast0_05.csv'));
writetable(CSFtable_ContrastsAll,fullfile(expID,'ExtractedResults', 'MinContrastPerSF_allCells0_05.csv'));
writetable(CSFtable_ContrastsAllCellsFormattedForR, fullfile(expID,'ExtractedResults', 'MinContrastPerSF_allCells_FormattedForR0_05.csv'));

%VisRespCellPlot(cutSignalsNeural,NvisRespCells,VisResp2Indices,expID, exp);


%%%get Data for plot

end 

end

