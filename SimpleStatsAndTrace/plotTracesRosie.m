function plotTracesRosie() 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% written by RC 2024 %%%%%%
%%%% largely based on %%%%%%%%
%%%% code written by %%%%%%%%%
%%%% Asta Vasalauskaite %%%%%%
%%%% and Adam Ranson 2014 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%select the animalID
% animaldir=uigetdir; 
% 
% cd (animaldir); 
% expfolds=dir(animaldir);

%%get animalID
[aniID]=GetAnimalID2();

% localReposPath= "E:\\Local_Repository"; 
%localReposPath= "G:\\FfionData_toFindAnalysis";
%localReposPath= "E:\\FfionData_toFindAnalysis";
%%%change path here
localReposPath="E:\\MyData";
%%if local copy of folder exists, cd to there, otherwise stop and print error
if isfolder(fullfile(localReposPath, aniID))==1
	cd(fullfile(localReposPath,aniID)); 
	
else 
	disp("Error: Local copy of Animal data not found")
	clearvars;
	return;
end 

%%get Exps
[Exps]=WhichExps5(localReposPath, aniID);


%%for all experiments in the folder
for i=1:length(Exps)
%%find exp folder and store to memory
    exp_ID=string(Exps(i));
    current_exp=fullfile(localReposPath,aniID,exp_ID);
%%load processed data
	processed_file_path=fullfile(current_exp,'Processed'); 

    load(fullfile(processed_file_path, 'Stim onset trials neural all 1.mat')); 
    %%make folder for traces
	cd(current_exp);
    
	mkdir Traces; 
    cd Traces;


    
    
    %first_rec=-(cutSignalsNeural.trialProperties(1, 29)); 
    first_rec=-1;
	%first_rec=-3;
    %%Rosie changed 25/07/2023
    first_recording_used=first_rec;
    
    stim_START=0; 
    
    stim_START_used=0.2; 
    last_rec=(cutSignalsNeural.trialProperties(1,21));
    last_rec_used=last_rec;
    %#############################################
	cellN=length(cutSignalsNeural.signals); 
    %%assign field types for values to hold in 'Timings' struct
    field1='Uncut';
    field2='cutting';
    %%find positions where the below are true
	Uncut=find(cutSignalsNeural.timeVector>=first_rec & cutSignalsNeural.timeVector<=last_rec);
    cutting=find(cutSignalsNeural.timeVector>=first_recording_used & cutSignalsNeural.timeVector<=last_rec_used);
    %%write list of positions into struct, using field allocations
    %%previously specified
	

	
    %%for each stimulus type, find trials
    for k= 1:max(cutSignalsNeural.trialProperties(:,28))
		visstim(k,:)=find(cutSignalsNeural.trialProperties(:,28) ==k); 
    end

	stimONline=stim_START;
	stimRESPONSEline=stim_START_used;
	figure1=figure('visible','off');
    
    nStims=length(visstim);

    indicesForFirstRepeat=visstim(:,1)';
    visstimparams=cutSignalsNeural.trialProperties(indicesForFirstRepeat,:);
	
    %%plot for each cell
	for cEll= 1:cellN
        %%plot for each stimulus
		for k=1:max(cutSignalsNeural.trialProperties(:,28))
            %%stimulus title
            
            %%change subtitles
            %%here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ie between the %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%             xposypos=strcat('X position: ',string(cutSignalsNeural.trialProperties(k,11)), ', Y position: ', string(cutSignalsNeural.trialProperties(k,12))); 
%             oriVal= strcat(' Orientation: ',string(cutSignalsNeural.trialProperties(k,16)));
%            contrast=strcat('C: ', string(cutSignalsNeural.trialProperties(k,15))); 
%            spatialF= strcat('SF: ', string(cutSignalsNeural.trialProperties(k,17)));

            contrast=strcat('C: ', string(visstimparams(k,15))); 
            spatialF= strcat('SF: ', string(visstimparams(k,17)));
%             tt{k}= sprintf('%s%s', xposypos, oriVal); 
            tt{k}= sprintf('%s%s', contrast, spatialF); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x=cutSignalsNeural.timeVector(1,cutting);
            y=cutSignalsNeural.signals{1,cEll}(visstim(k,:),cutting);
            meany=mean(cutSignalsNeural.signals{1,cEll}(visstim(k,:),cutting));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%change here to change the subplot layout%%%%%%%
            %%%% eg subplot(4,1,k) would give you 4 plots side-by-side%%%%%
			subplot(6,7,k);
			plot(x, y, 'Color', [0.5 0.5 0.5]); 
			hold on; 
			plot(x, meany, 'Color', 'r');
			ylim ([-1 3]);
            xlim ([first_recording_used last_rec_used]);
            ylabel('dF/F');
            xlabel('Time');
            set(gca, 'FontSize', 6);
            hold on;
			line ([stimRESPONSEline stimRESPONSEline], ylim, 'Color', 'b', 'LineStyle', '--');
			hold on; 
			line ([stimONline stimONline], ylim, 'Color', [0.7 0.7 0.7]);
            hold on; 
            t=title(tt);
            set(t, 'FontSize', 6);
            clearvars xposypos oriVal tt;
        end
        %%change main title here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T=sgtitle(['Cell', num2str(cEll), 'Experiment ',exp_ID]);
	set(T, 'FontSize', 6);
	%saveas (figure1, [exp_ID, num2str(cEll), 'VisResponses', '.bmp']);
    
    %%%%%%change filetype here. Change to .bmpfor higher quality%%%%%%%%%%
    filename=strcat(exp_ID, '_Cell', num2str(cEll), 'Traces.jpg');
    	saveas (figure1, filename);
        %%clear figure window
	clf;
	end
	close all; 
	%%?
	cd (current_exp);
	disp(['Plot for ', num2str(i), '/', num2str(length(Exps)) ' complete']); 
	clearvars -except i aniID Exps localReposPath; 
end 

disp 'All plots complete'; 
cd(fullfile(localReposPath,aniID));



end

	
	
	
	
