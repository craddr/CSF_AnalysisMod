%% RC 2024

function [aniID]= GetAnimalID2()

prompt= {'Provide AnimalID'}; 

dlgtitle='Access Local Repository'; 

dims=[1 20]; 

definput={'CFAA###'};

IDanswer=inputdlg(prompt,dlgtitle,dims,definput);

aniID= IDanswer{1,1}; 


end