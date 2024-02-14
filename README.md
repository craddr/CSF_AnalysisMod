# CSF_AnalysisMod
Codes required for the analysis of imaging data and visual stimulus data from Rosie Craddock Contrast Sensitivity Function Experiment 2024.

PC Setup for analysis:
Install MATLAB 2018a or later 
Install Suite2P as per Github instructions https://suite2p.readthedocs.io/en/latest/
Download all code in reops and add folder and subdirs to MATLAB path. 

Running Analysis: 
create cell masks from 2p GCaMP6f imaging data using suite2p and the ops file contained in the main repo.

in MATLAB run: 
collateDataBatchRosie 
makeEventFilePsychStimRosie
pullStimOnsetPsychStim_s2pRosie
NcellsAndRespSize

Further data analysis can be completed in R, codes are available on request from Rosie Craddock (craddock.rosie@gmail.com)

Most codes are written by Rosie Craddock (2024). Many codes are based on those of Adam Ranson and Asta Vasalauskaite during their time at Cardiff University.
Author and contributions are indicated in the first few lines of each code.


