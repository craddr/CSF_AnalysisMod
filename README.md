# CSF_AnalysisMod
Codes required for the analysis of imaging data and visual stimulus data from Rosie Craddock Contrast Sensitivity Function Experiment 2024.

PC Setup for analysis:

Install MATLAB 2018a or later 

Install Suite2P as per Github instructions https://suite2p.readthedocs.io/en/latest/

Download all code in reop, save to C:\\Analysis. 

Save Analysis directory and all subdirectoties to MATLAB path.

Change makeEventFilePsychStimRosie.m on line 15 so localReposPath points to the directory in which your experimental data is found.

Change collateDataBatchRosie.m on line 35 so localReposPath points to the directory in which your experimental data is found.

Change pullStimOnsetPsychStim_s2pRosie on line 31 so localReposPath points to the directory in which your experimental data is found.

Change NcellsAndRespSize4.m on line 8 so localReposPath points to the directory in which your experimental data is found.

Running Analysis: 

create cell masks from 2p GCaMP6f imaging data using suite2p and the ops file contained in the main repo.

in MATLAB run: 

collateDataBatchRosie 

makeEventFilePsychStimRosie

pullStimOnsetPsychStim_s2pRosie

NcellsAndRespSize4

Further data analysis can be completed in R, codes are available on request from Rosie Craddock (craddock.rosie@gmail.com)

Most codes are written by Rosie Craddock (2024). Many codes are based on those of Adam Ranson and Asta Vasalauskaite during their time at Cardiff University.
Author and contributions are indicated in the first few lines of each code.


Suite2p is used in this analysis pipeline. Please cite suite2p as required: 

Suite2p: beyond 10,000 neurons with standard two-photon microscopy
Marius Pachitariu, Carsen Stringer, Mario Dipoppa, Sylvia Schröder, L. Federico Rossi, Henry Dalgleish, Matteo Carandini, Kenneth D. Harris
bioRxiv 061507; doi: https://doi.org/10.1101/061507


