# El Tabbal et al. 2020 "Removal of the extracellular matrix biases auditory cortical layer dynamics towards supragranular frequency integration"
 
A read me file for the code and data to reproduce the analysis done in the above mentioned paper
 
## Dowloading the data required to reproduce the figures: 

For running the code you will need to download the following data from the following link 
https://cloud.lin-magdeburg.de/s/8Ee2nFETjq3ZaCb
Within the folder, you will find the following data files:
* all_animals_BF_singlesink_hyase_10thFeb2020.mat
* ControlECMDATA.mat
* DataforGroupanalysisHYASE.mat
* EVOKEDPSD_10thFeb2020.mat
* SpontaneousECMDATA.mat
* SpontaneousPSD_10thFEB202.mat

The folder further contains the used code for n-way repeated measures ANOVAs within a zip-folder. After downloading all the files you need to download the MVGC_toolbox for granger causality (https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html) along with the Chronux toolbox for oscillation analysis (http://chronux.org/). The available data are the single-trial traces of the current-source density (CSD) calculated by the laminar local field potential (cf. Happel et al., 2010; J Neurosci). Raw recording files (plx and mat format) can be obtained by the authors by request.  

## Usage of the code and reproduction of figures:
the following functions will help reproduce their corresponding figures
1) Figure 2 and 3 can be reproduced using the function "Tabbal_et_al_CSD_group_Hyaseanalysis"  
2) Figure 4 can be reproduced using the function "Tabbal_et_al_OscillationHyase_at_BF_function_chronux"  
3) Figure 5 can be reproduced using the function "Tabbal_et_al_singletrial_granger"  
4) Figure 6 can be reproduced using the function "Tabbal_et_al_CSD_Group_SpontAnalysis_function_00"  
5) Figure Supplementary can be reproduced using the function "Tabbal_et_al_Control_CSD_group_analysis"  
