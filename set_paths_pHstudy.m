folder=cd;

tempfolder = [folder,'/requirements']; % assays
addpath(genpath(tempfolder)); 

tempfolder = [folder,'/parameter_estimation']; % parameterEstimation
addpath(genpath(tempfolder)); 

tempfolder = [folder,'/data']; % experimentalData
addpath(genpath(tempfolder)); 

tempfolder = [folder,'/models'];
addpath(genpath(tempfolder)); 

tempfolder = [folder,'/manuscript']; % results
addpath(genpath(tempfolder)); 

% temporarily created
tempfolder = [folder,'/safecopy_late'];
addpath(genpath(tempfolder)); 
