% % MAIN_IMPORT_DATA.m
% This code imports the raw experimental data and reorganized in a way that
% it can be used for further processing.
% The following scripts are used
    % (root folder) set_paths_pHstudy.m
    % (requirements\experimental_data_import) concatenateImpData.m, impData_case.m, selectSetup_pH.m
    % (data\raw_data) original data, found in .xlsx files
% Processed datasets are saved in 'data\processed_data'
    
% Section structure
    % Section 1. Import data 
    % Section 2. Check correct import
    % Section 3. Save processed data

%% Section 1 Import data
% Performed for each progression curve
dbstop if error
set_paths_pHstudy;
setup.importAll = 1; 
    % 0: only selected dataset imported (lines below)
    % 1: import all data sets
setup.plotOutput = 0;
    % 0: do not show any plot
    % 1: show plots of the imported data
setup.saveOutput = 0;
    % 0: do not save plots
    % 1: asve plots in 'data/processed_data' folder ('setup.plotOutput' must be equal to 1).
    
if setup.importAll == 0 % when we want to select a specific dataset to import
    setup.caseALD = 1;
    setup.caseENO = 0;
    setup.caseGAPDH = 0;
    setup.caseGAPDHr = 0;
    setup.caseHXK = 0;
    setup.casePDC = 0;
    setup.casePFK = 0;
    setup.casePGI = 0;
    setup.casePGM = 0;
    setup.casePYK = 0;
    setup.caseTPI = 0;
elseif setup.importAll == 1 % when all the data are to be imported
    setup.caseALD = 1;
    setup.caseENO = 1;
    setup.caseGAPDH = 1;
    setup.caseGAPDHr = 1;
    setup.caseHXK = 1;
    setup.casePDC = 1;
    setup.casePFK = 1;
    setup.casePGI = 1;
    setup.casePGM = 1;
    setup.casePYK = 1;
    setup.caseTPI = 1;    
end

% sequence, one-by-one
if setup.caseALD == 1
    % only one case active
    setup.caseStudyALD = 1;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);   
    data.ALD = tempData;
end
if setup.caseENO == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 1;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.ENO = tempData;
end
if setup.caseGAPDH == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 1;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.GAPDH = tempData;
end
if setup.caseGAPDHr == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 1;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.GAPDHr = tempData;
end
if setup.caseHXK == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 1;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.HXK = tempData;
end
if setup.casePDC == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 1;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.PDC = tempData;
end
if setup.casePFK == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 1;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.PFK = tempData;
end
if setup.casePGI == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 1;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.PGI = tempData;
end
if setup.casePGM == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 1;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.PGM = tempData;
end
if setup.casePYK == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 1;
    setup.caseStudyTPI = 0;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.PYK = tempData;
end
if setup.caseTPI == 1
    % only one case active
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 1;
    
    selectSetup_pH;
    [tempData] = impData_case(setup);    
    data.TPI = tempData;
end

%% Section 2. Check correct import
control_data_import

%% Section 3. Save processed data
[expData] = concatenateImpData(data); %by now all the inputs were created and save in .mat files

    
% disp('stop here'); save('dataPYK.mat','tempData');
% % saving imported data (when all can be simultaneously run)
% % saveName = 'results\expData_pH.mat';
% saveName = 'expData_pH.mat';
% save(saveName,'data');
% %%
% load('expData.mat');
% expData.eno = tempData;
% save('expData.mat','expData');


