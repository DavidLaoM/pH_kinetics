% In this script all the datasets provided by Laura are imported. The
% pipeline is run in a case by case basis.

% (0) Select which data is imported.
% (1) Get required setup ('setup_(...).m')
% (2) Import data for the chosen case ('impData_case.m')
% (3) Plot experimental data for the chosen case ('plotExpData.m')
clear
dbstop if error

%%0. Which datasets are being imported? (1 = YES, 0 = NO, 2 = NO (already imported))
set_paths_pHstudy;

setup.importAll = 0;
if setup.importAll == 0 % when we want to select a specific dataset to import
    setup.caseALD = 0;
    setup.caseENO = 0;
    setup.caseGAPDH = 0;
    setup.caseGAPDHr = 0;
    setup.caseHXK = 0;
    setup.casePDC = 0;
    setup.casePFK = 0;
    setup.casePGI = 0;
    setup.casePGM = 0;
    setup.casePYK = 0;
    setup.caseTPI = 1;
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

%%import datasets (step 1, 2 and 3 for each enzyme)
% caseStudy... will avoid troubles in the scripts that are called.
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

if setup.importAll == 1
    [expData] = concatenateImpData; %by now all the inputs were created and save in .mat files
end
    
% disp('stop here'); save('dataPYK.mat','tempData');
% % saving imported data (when all can be simultaneously run)
% % saveName = 'results\expData_pH.mat';
% saveName = 'expData_pH.mat';
% save(saveName,'data');
% %%
% load('expData.mat');
% expData.eno = tempData;
% save('expData.mat','expData');


