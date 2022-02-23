% % MAIL_PARAMETER_ESTIMATION.m
% In this file the parameter estimation assays are run for each progression
% curve in glycolysis
% The stages are
% Stage 1: HXK
% Stage 2: PGI
% Stage 3: PFK
% Stage 4: ALD
% Stage 5: TPI
% Stage 6: GAPDH
% Stage 7: GAPDHR
% Stage 8: PGM
% Stage 9: ENO
% Stage 10: PYK
% Stage 11: PDC
% Stage 12: control 

close all
dbstop if error
set_paths_pHstudy;
setup.plotOutput = 0;
    % 0: do not show any plot
    % 1: show plots of the imported data
setup.saveOutput = 0;
    % 0: do not save plots
    % 1: asve plots in 'data/processed_data' folder ('setup.plotOutput' must be equal to 1).
% select which cases to estimate parameters
setup.estimation.caseALD = 1;
setup.estimation.caseENO = 1;
setup.estimation.caseGAPDH = 1;
setup.estimation.caseGAPDHR = 1;
setup.estimation.caseHXK = 1;
setup.estimation.casePDC = 1;
setup.estimation.casePFK = 1;
setup.estimation.casePGI = 1;
setup.estimation.casePGM = 1;
setup.estimation.casePYK = 1;
setup.estimation.caseTPI = 1;


%% STAGE 1 HXK
if setup.estimation.caseHXK == 1
    pEst_hxk
    close all
end
%% STAGE 2 PGI
if setup.estimation.casePGI == 1
    pEst_pgi
    close all
end
%% STAGE 3 PFK
if setup.estimation.casePFK == 1
    pEst_pfk
    close all
end
%% STAGE 4 ALD
if setup.estimation.caseALD == 1
    pEst_ald
    close all
end
%% STAGE 5 TPI
if setup.estimation.caseTPI == 1
    pEst_tpi
    close all
end
%% STAGE 6 GAPDH
if setup.estimation.caseGAPDH == 1
    pEst_gapdh
    close all
end
%% STAGE 7 GAPDHR
if setup.estimation.caseGAPDHR == 1
    pEst_gapdhr
    close all
end
%% STAGE 8 PGM
if setup.estimation.casePGM == 1
    pEst_pgm
    close all
end
%% STAGE 9 ENO
if setup.estimation.caseENO == 1
    pEst_eno
    close all
end
%% STAGE 10 PYK
if setup.estimation.casePYK == 1
    pEst_pyk
    close all
end
%% STAGE 11 PDC
if setup.estimation.casePDC == 1
    pEst_pdc
    close all
end


%% STAGE 12 control
if setup.estimation.caseHXK == 1
    clearvars -except setup
    setup.enzymeName = 'hxk';
    control_parameter_estimation
end
if setup.estimation.casePGI == 1
    clearvars -except setup
    setup.enzymeName = 'pgi';
    control_parameter_estimation
end
if setup.estimation.casePFK == 1
    clearvars -except setup
    setup.enzymeName = 'pfk';
    control_parameter_estimation
end
if setup.estimation.caseALD == 1
    clearvars -except setup
    setup.enzymeName = 'ald';
    control_parameter_estimation
end
if setup.estimation.caseTPI == 1
    clearvars -except setup
    setup.enzymeName = 'tpi';
    control_parameter_estimation
end
if setup.estimation.caseGAPDH == 1
    clearvars -except setup
    setup.enzymeName = 'gapdh';
    control_parameter_estimation
end
if setup.estimation.caseGAPDHR == 1
    clearvars -except setup
    setup.enzymeName = 'gapdhr';
    control_parameter_estimation
end
if setup.estimation.casePGM == 1
    clearvars -except setup
    setup.enzymeName = 'pgm';
    control_parameter_estimation
end
if setup.estimation.caseENO == 1
    clearvars -except setup
    setup.enzymeName = 'eno';
    control_parameter_estimation
end
if setup.estimation.casePYK == 1
    clearvars -except setup
    setup.enzymeName = 'pyk';
    control_parameter_estimation
end
if setup.estimation.casePDC == 1
    clearvars -except setup
    setup.enzymeName = 'pdc';
    control_parameter_estimation
end

