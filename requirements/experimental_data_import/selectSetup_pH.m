
if setup.caseStudyALD == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20190816_ald_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20190816_ald_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20190816_ald_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'ald';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = [0 0 1 0 0 1 1 1 1 1 1 1];
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1/2;
    % Protein data
    setup.concProtein = 1.778064829665; % needs updating 2.0106363826145692;
    setup.unitsProtein = 'mg mL^{-1}'; % needs updating 
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1; % needs updating 
    setup.constantVm = 0;
    % Reaction specifics
    setup.PSAmets = {'FBPo';'DHAPo';'GAPo';'G3Po';'NADHo';'NADo'};
    setup.observableMetabolite = 5;
    setup.params = {'K_{fbp} [mM]'; 'K_{gap} [mM]'; 'K_{dhap} [mM]'; 'v_{m} [mM s^{-1}]'; 'K_{eq} []'};
    pH_vals_real= [6.32 6.81 7.06 7.29 7.51 7.68 7.81 7.90];
    setup.Keq_FBA = [1.0E-3 7.9E-4 7.3E-4 7E-4 6.8E-4 6.7E-4 6.6E-4 6.5E-4];  %dir+
    setup.Keq_TPI = [1/(8.31) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38) 1/(9.39)];  %dir-
    setup.Keq_GPD = [1/(4.2E-6) 1/(1.5E-5) 1/(2.7E-5) 1/(4.7E-5) 1/(7.9E-5) 1/(1.2E-4) 1/(1.6E-4) 1/(2.0E-4) ]; %dir-
    % 3 + 8
    setup.parameterNames = {'K_{fbp} [mM]'; 'K_{gap} [mM]'; 'K_{dhap} [mM]';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';...
        };
end
if setup.caseStudyENO == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20191004_eno_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20191004_eno_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20191004_eno_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'eno';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/1.2223; % needs updating IS THIS RIGHT?
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = [1 1 1 1 1 1 1 1 1 1 1 1];
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
%     setup.keq = [5.22 5.22 5.21 5.21 5.2 5.2 5.2 5.19 5.19 5.19];
    setup.keq = [5.22 5.22 5.22 5.21 5.21 5.2 5.2 5.2 5.2 5.19 5.19 5.19];
    setup.branchFactor = 1;
    % Protein data
    setup.concProtein = 2.03139676907012;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60;
    % Reaction specifics
    setup.PSAmets = {'P2G';'PEP'};
    setup.observableMetabolite = 2;
    setup.params = {'v_{m} [mM s^{-1}]'; 'K_{2pg} [mM]'; 'K_{3pg} [mM]'; 'K_{eq} []'};
    % added to fit the other profiles
    setup.constantVm = 0;
    % 2 + 12
    setup.parameterNames = {'K_{2pg} [mM]';'K_{3pg} [mM]';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';'v_{m.pH11} [mM s^{-1}]';'v_{m.pH12} [mM s^{-1}]';...
        };
end
if setup.caseStudyGAPDH == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20190830_gapdh_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20190830_gapdh_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20190830_gapdh_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 5;
    % Data-Enzyme specifics
    setup.enzymeName = 'gapdh';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = ones(1,12);
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1;
    % Protein data
    setup.concProtein = 2.19111549413548;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    % Reaction specifics
    setup.DFactorsTotal = 4;
%     setup.PSAmets = {'P3Go';'ATPo';'BPGo';'ADPo';'NADo';'GAPo';'Pio';'NADH'};
%     setup.PSAobservableMetabolite = 8;
%     setup.params = {'K_{eq} []'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'k_{cat} [mM s^{-1}]';'p_{LinkReaction}'};
    % Reaction specifics
    setup.PSAmets = {'P3G';'ATP';'BPG';'ADP';'NAD';'GAP';'PHOS';'NADH'};
    setup.observableMetabolite = 8;
%     setup.params = {'K_{eq} []'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'v_{max} [mM s^{-1}]'; 'p_{LinkReaction}'};
    setup.params = {'v_{maxFWD} [mM s^{-1}]'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'v_{maxREV} [mM s^{-1}]'};%; 'kPGK_{eq} []'};
    setup.paramRefVals = [0.0056, 2.48, 1.18, 2.92, 0.022, 1, 1];
    setup.costfun2 = 0;
    setup.constantVm = 0;
    setup.selectedLambda = 0;
    % pH for the Keq
% % % %     setup.pH_vals = [6.1900    6.2600    6.4100    6.6000    6.8100    7.0600    7.2900    7.5100    7.6800    7.8100];
    setup.pH_vals = [6.1900    6.2600    6.3200    6.4100    6.6000    6.8100    7.0600    7.2900    7.5100    7.6800    7.8100    7.9000];
    % Keq_gapdh obtained via direct calculation (previous was made worng
    % (reverse way))
%     setup.pH_Keq_gapdh = [0.009361986321002,   0.007968342670936,   0.005641154602417,   0.003642235325689,   0.002245784097423,   0.001262897205436,   0.000743649006787,   0.000448092837714,   0.000302947938974,   0.000224578409742];
    setup.pH_Keq_gapdh = 0.029 * 0.05 .* 10 .^ (-7 + setup.pH_vals);
    % Keq_gapdh obtained with eQuilibrator
% % % %     eQvals = [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1];
% % % %     setup.pH_Keq_gapdh_eQ = 0.05 * (eQvals);
    eQvals = [1.7E-3, 2.3E-3, 2.9E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1, 3.05E-1];
    setup.pH_Keq_gapdh_eQ_fwd = 0.05 * (eQvals);
    % Keq_pgk obtained with eQuilibrator
% % % %     setup.pH_Keq_pgk = [1/(7.4E-4),    1/(7.2E-4),    1/(6.9E-4),    1/(6.5E-4),    1/(6.1E-4),    1/(5.7E-4),    1/(5.5E-4),    1/(5.3E-4),    1/(5.2E-4),   1/(5.2E-4)];
    setup.pH_Keq_pgk_fwd = [1/(7.4E-4),    1/(7.2E-4),    1/(7.1E-4),    1/(6.9E-4),    1/(6.5E-4),    1/(6.1E-4),    1/(5.7E-4),    1/(5.5E-4),    1/(5.3E-4),    1/(5.2E-4),   1/(5.2E-4), 1/(5.1E-4)];
    setup.ode = 'vanHeerden2014';
    setup.exp_vmax_gapdhr = [0.001672560425949,   0.001708335558995,   0.001863340099700,   0.002134052317366,   0.001625999420026,   0.001392974025436,   0.001102411496380,   0.001004022996442,   0.000986761856102,   0.000877728986288];
    setup.exp_vmax_gapdhf = [1.878012068989313e-04, 8.699614731347554e-05, 1.222088736070264e-04, 1.498266981509864e-04, 2.306088349420712e-04, 4.253144979769960e-04, 5.882596627863640e-04, 8.174876065012354e-04, 9.314111327450747e-04, 0.001119902785258];
    setup.sourceVm = 'experimentalSlopes';
%     setup.sourceVm = 'literature'; 
    setup.simAllProfiles = 0;
    % 4 + 12
    setup.parameterNames = {'K_{gap} [mM]';'K_{bpg} [mM]';'K_{nad} [mM]';'K_{nadh} [mM]';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';'v_{m.pH11} [mM s^{-1}]';'v_{m.pH12} [mM s^{-1}]';...
        };
end
if setup.caseStudyGAPDHr == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20190814_gapdhr_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20190814_gapdhr_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20190814_gapdhr_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'gapdhr';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = [1 1 0 1 1 1 1 1 1 1 1 0];
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1;
    % Protein data
    setup.concProtein = 2.4635559199658;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    % Reaction specifics
    setup.PSAmets = {'P3G';'ATP';'BPG';'ADP';'NAD';'GAP';'PHOS';'NADH'};
    setup.observableMetabolite = 8;
%     setup.params = {'K_{eq} []'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'v_{max} [mM s^{-1}]'; 'p_{LinkReaction}'};
    setup.params = {'v_{maxFWD} [mM s^{-1}]'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'v_{maxREV} [mM s^{-1}]'};%; 'kPGK_{eq} []'};
    setup.paramRefVals = [0.0056, 2.48, 1.18, 2.92, 0.022, 1, 1];
    setup.costfun2 = 0;
    setup.constantVm = 0;
    setup.selectedLambda = 0;
    % pH for the Keq
    setup.pH_vals = [6.1900    6.2600    6.4100    6.6000    6.8100    7.0600    7.2900    7.5100    7.6800    7.8100];
    % Keq_gapdh obtained via direct calculation (previous was made worng
    % (reverse way))
%     setup.pH_Keq_gapdh = [0.009361986321002,   0.007968342670936,   0.005641154602417,   0.003642235325689,   0.002245784097423,   0.001262897205436,   0.000743649006787,   0.000448092837714,   0.000302947938974,   0.000224578409742];
    setup.pH_Keq_gapdh = 0.029 * 0.05 .* 10 .^ (-7 + setup.pH_vals);
    % Keq_gapdh obtained with eQuilibrator
    eQvals = [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1];
    setup.pH_Keq_gapdh_eQ = 0.05 * (eQvals);
    % Keq_pgk obtained with eQuilibrator
    setup.pH_Keq_pgk = [1/(7.4E-4),    1/(7.2E-4),    1/(6.9E-4),    1/(6.5E-4),    1/(6.1E-4),    1/(5.7E-4),    1/(5.5E-4),    1/(5.3E-4),    1/(5.2E-4),   1/(5.2E-4)];
    setup.ode = 'vanHeerden2014';
    setup.exp_vmax_gapdhr = [0.001672560425949,   0.001708335558995,   0.001863340099700,   0.002134052317366,   0.001625999420026,   0.001392974025436,   0.001102411496380,   0.001004022996442,   0.000986761856102,   0.000877728986288];
    setup.exp_vmax_gapdhf = [1.878012068989313e-04, 8.699614731347554e-05, 1.222088736070264e-04, 1.498266981509864e-04, 2.306088349420712e-04, 4.253144979769960e-04, 5.882596627863640e-04, 8.174876065012354e-04, 9.314111327450747e-04, 0.001119902785258];
    setup.sourceVm = 'experimentalSlopes';
%     setup.sourceVm = 'literature'; 
    setup.simAllProfiles = 0;
    % 4 + 12
    setup.parameterNames = {'K_{gap} [mM]';'K_{bpg} [mM]';'K_{nad} [mM]';'K_{nadh} [mM]';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';...
        };
end
if setup.caseStudyHXK == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20190820_hxk_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20190820_hxk_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20190820_hxk_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'hxk';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = ones(1,12);
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1;
    % Protein data
    setup.concProtein = 1.90841097484276;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    % Reaction specifics
    setup.Keq_G6PDH = [0.402 0.471 0.54 0.663 1.02 1.65 2.92 4.95 8.2 12.1 16.3 16.3];
    setup.Keq_HXK = [2.74E2 3.03E2 3.31E2 3.79E2 5.16E2 7.45E2 1.2E3 1.9E3 3E3 4.3E3 5.7E3 7.0E3];  %dir+
    setup.PSAmets = {'Glc';'G6P';'6PG';'ATP';'ADP';'NADP';'NADPH'};
    setup.observableMetabolite = 7;
    setup.params = {'K_{adp}'; 'K_{atp}'; 'K_{g6p}'; 'K_{glc}'; 'v_{m}'};
    % added to fit the other profiles
    setup.constantVm = 0;
    % 4 + 12
    setup.parameterNames = {'K_{adp} [mM]';'K_{atp} [mM]';'K_{g6p} [mM]';'K_{glc} [mM]';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';'v_{m.pH11} [mM s^{-1}]';'v_{m.pH12} [mM s^{-1}]';...
        };
end
if setup.caseStudyPDC == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20191002_pdc_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20191002_pdc_kinetics_measurement.xlsx';
    setup.fileNameBackground = '20191002_pdc_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'pdc';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = ones(1,12);
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1;
    % Protein data
    setup.concProtein = 1.87306223194818;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    % kinetics
    setup.Keq_ADH = [1/(2.4E-5) 1/(2.8E-5) 1/(3.3E-5) 1/(4.0E-5) 1/(6.3E-5) 1/(1.0E-4) 1/(1.8E-4) 1/(3.1E-4) 1/(5.1E-4) 1/(7.6E-4) 1/(1.0E-3) 1/(1.3E-3)];
    setup.PSAmets = {'PYR';'CO2';'AcAld';'ETOH';'NADH';'NAD'};
    setup.observableMetabolite = 5;
    setup.params = {'K_{pyr}'; 'n_{hill}'; 'v_{m}'};
    setup.constantVm = 0;
    % 2 + 12
    setup.parameterNames = {'K_{pyk} [mM]';'hill []';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';'v_{m.pH11} [mM s^{-1}]';'v_{m.pH12} [mM s^{-1}]';...
        };
end
if setup.caseStudyPFK == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20191203_pfk_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20200225_pfk_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20200225_pfk_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 700; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'pfk';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = ones(1,12);
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1/2;
    % Protein data
    setup.concProtein = 1.90533229559748;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    % kinetics
    setup.Keq_FBA = [1.2E-3 1.1E-3 1.0E-3 9.7E-4 8.6E-4 7.9E-4 7.3E-4 7E-4 6.8E-4 6.7E-4 6.6E-4 6.5E-4];  %dir+ [EC 4.1.2.13]
    setup.Keq_TPI = [1/(8.07) 1/(8.21) 1/(8.30) 1/(8.46) 1/(8.74) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38) 1/(9.38)];  %dir-
    setup.Keq_GPD = [1/(2.9E-6) 1/(3.5E-6) 1/(4.2E-6) 1/(5.3E-6) 1/(8.7E-6) 1/(1.5E-5) 1/(2.7E-5) 1/(4.7E-5) 1/(7.9E-5) 1/(1.2E-4) 1/(1.6E-4) 1/(2.0E-4) ]; %dir-
    setup.Keq_PFK = [2.56E1 2.96E1 3.34E1 4.03E1 5.95E1 9.15E1 1.54E2 2.5E2 4.05E2 5.91E2 7.95E2 9.77E2];  %dir+, E.C.num.: EC 2.7.1.56

    setup.PSAmets = {'FBP';'DHAP';'GAP';'G3P';'NADH';'NAD';'F6P';'ATP';'ADP'};
    setup.observableMetabolite = 5;
    setup.params = {'g_{R}'; 'K_{f6p}'; 'K_{atp}'; 'L'; 'ci_{atp}'; 'Ki_{atp}'; 'c_{amp}'; 'K_{amp}'; 'c_{f26bp}'; 'K_{f26bp}'; 'c_{f16bp}'; 'K_{f16bp}'; 'c_{atp}'; 'v_{m}'};
    setup.constantVm = 0;
    % 13 + 12
    setup.parameterNames = {'g_{R}'; 'K_{f6p}'; 'K_{atp}'; 'L'; 'ci_{atp}'; 'Ki_{atp}'; 'c_{amp}'; 'K_{amp}'; 'c_{f26bp}'; 'K_{f26bp}'; 'c_{f16bp}'; 'K_{f16bp}'; 'c_{atp}';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';'v_{m.pH11} [mM s^{-1}]';'v_{m.pH12} [mM s^{-1}]';...
        };
end
if setup.caseStudyPGI == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20191213_pgi_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20191213_pgi_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20191213_pgi_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'pgi';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = ones(1,12);
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1/2;
    % Protein data
    setup.concProtein = 1.9720450050505;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    setup.constantVm = 0;
    % Equilibrium constants involved
    setup.Keq_PGI = [3.6E-1 3.6E-1 3.6E-1 3.6E-1 3.6E-1 3.61E-1 3.61E-1 3.62E-1 3.63E-1 3.64E-1 3.65E-1 3.66E-1];  %dir+
    setup.Keq_PFK = [2.56E1 2.96E1 3.34E1 4.03E1 5.95E1 9.15E1 1.54E2 2.5E2 4.05E2 5.91E2 7.95E2 9.77E2];  %dir+
    setup.Keq_FBA = [1.2E-3 1.1E-3 1.0E-3 9.7E-4 8.6E-4 7.9E-4 7.3E-4 7E-4 6.8E-4 6.7E-4 6.6E-4 6.5E-4];  %dir+
    setup.Keq_TPI = [1/(8.07) 1/(8.21) 1/(8.31) 1/(8.46) 1/(8.74) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38) 1/(9.39)];  %dir-
    setup.Keq_GPD = [1/(2.9E-6) 1/(3.5E-6) 1/(4.2E-6) 1/(5.3E-6) 1/(8.7E-6) 1/(1.5E-5) 1/(2.7E-5) 1/(4.7E-5) 1/(7.9E-5) 1/(1.2E-4) 1/(1.6E-4) 1/(2.0E-4) ]; %dir-
    % Reaction specifics
    setup.PSAmets = {'G6Po';'F6Po';'ATPo';'FBPo';'ADPo';'DHAPo';'GAPo';'NADHo';'G3Po';'NADo'};
    setup.observableMetabolite = 8;
    setup.params = {'K_{eq} []'; 'K_{g6p} [mM]'; 'K_{f6p} [mM]'; 'v_{m} [mM s^{-1}]'};
    % 2 + 12
    setup.parameterNames = {'K_{g6p} [mM]';'K_{f6p} [mM]';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';'v_{m.pH11} [mM s^{-1}]';'v_{m.pH12} [mM s^{-1}]';...
        };
end
if setup.caseStudyPGM == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20190817_pgm_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20190817_pgm_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20190817_pgm_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'pgm';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = ones(1,12);
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1;
    % Protein data
    setup.concProtein = 2.56229648411088;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    % kinetics
    setup.Keq_PGM = [1/6.62 1/6.44 1/6.29 1/6.12 1/5.83 1/5.62 1/5.46 1/5.38 1/5.33 1/5.3 1/5.29 1/5.29]; %dir-
    setup.Keq_ENO = [5.22 5.22 5.22 5.21 5.21 5.2 5.2 5.2 5.19 5.19 5.19 5.19];  %dir ?
    setup.Keq_PYK = [1/(3.1E-6) 1/(3.5E-6) 1/(3.9E-6) 1/(4.5E-6) 1/(6.5E-6) 1/(9.7E-6) 1/(1.6E-5) 1/(2.5E-5) 1/(4.1E-5) 1/(5.9E-5) 1/(7.9E-5) 1/(9.6E-5)];  %dir-
    setup.Keq_LDH = [1/(2.3E-6)	1/(2.7E-6) 1/(3.2E-6) 1/(3.9E-6) 1/(6.1E-6)	1/(9.9E-6) 1/(1.8E-5) 1/(3.0E-5) 1/(5.0E-5) 1/(7.3E-5) 1/(9.9E-5) 1/(1.2E-4)];
    setup.PSAmets = {'P3G';'P2G';'PEP';'PYR';'LAC';'ADP';'ATP';'NADH';'NAD'};
    setup.observableMetabolite = 8;
    setup.params = {'K_{2pg}'; 'K_{3pg}'; 'v_{m}'};
    setup.constantVm = 0;
    % 2 + 12
    setup.parameterNames = {'K_{2pg} [mM]';'K_{3pg} [mM]';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';'v_{m.pH11} [mM s^{-1}]';'v_{m.pH12} [mM s^{-1}]';...
        };
end
if setup.caseStudyPYK == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20191203_pyk_layout_kinetics_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20190203_pyk_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20190203_pyk_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'pyk';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = ones(1,12);
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6	6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1;
    % Protein data
    setup.concProtein = 2.22158843695271;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    setup.constantVm = 0;
%     % Reaction specifics
    setup.PSAmets = {'ADP';'NADH';'FBP';'PEP';'ATP';'NAD';'PYR';'LAC'};
    setup.observableMetabolite = 2;
    setup.params = {'K_{pep} [mM]';'K_{adp} [mM]';'K_{atp} [mM]';'K_{fbp} [mM]';'L []';'v_{m} [mM s^{-1}]';};
    setup.Keq_PYK = [1/(3.1E-6) 1/(3.5E-6) 1/(3.9E-6) 1/(4.5E-6) 1/(6.5E-6) 1/(9.7E-6) 1/(1.6E-5) 1/(2.5E-5) 1/(4.1E-5) 1/(5.9E-5) 1/(7.9E-5) 1/(9.6E-5)];  %dir-
    setup.Keq_LDH = [1/(2.3E-6)	1/(2.7E-6) 1/(3.2E-6) 1/(3.9E-6) 1/(6.1E-6)	1/(9.9E-6) 1/(1.8E-5) 1/(3.0E-5) 1/(5.0E-5) 1/(7.3E-5) 1/(9.9E-5) 1/(1.2E-4)];
    % 6 + 12
    setup.parameterNames = {'K_{adp} [mM]';'K_{atp} [mM]';'K_{f16bp} [mM]';'K_{pep} [mM]';'L []';'hill []';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';'v_{m.pH8} [mM s^{-1}]';'v_{m.pH9} [mM s^{-1}]';...
        'v_{m.pH10} [mM s^{-1}]';'v_{m.pH11} [mM s^{-1}]';'v_{m.pH12} [mM s^{-1}]';...
        };
end
if setup.caseStudyTPI == 1
%     setup.saveOutput = 1;
%     setup.plotOutput = 1;
    setup.deleteWrongDatapoints = 1;
    % Layout data file
    setup.layoutFileName = '20191212_tpi_layout_kinetics_baseline_experiment.xlsx';
    setup.layoutSheetNumberMax = 5; % old setup.layoutSheetNumber
    setup.sizePlates = 'B2:M9';
    % Data and background file
    setup.FileName = '20191212_tpi_kinetics_experiment.xlsx';
    setup.fileNameBackground = '20191212_tpi_kinetics_baseline.xlsx';
    setup.dataSheetNumberMax = 12; % old setup.dataSheetNumber
    setup.dataPointsMax = 250; % old setup.dataPoints = 50;
    setup.lettermax = 5;
    setup.numberReplicates = 3;
    setup.DFactorsTotal = 4;
    % Data-Enzyme specifics
    setup.enzymeName = 'tpi';
    setup.casesInColums = 2; % review when the line
    setup.extinction_coefficient = 1/4.8278;
    setup.EE_untis = 'UA mM^{-1}';    
    setup.pHtested = [1 0 1 0 1 1 0 1 0 1 0 1];
    setup.fullpHarray = [6.19	6.26	6.32	6.41	6.6     6.81	7.06	7.29	7.51	7.68	7.81	7.9];
    setup.branchFactor = 1/2; %1;
    % Protein data
    setup.concProtein = 1.82951772639691;
    setup.unitsProtein = 'mg mL^{-1}';
    setup.extraDF = 60; setup.DFstudy = 4; setup.costfun = 1;
    % kinetics
    setup.Keq_TPI = [1/(8.07) 1/(8.21) 1/(8.30) 1/(8.46) 1/(8.74) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38) 1/(9.38)];  %dir-
    setup.Keq_GPD = [1/(2.9E-6) 1/(3.5E-6) 1/(4.2E-6) 1/(5.3E-6) 1/(8.7E-6) 1/(1.5E-5) 1/(2.7E-5) 1/(4.7E-5) 1/(7.9E-5) 1/(1.2E-4) 1/(1.6E-4) 1/(2.0E-4) ]; %dir-
    setup.PSAmets = {'FBP';'DHAP';'GAP';'G3P';'NADH';'NAD'};
    setup.observableMetabolite = 5;
    setup.params = {'K_{dhap}'; 'K_{gap}'; 'v_{m}'};
    setup.constantVm = 0;
    % 2 + 7
    setup.parameterNames = {'K_{dhap} [mM]';'K_{gap} [mM]';...
        'v_{m.pH1} [mM s^{-1}]';'v_{m.pH2} [mM s^{-1}]';'v_{m.pH3} [mM s^{-1}]';...
        'v_{m.pH4} [mM s^{-1}]';'v_{m.pH5} [mM s^{-1}]';'v_{m.pH6} [mM s^{-1}]';...
        'v_{m.pH7} [mM s^{-1}]';...
        };
end


% % % % Background from previous application
% %     %     %
% %     setup.layoutFileName = '20190814_layout_gapdh_r.xlsx';
% %     setup.layoutSheetNumber = 2;
% %     setup.sizePlates = 'B2:M9';
% %     setup.FileName = '20190814_GAPDH_R_py.xlsx';
% %     setup.dataSheetNumber = 5;
% %     setup.dataPoints = 50;
% %     setup.lettermax = 3;
% %     setup.fileNameBackground = '20190814_GAPDH_R_background.xlsx';
% %     setup.numberReplicates = 3;
% %     setup.casesInColums = 2;
% %     setup.extinction_coefficient = 1/4.3992;
% %     setup.EE_untis = 'UA mM^{-1}';
% %     setup.numpHtested = 10;
% %     setup.enzymeName = 'GAPDH_{reverse}';
% %     setup.saveOutput = 1;
% %     setup.deleteWrongExperiments = 1;
% %     setup.concProtein = 2.0106363826145692;
% %     setup.unitsProtein = 'mg mL^{-1}';
% %     setup.extraDF = 60;
% %     setup.DFactorsTotal = 4;
% %     setup.PSAmets = {'P3Go';'ATPo';'BPGo';'ADPo';'NADo';'GAPo';'Pio';'NADH'};
% %     setup.PSAobservableMetabolite = 8;
% %     setup.params = {'K_{eq} []'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'k_{cat} [mM s^{-1}]';'LinkReaction'};
% %     %     %