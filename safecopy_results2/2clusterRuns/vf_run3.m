% GAPDH Parameter estimation after the changes discussed with Bas an Laura:
% - Additional assays to test robutness of regularization results
% - Keq is no more a parameter to be estimated, but a pH-dependent constant
% - Haldane relationship used to calculate Vr and Vf
% - Kinetics from van Heerden 2014 are used.

% Files changed from previous sections
% 1 - 'odeGAPDHr.m' to 'odeGAPDH_vHeerden.m'
% 2 - (different ode to be called)

% Section order in this file
% (0) Setup and data load
% % % (1) System simulation and PSA


%% (0) Setup and data load
clear
set_paths_pHstudy;
dbstop if error

% select specific case and recall data
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
% added
setup.saveOutput = 0;

load('expData.mat','expData');
import_gapdhR = expData.gapdhr;

DFs = setup.DFactorsTotal;
pHtested = setup.pHtested;
numpHtested = nnz(pHtested);
pHs = numpHtested;
blank = zeros(pHs,DFs);
blankCell = cell(pHs,DFs);

% data reorganization
pHTemp = blank';
DFTemp = blank';
abs_meanTemp = blankCell';
abs_stdTemp = blankCell';
conc_meanTemp = blankCell';
conc_stdTemp = blankCell';
timeTemp = blankCell';
RRsTemp = blankCell';

pHarray = unique(import_gapdhR.treatedData.pH_corrected);
for i = 1:numpHtested
    pHval = pHarray(i);
    tempID = find(import_gapdhR.treatedData.pH_corrected==pHval);
    pHTemp(:,i) = import_gapdhR.treatedData.pH_corrected(tempID);
    DFTemp(:,i) = import_gapdhR.treatedData.dilution_corrected(tempID);
    for j = 1:4
        abs_meanTemp{j,i} = import_gapdhR.treatedData.absorbance_mean{tempID(j)};
        abs_stdTemp{j,i} = import_gapdhR.treatedData.absorbance_std{tempID(j)};
        conc_meanTemp{j,i} = import_gapdhR.treatedData.concentration_mean{tempID(j)};
        conc_stdTemp{j,i} = import_gapdhR.treatedData.concentration_std{tempID(j)};
        timeTemp{j,i} = import_gapdhR.treatedData.time{tempID(j)};
        RRsTemp{j,i} = import_gapdhR.treatedData.reaction_rate{tempID(j)};
    end
end

pH = pHTemp';
DF = DFTemp';
abs_mean = abs_meanTemp';
abs_std = abs_stdTemp';

conc_mean = conc_meanTemp';
conc_std = conc_stdTemp';
time = timeTemp';
RRs = RRsTemp';
Vmax = blank';

NADH = blankCell;
Vmax = blank;
for i = 1:(DFs*numpHtested)
    tempDiff = conc_mean{i} - conc_mean{i}(1); % all stoichiometries are 1-to-1.
    NADH{i} = conc_mean{i};
%     RRs2 = RRs';
    Vmax(i) = max(abs(RRs{i}));
end

pH = pHTemp';
DF = DFTemp';
abs_mean = abs_meanTemp';
abs_std = abs_stdTemp';
conc_mean = conc_meanTemp';
conc_std = conc_stdTemp';
time = timeTemp';
RRs = RRsTemp';
clear pHTemp DFTemp abs_meanTemp abs_stdTemp conc_meanTemp conc_stdTemp timeTemp RRsTemp

% save in data
data.pH = pH;
data.DF = DF;
data.abs_mean = abs_mean;
data.abs_std = abs_std;
data.conc_mean = conc_mean;
data.conc_std = conc_std;
data.time = time;
data.RRs = RRs;
data.Vmax = Vmax;
data.chosenVmax = max(max(Vmax));
data.chosenNADini = 0.15;
temp1 = import_gapdhR.rawData.absorbance_corrected{4,4};
temp2 = import_gapdhR.rawData.absorbance_corrected{5,4};
temp3 = import_gapdhR.rawData.absorbance_corrected{6,4};
data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
data.raw.time = import_gapdhR.rawData.time{1};
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Directly changing the concentration here, sicne the extinction
% coefficient did not change.
dps = length(NADH{1,1});
for i = 1:DFs
    for j = 1:numpHtested
        for k = 1:dps
            switch j
                case {1,2,3,4,5,6,7}
                    NADH{j,i}(k) = NADH{j,i}(k) - NADH{j,DFs}(dps);
                case {8,9,10,11,12}
                    NADH{j,i}(k) = NADH{j,i}(k) - 0.0453;
%                     NADH{j,i}(k) = NADH{j,i}(k) - NADH{6,DFs}(dps);
            end
        end
    end
end
data.conc_mean = NADH;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

pHvals = unique(import_gapdhR.treatedData.pH_corrected);
% % visualize: check calculations made
% figure('units','normalized','outerposition',[0 0 1 1])
% for i = 1:numpHtested
%     subplot(3,4,i)
%     for j = 1:DFs
%         plot(time{i,j},NADH{i,j},'.-')
%         hold on
%     end
%     title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
%     if i == numpHtested
%         if setup.caseStudyGAPDHr == 1
%             legend('DF 8','DF 4','DF 2','DF 1')
%         end
%     end
% end
% suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
% suptitle(suptitleName);


%% Full estimationg: variable parameters + vm constant + Km cosntant
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode = 'vanHeerden2014';
setup.ode_pH = 'on';
ode_pH = setup.ode_pH;
setup.numpHtested = numpHtested;

setup.DFstudy = 4;
setup.costfun = 1;

setup.plotResults = 0;
setup.plotEachSim = 0;
setup.plotEachSimCF = 0;

% % All variable
% optfun = @costfun_allVars;
% plength = 60; 
% x_temp = zeros(1,plength);
% ub = 3*ones(1,plength);
% lb = -3*ones(1,plength);
% ub([1,7,13,19,25,31,37,43,49,55]) = 6; %p1
% ub([6,12,18,24,30,36,42,48,54,60]) = 6; % p6
% lb([1,7,13,19,25,31,37,43,49,55]) = -6; %p1
% lb([6,12,18,24,30,36,42,48,54,60]) = -6; % p6
% options = optimset('Display','iter');

% Km variable, Vm fixed
optfun = @costfun_Vmfixed;
plength = 42; % Kms (4) + Vms (2) * numpH (10)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
ub([1,2]) = 6;
lb([1,2]) = -6;
options = optimset('Display','iter');
% weightTest = [1E-8 1E-7 1E-6 1E-5 1E-4 3E-4 1E-3 3E-3 1E-2 3E-2 1E-1 3E-1 1E0 3E0 1E1 3E1 1E2 3E2 1E3 3E3 1E4 3E4 1E5 1E6 1E7 1E8 1E9 1E10];
% weightTest = [1E-8 1E-7 1E-6 1E-5 1E-4 3E-4]; %#1
% weightTest = [1E-3 3E-3 1E-2 3E-2 1E-1 3E-1]; %#2
weightTest = [1E0 1E1 1E2 1E3 1E4 1E5]; %#3
% weightTest = [3E5 1E6 3E6 1E7 3E7 1E8]; %#4
% weightTest = [3E8 1E9 3E9 1E10 3E10 1E11]; %#5
% weightTest = [3E11 1E12 3E12 1E13 3E13 1E14]; %#6

% % Vm variable, Km fixed
% optfun = @costfun_Kmfixed;
% plength = 24; % Kms (4) + Vms (2) * numpH (10)
% x_temp = zeros(1,plength);
% ub = 3*ones(1,plength);
% lb = -3*ones(1,plength);
% ub(5:end) = 6;
% lb(5:end) = -6;
% options = optimset('Display','iter','MaxIter',1);
% % weightTest = [1E-6 3E-6 1E-5 3E-5 1E-4 3E-4 1E-3 3E-3 1E-2 3E-2 1E-1 3E-1 1E0 3E0 1E1 3E1 1E2 3E2 1E3 3E3 1E4 3E4 1E5 3E5 1E6 3E6 1E8]; 
% weightTest = [1E-4 3E-4 1E-3 3E-3 1E-2 3E-2]; %#1
% % weightTest = [1E-1 3E-1 1E0 3E0 1E1 3E1]; %#2
% % weightTest = [1E2 3E2 1E3 3E3 1E4 3E4]; %#3
% % weightTest = [1E5 3E5 1E6 3E6 1E7 3E7]; %#4
% % weightTest = [1E8 3E8 1E9 3E9 1E10 3E10]; %#5
% % weightTest = [1E11 3E11 1E12 3E12 1E13 3E13]; %#6

% testing several regularization factors
% weightTest = [0 2E-4 1E-3 2E-3 1E-2 2E-2 1E-1 2E-1 1E0 2E0];
% weightTest = [0 1E-2 2E-2 1E-1 2E-1 1E0 2E0 1E1 2E1 3E1 3E2];
% weightTest = [0 3E2];
% weightTest = [1E-4 1E-2 1E0 1E2 1E4];
xres_cell = cell(length(weightTest),1);
errorData_cell = cell(length(weightTest),1);
errorHaldane_cell = cell(length(weightTest),1);
errorRegpars_cell = cell(length(weightTest),1);

for o = 1:length(weightTest)
    % early steps
    setup.weightData = 1;
    setup.weightDataEsp = ones(1,setup.numpHtested);
    setup.weightHaldane = weightTest(o);
    setup.selectedLambda = 0;
    
    % parameter estiamtion
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;
    
    % calculating error
    temp_wD = setup.weightData;
    temp_wH = setup.weightHaldane;
    temp_wL = setup.selectedLambda;
    setup.weightData = 1;
    setup.weightHaldane = 1;
    setup.selectedLambda = 0;
    [error] = optfun(xres,data,setup);
    setup.weightData = temp_wD;
    setup.weightHaldane = temp_wH;
    setup.selectedLambda = temp_wL;
    
    % closing, placing thigns in matrix and change reset (just in case for later)
    xres_cell{o} = xres;
    errorData_cell{o} = error(1:260); %the first 260
    errorHaldane_cell{o} = error(261:270); % residual(261:270)
    errorRegpars_cell{o} = error(271:end); % not needed
    setup.selectedLambda = 0;
end

% disp('Here with Km fixed');
% %% Visualization

eData = zeros(length(weightTest),1);
eHaldane = zeros(length(weightTest),1);
for i = 1:length(weightTest)
    eData(i) = sum(abs(errorData_cell{i}));
    eHaldane(i) = sum(abs(errorHaldane_cell{i}));
end

output.xres = xres_cell;
output.weightTest = weightTest;
output.errorData = errorData_cell;
output.errorHaldane = errorHaldane_cell;
output.errorRegpars = errorRegpars_cell;
output.eData = eData;
output.eHaldane = eHaldane;

% save('output_vf_run1.mat','output') %#1
% save('output_vf_run2.mat','output') %#2
save('output_vf_run3.mat','output') %#3
% save('output_vf_run4.mat','output') %#4
% save('output_vf_run5.mat','output') %#5
% save('output_vf_run6.mat','output') %#6

% save('output_kf_run1.mat','output') %#1
% save('output_kf_run2.mat','output') %#2
% save('output_kf_run3.mat','output') %#3
% save('output_kf_run4.mat','output') %#4
% save('output_kf_run5.mat','output') %#5
% save('output_kf_run6.mat','output') %#6
