function [error] = costfun_gapdhFWDREV(x_temp,data,setup)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%     x_temp(1) = Kgap
%     x_temp(2) = Kbpg
%     x_temp(3) = Knad
%     x_temp(4) = Knadh
%     x_temp([5:6]) = {Vmf, Vmr}, pH#1
%     x_temp([7:8]) = {Vmf, Vmr}, pH#2
%     x_temp([9:10]) = {Vmf, Vmr}, pH#3
%     x_temp([11:12]) = {Vmf, Vmr}, pH#4
%     x_temp([13:14]) = {Vmf, Vmr}, pH#5
%     x_temp([15:16]) = {Vmf, Vmr}, pH#6
%     x_temp([17:18]) = {Vmf, Vmr}, pH#7
%     x_temp([19:20]) = {Vmf, Vmr}, pH#8
%     x_temp([21:22]) = {Vmf, Vmr}, pH#9
%     x_temp([23:24]) = {Vmf, Vmr}, pH#10

% Early setup cases converted to FWD, could also be REV instead.
% enzyme = setup.enzymeName;
DFs = setup.FWD.DFactorsTotal;
DFstudy = setup.DFstudy; % default is the lowest dilution factor (usually DF1, location 4)
obsMet = setup.FWD.observableMetabolite;
costfun = setup.costfun; % default value is 1. No regularization
lambda = setup.selectedLambda;
numpH = setup.numpHtested;
wD = setup.weightData;
wH = setup.weightHaldane;
wL = setup.selectedLambda;
plotEachSimCF = setup.plotEachSimCF;
simAllProfiles = setup.simAllProfiles;

% reaction-specific setup
wDesp_REV = setup.REV.weightDataEsp;
sourceVm_REV = setup.REV.sourceVm;
ode_pH_REV = setup.REV.ode_pH;
wDesp_FWD = setup.FWD.weightDataEsp;
sourceVm_FWD = setup.FWD.sourceVm;
ode_pH_FWD = setup.FWD.ode_pH;

caseKm = setup.caseKm;

% REVERSE DIRECTON: simulaton
setup.enzymeName = 'gapdhr';
numpH_REV = length(setup.REV.pH_vals);
simNADH_REV = cell(DFs,numpH_REV);
expNADH_REV = cell(DFs,numpH_REV);
simGAPDHr_REV = cell(DFs,numpH_REV);
expGAPDHr_REV = cell(DFs,numpH_REV);
temp_simResult_REV = cell(DFs,numpH_REV);
    % emptying 2 holes for pH values between rev and fwd
    tempLen = length(setup.REV.pH_vals);
    idxsTranspose = zeros(1,tempLen); 
    for k = 1:tempLen
        tempNum = setup.REV.pH_vals(k);
        idxsTranspose(k) = find(setup.FWD.pH_vals == tempNum);
    end
    
for j = 1:numpH_REV
    % select required data
    data.REV.Vmaxs = data.REV.Vmax(j,:);
    data.REV.NADH = data.REV.conc_mean(j,:);
    data.REV.Vprofs = data.REV.RRs(j,:);
    data.REV.tempTime = data.REV.time(j,:);
    % inputs to be selected
    data.REV.KeqGAPDH = setup.REV.pH_Keq_gapdh_eQ; %setup.pH_Keq_gapdh(i);
    data.REV.KeqPGK = setup.REV.pH_Keq_pgk;
    data.REV.chosenKeqGAPDH = data.REV.KeqGAPDH(j);
    data.REV.chosenKeqPGK = data.REV.KeqPGK(j);
    % parameter setup
    switch caseKm
        case 'pH_independent'
            % Mode 1. still keeping kms, but fixed
            idx5 = 4 + j;
            xassay = zeros(1,5);
            xassay(1) = x_temp(1);
            xassay(2) = x_temp(2);
            xassay(3) = x_temp(3);
            xassay(4) = x_temp(4);
            xassay(5) = x_temp(idx5);
        case 'pH_dependent'
            idx1 = j + numpH * 0;
            idx2 = j + numpH * 1;
            idx3 = j + numpH * 2;
            idx4 = j + numpH * 3;
            idx5 = j + numpH * 4;
            xassay = zeros(1,5);
            xassay(1) = x_temp(idx1);
            xassay(2) = x_temp(idx2);
            xassay(3) = x_temp(idx3);
            xassay(4) = x_temp(idx4);
            xassay(5) = x_temp(idx5);
        otherwise
            disp('Warning: no specification has been made on Km being pH dependent or independent');
    end
    % simulations
    for i = DFstudy
        % recall vmax for the specific value and simulate
        data.REV.chosenVmax = data.REV.Vmaxs(1,4)/data.REV.DF(1,i); % vmax from the highest DF is taken and then divided
        data.REV.chosenLink = data.REV.DF(1,i);
        data.REV.chosenNADini = data.REV.NADH{i}(1);
        data.REV.chosenDF = data.REV.DF(j,i);
        setup.REV.excessPGK = 1;

        data.REV.NADH = data.REV.conc_mean(j,:);
        data.REV.Vprofs = data.REV.RRs(j,:);
        data.REV.tempTime = data.REV.time(j,:);                
        data.REV.i = j;

        % simulate metabolites
            setup.REV.plotEachSim = setup.plotEachSim;
        [simResult] = simSys(xassay,data.REV,setup.REV);
        % calculation of reaction rates
        [vObs,~] = calcRates(xassay,simResult,data.REV,setup.REV);   
        % cost function (NADHexp + vGAPDHr)
        simTime = simResult.t;
        simMet = simResult.y(:,obsMet);
        simRate = vObs;

        simNADH_REV{i,j} = interp1(simTime,simMet,data.REV.tempTime{i},'pchip');
        simGAPDHr_REV{i,j} = interp1(simTime,simRate,data.REV.tempTime{i},'pchip');
        expNADH_REV{i,j} = data.REV.NADH{i};
        expGAPDHr_REV{i,j} = -data.REV.Vprofs{i};
        simResult.vObs = vObs;
        temp_simResult_REV{i,j} = simResult;
    end
end

% FORWARD DIRECTON: simulaton
setup.enzymeName = 'gapdh';
numpH_FWD = length(setup.FWD.pH_vals);
simNADH_FWD = cell(DFs,numpH);
expNADH_FWD = cell(DFs,numpH);
simGAPDHr_FWD = cell(DFs,numpH);
expGAPDHr_FWD = cell(DFs,numpH);
temp_simResult_FWD = cell(DFs,numpH);

for j = 1:numpH_FWD
    % select required data
    data.FWD.Vmaxs = data.FWD.Vmax(j,:);
    data.FWD.NADH = data.FWD.conc_mean(j,:);
    data.FWD.Vprofs = data.FWD.RRs(j,:);
    data.FWD.tempTime = data.FWD.time(j,:);
    % inputs to be selected
    data.FWD.KeqGAPDH = setup.FWD.pH_Keq_gapdh_eQ_fwd; %setup.pH_Keq_gapdh(i);
    data.FWD.KeqPGK = setup.FWD.pH_Keq_pgk_fwd;
    data.FWD.chosenKeqGAPDH = data.FWD.KeqGAPDH(j);
    data.FWD.chosenKeqPGK = data.FWD.KeqPGK(j);

    % selecting the right parameters
    switch caseKm
        case 'pH_independent'
            % Mode 1. still keeping kms, but fixed
            idx5 = 4 + j;
            xassay = zeros(1,5);
            xassay(1) = x_temp(1);
            xassay(2) = x_temp(2);
            xassay(3) = x_temp(3);
            xassay(4) = x_temp(4);
            xassay(5) = x_temp(idx5);
        case 'pH_dependent'
            idx1 = j + numpH * 0;
            idx2 = j + numpH * 1;
            idx3 = j + numpH * 2;
            idx4 = j + numpH * 3;
            idx5 = j + numpH * 4;
            xassay = zeros(1,5);
            xassay(1) = x_temp(idx1);
            xassay(2) = x_temp(idx2);
            xassay(3) = x_temp(idx3);
            xassay(4) = x_temp(idx4);
            xassay(5) = x_temp(idx5);
        otherwise
            disp('Warning: no specification has been made on Km being pH dependent or independent');
    end

    % simulations
    for i = DFstudy
        % recall vmax for the specific value and simulate
%                 data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,i); % vmax from the highest DF is taken and then divided
        data.FWD.chosenVmax = data.FWD.Vmax(j,4)/data.FWD.DF(j,i); % vmax from the highest DF is taken and then divided
        data.FWD.chosenLink = data.FWD.DF(1,i);
        data.FWD.chosenNADini = data.FWD.NADH{i}(1);
        data.FWD.chosenDF = data.FWD.DF(j,i);
        setup.FWD.excessPGK = 1;

        data.FWD.NADH = data.FWD.conc_mean(j,:);
        data.FWD.Vprofs = data.FWD.RRs(j,:);
        data.FWD.tempTime = data.FWD.time(j,:);                
        data.FWD.i = j;

        % simulate metabolites
            setup.FWD.plotEachSim = setup.plotEachSim;
        [simResult] = simSys(xassay,data.FWD,setup.FWD);
        % calculation of reaction rates
        [vObs,~] = calcRates(xassay,simResult,data.FWD,setup.FWD);   
        % cost function (NADHexp + vGAPDHr)
        simTime = simResult.t;
        simMet = simResult.y(:,obsMet);
        simRate = vObs;

        simNADH_FWD{i,j} = interp1(simTime,simMet,data.FWD.tempTime{i},'pchip');
        simGAPDHr_FWD{i,j} = interp1(simTime,simRate,data.FWD.tempTime{i},'pchip');
        expNADH_FWD{i,j} = data.FWD.NADH{i};
        expGAPDHr_FWD{i,j} = data.FWD.Vprofs{i};
        simResult.vObs = vObs;
        temp_simResult_FWD{i,j} = simResult;
    end
end



% OVERAL DIRECTON: plotting
for j = 1:numpH
    if plotEachSimCF == 1
        simulationVisualization_gapdhFWDREV;
    end
end

% OVERALL costfun 
switch costfun
    case 3
        % Data error (FORWARD DIREACTION)        
        % DF1
        errorNADH1_df1_FWD = simNADH_FWD{4,1} - expNADH_FWD{4,1};
        errorNADH2_df1_FWD = simNADH_FWD{4,2} - expNADH_FWD{4,2};
        errorNADH3_df1_FWD = simNADH_FWD{4,3} - expNADH_FWD{4,3};
        errorNADH4_df1_FWD = simNADH_FWD{4,4} - expNADH_FWD{4,4};
        errorNADH5_df1_FWD = simNADH_FWD{4,5} - expNADH_FWD{4,5};
        errorNADH6_df1_FWD = simNADH_FWD{4,6} - expNADH_FWD{4,6};
        errorNADH7_df1_FWD = simNADH_FWD{4,7} - expNADH_FWD{4,7};
        errorNADH8_df1_FWD = simNADH_FWD{4,8} - expNADH_FWD{4,8};
        errorNADH9_df1_FWD = simNADH_FWD{4,9} - expNADH_FWD{4,9};
        errorNADH10_df1_FWD = simNADH_FWD{4,10} - expNADH_FWD{4,10};
        errorNADH11_df1_FWD = simNADH_FWD{4,11} - expNADH_FWD{4,11};
        errorNADH12_df1_FWD = simNADH_FWD{4,12} - expNADH_FWD{4,12};
        % DF2
        errorNADH1_df2_FWD = simNADH_FWD{3,1} - expNADH_FWD{3,1};
        errorNADH2_df2_FWD = simNADH_FWD{3,2} - expNADH_FWD{3,2};
        errorNADH3_df2_FWD = simNADH_FWD{3,3} - expNADH_FWD{3,3};
        errorNADH4_df2_FWD = simNADH_FWD{3,4} - expNADH_FWD{3,4};
        errorNADH5_df2_FWD = simNADH_FWD{3,5} - expNADH_FWD{3,5};
        errorNADH6_df2_FWD = simNADH_FWD{3,6} - expNADH_FWD{3,6};
        errorNADH7_df2_FWD = simNADH_FWD{3,7} - expNADH_FWD{3,7};
        errorNADH8_df2_FWD = simNADH_FWD{3,8} - expNADH_FWD{3,8};
        errorNADH9_df2_FWD = simNADH_FWD{3,9} - expNADH_FWD{3,9};
        errorNADH10_df2_FWD = simNADH_FWD{3,10} - expNADH_FWD{3,10};
        errorNADH11_df2_FWD = simNADH_FWD{3,11} - expNADH_FWD{3,11};
        errorNADH12_df2_FWD = simNADH_FWD{3,12} - expNADH_FWD{3,12};
        % DF4
        errorNADH1_df4_FWD = simNADH_FWD{2,1} - expNADH_FWD{2,1};
        errorNADH2_df4_FWD = simNADH_FWD{2,2} - expNADH_FWD{2,2};
        errorNADH3_df4_FWD = simNADH_FWD{2,3} - expNADH_FWD{2,3};
        errorNADH4_df4_FWD = simNADH_FWD{2,4} - expNADH_FWD{2,4};
        errorNADH5_df4_FWD = simNADH_FWD{2,5} - expNADH_FWD{2,5};
        errorNADH6_df4_FWD = simNADH_FWD{2,6} - expNADH_FWD{2,6};
        errorNADH7_df4_FWD = simNADH_FWD{2,7} - expNADH_FWD{2,7};
        errorNADH8_df4_FWD = simNADH_FWD{2,8} - expNADH_FWD{2,8};
        errorNADH9_df4_FWD = simNADH_FWD{2,9} - expNADH_FWD{2,9};
        errorNADH10_df4_FWD = simNADH_FWD{2,10} - expNADH_FWD{2,10};
        errorNADH11_df4_FWD = simNADH_FWD{2,11} - expNADH_FWD{2,11};
        errorNADH12_df4_FWD = simNADH_FWD{2,12} - expNADH_FWD{2,12};
        % DF8
        errorNADH1_df8_FWD = simNADH_FWD{1,1} - expNADH_FWD{1,1};
        errorNADH2_df8_FWD = simNADH_FWD{1,2} - expNADH_FWD{1,2};
        errorNADH3_df8_FWD = simNADH_FWD{1,3} - expNADH_FWD{1,3};
        errorNADH4_df8_FWD = simNADH_FWD{1,4} - expNADH_FWD{1,4};
        errorNADH5_df8_FWD = simNADH_FWD{1,5} - expNADH_FWD{1,5};
        errorNADH6_df8_FWD = simNADH_FWD{1,6} - expNADH_FWD{1,6};
        errorNADH7_df8_FWD = simNADH_FWD{1,7} - expNADH_FWD{1,7};
        errorNADH8_df8_FWD = simNADH_FWD{1,8} - expNADH_FWD{1,8};
        errorNADH9_df8_FWD = simNADH_FWD{1,9} - expNADH_FWD{1,9};
        errorNADH10_df8_FWD = simNADH_FWD{1,10} - expNADH_FWD{1,10};
        errorNADH11_df8_FWD = simNADH_FWD{1,11} - expNADH_FWD{1,11};
        errorNADH12_df8_FWD = simNADH_FWD{1,12} - expNADH_FWD{1,12};
        % all together
        wDesp_t_FWD = wDesp_FWD';
        errorNADH_FWD = [...
            wDesp_t_FWD(1,1)*errorNADH1_df8_FWD; wDesp_t_FWD(2,1)*errorNADH1_df4_FWD; wDesp_t_FWD(3,1)*errorNADH1_df2_FWD; wDesp_t_FWD(4,1)*errorNADH1_df1_FWD; 
            wDesp_t_FWD(1,2)*errorNADH2_df8_FWD; wDesp_t_FWD(2,2)*errorNADH2_df4_FWD; wDesp_t_FWD(3,2)*errorNADH2_df2_FWD; wDesp_t_FWD(4,2)*errorNADH2_df1_FWD;
            wDesp_t_FWD(1,3)*errorNADH3_df8_FWD; wDesp_t_FWD(2,3)*errorNADH3_df4_FWD; wDesp_t_FWD(3,3)*errorNADH3_df2_FWD; wDesp_t_FWD(4,3)*errorNADH3_df1_FWD;
            wDesp_t_FWD(1,4)*errorNADH4_df8_FWD; wDesp_t_FWD(2,4)*errorNADH4_df4_FWD; wDesp_t_FWD(3,4)*errorNADH4_df2_FWD; wDesp_t_FWD(4,4)*errorNADH4_df1_FWD;
            wDesp_t_FWD(1,5)*errorNADH5_df8_FWD; wDesp_t_FWD(2,5)*errorNADH5_df4_FWD; wDesp_t_FWD(3,5)*errorNADH5_df2_FWD; wDesp_t_FWD(4,5)*errorNADH5_df1_FWD;
            wDesp_t_FWD(1,6)*errorNADH6_df8_FWD; wDesp_t_FWD(2,6)*errorNADH6_df4_FWD; wDesp_t_FWD(3,6)*errorNADH6_df2_FWD; wDesp_t_FWD(4,6)*errorNADH6_df1_FWD;
            wDesp_t_FWD(1,7)*errorNADH7_df8_FWD; wDesp_t_FWD(2,7)*errorNADH7_df4_FWD; wDesp_t_FWD(3,7)*errorNADH7_df2_FWD; wDesp_t_FWD(4,7)*errorNADH7_df1_FWD;
            wDesp_t_FWD(1,8)*errorNADH8_df8_FWD; wDesp_t_FWD(2,8)*errorNADH8_df4_FWD; wDesp_t_FWD(3,8)*errorNADH8_df2_FWD; wDesp_t_FWD(4,8)*errorNADH8_df1_FWD;
            wDesp_t_FWD(1,9)*errorNADH9_df8_FWD; wDesp_t_FWD(2,9)*errorNADH9_df4_FWD; wDesp_t_FWD(3,9)*errorNADH9_df2_FWD; wDesp_t_FWD(4,9)*errorNADH9_df1_FWD;
            wDesp_t_FWD(1,10)*errorNADH10_df8_FWD; wDesp_t_FWD(2,10)*errorNADH10_df4_FWD; wDesp_t_FWD(3,10)*errorNADH10_df2_FWD; wDesp_t_FWD(4,10)*errorNADH10_df1_FWD;
            wDesp_t_FWD(1,11)*errorNADH11_df8_FWD; wDesp_t_FWD(2,11)*errorNADH11_df4_FWD; wDesp_t_FWD(3,11)*errorNADH11_df2_FWD; wDesp_t_FWD(4,11)*errorNADH11_df1_FWD;
            wDesp_t_FWD(1,12)*errorNADH12_df8_FWD; wDesp_t_FWD(2,12)*errorNADH12_df4_FWD; wDesp_t_FWD(3,12)*errorNADH12_df2_FWD; wDesp_t_FWD(4,12)*errorNADH12_df1_FWD];
        % Data error (REVERSIBLE DIREACTION)
        % DF1
        errorNADH1_df1_REV = simNADH_REV{4,1} - expNADH_REV{4,1};
        errorNADH2_df1_REV = simNADH_REV{4,2} - expNADH_REV{4,2};
        errorNADH3_df1_REV = simNADH_REV{4,3} - expNADH_REV{4,3};
        errorNADH4_df1_REV = simNADH_REV{4,4} - expNADH_REV{4,4};
        errorNADH5_df1_REV = simNADH_REV{4,5} - expNADH_REV{4,5};
        errorNADH6_df1_REV = simNADH_REV{4,6} - expNADH_REV{4,6};
        errorNADH7_df1_REV = simNADH_REV{4,7} - expNADH_REV{4,7};
        errorNADH8_df1_REV = simNADH_REV{4,8} - expNADH_REV{4,8};
        errorNADH9_df1_REV = simNADH_REV{4,9} - expNADH_REV{4,9};
        errorNADH10_df1_REV = simNADH_REV{4,10} - expNADH_REV{4,10};
        % DF2
        errorNADH1_df2_REV = simNADH_REV{3,1} - expNADH_REV{3,1};
        errorNADH2_df2_REV = simNADH_REV{3,2} - expNADH_REV{3,2};
        errorNADH3_df2_REV = simNADH_REV{3,3} - expNADH_REV{3,3};
        errorNADH4_df2_REV = simNADH_REV{3,4} - expNADH_REV{3,4};
        errorNADH5_df2_REV = simNADH_REV{3,5} - expNADH_REV{3,5};
        errorNADH6_df2_REV = simNADH_REV{3,6} - expNADH_REV{3,6};
        errorNADH7_df2_REV = simNADH_REV{3,7} - expNADH_REV{3,7};
        errorNADH8_df2_REV = simNADH_REV{3,8} - expNADH_REV{3,8};
        errorNADH9_df2_REV = simNADH_REV{3,9} - expNADH_REV{3,9};
        errorNADH10_df2_REV = simNADH_REV{3,10} - expNADH_REV{3,10};
        % DF4
        errorNADH1_df4_REV = simNADH_REV{2,1} - expNADH_REV{2,1};
        errorNADH2_df4_REV = simNADH_REV{2,2} - expNADH_REV{2,2};
        errorNADH3_df4_REV = simNADH_REV{2,3} - expNADH_REV{2,3};
        errorNADH4_df4_REV = simNADH_REV{2,4} - expNADH_REV{2,4};
        errorNADH5_df4_REV = simNADH_REV{2,5} - expNADH_REV{2,5};
        errorNADH6_df4_REV = simNADH_REV{2,6} - expNADH_REV{2,6};
        errorNADH7_df4_REV = simNADH_REV{2,7} - expNADH_REV{2,7};
        errorNADH8_df4_REV = simNADH_REV{2,8} - expNADH_REV{2,8};
        errorNADH9_df4_REV = simNADH_REV{2,9} - expNADH_REV{2,9};
        errorNADH10_df4_REV = simNADH_REV{2,10} - expNADH_REV{2,10};
        % DF8
        errorNADH1_df8_REV = simNADH_REV{1,1} - expNADH_REV{1,1};
        errorNADH2_df8_REV = simNADH_REV{1,2} - expNADH_REV{1,2};
        errorNADH3_df8_REV = simNADH_REV{1,3} - expNADH_REV{1,3};
        errorNADH4_df8_REV = simNADH_REV{1,4} - expNADH_REV{1,4};
        errorNADH5_df8_REV = simNADH_REV{1,5} - expNADH_REV{1,5};
        errorNADH6_df8_REV = simNADH_REV{1,6} - expNADH_REV{1,6};
        errorNADH7_df8_REV = simNADH_REV{1,7} - expNADH_REV{1,7};
        errorNADH8_df8_REV = simNADH_REV{1,8} - expNADH_REV{1,8};
        errorNADH9_df8_REV = simNADH_REV{1,9} - expNADH_REV{1,9};
        errorNADH10_df8_REV = simNADH_REV{1,10} - expNADH_REV{1,10};
        % all together
        wDesp_t_REV = wDesp_REV';
        errorNADH_REV = [...
            wDesp_t_REV(1,1)*errorNADH1_df8_REV; wDesp_t_REV(2,1)*errorNADH1_df4_REV; wDesp_t_REV(3,1)*errorNADH1_df2_REV; wDesp_t_REV(4,1)*errorNADH1_df1_REV; 
            wDesp_t_REV(1,2)*errorNADH2_df8_REV; wDesp_t_REV(2,2)*errorNADH2_df4_REV; wDesp_t_REV(3,2)*errorNADH2_df2_REV; wDesp_t_REV(4,2)*errorNADH2_df1_REV;
            wDesp_t_REV(1,3)*errorNADH3_df8_REV; wDesp_t_REV(2,3)*errorNADH3_df4_REV; wDesp_t_REV(3,3)*errorNADH3_df2_REV; wDesp_t_REV(4,3)*errorNADH3_df1_REV;
            wDesp_t_REV(1,4)*errorNADH4_df8_REV; wDesp_t_REV(2,4)*errorNADH4_df4_REV; wDesp_t_REV(3,4)*errorNADH4_df2_REV; wDesp_t_REV(4,4)*errorNADH4_df1_REV;
            wDesp_t_REV(1,5)*errorNADH5_df8_REV; wDesp_t_REV(2,5)*errorNADH5_df4_REV; wDesp_t_REV(3,5)*errorNADH5_df2_REV; wDesp_t_REV(4,5)*errorNADH5_df1_REV;
            wDesp_t_REV(1,6)*errorNADH6_df8_REV; wDesp_t_REV(2,6)*errorNADH6_df4_REV; wDesp_t_REV(3,6)*errorNADH6_df2_REV; wDesp_t_REV(4,6)*errorNADH6_df1_REV;
            wDesp_t_REV(1,7)*errorNADH7_df8_REV; wDesp_t_REV(2,7)*errorNADH7_df4_REV; wDesp_t_REV(3,7)*errorNADH7_df2_REV; wDesp_t_REV(4,7)*errorNADH7_df1_REV;
            wDesp_t_REV(1,8)*errorNADH8_df8_REV; wDesp_t_REV(2,8)*errorNADH8_df4_REV; wDesp_t_REV(3,8)*errorNADH8_df2_REV; wDesp_t_REV(4,8)*errorNADH8_df1_REV;
            wDesp_t_REV(1,9)*errorNADH9_df8_REV; wDesp_t_REV(2,9)*errorNADH9_df4_REV; wDesp_t_REV(3,9)*errorNADH9_df2_REV; wDesp_t_REV(4,9)*errorNADH9_df1_REV;
            wDesp_t_REV(1,10)*errorNADH10_df8_REV; wDesp_t_REV(2,10)*errorNADH10_df4_REV; wDesp_t_REV(3,10)*errorNADH10_df2_REV; wDesp_t_REV(4,10)*errorNADH10_df1_REV];
        % Haldande error
        errorHaldane = 0;    
        % Regularization error
        errorReg = lambda * x_temp(1:48)';  
    otherwise
        disp('Warning: No cost function has been selected.');
end
% Overal error (+ weighting for sample size)
error = [
    wD * errorNADH_FWD .* length(errorNADH_REV) ./ length(errorNADH_FWD); % weight for sample size
    wD * errorNADH_REV;
    wH * errorHaldane;
    wL * errorReg;
    ];

end

% %% memoryDump
% figure
% subplot(141), plot(errorNADH_FWD .* length(errorNADH_REV) ./ length(errorNADH_FWD),'.-')
% subplot(142), plot(errorNADH_REV,'.-')
% subplot(143), plot(errorHaldane,'.-')
% subplot(144), plot(errorReg,'.-')
