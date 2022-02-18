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


%% PART 1. ESTIMATION USING REGULARIZATION AND OTHER CONSTRAINTS


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
% visualize: check calculations made
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:numpHtested
    subplot(3,4,i)
    for j = 1:DFs
        plot(time{i,j},NADH{i,j},'.-')
        hold on
    end
    title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
    if i == numpHtested
        if setup.caseStudyGAPDHr == 1
            legend('DF 8','DF 4','DF 2','DF 1')
        end
    end
end
suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
suptitle(suptitleName);


%% (1) System simulation (@ph7) and PSA
setup.ode = 'vanHeerden2014';
% setup.ode = 'gapdhr_s_revMM';
setup.sourceVm = 'experimentalSlopes';
% setup.sourceVm = 'literature'; 
setup.nLinSpace = 21;
setup.PSArefval = 11;
setup.plotResults = 1;
setup.plotEachSim = 0;

%% (2.1) (pre-)Parameter estimation: simulation of all cases using the cost funtion
% In a way, this is a check that keq increases with time right now
setup.DFstudy = 4;
setup.costfun = 1;

setup.ode = 'vanHeerden2014';
setup.ode_pH = 'off';

setup.plotResults = 0;
setup.plotEachSim = 0;
setup.plotEachSimCF = 1;
setup.weightData = 0;
setup.weightHaldane = 0;
setup.selectedLambda = 0;

numpH = numpHtested;
for i = 1:numpH    
    % inputs to be selected
    x_temp = zeros(1,6);
    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(i);
    setup.excessPGK = 1;
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;
        
    [error] = costfun_pH_new(x_temp,data,setup);
end


%% (2.2) Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries
% setup.sourceVm = 'experimentalSlopes';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode = 'vanHeerden2014';
setup.ode_pH = 'on';

setup.DFstudy = 4;
setup.costfun = 1;

setup.plotResults = 0;
setup.plotEachSim = 0;
setup.plotEachSimCF = 0;
setup.weightData = 1; %0.1;
setup.weightHaldane = 1E2; %1;
setup.selectedLambda = 1E-2;
    temp_wD = setup.weightData;
    temp_wH = setup.weightHaldane;

ode_pH = setup.ode_pH;
numpH = numpHtested;
plength = length(setup.params);
pvals = zeros(numpH,plength);
pcis = zeros(numpH,plength);
x_temp = zeros(1,plength);
ub = 3*ones(1,plength); ub([1,6]) = 6;
lb = -3*ones(1,plength); lb([1,6]) = -6;
% ub = 6*ones(1,plength);
% lb = -6*ones(1,plength);
options = optimset('Display','iter');

errorData = zeros(numpH,1);
errorHaldane = zeros(numpH,1);
errorRegpars = zeros(numpH,1);
for i = 1:numpH 
% for i = 1
    % inputs to be selected
    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(i);
    setup.excessPGK = 1;
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;
    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_new,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    % confidence intervals estimated from covar./FIM. Only experimental
    % datapoins are considered for total N, and not regularization.
    lN = length(setup.DFstudy);
    switch lN
        case 1
            N = length(data.NADH{4});
        case 2
            N = length(data.NADH{4}) + length(data.NADH{3});
        case 4
            N = length(data.NADH{4}) + length(data.NADH{3}) + length(data.NADH{2}) + length(data.NADH{1});
        otherwise
            disp('No N has been selected');
    end
    Jacobian = full(Jacobian);  
    varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
    stdp = sqrt(diag(varp));
    pcis(i,:) = stdp; % confidence intervals

    % simulate and plot results
%     setup.weightData = 1;
%     setup.weightHaldane = 1; %0.1; %1;
%     setup.selectedLambda = 0;
    setup.plotEachSimCF = 1; %1;
    setup.simAllProfiles = 0; %1;
    setup.weightData = 1;
    setup.weightHaldane = 1;
    [error] = costfun_pH_new(xres,data,setup);
    setup.weightData = temp_wD;
    setup.weightHaldane = temp_wH;
    setup.simAllProfiles = 0;
    setup.plotEachSimCF = 0;
    % calculating errors
    errorData(i) = sum(abs(error(1:end-7)));
    errorHaldane(i) = sum(abs(error(end-6)));
    errorRegpars(i) = sum(abs(error(end-5,end)));
end

% %% parameter visualization + add confidence intervals
sourceVm = setup.sourceVm;
figure
for i = 1:plength
    % plot parameter values
    subplot(3,3,i)
%     plot([pHvals(1:2);pHvals(4:end)],[pvals(1:2,i);pvals(4:end,i)],'.-')
%     plot(pHvals,pvals(:,i),'.-')
    errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    titleName = setup.params{i};
    title(titleName);
    % plot errors
    if i == plength
        subplot(3,3,i+1)
        plot(pHvals,errorData,'.-')
        hold on
        plot(pHvals,errorHaldane,'.-')
        hold on
        plot(pHvals,errorRegpars,'.-')
        legend('error_{Data}','error_{Haldane}','error_{Regpars}','location','southoutside','orientation','horizontal')
        ylim([0 0.1])
    end
    % plot haldaner relationship
    if i == plength
        Keq_haldane_estimated = zeros(1,numpH);
        for j = 1:numpH
            data.i = j;
            switch sourceVm
                case 'literature'
                    vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                    vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                case 'experimentalSlopes'
                    vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                    vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                case 'experimentalSlopesFixed'
                    vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                    vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                otherwise
                    disp('No source for vmax has been selected');
            end
            ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
            ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
            kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
            kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
            switch ode_pH
                case 'on'
                    H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
%                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                otherwise
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
            end
        end
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        subplot(3,3,i+3)
        semilogy(pHvals,Keq_haldane_estimated)
        hold on
        semilogy(pHvals,Keq_haldane_theory,'k+')
        legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')


%% (2.3) Increasing weight Haldane. WeightData = 1, WeightLambda = 0.01
% setup.sourceVm = 'experimentalSlopes';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode = 'vanHeerden2014';
setup.ode_pH = 'on';
ode_pH = setup.ode_pH;

setup.DFstudy = 4;
setup.costfun = 1;

pvals_cell = cell(10,1);
pcis_cell = cell(10,1);
errorData_cell = cell(10,1);
errorHaldane_cell = cell(10,1);
errorRegpars_cell = cell(10,1);
% weightTest = [0 1E-2 2E-2 1E-1 2E-1 1E0 2E0 1E1 2E1 3E1 3E2];
weightTest = [0 1E-2 1E-1 1E0  1E1  1E2 1E3 1E4 1E5 1E6];
for o = 1:10
% for o = 1
    setup.plotResults = 0;
    setup.plotEachSim = 0;
    setup.plotEachSimCF = 0;
    % % % % setup.weightData = 0.1;%0; %1;
    % % % % setup.weightHaldane = 1; %0.1; %1;
    setup.weightData = 1;%0; %1;
    setup.weightHaldane = weightTest(o); %0.1; %1;
%     setup.selectedLambda = 0;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     setup.weightData = weightTest(o);%0; %1;
%     setup.weightHaldane = 1; %0.1; %1;
    setup.selectedLambda = 0.01;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
        temp_wD = setup.weightData;
        temp_wH = setup.weightHaldane;

    numpH = numpHtested;
    plength = length(setup.params);
    pvals = zeros(numpH,plength);
    pcis = zeros(numpH,plength);
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength); ub([1,6]) = 6;
    lb = -3*ones(1,plength); lb([1,6]) = -6;
    % ub = 6*ones(1,plength);
    % lb = -6*ones(1,plength);
    options = optimset('Display','iter');

    errorData = zeros(numpH,1);
    errorHaldane = zeros(numpH,1);
    errorRegpars = zeros(numpH,1);
    for i = 1:numpH 
    % for i = 1
        % inputs to be selected
        data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
        data.KeqPGK = setup.pH_Keq_pgk(i);
        setup.excessPGK = 1;
        data.NADH = data.conc_mean(i,:);
        data.Vprofs = data.RRs(i,:);
        data.tempTime = data.time(i,:);
        data.i = i;
        % parameter estimation
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_new,x_temp,lb,ub,options,data,setup);
        t = toc;
        pvals(i,:) = xres;
        fprintf('Pest finished for pH #%d, time %d s\n',i,t);
        % confidence intervals estimated from covar./FIM. Only experimental
        % datapoins are considered for total N, and not regularization.
        lN = length(setup.DFstudy);
        switch lN
            case 1
                N = length(data.NADH{4});
            case 2
                N = length(data.NADH{4}) + length(data.NADH{3});
            case 4
                N = length(data.NADH{4}) + length(data.NADH{3}) + length(data.NADH{2}) + length(data.NADH{1});
            otherwise
                disp('No N has been selected');
        end
        Jacobian = full(Jacobian);  
        varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
        stdp = sqrt(diag(varp));
        pcis(i,:) = stdp; % confidence intervals

        % simulate and plot results
        setup.plotEachSimCF = 0; %1;
        setup.simAllProfiles = 0; %1;
        % %     
        setup.weightData = 1;
        setup.weightHaldane = 1;
        % % 
        [error] = costfun_pH_new(xres,data,setup);
        % %     
        setup.weightData = temp_wD;
        setup.weightHaldane = temp_wH;
        % %     
        setup.simAllProfiles = 0;
        setup.plotEachSimCF = 0;   

        % calculating errors
        errorData(i) = sum(abs(error(1:end-7)));
        errorHaldane(i) = sum(abs(error(end-6)));
        errorRegpars(i) = sum(abs(error(end-5,end)));
    end
    pvals_cell{o} = pvals;
    pcis_cell{o} = pcis;
    errorData_cell{o} = errorData;
    errorHaldane_cell{o} = errorHaldane;
    errorRegpars_cell{o} = errorRegpars;

    % %% parameter visualization + add confidence intervals
    sourceVm = setup.sourceVm;
    figure
    for i = 1:plength
        % plot parameter values
        subplot(3,3,i)
        plot(pHvals,pvals(:,i),'.-')
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
        titleName = setup.params{i};
        title(titleName);
        % plot errors
        if i == plength
            subplot(3,3,i+1)
            plot(pHvals,errorData,'.-')
            hold on
            plot(pHvals,errorHaldane,'.-')
            hold on
            plot(pHvals,errorRegpars,'.-')
            legend('error_{Data}','error_{Haldane}','error_{Regpars}','location','southoutside','orientation','horizontal')
%             ylim([0 0.1])
        end
        % plot haldaner relationship
        if i == plength
            Keq_haldane = zeros(1,numpH);
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
%                 Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                switch ode_pH
                    case 'on'
                        H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
    %                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                    otherwise
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
                
            end
            Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated)
            hold on
            semilogy(pHvals,Keq_haldane_theory,'k+')
            legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
        end

    end
    suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')

end


% %% (2.3b) [Visualization all together] Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries

% parameter visualization + add confidence intervals
sourceVm = setup.sourceVm;
c = cool(10);
figure
for i = 1:plength
    
    % plot parameter values
    subplot(3,3,i)
    for o = 1:10
        pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    titleName = setup.params{i};
    title(titleName);
    
    % plot error Data
    if i == plength
        subplot(3,3,i+1)
        for o = 1:10
            errorData = errorData_cell{o};
            plot(pHvals,errorData,'.-','color',c(o,:))
            hold on
        end
        title('error_{Data}')
    end
    
    % plot haldane relationship
    if i == plength
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        Keq_haldane_estimated = zeros(1,numpH);
        for o = 1:10
            pvals = pvals_cell{o};
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
% % % %                 Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                switch ode_pH
                    case 'on'
                        H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
    %                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                    otherwise
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated,'.-','color',c(o,:))
            hold on
            if o == 10
                semilogy(pHvals,Keq_haldane_theory,'k.','MarkerSize',10)
            end
        end
%         legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('Fitting the experimental data or the haldane relationship (k_{eq})')


% %% (2.3c) [Pareto front visualization] Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries
sumErrorHaldane = zeros(10,1);
sumErrorData = zeros(10,1);
for o = 1:10
    sumErrorHaldane(o) = sum(abs(errorHaldane_cell{o}));
    sumErrorData(o) = sum(abs(errorData_cell{o}));
end

figure
plot(sumErrorHaldane,sumErrorData,'k.-')
xlabel('Error Haldane (k_{eq})')
ylabel('Error Data (NADH)')
% adding labels (locations 3 to 6)
for i = 3:6
    hold on
    txt = erase(sprintf('w_{H}=%d, w_{D}=%d',weightTest(i),1),".000000");
%     txt = sprintf('w_{H}=%d, w_{D}=%d',weightTest(i),1);
    if i == 3, sumErrorHaldane(i) = sumErrorHaldane(i)-0.4; end
    text(sumErrorHaldane(i),sumErrorData(i),txt,'color','blue')
end


%% (2.4) Increasing weight Regulatization. WeightData = 1, WeightHaldane = 1E2
% setup.sourceVm = 'experimentalSlopes';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode = 'vanHeerden2014';
setup.ode_pH = 'on';
ode_pH = setup.ode_pH;

setup.DFstudy = 4;
setup.costfun = 1;

pvals_cell = cell(10,1);
pcis_cell = cell(10,1);
errorData_cell = cell(10,1);
errorHaldane_cell = cell(10,1);
errorRegpars_cell = cell(10,1);
% weightTest = [0 1E-2 2E-2 1E-1 2E-1 1E0 2E0 1E1 2E1 3E1 3E2];
% weightTest = [0 1E-2 1E-1 1E0  1E1  1E2 1E3 1E4 1E5 1E6];
weightTest = [0 2E-4 1E-3 2E-3 1E-2 2E-2 1E-1 2E-1 1E0 2E0];
for o = 1:10
% for o = 1
    setup.plotResults = 0;
    setup.plotEachSim = 0;
    setup.plotEachSimCF = 0;
    % % % % setup.weightData = 0.1;%0; %1;
    % % % % setup.weightHaldane = 1; %0.1; %1;
    setup.weightData = 1;%0; %1;
    setup.weightHaldane = 1E2; %0.1; %1;
% % % %     setup.weightHaldane = 1E4; %0.1; %1;
%     setup.selectedLambda = 0;
% % % % % % % % % % % % % % % % % % % % % % %x % % % % % % % % % % % % % % %
%     setup.weightData = weightTest(o);%0; %1;
%     setup.weightHaldane = 1; %0.1; %1;
    setup.selectedLambda = weightTest(o);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
        temp_wD = setup.weightData;
        temp_wH = setup.weightHaldane;

    numpH = numpHtested;
    plength = length(setup.params);
    pvals = zeros(numpH,plength);
    pcis = zeros(numpH,plength);
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength); ub([1,6]) = 6;
    lb = -3*ones(1,plength); lb([1,6]) = -6;
    % ub = 6*ones(1,plength);
    % lb = -6*ones(1,plength);
    options = optimset('Display','iter');

    errorData = zeros(numpH,1);
    errorHaldane = zeros(numpH,1);
    errorRegpars = zeros(numpH,1);
    for i = 1:numpH 
    % for i = 1
        % inputs to be selected
        data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
        data.KeqPGK = setup.pH_Keq_pgk(i);
        setup.excessPGK = 1;
        data.NADH = data.conc_mean(i,:);
        data.Vprofs = data.RRs(i,:);
        data.tempTime = data.time(i,:);
        data.i = i;
        % parameter estimation
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_new,x_temp,lb,ub,options,data,setup);
        t = toc;
        pvals(i,:) = xres;
        fprintf('Pest finished for pH #%d, time %d s\n',i,t);
        % confidence intervals estimated from covar./FIM. Only experimental
        % datapoins are considered for total N, and not regularization.
        lN = length(setup.DFstudy);
        switch lN
            case 1
                N = length(data.NADH{4});
            case 2
                N = length(data.NADH{4}) + length(data.NADH{3});
            case 4
                N = length(data.NADH{4}) + length(data.NADH{3}) + length(data.NADH{2}) + length(data.NADH{1});
            otherwise
                disp('No N has been selected');
        end
        Jacobian = full(Jacobian);  
        varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
        stdp = sqrt(diag(varp));
        pcis(i,:) = stdp; % confidence intervals

        % simulate and plot results
        setup.plotEachSimCF = 0; %1;
        setup.simAllProfiles = 0; %1;
        % %     
        setup.weightData = 1;
        setup.weightHaldane = 1;
        % % 
        [error] = costfun_pH_new(xres,data,setup);
        % %     
        setup.weightData = temp_wD;
        setup.weightHaldane = temp_wH;
        % %     
        setup.simAllProfiles = 0;
        setup.plotEachSimCF = 0;   

        % calculating errors
        errorData(i) = sum(abs(error(1:end-7)));
        errorHaldane(i) = sum(abs(error(end-6)));
        errorRegpars(i) = sum(abs(error(end-5,end)));
    end
    pvals_cell{o} = pvals;
    pcis_cell{o} = pcis;
    errorData_cell{o} = errorData;
    errorHaldane_cell{o} = errorHaldane;
    errorRegpars_cell{o} = errorRegpars;

    % %% parameter visualization + add confidence intervals
    sourceVm = setup.sourceVm;
    figure
    for i = 1:plength
        % plot parameter values
        subplot(3,3,i)
%         plot(pHvals,pvals(:,i),'.-')
        errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
        titleName = setup.params{i};
        title(titleName);
        % plot errors
        if i == plength
            subplot(3,3,i+1)
            plot(pHvals,errorData,'.-')
            hold on
            plot(pHvals,errorHaldane,'.-')
            hold on
            plot(pHvals,errorRegpars,'.-')
            legend('error_{Data}','error_{Haldane}','error_{Regpars}','location','southoutside','orientation','horizontal')
%             ylim([0 0.1])
        end
        % plot haldaner relationship
        if i == plength
            Keq_haldane = zeros(1,numpH);
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
%                 Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                switch ode_pH
                    case 'on'
                        H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
    %                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                    otherwise
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
                
            end
            Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated)
            hold on
            semilogy(pHvals,Keq_haldane_theory,'k+')
            legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
        end

    end
    suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')

end

% %% (2.3b) [Visualization all together] Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries

% parameter visualization + add confidence intervals
sourceVm = setup.sourceVm;
c = cool(10);
figure
for i = 1:plength
    
    % plot parameter values
    subplot(3,3,i)
    for o = 1:10
        pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    titleName = setup.params{i};
    title(titleName);
    
    % plot error Data
    if i == plength
        subplot(3,3,i+1)
        for o = 1:10
            errorData = errorData_cell{o};
            plot(pHvals,errorData,'.-','color',c(o,:))
            hold on
        end
        title('error_{Data}')
    end
    
    % plot haldane relationship
    if i == plength
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        Keq_haldane_estimated = zeros(1,numpH);
        for o = 1:10
            pvals = pvals_cell{o};
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
% % % %                 Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                switch ode_pH
                    case 'on'
                        H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
    %                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                    otherwise
                        Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated,'.-','color',c(o,:))
            hold on
            if o == 10
                semilogy(pHvals,Keq_haldane_theory,'k.','MarkerSize',10)
            end
        end
%         legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('Fitting the experimental data or the haldane relationship (k_{eq})')


% %% (2.3c) [Pareto front visualization] Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries
sumErrorHaldane = zeros(10,1);
sumErrorData = zeros(10,1);
for o = 1:10
    sumErrorHaldane(o) = sum(abs(errorHaldane_cell{o}));
    sumErrorData(o) = sum(abs(errorData_cell{o}));
end

figure
plot(sumErrorHaldane,sumErrorData,'k.-')
xlabel('Error Haldane (k_{eq})')
ylabel('Error Data (NADH)')
% adding labels (locations 3 to 6)
for i = 3:6
    hold on
    txt = erase(sprintf('w_{H}=%d, w_{D}=%d',weightTest(i),1),".000000");
%     txt = sprintf('w_{H}=%d, w_{D}=%d',weightTest(i),1);
    if i == 3, sumErrorHaldane(i) = sumErrorHaldane(i)-0.4; end
    text(sumErrorHaldane(i),sumErrorData(i),txt,'color','blue')
end


%% taking the values in figure 5, recalculating parameter values to compare with literature
% wData = 1, wHaldane = 1E2, wLambda = 1E-2.	
pvals_mMs = zeros(size(pvals_cell{5}));
pvals_mMs(:,1) = 10 .^ pvals_cell{5}(:,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1} %vmf
pvals_mMs(:,2) = 10 .^ pvals_cell{5}(:,2) .* 2.48; % mM %ks1, k_gap
pvals_mMs(:,4) = 10 .^ pvals_cell{5}(:,4) .* 2.92; % mM %ks2, k_nad
pvals_mMs(:,3) = 10 .^ pvals_cell{5}(:,3) .* 1.18; % mM %kp1, k_bpg
pvals_mMs(:,5) = 10 .^ pvals_cell{5}(:,5) .* 0.022; % mM %kp2, k_nadh
pvals_mMs(:,6) = 10 .^ pvals_cell{5}(:,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1} %vmr

figure
for i = 1:plength
    
    % plot parameter values
    subplot(2,3,i)
    for o = 1:numpHtested
        plot(pHvals,pvals_mMs(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    titleName = setup.params{i};
    title(titleName);
end

figure
for i = 1:plength
    
    % plot parameter values
    subplot(2,3,i)
    for o = 1:numpHtested
        switch i
            case {2,3,4,5}
                plot(pHvals,pvals_mMs(:,i),'.-','color',c(o,:))
            case {1,6}
                plot(pHvals,pvals_mMs(:,i)*60,'.-','color',c(o,:))
        end
        hold on
        %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    switch i
        case {2,3,4,5}
            titleName = setup.params{i};
            title(titleName);
        case 1
            titleName = 'v_{maxFWD} [mM min^{-1}]';
            title(titleName);
        case 6
            titleName = 'v_{maxREV} [mM min^{-1}]';
            title(titleName);
    end
end


%% PART 2. ESTIMATION USING VM constant and Keq constant idea
% Use of 'costfun_pH_Kmconstant' and function for vm. Based on structure in
% 'pEst_gapdhr.m'.

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
    
setup.weightData = 1;
setup.weightDataEsp = ones(1,10);
% setup.weightDataEsp = zeros(1,10); setup.weightDataEsp(8) = 1;
setup.weightHaldane = 0; %(change here)
setup.selectedLambda = 0;

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

% % Vm fixed
% optfun = @costfun_Vmfixed;
% plength = 42; % Kms (4) + Vms (2) * numpH (10)
% x_temp = zeros(1,plength);
% ub = 3*ones(1,plength);
% lb = -3*ones(1,plength);
% ub([1,2]) = 6;
% lb([1,2]) = -6;
% options = optimset('Display','iter');

% Km fixed
optfun = @costfun_Kmfixed;
plength = 24; % Kms (4) + Vms (2) * numpH (10)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
ub(5:end) = 6;
lb(5:end) = -6;
options = optimset('Display','iter');

tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t = toc;

%%
% simulate and plot results
setup.plotEachSimCF = 1; %1;
setup.simAllProfiles = 1; %1;
setup.DFstudy = 1:4;
    
temp_wD = setup.weightData;
temp_wH = setup.weightHaldane;
temp_wL = setup.selectedLambda;
setup.weightData = 1;
setup.weightHaldane = 1;
setup.selectedLambda = 1;

[error] = optfun(xres,data,setup);
    
setup.weightData = temp_wD;
setup.weightHaldane = temp_wH;
setup.selectedLambda = temp_wL;
    
setup.simAllProfiles = 0;
setup.plotEachSimCF = 0;  
        
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
weightTest = [1E-8 1E-7 1E-6 1E-5 1E-4 3E-4 1E-3 3E-3 1E-2 3E-2 1E-1 3E-1 1E0 3E0 1E1 3E1 1E2 3E2 1E3 3E3 1E4 3E4 1E5 1E6 1E7 1E8 1E9 1E10];
% weightTest = [1E-8 1E-7 1E-6 1E-5 1E-4 3E-4]; %#1
% weightTest = [1E-3 3E-3 1E-2 3E-2 1E-1 3E-1]; %#2
% weightTest = [1E0 1E1 1E2 1E3 1E4 1E5]; %#3
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
% options = optimset('Display','iter');
% weightTest = [1E-6 3E-6 1E-5 3E-5 1E-4 3E-4 1E-3 3E-3 1E-2 3E-2 1E-1 3E-1 1E0 3E0 1E1 3E1 1E2 3E2 1E3 3E3 1E4 3E4 1E5 3E5 1E6 3E6 1E8]; 
% % weightTest = [1E-4 3E-4 1E-3 3E-3 1E-2 3E-2]; %#1
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
    % errorRegpars_cell(o) = 1; % not needed
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

figure
yyaxis left
semilogx(weightTest,eHaldane,'o-')
% ax1 = gca;
hold on
yyaxis right
semilogx(weightTest,eData,'o-')
% ax2 = gca;
legend('e_{Haldane}','e_{Data}')
% suptitle('All variable: Haldane Regularization')
suptitle('Km variable, Vm fixed: Haldane Regularization')
% suptitle('Vm variable, Km fixed: Haldane Regularization')

% mins = [ax1.YLim(1), ax2.YLim(1)];
% maxs = [ax1.YLim(2), ax2.YLim(2)];
% line([0.01 0.01],[min(mins) max(maxs)],'Color','black','LineStyle','--')
% if i == plength
%     hL = legend('errorData','errorParameters','lambda = 0.01');
%     newPosition = [0.7 0.1 0.2 0.2];
%     newUnits = 'normalized';
%     set(hL,'Position', newPosition,'Units', newUnits);
%     hold off
% end

% figure,
% bar(errorData_cell{1})
% hold on
% bar(errorData_cell{10})
% legend('1','10')



