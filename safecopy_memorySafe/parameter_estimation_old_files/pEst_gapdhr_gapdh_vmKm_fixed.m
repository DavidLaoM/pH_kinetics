% % PEST_GAPDHR_GAPDH.M 
% Parameter estimation for GAPDH was tested with the gapdh_r dataset in the
% code 'pEst_gapdhr_updatedTopology'. In this code, the exercise is 
% supplemented with data for the GAPDH (forward sense) reaction.

% The idea is to see if the parameters can be well-quantified with both
% data sets at the same time and see if this adds information to the
% problem. In principle, the Haldane constraint will not be considered (the
% weight will be set to zero) since the Keq approach was not solved in the
% previous assay, and needs to be further discussed.

% Sections
% (0) Setup: The data sets to be used are loaded.
% (1) Test for simulation file
% (2) Test for parameter the cost function
% (3) Parameter estimation: no regularization
% (4) Parameter estimation: added regularization


%% (0) Setup
clear
set_paths_pHstudy;
dbstop if error

%% (0.1) Setup and data load: gapdh_forward
for reducingVisualLines = 1
    % select specific case and recall data
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
    setup.saveOutput = 0;

    load('expData.mat','expData');
    import_gapdhF_temp = expData.gapdh;
    % avoid extra 2 pH values
    import_gapdhF.rawData = import_gapdhF_temp.rawData;
    
    import_gapdhF.treatedData.excelSheet_corrected = import_gapdhF_temp.treatedData.excelSheet_corrected(:,1:20);
    import_gapdhF.treatedData.pH_corrected = import_gapdhF_temp.treatedData.pH_corrected(:,1:20);
    import_gapdhF.treatedData.dilution_corrected = import_gapdhF_temp.treatedData.dilution_corrected(:,1:20);
    import_gapdhF.treatedData.absorbance_mean = import_gapdhF_temp.treatedData.absorbance_mean(:,1:20);
    import_gapdhF.treatedData.absorbance_samples = import_gapdhF_temp.treatedData.absorbance_samples(:,1:20);
    import_gapdhF.treatedData.absorbance_std = import_gapdhF_temp.treatedData.absorbance_std(:,1:20);
    import_gapdhF.treatedData.concentration_mean = import_gapdhF_temp.treatedData.concentration_mean(:,1:20);
    import_gapdhF.treatedData.concentration_std = import_gapdhF_temp.treatedData.concentration_std(:,1:20);
    import_gapdhF.treatedData.time = import_gapdhF_temp.treatedData.time(:,1:20);
    import_gapdhF.treatedData.reaction_rate = import_gapdhF_temp.treatedData.reaction_rate(:,1:20);
    import_gapdhF.treatedData.unitsTime = import_gapdhF_temp.treatedData.unitsTime;
    import_gapdhF.treatedData.unitsAbsorbance = import_gapdhF_temp.treatedData.unitsAbsorbance;
    import_gapdhF.treatedData.unitsConcentration = import_gapdhF_temp.treatedData.unitsConcentration;
    import_gapdhF.treatedData.unitsRates = import_gapdhF_temp.treatedData.unitsRates;
    import_gapdhF.treatedData.concProtein = import_gapdhF_temp.treatedData.concProtein;
    import_gapdhF.treatedData.unitsProtein = import_gapdhF_temp.treatedData.unitsProtein;
    import_gapdhF.treatedData.protDF = import_gapdhF_temp.treatedData.protDF;
    
%     DFs = setup.DFactorsTotal;
    DFs = 4;
%     pHtested = setup.pHtested;
    pHtested = [1 1 0 1 1 1 1 1 1 1 1 0];
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

    pHarray = unique(import_gapdhF.treatedData.pH_corrected);
    for i = 1:numpHtested
%         disp(i);
        pHval = pHarray(i);
        tempID = find(import_gapdhF.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_gapdhF.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_gapdhF.treatedData.dilution_corrected(tempID);
        for j = 1:4
%             disp(j);
            abs_meanTemp{j,i} = import_gapdhF.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_gapdhF.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_gapdhF.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_gapdhF.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_gapdhF.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_gapdhF.treatedData.reaction_rate{tempID(j)};
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
    data_gapdhF.pH = pH;
    data_gapdhF.DF = DF;
    data_gapdhF.abs_mean = abs_mean;
    data_gapdhF.abs_std = abs_std;
        conc_mean{1,1} = conc_mean{1,1} * 2;
        conc_mean{1,2} = conc_mean{1,2} * 2;
        conc_mean{1,3} = conc_mean{1,3} * 2;
        conc_mean{1,4} = conc_mean{1,4} * 2;
        conc_mean{2,1} = conc_mean{2,1} * 2;
        conc_mean{2,2} = conc_mean{2,2} * 2;
        conc_mean{2,3} = conc_mean{2,3} * 2;
        conc_mean{2,4} = conc_mean{2,4} * 2;
    data_gapdhF.conc_mean = conc_mean;              % 2
        conc_std{1,1} = conc_std{1,1} * 2;
        conc_std{1,2} = conc_std{1,2} * 2;
        conc_std{1,3} = conc_std{1,3} * 2;
        conc_std{1,4} = conc_std{1,4} * 2;
        conc_std{2,1} = conc_std{2,1} * 2;
        conc_std{2,2} = conc_std{2,2} * 2;
        conc_std{2,3} = conc_std{2,3} * 2;
        conc_std{2,4} = conc_std{2,4} * 2;
    data_gapdhF.conc_std = conc_std;                % 2
    data_gapdhF.time = time;
        RRs{1,1} = RRs{1,1} * 2;
        RRs{1,2} = RRs{1,2} * 2;
        RRs{1,3} = RRs{1,3} * 2;
        RRs{1,4} = RRs{1,4} * 2;
        RRs{2,1} = RRs{2,1} * 2;
        RRs{2,2} = RRs{2,2} * 2;
        RRs{2,3} = RRs{2,3} * 2;
        RRs{2,4} = RRs{2,4} * 2;
    data_gapdhF.RRs = RRs;                          % 2
        Vmax(1,1) = Vmax(1,1) * 2;
        Vmax(1,2) = Vmax(1,2) * 2;
        Vmax(1,3) = Vmax(1,3) * 2;
        Vmax(1,4) = Vmax(1,4) * 2;
        Vmax(2,1) = Vmax(2,1) * 2;
        Vmax(2,2) = Vmax(2,2) * 2;
        Vmax(2,3) = Vmax(2,3) * 2;
        Vmax(2,4) = Vmax(2,4) * 2;
    data_gapdhF.Vmax = Vmax;                        % 2
    data_gapdhF.chosenVmax = max(max(Vmax));
    data_gapdhF.chosenNADini = 0.05;
    temp1 = import_gapdhF.rawData.absorbance_corrected{4,4};
    temp2 = import_gapdhF.rawData.absorbance_corrected{5,4};
    temp3 = import_gapdhF.rawData.absorbance_corrected{6,4};
    data_gapdhF.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data_gapdhF.raw.time = import_gapdhF.rawData.time{1};

    pHvals = unique(import_gapdhF.treatedData.pH_corrected);
    % visualize: check calculations made
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numpHtested
        subplot(3,4,i)
        for j = 1:DFs
            plot(time{i,j},NADH{i,j},'.-')
            hold on
        end
        title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
        if((i == 1)||(i == 2))
            text(200,0.039,'DF.2,4,8,16');
        end
        xlim([0 600])
        ylim([0 0.2])
        if i == numpHtested
            if setup.caseStudyGAPDHr == 1
                legend('DF 8','DF 4','DF 2','DF 1')
            end
        end
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
    suptitle(suptitleName);
end

%% (0.2) Setup and data load: gapdh_reverse
for reducingVisualLines = 1
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
    data_gapdhR.pH = pH;
    data_gapdhR.DF = DF;
    data_gapdhR.abs_mean = abs_mean;
    data_gapdhR.abs_std = abs_std;
    data_gapdhR.conc_mean = conc_mean;
    data_gapdhR.conc_std = conc_std;
    data_gapdhR.time = time;
    data_gapdhR.RRs = RRs;
    data_gapdhR.Vmax = Vmax;
    data_gapdhR.chosenVmax = max(max(Vmax));
    data_gapdhR.chosenNADini = 0.15;
    temp1 = import_gapdhR.rawData.absorbance_corrected{4,4};
    temp2 = import_gapdhR.rawData.absorbance_corrected{5,4};
    temp3 = import_gapdhR.rawData.absorbance_corrected{6,4};
    data_gapdhR.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data_gapdhR.raw.time = import_gapdhR.rawData.time{1};

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
end

%% (10a) Preliminary setup
data.gapdhF = data_gapdhF;
data.gapdhR = data_gapdhR;
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopesFixed';
% % % % setup.ode_pH = 'on';
% setup.ode_pH = 'on_vmf_vmr';
% setup.ode_pH = 'on_vm_keq';
setup.ode_pH = 'on_revMM_bitri';
setup.typeVm = 'specific';
typeVm = setup.typeVm;
setup.params =  {'v_{maxFWD}^{REVdata} [mM s^{-1}]'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'v_{maxREV}^{REVdata} [mM s^{-1}]'; 'v_{maxFWD}^{FWDdata} [mM s^{-1}]'; 'v_{maxREV}^{FWDdata} [mM s^{-1}]'};
% setup.DFstudy = 4;
setup.DFstudy = [3 4];
setup.costfun = 1;

% no intermediate simulations, only reverse experiment data ued for
% optimization
setup.plotResults = 0;
setup.plotEachSim = 0;
setup.plotEachSimCF = 0;
setup.weightDataR = 1;
setup.weightDataF = 0.1; % 0.1
setup.weightHaldane = 0;
setup.selectedLambda = 0.001; %0;
    temp_wDR = setup.weightDataR;
    temp_wDF = setup.weightDataF;
    temp_wH = setup.weightHaldane;
    temp_wL = setup.selectedLambda;

o = 10;
ode_pH = setup.ode_pH;
numpH = numpHtested;
plength = length(setup.params);
pvals = zeros(numpH,plength);
pcis = zeros(numpH,plength);
x_temp = zeros(1,plength);
% % % % ub = 1*ones(1,plength); ub([1,6,7,8]) = 3; % ub(6) = 6;
% % % % lb = -1*ones(1,plength); lb([1,6,7,8]) = -3;
ub = 1*ones(1,plength); ub([1,6,7,8]) = 6;
lb = -1*ones(1,plength); lb([1,6,7,8]) = -6;
options = optimoptions('lsqnonlin','Display','iter');
% % % % % ub = 3*ones(1,plength); ub([1,6]) = 6;
% % % % % lb = -3*ones(1,plength); lb([1,6]) = -6;
% % % % % % ub = 1*ones(1,plength); ub([1,6]) = 10;
% % % % % % lb = -1*ones(1,plength); lb([1,6]) = -10;
% ub = 1*ones(1,plength); ub(1) = 6; ub(6) = 15;
% lb = -1*ones(1,plength); lb(1) = 0; lb(6) = 0;
% % % % ub = 1*ones(1,plength); ub(1) = 6; ub(6) = 15;
% % % % lb = -1*ones(1,plength); lb(1) = 2; lb(6) = 2;

%% quick simulations
xtemp = zeros(8,1);
data.chosenDF = 1;
data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(6);
data.chosenKeqPGK = setup.pH_Keq_pgk(6);
data.i = 6;
setup.excessPGK = 1;
setup.plotEachSim = 1;
    xtemp(1) = -3;
[testSim] = simSys_FandR(xtemp,data,setup);

%% (10a) Specific Vm. Estimation all parameters variable
% setup.plotEachSimCF = 1;
tic
[xres_temp,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_FandR_fixed,x_temp,lb,ub,options,data,setup);
t = toc;
disp(xres_temp);
% setup.plotEachSimCF = 0;
% %%
setup.plotEachSimCF = 1;
[error] = costfun_pH_FandR_fixed(xres_temp,data,setup);
setup.plotEachSimCF = 0;
% disp(error(1:26));

%% (10b) Specific Vm. Vms constant value estimated + (10c) Specific Vm. KMs constant value estimated
% x_temp = [3.4869   -0.5192    0.1867   -0.5119   -0.1666    3.4291         0         0];
% x_temp = [1.9694   -0.9122    0.0510   -0.8527    0.9665    2.3869   -1.1284   -1.1742];
% x_temp = [1.9684   -0.9993    0.3829   -0.9773    0.9914    2.7449   -1.1372   -2.7234];
x_temp = [2.7481   -0.4878   -0.2746   -0.5066    0.6720    1.9146   -0.2367    2.8005];

errorDataReverse = zeros(numpH,1);
errorDataForward = zeros(numpH,1);
errorHaldane = zeros(numpH,1);
errorRegpars = zeros(numpH,1);
for i = 1:numpH 
    % inputs to be selected
    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(i);
    setup.excessPGK = 1;
    data.gapdhR.NADH = data.gapdhR.conc_mean(i,:);
    data.gapdhR.Vprofs = data.gapdhR.RRs(i,:);
    data.gapdhR.tempTime = data.gapdhR.time(i,:);
    data.gapdhF.NADH = data.gapdhF.conc_mean(i,:);
    data.gapdhF.Vprofs = data.gapdhF.RRs(i,:);
    data.gapdhF.tempTime = data.gapdhF.time(i,:);
    data.i = i;
    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_FandR,x_temp,lb,ub,options,data,setup);
    t = toc;
%     xres = x_temp;
    pvals(i,:) = xres;
% % % %     xres = pvals(i,:);
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    % confidence intervals estimated from covar./FIM. Only experimental
    % datapoins are considered for total N, and not regularization.
    lN = length(setup.DFstudy);
    switch lN
        case 1
            N = length(data.gapdhR.NADH{4});
        case 2
            N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3});
        case 4
            N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3}) + length(data.gapdhR.NADH{2}) + length(data.gapdhR.NADH{1});
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
    setup.weightDataR = 1;
    setup.weightDataF = 1;
    setup.weightHaldane = 1;
    [error] = costfun_pH_FandR(xres,data,setup);
%     [error] = costfun_pH_FandR(xres_temp,data,setup);
    setup.weightDataR = temp_wDR;
    setup.weightDataF = temp_wDF;
    setup.weightHaldane = temp_wH;
    setup.simAllProfiles = 0;
    setup.plotEachSimCF = 0;
    % calculating errors
    errorData(i) = sum(abs(error(1:end-7)));
    errorDataReverse(i) = sum(abs(error(1:26)));
% % % %     errorDataReverse(i) = sum(abs(error(1:52)));
    errorDataForward(i) = sum(abs(error(27:end-7)));
    errorHaldane(i) = sum(abs(error(end-6)));
    errorRegpars(i) = sum(abs(error(end-5,end)));
end

tempError = sum(abs(errorDataReverse)); disp(tempError);
% pvals(:,1) = 1.8620; pvals(:,6) = 1.9253; pvals(:,7) = -0.6412; pvals(:,8) = -0.0939; pvals_VMconst = pvals;
% pvals(:,2) = -0.9592; pvals(:,3) = -0.4083; pvals(:,4) = -0.9628; pvals(:,5) = 0.9151; pvals(:,7) = -0.6412; pvals(:,8) = -0.0939; pvals_KMconst = pvals;

% %% parameter visualization + add confidence intervals
sourceVm = setup.sourceVm;
figure
for i = 1:plength
    % plot parameter values
    subplot(4,3,i)
    plot(pHvals,pvals(:,i),'.-')
%         errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    titleName = setup.params{i};
    title(titleName);
    % plot errors
    if i == plength
        subplot(4,3,i+1)
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
            if setup.ode_pH == 'on'
                H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
            else
                Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
            end
        end
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        subplot(4,3,i+3)
        semilogy(pHvals,Keq_haldane_estimated)
        hold on
        semilogy(pHvals,Keq_haldane_theory,'k+')
        legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end

end
suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')

figure
for i = 1:plength
    % plot parameter values
    subplot(4,3,i)
    errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    titleName = setup.params{i};
    title(titleName);
    ylim([-3 3])
    line([6 8],[1 1],'Color','black','LineStyle','--')
    line([6 8],[-1 -1],'Color','black','LineStyle','--')
end

%%

% % REV data optimization
% errorNonReg = [0.0303, 0.0337, 0.0270, 0.0215, 0.0038, 0.0048, 0.0040, 0.0090, 0.0048, 0.0029];
% errorReg = [0.0340, 0.0351, 0.0312, 0.0216, 0.0108, 0.0063, 0.0026, 0.0043, 0.0033, 0.0034];
% errorKmconst = [0.0546, 0.0616, 0.0557, 0.0410, 0.0513, 0.0479, 0.0351, 0.0269, 0.0178, 0.0122];
% errorVmconst = [0.0485, 0.0563, 0.0559, 0.0344, 0.0214, 0.0107, 0.0099, 0.0084, 0.0015, 0.0032];
% 
% REF & FWD optimization
errorNonReg = [0.0366, 0.0311, 0.0292, 0.0129, 0.0027, 0.0040, 0.0029, 0.0067, 0.0149, 0.0144];
errorReg = [0.0341, 0.0358, 0.0320, 0.0186, 0.0075, 0.0049, 0.0037, 0.0043, 0.0035, 0.0040];
errorKmconst = [0.0449, 0.0505, 0.0447, 0.0321, 0.0309, 0.0242, 0.0132, 0.0065, 0.0023, 0.0049];
errorVmconst = [0.0404, 0.0429, 0.0329, 0.0212, 0.0151, 0.0366, 0.0135, 0.0213, 0.0277, 0.0314];

figure,
plot(pHvals,errorNonReg,'.-')
hold on
plot(pHvals,errorReg,'.-')
hold on
plot(pHvals,errorKmconst,'.-')
hold on
plot(pHvals,errorVmconst,'.-')
legend('errorNonReg','errorReg','errorKmconst','errorVmconst')
xlabel('pH value')
ylabel('error_{NADH}')

% % comparing models
% errorNonReg = [0.0303, 0.0337, 0.0270, 0.0215, 0.0038, 0.0048, 0.0040, 0.0090, 0.0048, 0.0029];
% errorSimpleModel1 = [0.0528, 0.0592, 0.0533, 0.0392, 0.0471, 0.0429, 0.0303, 1.7444, 0.0132, 0.0085];
% errorSimpleModel2 = [0.1940, 0.2159, 0.1952, 0.1605, 0.0471, 0.4409, 0.4815, 0.0223, 0.0131, 0.0085];
% errorSimpleModel3 = [0.0287, 0.0331, 0.0259, 0.0192, 1.3461, 0.0127, 0.0078, 0.0022, 0.0146, 0.0063];
% 
% figure,
% plot(pHvals,errorNonReg,'.-','Linewidth',1.2,'MarkerSize',12)
% hold on
% plot([pHvals(1:7); pHvals(9:10)],[errorSimpleModel1(1:7) errorSimpleModel1(9:10)],'.-','Linewidth',1.2,'MarkerSize',12)
% hold on
% plot(pHvals,errorSimpleModel2,'.-','Linewidth',1.2,'MarkerSize',12)
% hold on
% plot([pHvals(1:4); pHvals(6:10)],[errorSimpleModel3(1:4) errorSimpleModel3(6:10)],'.-','Linewidth',1.2,'MarkerSize',12)
% legend('errorNonReg','error_{Model 1}','error_{Model 2}','error_{Model 3}','location','SouthOutside','Orientation','Horizontal')
% xlabel('pH value')
% ylabel('error_{NADH}')
% ylim([0 0.25])

% tempError2 = errorDataReverse + errorDataForward;
% pvals_KMconst = pvals;
% error_KMconst78 = tempError2;
% pvals_VMconst = pvals;
% error_VMconst78 = tempError2;
% % pvals_Allconst = xres_temp;
% % error_Allconst = error;
% save('resultsFixed_R&F_pvals_error.mat','pvals_Allconst','pvals_VMconst','pvals_KMconst','error_Allconst','error_VMconst78','error_KMconst78');

%% plotting
% load('resultsFixed_pvals_error.mat');
load('resultsFixed_R&F_pvals_error.mat');

figure
for i = 1:8
    subplot(3,3,i)  
    plot(pHvals, pvals_Allconst(i)*ones(size(pHvals)),'.-')
    hold on
    plot(pHvals, pvals_VMconst(:,i),'.-')
    hold on
    plot(pHvals, pvals_KMconst(:,i),'.-')
    title(setup.params{i})
end
subplot(3,3,9)
plot(pHvals, error_Allconst,'.-')
hold on
plot(pHvals, error_VMconst78,'.-')
hold on
plot(pHvals, error_KMconst78,'.-')
legend('All constant','vm constant','kM constant')
title('Error profiles')


suptitle('resultsFixed^{pvals}_{error}');


%% test

errArr = zeros(numpH,1);
for i = 1:numpH 
    % inputs to be selected
    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(i);
    setup.excessPGK = 1;
    data.gapdhR.NADH = data.gapdhR.conc_mean(i,:);
    data.gapdhR.Vprofs = data.gapdhR.RRs(i,:);
    data.gapdhR.tempTime = data.gapdhR.time(i,:);
    data.gapdhF.NADH = data.gapdhF.conc_mean(i,:);
    data.gapdhF.Vprofs = data.gapdhF.RRs(i,:);
    data.gapdhF.tempTime = data.gapdhF.time(i,:);
    data.i = i;

    setup.plotEachSimCF = 1; %1;
    setup.simAllProfiles = 1; %1;
    setup.weightDataR = 1;
    setup.weightDataF = 1;
    setup.weightHaldane = 1;
    [error] = costfun_pH_FandR(xres2,data,setup);
%     [error] = costfun_pH_FandR(xres_temp,data,setup);
    setup.weightDataR = temp_wDR;
    setup.weightDataF = temp_wDF;
    setup.weightHaldane = temp_wH;
    setup.simAllProfiles = 0;
    setup.plotEachSimCF = 0;
 
    errArr(i) = sum(abs(error(1:26)));

end
disp(sum(errArr));
