% % % % function x=dev_pEst_gapdh_fwd
% % PEST_ALD.m
% Parameter estimation for the data in the Enolase assay.
% Vm are estimated changing with pH
% Kms are assumed constant
% Keq from the eQuilibrator


%% (0) Setup and data load
clear, close all
set_paths_pHstudy;
dbstop if error
for step0 = 1
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
    setup.caseStudyENO = 0;
    selectSetup_pH;
    % added
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_gapdh = expData.gapdh;
    % changes in the dataset for gapdh_fwd
        %         excelSheet_corrected: [2×23 double]
        import_gapdh.treatedData.excelSheet_corrected(1,24) = import_gapdh.treatedData.excelSheet_corrected(1,23);
        import_gapdh.treatedData.excelSheet_corrected(2,24) = import_gapdh.treatedData.excelSheet_corrected(2,23);
        %             pH_corrected: [2×23 double]
        import_gapdh.treatedData.pH_corrected(1,24) = import_gapdh.treatedData.pH_corrected(1,23);
        import_gapdh.treatedData.pH_corrected(2,24) = import_gapdh.treatedData.pH_corrected(2,23);
        %       dilution_corrected: [2×23 double]
        import_gapdh.treatedData.dilution_corrected(1,24) = import_gapdh.treatedData.dilution_corrected(1,23);
        import_gapdh.treatedData.dilution_corrected(2,24) = import_gapdh.treatedData.dilution_corrected(2,23);
        %          absorbance_mean: {2×23 cell}
        import_gapdh.treatedData.absorbance_mean{1,24} = import_gapdh.treatedData.absorbance_mean{1,23};
        import_gapdh.treatedData.absorbance_mean{2,24} = import_gapdh.treatedData.absorbance_mean{2,23};
        %       absorbance_samples: {2×23 cell}
        import_gapdh.treatedData.absorbance_samples{1,24} = import_gapdh.treatedData.absorbance_samples{1,23};
        import_gapdh.treatedData.absorbance_samples{2,24} = import_gapdh.treatedData.absorbance_samples{2,23};
        %           absorbance_std: {2×23 cell}
        import_gapdh.treatedData.absorbance_std{1,24} = import_gapdh.treatedData.absorbance_std{1,23};
        import_gapdh.treatedData.absorbance_std{2,24} = import_gapdh.treatedData.absorbance_std{2,23};
        %       concentration_mean: {2×23 cell}
        import_gapdh.treatedData.concentration_mean{1,24} = import_gapdh.treatedData.concentration_mean{1,23};
        import_gapdh.treatedData.concentration_mean{2,24} = import_gapdh.treatedData.concentration_mean{2,23};
        %        concentration_std: {2×23 cell}
        import_gapdh.treatedData.concentration_std{1,24} = import_gapdh.treatedData.concentration_std{1,23};
        import_gapdh.treatedData.concentration_std{2,24} = import_gapdh.treatedData.concentration_std{2,23};
        %                     time: {2×23 cell}
        import_gapdh.treatedData.time{1,24} = import_gapdh.treatedData.time{1,23};
        import_gapdh.treatedData.time{2,24} = import_gapdh.treatedData.time{2,23};
        %            reaction_rate: {2×23 cell}
        import_gapdh.treatedData.reaction_rate{1,24} = import_gapdh.treatedData.reaction_rate{1,23};
        import_gapdh.treatedData.reaction_rate{2,24} = import_gapdh.treatedData.reaction_rate{2,23};
    
    
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
    
    pHarray = unique(import_gapdh.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_gapdh.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_gapdh.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_gapdh.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_gapdh.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_gapdh.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_gapdh.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_gapdh.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_gapdh.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_gapdh.treatedData.reaction_rate{tempID(j)};
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
        % Option 1. Vmax from the values obtained
        Vmax(i) = max(abs(RRs{i}));
%         % Option 2. Vmax naive approach. First datapoints
%         Vmax(i) = (conc_mean{i}(end) - conc_mean{i}(1)) ./ (time{i}(end) - time{i}(1)); 
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
    data.chosenNADHini = 0.15;
    temp1 = import_gapdh.rawData.absorbance_corrected{4,4};
    temp2 = import_gapdh.rawData.absorbance_corrected{5,4};
    temp3 = import_gapdh.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_gapdh.rawData.time{1};
    
        % (1) Correct for minimum value
        % (2) Bring the minimum to zero (apply to all)
        % (3) In principle, use the 3 first dilution rates
        % (4) Watch out with the dilution factors (first 2 cases are
        % reversed)
        
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Directly changing the concentration here, sicne the extinction
    % coefficient did not change.
    dps = length(NADH{1,1});
    for i = 1:DFs
        for j = 1:numpHtested
            dps = length(NADH{j,i});
            k2 = 1:dps;
            for k = 1:dps
                k3 = k2(end-k+1);
%                 switch j
%                     case {1,2,3,4,5,6,7}
%                         NADH{j,i}(k) = NADH{j,i}(k) - NADH{j,DFs}(dps);
%                         NADH{j,i}(k) = NADH{j,i}(k) - NADH{j,DFs}(1);
                        NADH{j,i}(k3) = NADH{j,i}(k3) - NADH{j,i}(1);
%                     case {8,9,10,11,12}
%                         NADH{j,i}(k) = NADH{j,i}(k) - 0.0453;
%     %                     NADH{j,i}(k) = NADH{j,i}(k) - NADH{6,DFs}(dps);
%                 end
            end
        end
    end
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    pHvals = unique(import_gapdh.treatedData.pH_corrected);
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
        if setup.caseStudyENO == 1
            ylim([0 1.5])
        end
        if setup.caseStudyHXK == 1
            ylim([0 0.15])
        end
        if setup.caseStudyALD == 1
            ylim([0 0.15])
        end
        if setup.caseStudyGAPDHr == 1
            ylim([0 0.15])
        end
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
    suptitle(suptitleName);
    
%     figure
%     plot(pHvals, Vmax(:,4),'.-')
%     title('Starting estimate: Vmax [mM s-1] vs pH')    
end


%% (0.1) Experimental Vmax determination

% %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 6; % minimum size of the window
    limRates = [0 1E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.05];  %Ylims plot conc
% select start point (this needs manual selection deppending on previous plots)
dp_start = ones(size(data.conc_mean));
% blank cell total length
total_len = zeros(size(dp_start));
% DFs considered
DFarray = [1/8 1/4 1/2 1/1];
% idxs2consider
idxs2consider = [1 1 1 1;
                1 1 1 1;
                1 1 1 0;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 0];

% Experimental rates determination and plotting
expRatesDetermination;


%% (1.1) Simple parameter fit. Parameter estimation
setup.ode = 'gapdh_fwd_simplified';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';

setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
setup.costfun = 3;

setup.weightData = 1;
idxs2consider2 = ones(size(idxs2consider));

setup.weightDataEsp = idxs2consider2;

setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 1E-5; % by now just testing

% Km fixed
optfun = @costfun_Kmfixed;
plength = 16; % Kms (2) + Vms (1) * numpH (12)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');
setup.Kmvariable = 0;

% %% simple optimization
% % parameter estimation
% tic
% [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t = toc;
% 
% 
% % testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(xres,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;


%% PSA on the effect of Keq for gapdh_reverse
setup.Kmvariable = 0;
setup.KeqPSA = 1;
x_temp = [-0.2096   -2.2655   -0.5194   -0.3304    0.4115    0.3176...
           0.3641    0.3567    0.3358    0.3991    0.4211    0.4853...
           0.5134    0.5195    0.5572    0.5481];
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.pH_Keq_gapdh);
simResults = cell(1,nPSAvals);
for i = 1:nPSAvals
    setup.nPSA = i;
    [error,simResult] = costfun_Kmfixed_simOutput(x_temp,data,setup);
    simResults{i} = simResult;
end

numpH = setup.numpHtested;
for j = 1:numpH   
    simulationVisualization_PSAkeq;
end
% %%
% % if setup.plotOutput == 1
% %     set(1,'color','white'), savefig(1,['results/',setup.enzymeName,'/',setup.enzymeName, '_tempName1.fig']);
% %     set(2,'color','white'), savefig(2,['results/',setup.enzymeName,'/',setup.enzymeName, '_tempName2.fig']);
% % end
% setup.Kmvariable = 0;
% setup.KeqPSA = 0;
% x_temp = [-0.2096   -2.2655   -0.5194   -0.3304    0.4115    0.3176...
%            0.3641    0.3567    0.3358    0.3991    0.4211    0.4853...
%            0.5134    0.5195    0.5572    0.5481];
% % % bypass on Keq
% % temp_keq.GAPDH = setup.pH_Keq_gapdh_eQ_fwd;
% % temp_keq.PGK = setup.pH_Keq_pgk_fwd;
% 
% for i = 1:12
%     setup.pH_Keq_gapdh_eQ_fwd(i) = temp_keq.GAPDH(12);
%     setup.pH_Keq_pgk_fwd(i) = temp_keq.PGK(12);
% end
% 
% % simulate
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error,simResult] = costfun_Kmfixed_simOutput(x_temp,data,setup);




