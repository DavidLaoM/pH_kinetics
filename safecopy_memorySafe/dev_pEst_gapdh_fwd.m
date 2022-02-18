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
setup.plotOutput = 0;
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

% %% (0.3) saving
if setup.plotOutput == 1
    save(['results/',setup.enzymeName,'/',setup.enzymeName, '_initial_variables.mat'],'Vmax_mw_opt_corr','idxs2consider','DF','pH');
    set(1,'color','white'), savefig(1,['results/',setup.enzymeName,'/',setup.enzymeName, '_concentrations_basezero.fig']);
    set(2,'color','white'), savefig(2,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_df.fig']);
    set(3,'color','white'), savefig(3,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_movingWindow.fig']);
    set(4,'color','white'), savefig(4,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_R2_vs_movingWindow.fig']);
    set(5,'color','white'), savefig(5,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_iniPoint_vs_movingWindow.fig']);
    set(6,'color','white'), savefig(6,['results/',setup.enzymeName,'/',setup.enzymeName, '_experimental_vmax_vs_pH.fig']);
end


%% (1.1) Simple parameter fit. Parameter estimation
setup.ode = 'gapdh_fwd_simplified';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';
setup.caseKm = 'pH_independent';
% setup.caseKm = 'pH_dependent';

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
caseKm = setup.caseKm;
switch caseKm
    case 'pH_independent'
        optfun = @costfun_Kmfixed;
        plength = 16; % Kms (4) + Vms (1) * numpH (12)
        x_temp = zeros(1,plength);
        ub = 3*ones(1,plength);
        lb = -3*ones(1,plength);
        options = optimset('Display','iter');
    case 'pH_dependent'
        optfun = @costfun_Kmfixed;
        plength = 60; % Kms (4) * numpH (12) + Vms (1) * numpH (12)
        x_temp = zeros(1,plength);
        ub = 3*ones(1,plength);
        lb = -3*ones(1,plength);
        options = optimset('Display','iter');
    otherwise
        disp('Warning: no specification has been made on Km being pH dependent or independent');
end
        
%%
% % changing constants
    % original
%     setup.pH_Keq_gapdh_eQ_fwd = [0.000085000000000   0.000115000000000 ...
%        0.000145000000000   0.000200000000000...
%        0.000400000000000   0.000810000000000...
%        0.001730000000000   0.003280000000000...
%        0.005800000000000   0.008900000000000...
%        0.012250000000000   0.015250000000000];
%     setup.pH_Keq_pgk_fwd = 1.0e+03 * [1.351351351351352   1.388888888888889...
%        1.408450704225352   1.449275362318841...
%        1.538461538461539   1.639344262295082...
%        1.754385964912281   1.818181818181818...
%        1.886792452830189   1.923076923076923...
%        1.923076923076923   1.960784313725490];
    % lower equilibrium constants -> slower reactions
        % 1E-1 no change
        % 1E-2 no change
        % 1E-3 no change (already some change in x_temp2)
        % 1E-4 change starts
        % 1E-5 really slow
    % higher equilibrium constants -> faster reactions (it's already at the
    % fastest)
        % 1E3 no change
        % 1E5 no change
        % 1E10 no change
%     setup.pH_Keq_gapdh_eQ_fwd = 1E0 * [0.000085000000000   0.000115000000000 ...
%        0.000145000000000   0.000200000000000...
%        0.000400000000000   0.000810000000000...
%        0.001730000000000   0.003280000000000...
%        0.005800000000000   0.008900000000000...
%        0.012250000000000   0.015250000000000];
%     setup.pH_Keq_pgk_fwd = 1E0 * 1.0e+03 * [1.351351351351352   1.388888888888889...
%        1.408450704225352   1.449275362318841...
%        1.538461538461539   1.639344262295082...
%        1.754385964912281   1.818181818181818...
%        1.886792452830189   1.923076923076923...
%        1.923076923076923   1.960784313725490];
% fixed constans
% setup.pH_Keq_gapdh_eQ_fwd = setup.pH_Keq_gapdh_eQ_fwd(6) .* ones(size(setup.pH_Keq_gapdh_eQ_fwd));
% setup.pH_Keq_pgk_fwd = setup.pH_Keq_pgk_fwd(6) .* ones(size(setup.pH_Keq_pgk_fwd));
% x_temp2 = [-0.2096   -2.2655   -0.5194   -0.3304    0.4115    0.3176...
%            0.3641    0.3567    0.3358    0.3991    0.4211    0.4853...
%            0.5134    0.5195    0.5572    0.5481];
% % testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(x_temp,data,setup);
% % [error] = optfun(x_temp2,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;

%% Parameter fit: all parameters (vmax and Km) variable
% % % % % check ifthe fit can improve when Km can change
% % % % 
% % % % % parameter estimation setup
% % % % setup.Kmvariable = 1;
% % % % plength = 60;
% % % % x_temp = zeros(1,plength);
% % % % ub = 3*ones(1,plength);
% % % % lb = -3*ones(1,plength);
% % % %     
% % % % % parameter estimation
% % % % tic
% % % % [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% % % % t = toc;
% % % % 
% % % % % simulation
% % % % setup.plotEachSimCF = 1;
% % % % setup.plotEachSim = 0;
% % % % setup.simAllProfiles = 1;
% % % % [error] = optfun(xres,data,setup);
% % % % setup.plotEachSimCF = 0;
% % % % setup.plotEachSim = 0;
% % % % setup.simAllProfiles = 0;


%% Multi Start parameter fit: Vm fixed, Km variable
% % % % % parameter estimation: with MultiStart
% % % %     % --> cannot use the 'setup' as input. Create own MultiStart
% % % %     % problem = createOptimProblem('lsqnonlin','x0',x_temp,'objective',optfun,...
% % % %     %     'lb',lb,'ub',ub,'xdata',data,'ydata',data,'setup',setup);
% % % %     % ms = MultiStart;
% % % %     % [xmulti,errormulti] = run(ms,problem,50);
% % % % plength = 60;
% % % % nsamples = 50;
% % % % fullError_Cell = cell(1,nsamples);
% % % % error_Cell = cell(1,nsamples);
% % % % time_Cell = cell(1,nsamples);
% % % % xres_Cell = cell(1,nsamples);
% % % % for i = 1:nsamples
% % % %     % call specific starting point
% % % %     rng(i);
% % % %     x_temp = -1 + (1+1) * rand(1,plength);
% % % %     % common setup
% % % %     ub = 3*ones(1,plength);
% % % %     lb = -3*ones(1,plength);
% % % %     x_temp([5 10 15 20 25 30 35 40 45 50 55 60]) = zeros;
% % % %     ub([5 10 15 20 25 30 35 40 45 50 55 60]) = zeros;
% % % %     lb([5 10 15 20 25 30 35 40 45 50 55 60]) = zeros;
% % % %     % parameter estimation
% % % %     tic
% % % %     [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% % % %     xres_Cell{i} = xres;
% % % %     t = toc;
% % % %     time_Cell{i} = t;
% % % %     % error calculation
% % % %     setup.plotEachSimCF = 0;
% % % %     setup.plotEachSim = 0;
% % % %     setup.simAllProfiles = 0;
% % % %     [error] = optfun(xres,data,setup);
% % % %     fullError_Cell{i} = error;
% % % %     error_Cell{i} = sum(abs(error));
% % % % end
% % % % % pretreatment for plot
% % % % A = cell2mat(error_Cell);
% % % % B = cell2mat(time_Cell);
% % % % % plot
% % % % figure,
% % % % subplot(121), plot(A,'k.-'),title('overall error per iteration')
% % % % subplot(122), plot(B,'k.-'),title('time per iteration')
% % % % 
% % % % % selecting and simulating the lowest error case 
% % % % [~,idx] = min(A);
% % % % setup.plotEachSimCF = 1;
% % % % setup.plotEachSim = 0;
% % % % setup.simAllProfiles = 1;
% % % % [error] = optfun(xres_Cell{idx},data,setup);
% % % % setup.plotEachSimCF = 0;
% % % % setup.plotEachSim = 0;
% % % % setup.simAllProfiles = 0;
% % % % 
% % % % % save
% % % % saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_MultiStart_Results.mat'];
% % % % save(saveName,'xres_Cell','time_Cell','fullError_Cell','error_Cell')
% % % % set(101,'color','white'), savefig(101,['results/',setup.enzymeName,'/',setup.enzymeName, '_MultiStart_minError_metabolites_reg.fig']);
% % % % set(102,'color','white'), savefig(102,['results/',setup.enzymeName,'/',setup.enzymeName, '_MultiStart_minError_fluxes_reg.fig']);
% % % % set(1,'color','white'), savefig(1,['results/',setup.enzymeName,'/',setup.enzymeName, '_MultiStart_minError.fig']);


%% (2.1) Parameter estimation with regularization
% parameter estimation
% lambdalist = 1E-5;
lambdalist = [...
    1E-5, 2E-5, 5E-5,...
    1E-4, 2E-4, 5E-4,...
    1E-3, 2E-3, 5E-3,...
    1E-2, 2E-2, 5E-2,...
    1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 7E-1,... %area of change
    1E0, 2E0, 3E0, 4E0, 5E0, 7E0,...%area of change
    1E1, 2E1, 5E1,...
    1E2, 2E2, 5E2,...
    1E3, 2E3, 5E3,...
    1E4, 2E4, 5E4,...
    1E5, 2E5, 5E5,...
    ];
parameterEstimation_lambdalist;


%% (2.2) Regularization. Results Visualization
selLambdaPos = 1; %1; %15;
regularizationVisualization;


% % % % saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_temp_regRes.mat'];
% % % % save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected','array_resnorm','array_residual','array_Jacobian','array_raw_error')
% %% (2.3) Saving data after regularization
if setup.plotOutput == 1
    saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_regularizationResults.mat'];
%     save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected')
        save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected','array_resnorm','array_residual','array_Jacobian','array_raw_error')
    set(101,'color','white'), savefig(101,['results/',setup.enzymeName,'/',setup.enzymeName, '_trainData_fit_metabolites_reg.fig']);
    set(102,'color','white'), savefig(102,['results/',setup.enzymeName,'/',setup.enzymeName, '_testData_fit_fluxes_reg.fig']);
    set(103,'color','white'), savefig(103,['results/',setup.enzymeName,'/',setup.enzymeName, '_regularization.fig']);
end


%% (2.3) Studying confidence intervals vs regularization
plot_ParsCIs_lambdalist;
if setup.plotOutput == 1
    savefig(104,['results/',setup.enzymeName,'/',setup.enzymeName, '_ParsCIs_lambdalist.fig']);
end


%% (3.1) Study on paarameter values: estimation with the lambda value
% % select lambda
% setup.selectedLambda = lambdalist(selLambdaPos);
% % estimate parameters
% tic
% [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t = toc;
% % simulate for error
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;
% [error] = optfun(xres,data,setup);
%  % confidence intervals
% N = length(error);
% Jacobian = full(Jacobian);  
% varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
% stdp = sqrt(diag(varp));

% recall specific values
setup.selectedLambda = lambdalist(selLambdaPos);
error = array_raw_error{selLambdaPos};
N = cell_N{selLambdaPos};
Jacobian = cell_Jacobian{selLambdaPos};
varp = cell_varp{selLambdaPos};
stdp = cell_stdp{selLambdaPos};


%% ENZYME-SPECIFIC
%% (3.2) Study on paarameter values: recalculation
% parameter values
% xres_selected = xres; %lambdalist based in 'ones', lam=0.1, loc=5.
% xres_selected = array_xres{15}; %lambdalist based in 'ones', lam=0.1, loc=5.
xres_selected = array_xres{selLambdaPos};

kgap = ones(numpHtested,1);
kbpg = ones(numpHtested,1);
knad = ones(numpHtested,1);
knadh = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    kgap(i) = 2.48 * 10.^xres_selected(1); %mM
    kbpg(i) = 1.18 * 10.^xres_selected(2); %mM
    knad(i) = 2.92 * 10.^xres_selected(3); %mM
    knadh(i) = 0.022 * 10.^xres_selected(4); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+4);
end
keq_gapdh = setup.pH_Keq_gapdh_eQ_fwd;
keq_pgk = setup.pH_Keq_pgk_fwd;

% limits
kgap_up = ones(numpHtested,1);
kbpg_up = ones(numpHtested,1);
knad_up = ones(numpHtested,1);
knadh_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);
kgap_down = ones(numpHtested,1);
kbpg_down = ones(numpHtested,1);
knad_down = ones(numpHtested,1);
knadh_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    %up
    kgap_up(i) = 2.48 * 10.^(xres_selected(1) + stdp(1)); %mM
    kbpg_up(i) = 1.18 * 10.^(xres_selected(2) + stdp(2)); %mM
    knad_up(i) = 2.92 * 10.^(xres_selected(3) + stdp(3)); %mM
    knadh_up(i) = 0.022 * 10.^(xres_selected(4) + stdp(4)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+4) + stdp(i+4));
    %down
    kgap_down(i) = 2.48 * 10.^(xres_selected(1) - stdp(1)); %mM
    kbpg_down(i) = 1.18 * 10.^(xres_selected(2) - stdp(2)); %mM
    knad_down(i) = 2.92 * 10.^(xres_selected(3) - stdp(3)); %mM
    knadh_down(i) = 0.022 * 10.^(xres_selected(4) - stdp(4)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+4) - stdp(i+4));
end

vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;
% % % % % Test resulting value
% % % % vm = data.Vmax(:,4) .* (10 .^ xres(4:11)') .* 60 .* 60 ./ 1.7781;
% % % % vm = mean([data.Vmax(:,1)*8 data.Vmax(:,2)*4],2) .* (10 .^ zeros(8,1)) .* 60 .* 60 ./ 1.7781;


%% (3.3) Study on paarameter values: plotting
figure(105)

subplot(331) % vm
    plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_uChange,'.-','color','black')
title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

subplot(332) % vm experimental
    errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
    title({'v_{m} experimental';'not normalized'})

subplot(333) % vm(log_difference)
    plot(pHarray,xres_selected(5:16)','.-','color','black')
    title('v_{m} difference []')

subplot(334) % kpyr
    plot(pHarray,kgap_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kgap_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kgap,'.-','color','black')
title('k_{gap} [mM]')

subplot(335) % kbpg
    plot(pHarray,kbpg_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kbpg_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kbpg,'.-','color','black')
title('k_{bpg} [mM]')

subplot(336) % knad
    plot(pHarray,knad_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,knad_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,knad,'.-','color','black')
title('k_{nad} [mM]')

subplot(337) % knadh
    plot(pHarray,knadh_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,knadh_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,knadh,'.-','color','black')
title('k_{nadh} [mM]')

subplot(338) % keq_adh
plot(pHarray,keq_gapdh,'.-','color','black')
title('k_{eq.GAPDH} [mM]')

subplot(339) % keq_adh
plot(pHarray,keq_pgk,'.-','color','black')
title('k_{eq.PGK} [mM]')

suptitle('GAPDH_{pgk}: parameter estimates vs pH')


%%
output_gapdh.xres_selected = xres_selected;

output_gapdh.pHarray = pHarray;

output_gapdh.kgap = kgap;% = ones(numpHtested,1);
output_gapdh.kbpg = kbpg;% = ones(numpHtested,1);
output_gapdh.knad = knad;% = ones(numpHtested,1);
output_gapdh.knadh = knadh;% = ones(numpHtested,1);
output_gapdh.vm = vm;% = zeros(numpHtested,1);
output_gapdh.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_gapdh.kgap_up = kgap_up;% = ones(numpHtested,1);
output_gapdh.kbpg_up = kbpg_up;% = ones(numpHtested,1);
output_gapdh.knad_up = knad_up;% = ones(numpHtested,1);
output_gapdh.knadh_up = knadh_up;% = ones(numpHtested,1);
output_gapdh.vm_up = vm_up;% = zeros(numpHtested,1);
output_gapdh.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_gapdh.kgap_down = kgap_down;% = ones(numpHtested,1);
output_gapdh.kbpg_down = kbpg_down;% = ones(numpHtested,1);
output_gapdh.knad_down = knad_down;% = ones(numpHtested,1);
output_gapdh.knadh_down = knadh_down;% = ones(numpHtested,1);
output_gapdh.vm_down = vm_down;% = zeros(numpHtested,1);
output_gapdh.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_gapdh.keq_gapdh = keq_gapdh;% = setup.Keq_ADH;
output_gapdh.keq_pgk = keq_pgk;% = setup.Keq_ADH;


%% saving output
output_gapdh.Vmax_experimental = Vmax_experimental;
output_gapdh.stDev_experimental = stDev_experimental;
saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_parEst.mat'];
save(saveName,'output_gapdh');
set(105,'color','white'), savefig(105,['results/',setup.enzymeName,'/',setup.enzymeName, '_parEstimates.fig']);


% % % % x = 1; end

