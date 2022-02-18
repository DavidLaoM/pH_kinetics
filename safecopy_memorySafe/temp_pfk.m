% % % % function x=dev_pEst_pfk
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
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 0;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 1;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    setup.caseStudyENO = 0;
    selectSetup_pH;
    % added
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_pfk = expData.pfk;
    
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
    
    pHarray = unique(import_pfk.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_pfk.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_pfk.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_pfk.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_pfk.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_pfk.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_pfk.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_pfk.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_pfk.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_pfk.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_pfk.rawData.absorbance_corrected{4,4};
    temp2 = import_pfk.rawData.absorbance_corrected{5,4};
    temp3 = import_pfk.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_pfk.rawData.time{1};
    
        % (1) Correct for minimum value
        % (2) Bring the minimum to zero (apply to all)
        % (3) In principle, use the 3 first dilution rates
        % (4) Watch out with the dilution factors (first 2 cases are
        % reversed)
        
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Adjusted to PGM
    % Directly changing the concentration here, since the extinction
    % coefficient did not change.
    dpsarray = zeros(numpHtested,DFs);
    for i = 1:DFs
        for j = 1:numpHtested
            dpsarray(j,i) = length(NADH{j,i});
        end
    end
    endPoint = zeros(numpHtested,DFs);
    endLocation = zeros(numpHtested,DFs);
    % locate the minimum
    for i = 1:DFs
        for j = 1:numpHtested
            [endPoint(j,i),endLocation(j,i)] = min(NADH{j,i});
        end
    end
    % bringing the minimum to zero
    for i = 1:DFs
        for j = 1:numpHtested
            dps = dpsarray(j,i);
            if ((i==3)&&(j>=3)&&(j<=9))
                for k = 1:dps
                    NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,3);
                end
            else
                for k = 1:dps
                    NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,4);
                end
            end
        end
    end
    % bring the late increase to zero
    for i = 1:DFs
        for j = 1:numpHtested
            % locate the minimum
            [tempval, tempidx] = min(NADH{j,i});
            % from the index to end(dps) make it the value of the minimum (zero)
            dps = dpsarray(j,i);
            for k = (tempidx+1):dps
                NADH{j,i}(k) = NADH{j,i}(tempidx);
            end
        end
    end
    
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    takenTime = [700, 700, 550, 500, 450, 450, 350, 350, 350, 350, 350, 350];
    cuttingPoints = [141, 141, 111, 101, 91, 91, 71, 71, 71, 71, 71, 71];
    for i = 1:DFs
        for j = 1:numpHtested
            startVal = cuttingPoints(j);
            data.abs_mean{j,i} = data.abs_mean{j,i}(startVal:end);
            data.abs_std{j,i} = data.abs_std{j,i}(startVal:end);
            data.conc_mean{j,i} = data.conc_mean{j,i}(startVal:end);
            data.conc_std{j,i} = data.conc_std{j,i}(startVal:end);
            data.time{j,i} = data.time{j,i}(startVal:end) - takenTime(j);
            data.RRs{j,i} = data.RRs{j,i}(startVal:end);
            time{j,i} = time{j,i}(startVal:end) - takenTime(j);
            NADH{j,i} = NADH{j,i}(startVal:end);
        end
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
    
    pHvals = unique(import_pfk.treatedData.pH_corrected);
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
        if setup.caseStudyPFK == 1
            ylim([0 0.12])
            xlim([0 2000])
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
    minwindow = 70; % minimum size of the window
    limRates = [0 1.5E-4]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
% select start point (this needs manual selection deppending on previous plots)
dp_start = 6 * ones(size(data.conc_mean));
% blank cell total length
total_len = zeros(size(dp_start));
% DFs considered
DFarray_pre = 1./data.DF;
DFarray = DFarray_pre(1,:);
% idxs2consider
idxs2consider = [0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1];

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
setup.ode = 'vanHeerden2014';
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
setup.weightDataEsp = idxs2consider;
% % % % idxs2consider2 = [0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1;
% % % %                 0 0 0 1];
% % % % setup.weightDataEsp = idxs2consider2;
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% three stypes of topologies and study
% % % % setup.problemStudy = 'onlyVmax';
setup.problemStudy = 'fullKinetics_paramsPartFixed';
% setup.problemStudy = 'fullKinetics_paramsAllFlexible';
% prepare setup based of 'problemStudy'
problemStudy = setup.problemStudy;
optfun = @costfun_Kmfixed;
switch problemStudy
    case 'onlyVmax'
        plength = 12; % only the vmax change with pH. Other params OFF
    case 'fullKinetics_paramsPartFixed'
        plength = 25; % only the vmax change with pH. Other params ON
    case 'fullKinetics_paramsAllFlexible'
        plength = 168; % All parameters change with pH. Other params ON
    otherwise
        disp('Warning: problem study (reaction kinetics) have not been selected');
end
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
%     ub(1:13) = zeros;
%     lb(1:13) = zeros; 
% options = optimset('Display','iter');
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4);
% % plength = 25; % vms (12) + other params (12)
% % x_temp = zeros(1,plength);
% plength = 12; % vms (12) + other params off
% % x_temp = -0.45*ones(1,plength);


% %%
% % testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% % [error] = optfun(x_temp,data,setup);
% [error] = optfun(xres,data,setup);
% % xres2 = xres; xres2(14:end) = zeros;
% % [error] = optfun(xres2,data,setup);
% % x_temp3 = x_temp; x_temp3(12:end) = 0.1;
% % [error] = optfun(x_temp3,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;
% %%
% % parameter estimation
% tic
% [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t = toc;
% 
% % testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(xres,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;

%% (2.1) Parameter estimation with regularization
% parameter estimation
% lambdalist = [...
%     1E-5, 2E-5, 5E-5,...
%     1E-4, 2E-4, 5E-4,...
%     1E-3, 2E-3, 5E-3,...
%     1E-2, 2E-2, 5E-2,...
%     1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 7E-1,... %area of change
%     1E0, 2E0, 3E0, 4E0, 5E0, 7E0,...%area of change
%     1E1, 2E1, 5E1,...
%     1E2, 2E2, 5E2,...
%     1E3, 2E3, 5E3,...
%     1E4, 2E4, 5E4,...
%     1E5, 2E5, 5E5,...
%     ];
lambdalist = 0;
parameterEstimation_lambdalist;


%% (2.2) Regularization. Results Visualization
selLambdaPos = 1;%1;%25; %8;
regularizationVisualization;

% % % % saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_temp_regRes.mat'];
% % % % save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected','array_resnorm','array_residual','array_Jacobian','array_raw_error')
% %% (2.3) Saving data after regularization
% if setup.plotOutput == 1
%     saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_regularizationResults.mat'];
% %     save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected')
%         save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected','array_resnorm','array_residual','array_Jacobian','array_raw_error')
%     set(101,'color','white'), savefig(101,['results/',setup.enzymeName,'/',setup.enzymeName, '_trainData_fit_metabolites_reg.fig']);
%     set(102,'color','white'), savefig(102,['results/',setup.enzymeName,'/',setup.enzymeName, '_testData_fit_fluxes_reg.fig']);
%     set(103,'color','white'), savefig(103,['results/',setup.enzymeName,'/',setup.enzymeName, '_regularization.fig']);
% end
% save('20201120_temp_xresPFKfull.mat','xres'); % Test for all parameters

%% (2.3) Studying confidence intervals vs regularization
% plot_ParsCIs_lambdalist;
% if setup.plotOutput == 1
%     savefig(104,['results/',setup.enzymeName,'/',setup.enzymeName, '_ParsCIs_lambdalist.fig']);
% end


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


vm = zeros(numpHtested,1);
for i = 1:numpHtested
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+13);
end
keq_fba = setup.Keq_FBA;
keq_gpd = setup.Keq_GPD;
keq_tpi = setup.Keq_TPI;
keq_pfk = setup.Keq_PFK;

% limits
vm_up = zeros(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    %up
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+13)+stdp(i+13));
    %down
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+13)-stdp(i+13));
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
    plot(pHarray,xres_selected(14:end)','.-','color','black')
    title('v_{m} difference []')

subplot(336) % keq_pfk
plot(pHarray,keq_pfk,'.-','color','black')
title('k_{eq.PFK} [mM]')

subplot(337) % keq_tpi
plot(pHarray,keq_tpi,'.-','color','black')
title('k_{eq.TPI} [mM]')

subplot(338) % keq_gpd
plot(pHarray,keq_gpd,'.-','color','black')
title('k_{eq.GPD} [mM]')

subplot(339) % keq_fba
plot(pHarray,keq_fba,'.-','color','black')
title('k_{eq.FBA} [mM]')

suptitle('PFK: parameter estimates vs pH')


%%
output_pfk.xres_selected = xres_selected;

output_pfk.pHarray = pHarray;
output_pfk.idxs2consider = idxs2consider;
output_pfk.DF = data.DF;

output_pfk.vm = vm;% = zeros(numpHtested,1);
output_pfk.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_pfk.vm_up = vm_up;% = zeros(numpHtested,1);
output_pfk.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_pfk.vm_down = vm_down;% = zeros(numpHtested,1);
output_pfk.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_pfk.keq_pfk = keq_pfk;% = setup.Keq_PFK;
output_pfk.keq_tpi = keq_tpi;% = setup.Keq_TPI;
output_pfk.keq_gpd = keq_gpd;% = setup.Keq_GPD;
output_pfk.keq_fba = keq_fba;% = setup.Keq_FBA;


%% saving output
output_pfk.Vmax_experimental = Vmax_experimental;
output_pfk.stDev_experimental = stDev_experimental;
saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_parEst.mat'];
save(saveName,'output_pfk');
set(105,'color','white'), savefig(105,['results/',setup.enzymeName,'/',setup.enzymeName, '_parEstimates.fig']);

% %% reverse to reload
% load('output_pfk.mat','output_pfk');
% 
% xres_selected = output_pfk.xres_selected;
% 
% pHarray = output_pfk.pHarray;
% 
% kfbp = output_pfk.kfbp;%
% k3pg = output_pfk.k3pg;%
% kdhap = output_pfk.kdhap;%
% vm = output_pfk.vm;%
% 
% kfbp_up = output_pfk.kfbp_up;%
% k3pg_up = output_pfk.k3pg_up;%
% kdhap_up = output_pfk.kdhap_up;%
% vm_up = output_pfk.vm_up;%
% 
% kfbp_down = output_pfk.kfbp_down;%
% k3pg_down = output_pfk.k3pg_down;%
% kdhap_down = output_pfk.kdhap_down;%
% vm_down = output_pfk.vm_down;%
% 
% keq_fba = output_pfk.keq_fba;%
% keq_tpi = output_pfk.keq_tpi;%
% keq_gpd = output_pfk.keq_gpd;%
% % % % x = 1; end


