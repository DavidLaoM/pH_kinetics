% % % % function x=dev_pEst_ald
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
    setup.caseStudyENO = 0;
    selectSetup_pH;
    % added
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_ald = expData.ald;
    
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
    
    pHarray = unique(import_ald.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_ald.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_ald.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_ald.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_ald.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_ald.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_ald.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_ald.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_ald.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_ald.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_ald.rawData.absorbance_corrected{4,4};
    temp2 = import_ald.rawData.absorbance_corrected{5,4};
    temp3 = import_ald.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_ald.rawData.time{1};
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Directly changing the concentration here, sicne the extinction
    % coefficient did not change.
    dps = length(NADH{1,1});
    endPoint = zeros(numpHtested,DFs);
    for i = 1:DFs
        for j = 1:numpHtested
%             endPoint(j,i) = min(NADH{j,i});
            endPoint(j,i) = min([NADH{j,1}' NADH{j,2}' NADH{j,3}' NADH{j,4}']);
        end
    end
    for i = 1:DFs
        for j = 1:numpHtested
            for k = 1:dps
                NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,i);
            end
        end
    end
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % % % % %     Addition ALDolase to delete the first part of the profile
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % NADH concentration
    NADH2 = cell(size(NADH));
    for i = 1:DFs
        for j = 1:numpHtested
            if j == 1
                NADH2{j,i} = NADH{j,i}(25:end);
            else
                NADH2{j,i} = NADH{j,i}(8:end);
            end
        end
    end
    % bring time to zero start
    for i = 1:DFs
        for j = 1:numpHtested
            for k = 1:dps
                if j == 1
                    time{j,i}(k) = time{j,i}(k) - 120;
                else
                    time{j,i}(k) = time{j,i}(k) - 35;
                end
            end
        end
    end
    % GPD reaction rate
    RRs2 = cell(size(RRs));
    for i = 1:DFs
        for j = 1:numpHtested
            if j == 1
                RRs2{j,i} = RRs{j,i}(25:end);
            else
                RRs2{j,i} = RRs{j,i}(8:end);
            end
        end
    end
    data.RRs = RRs2;
    % time
    time2 = cell(size(time));
    for i = 1:DFs
        for j = 1:numpHtested
            if j == 1
                time2{j,i} = time{j,i}(25:end);
            else
                time2{j,i} = time{j,i}(8:end);
            end
        end
    end
    clear NADH, NADH = NADH2; clear NADH2
    data.conc_mean = NADH;    
    clear time, time = time2; clear time2
    data.time = time;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
    pHvals = unique(import_ald.treatedData.pH_corrected);
    % visualize: check calculations made
    hs = figure('units','normalized','outerposition',[0 0 1 1]);
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
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
    suptitle(suptitleName);
%     if setup.plotOutput == 1
%         saveFigName = ['results\', setup.enzymeName,'\', setup.enzymeName, '_concentrations_basezero.fig'];
%         savefig(hs,saveFigName);
%     end
    
%     figure
%     plot(pHvals, Vmax(:,4),'.-')
%     title('Starting estimate: Vmax [mM s-1] vs pH')    
end


%% (0.1) Experimental Vmax determination
setup.plotOutput = 0;
% %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 6; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
% select start point (this needs manual selection deppending on previous plots)
dp_start = 2 * ones(size(data.conc_mean));
dp_start(:,1) = 10 * ones(size(dp_start(:,1)));
dp_start(:,2) = 5 * ones(size(dp_start(:,2)));
dp_start(:,3) = 3 * ones(size(dp_start(:,3)));
dp_start(2,4) = 3 * ones(size(dp_start(2,4)));
% blank cell total length
total_len = zeros(size(dp_start));
% DFs considered
DFarray = [1/8 1/4 1/2 1/1];
% idxs2consider
idxs2consider = [1 1 1 0;
                1 1 0 0;
                1 1 0 0;
                1 1 1 0;
                1 1 1 0;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1];

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
setup.plotEachSim = 0;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
setup.costfun = 3;

setup.weightData = 1;
setup.weightDataEsp = idxs2consider;
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% Km fixed
optfun = @costfun_Kmfixed;
% % % % plength = 13; % Kms (3) + Vms (1) * numpH (8) + vm linking (2) + (Keq is fixed to experimental data)
plength = 11; % Kms (3) + Vms (1) * numpH (8) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');
% %%
% % testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(x_temp,data,setup); % sum(abs(error)) = 9.5457
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;
% %%
% % parameter estimation
% tic
% [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t = toc;
% %%
% % testing the result
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(xres,data,setup); % sum(abs(error)) = 2.3872
% % [error] = optfun(array_xres{13},data,setup); % sum(abs(error)) = 2.3872
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;


%% (2.1) Parameter estimation with regularization
% parameter estimation
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
% % % % lambdalist = [...
% % % %     1E-5,...% 2E-5, 5E-5,...
% % % %     1E-4,...% 2E-4, 5E-4,...
% % % %     1E-3,...% 2E-3, 5E-3,...
% % % %     1E-2,...% 2E-2, 5E-2,...
% % % %     1E-1, 2E-1,...% 3E-1, 4E-1, 5E-1, 7E-1,... %area of change
% % % %     1E0, 2E0,...% 3E0, 4E0, 5E0, 7E0,...%area of change
% % % %     1E1,...% 2E1, 5E1,...
% % % %     1E2,...% 2E2, 5E2,...
% % % %     1E3,...% 2E3, 5E3,...
% % % %     1E4,...% 2E4, 5E4,...
% % % %     1E5,...% 2E5, 5E5,...
% % % %     ];
% % % % lambdalist = [1E-4 2E-1 1E5];
parameterEstimation_lambdalist;


%% (2.2) Regularization. Results Visualization
selLambdaPos = 1;%19;%1;%17;
regularizationVisualization;
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
% xres_selected = array_xres{17}; %lambdalist based in 'ones', lam=0.1, loc=5.
xres_selected = array_xres{selLambdaPos};

kfbp = ones(numpHtested,1);
kgap = ones(numpHtested,1);
kdhap = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    kfbp(i) = 0.451 * 10.^xres_selected(1); %mM
    kgap(i) = 2 * 10.^xres_selected(2); %mM
    kdhap(i) = 2.4 * 10.^xres_selected(3); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+3);
end
keq_fba = setup.Keq_FBA;% = [1.0E-3 7.9E-4 7.3E-4 7E-4 6.8E-4 6.7E-4 6.6E-4 6.5E-4];  %dir+
keq_tpi = setup.Keq_TPI;% = [1/(8.31) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38) 1/(9.39)];  %dir-
keq_gpd = setup.Keq_GPD;% = [1/(4.2E-6) 1/(1.5E-5) 1/(2.7E-5) 1/(4.7E-5) 1/(7.9E-5) 1/(1.2E-4) 1/(1.6E-4) 1/(2.0E-4) ]; %dir-

% limits
kfbp_up = ones(numpHtested,1);
kgap_up = ones(numpHtested,1);
kdhap_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);
kfbp_down = ones(numpHtested,1);
kgap_down = ones(numpHtested,1);
kdhap_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    %up
    kfbp_up(i) = 0.451 * 10.^(xres_selected(1)+stdp(1)); %mM
    kgap_up(i) = 2 * 10.^(xres_selected(2)+stdp(2)); %mM
    kdhap_up(i) = 2.4 * 10.^(xres_selected(3)+stdp(3)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+3)+stdp(i+3));
    %down
    kfbp_down(i) = 0.451 * 10.^(xres_selected(1)-stdp(1)); %mM
    kgap_down(i) = 2 * 10.^(xres_selected(2)-stdp(2)); %mM
    kdhap_down(i) = 2.4 * 10.^(xres_selected(3)-stdp(3)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+3)-stdp(i+3));
end


%% (3.3) Study on paarameter values: plotting
vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;

figure(105)
% figure(106)

subplot(331) % vm estimated
%     plot(pHarray,vm_up,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm_down,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm,'.-','color','black')
% title('v_{m} [mM]')
    plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_uChange,'.-','color','black')
title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

subplot(332) % vm experimental
    errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
    title({'v_{m} experimental';'not normalized'})

subplot(333) % vm(log_difference)
    plot(pHarray,xres_selected(4:11)','.-','color','black')
    title('v_{m} difference []')

subplot(334) % kfbp
    plot(pHarray,kfbp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kfbp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kfbp,'.-','color','black')
title('k_{fbp} [mM]')

subplot(335) % kgap
    plot(pHarray,kgap_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kgap_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kgap,'.-','color','black')
title('k_{gap} [mM]')

subplot(336) % kdhap
    plot(pHarray,kdhap_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kdhap_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kdhap,'.-','color','black')
title('k_{dhap} [mM]')

subplot(337) % keq_fba
plot(pHarray,keq_fba,'.-','color','black')
title('k_{eq.FBA} [mM]')

subplot(338) % keq_tpi
plot(pHarray,keq_tpi,'.-','color','black')
title('k_{eq.TPI} [mM]')

subplot(339) % keq_gpd
plot(pHarray,keq_gpd,'.-','color','black')
title('k_{eq.GPD} [mM]')

suptitle('ALD: parameter estimates vs pH')


%%
output_ald.xres_selected = xres_selected;

output_ald.pHarray = pHarray;
output_ald.idxs2consider = idxs2consider;
output_ald.DF = data.DF;

output_ald.kfbp = kfbp;% = ones(numpHtested,1);
output_ald.kgap = kgap;% = ones(numpHtested,1);
output_ald.kdhap = kdhap;% = ones(numpHtested,1);
output_ald.vm = vm;% = zeros(numpHtested,1);
output_ald.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_ald.kfbp_up = kfbp_up;% = ones(numpHtested,1);
output_ald.kgap_up = kgap_up;% = ones(numpHtested,1);
output_ald.kdhap_up = kdhap_up;% = ones(numpHtested,1);
output_ald.vm_up = vm_up;% = zeros(numpHtested,1);
output_ald.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_ald.kfbp_down = kfbp_down;% = ones(numpHtested,1);
output_ald.kgap_down = kgap_down;% = ones(numpHtested,1);
output_ald.kdhap_down = kdhap_down;% = ones(numpHtested,1);
output_ald.vm_down = vm_down;% = zeros(numpHtested,1);
output_ald.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_ald.keq_fba = keq_fba;% = setup.Keq_FBA;
output_ald.keq_tpi = keq_tpi;% = setup.Keq_TPI;
output_ald.keq_gpd = keq_gpd;% = setup.Keq_GPD;


%% saving output
output_ald.Vmax_experimental = Vmax_experimental;
output_ald.stDev_experimental = stDev_experimental;
saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_parEst.mat'];
save(saveName,'output_ald');
set(105,'color','white'), savefig(105,['results/',setup.enzymeName,'/',setup.enzymeName, '_parEstimates.fig']);

% save('output_ald.mat','output_ald');
% % % % %% reverse to reload
% % % % load('output_ald.mat','output_ald');
% % % % 
% % % % xres_selected = output_ald.xres_selected;
% % % % 
% % % % pHarray = output_ald.pHarray;
% % % % 
% % % % kfbp = output_ald.kfbp;%
% % % % kgap = output_ald.kgap;%
% % % % kdhap = output_ald.kdhap;%
% % % % vm = output_ald.vm;%
% % % % 
% % % % kfbp_up = output_ald.kfbp_up;%
% % % % kgap_up = output_ald.kgap_up;%
% % % % kdhap_up = output_ald.kdhap_up;%
% % % % vm_up = output_ald.vm_up;%
% % % % 
% % % % kfbp_down = output_ald.kfbp_down;%
% % % % kgap_down = output_ald.kgap_down;%
% % % % kdhap_down = output_ald.kdhap_down;%
% % % % vm_down = output_ald.vm_down;%
% % % % 
% % % % keq_fba = output_ald.keq_fba;%
% % % % keq_tpi = output_ald.keq_tpi;%
% % % % keq_gpd = output_ald.keq_gpd;%


%% final comparison case with and without regularization
% simulation display off
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
% no lambda selected
setup.selectedLambda = 0; % by now just testing
tic
[xres_nonreg,resnorm_nonreg,residual_nonreg,~,~,~,Jacobian_nonreg] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t_nonreg = toc;
    [error_nonreg] = optfun(xres_nonreg,data,setup);
    N_nonreg = length(error_nonreg);
    Jacobian_nonreg = full(Jacobian_nonreg);  
    varp_nonreg = resnorm_nonreg*inv(Jacobian_nonreg'*Jacobian_nonreg)/N_nonreg; % covariance matrix
    stdp_nonreg = sqrt(diag(varp_nonreg));

% selected lambda value
setup.selectedLambda = lambdalist(selLambdaPos); % by now just testing
tic
[xres_reg,resnorm_reg,residual_reg,~,~,~,Jacobian_reg] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t_reg = toc;
    [error_reg] = optfun(xres_reg,data,setup);
    N_reg = length(error_reg);
    Jacobian_reg = full(Jacobian_reg);  
    varp_reg = resnorm_reg*inv(Jacobian_reg'*Jacobian_reg)/N_reg; % covariance matrix
    stdp_reg = sqrt(diag(varp_reg));
    
% saving
saveName2 = ['results/',setup.enzymeName,'/',setup.enzymeName, '_regularizationONOFF.mat'];
save(saveName2,'xres_nonreg','resnorm_nonreg','residual_nonreg','Jacobian_nonreg','t_nonreg',...
    'xres_reg','resnorm_reg','residual_reg','Jacobian_reg','t_reg',...
    'error_nonreg','N_nonreg','Jacobian_nonreg','varp_nonreg','stdp_nonreg',...
    'error_reg','N_reg','Jacobian_reg','varp_reg','stdp_reg');

pNames = {'Kf16bp';...
    'Kglyceral3p';...
    'Kdhap';...
    'Vm_pH6_32';...
    'Vm_pH6_81';...
    'Vm_pH7_06';...
    'Vm_pH7_29';...
    'Vm_pH7_51';...
    'Vm_pH7_68';...
    'Vm_pH7_81';...
    'Vm_pH7_90';...
    'blank';...
    'blank'};
T = table(pNames,stdp_nonreg,stdp_reg);


