% % % % function x=dev_pEst_hxk
% % PEST_ENO.m
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
    setup.caseStudyHXK = 1;
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
    import_eno = expData.hxk;

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

    pHarray = unique(import_eno.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_eno.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_eno.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_eno.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_eno.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_eno.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_eno.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_eno.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_eno.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_eno.treatedData.reaction_rate{tempID(j)};
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

    NADPH = blankCell;
    Vmax = blank;
    for i = 1:(DFs*numpHtested)
        tempDiff = conc_mean{i} - conc_mean{i}(1); % all stoichiometries are 1-to-1.
        NADPH{i} = conc_mean{i};
        %     RRs2 = RRs';
%         % Option 1. Vmax from the values obtained
%         Vmax(i) = max(abs(RRs{i}));
        % Option 2. Vmax naive approach. Full profile
        Vmax(i) = (conc_mean{i}(end) - conc_mean{i}(1)) ./ (time{i}(end) - time{i}(1)); 
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
    data.chosenNADPHini = 0.05;
    temp1 = import_eno.rawData.absorbance_corrected{4,4};
    temp2 = import_eno.rawData.absorbance_corrected{5,4};
    temp3 = import_eno.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_eno.rawData.time{1};
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Directly changing the concentration here, sicne the extinction
    % coefficient did not change.
    dps = length(NADPH{1,1});
    iniPointDF4 = zeros(numpHtested,DFs);
    for i = 1:DFs
        for j = 1:numpHtested
            iniPointDF4(j,i) = NADPH{j,1}(1);
        end
    end
    for i = 1:DFs
        for j = 1:numpHtested
            for k = 1:dps
                NADPH{j,i}(k) = NADPH{j,i}(k) - iniPointDF4(j,1);
            end
        end
    end
    data.conc_mean = NADPH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    pHvals = unique(import_eno.treatedData.pH_corrected);
    % visualize: check calculations made
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numpHtested
        subplot(3,4,i)
        for j = 1:DFs
            plot(time{i,j},NADPH{i,j},'.-')
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
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': NADPH concentration profile'];
    suptitle(suptitleName);
end


%% (0.1) Experimental Vmax determination
setup.plotOutput = 0;
% %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 6; % minimum size of the window
    limRates = [0 5E-4]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.1];  %Ylims plot conc
    idxs_one_initial = ones(size(data.conc_mean));
    idxs_one_initial(1:2,:) = (20 - minwindow) * idxs_one_initial(1:2,:);
    idxs_one_initial(3:end,:) = (61 - minwindow) * idxs_one_initial(3:end,:);
% select start point (this needs manual selection deppending on previous plots)
% % % % dp_start = 2 * ones(size(data.conc_mean));
% % % % dp_start(:,1) = 10 * ones(size(dp_start(:,1)));
% % % % dp_start(:,2) = 5 * ones(size(dp_start(:,2)));
% % % % dp_start(:,3) = 3 * ones(size(dp_start(:,3)));
% % % % dp_start(2,4) = 3 * ones(size(dp_start(2,4)));
dp_start = ones(size(data.conc_mean));
dp_start(:,1) = 6 * ones(size(dp_start(:,1)));
% blank cell total length
total_len = zeros(size(dp_start));
% DFs considered
DFarray = [1/8 1/4 1/2 1/1];
% idxs2consider
idxs2consider = ones(size(DF));

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
plength = 16; % Kms (4) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');
% %%
% testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(x_temp,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;
% %%
% parameter estimation
% tic
% [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t = toc;
% % xres = [-1.4378    1.6545    1.3091   -0.0101   -2.9995   -3.0000,...
% %     -2.0691   -1.9189   -1.6789   -1.9170   -1.8883   -1.9901,...
% %     -2.0522   -2.0511   -2.0648   -2.0184];


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
parameterEstimation_lambdalist;


%% (2.2) Regularization. Results Visualization
selLambdaPos = 1;%11;%1;%14;
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
% xres_selected = array_xres{5}; %lambdalist based in 'ones', lam=0.1, loc=5.
xres_selected = array_xres{selLambdaPos};

kadp = ones(numpHtested,1);
katp = ones(numpHtested,1);
kg6p = ones(numpHtested,1);
kglc = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    kadp(i) = 0.23 * 10.^xres_selected(1); %mM
    katp(i) = 0.15 * 10.^xres_selected(2); %mM
    kg6p(i) = 30 * 10.^xres_selected(3); %mM
    kglc(i) = 0.08 * 10.^xres_selected(4); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+4);
end
keq_hxk = setup.Keq_HXK;
keq_g6pdh = setup.Keq_G6PDH;

% limits
kadp_up = ones(numpHtested,1);
katp_up = ones(numpHtested,1);
kg6p_up = ones(numpHtested,1);
kglc_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);
kadp_down = ones(numpHtested,1);
katp_down = ones(numpHtested,1);
kg6p_down = ones(numpHtested,1);
kglc_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    %up
    kadp_up(i) = 0.23 * 10.^(xres_selected(1)+stdp(1)); %mM
    katp_up(i) = 0.15 * 10.^(xres_selected(2)+stdp(2)); %mM
    kg6p_up(i) = 30 * 10.^(xres_selected(3)+stdp(3)); %mM
    kglc_up(i) = 0.08 * 10.^(xres_selected(4)+stdp(4)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+4)+stdp(i+4));
    %down
    kadp_down(i) = 0.23 * 10.^(xres_selected(1)-stdp(1)); %mM
    katp_down(i) = 0.15 * 10.^(xres_selected(2)-stdp(2)); %mM
    kg6p_down(i) = 30 * 10.^(xres_selected(3)-stdp(3)); %mM
    kglc_down(i) = 0.08 * 10.^(xres_selected(4)-stdp(4)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+4)-stdp(i+4));
end


%% (3.3) Study on paarameter values: plotting
vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;

figure(105)
% figure(106)

subplot(331) % vm
    plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_uChange,'.-','color','black')
    title({'v_{m} [umol_{NADPH} mg_{P}^{-1} min^{-1}]';'not normalized'})  

subplot(332) % vm experimental
    errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
    title({'v_{m} experimental';'not normalized'})

subplot(333) % keq_g6pdh
yyaxis left
    plot(pHarray,keq_hxk,'.-','color','black')
yyaxis right
    plot(pHarray,keq_g6pdh,'.-','color','black')
    legend('k_{eq.HXK} [mM]','k_{eq.G6PDH} [mM]')

subplot(335) % kadp
    plot(pHarray,kadp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kadp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kadp,'.-','color','black')
title('k_{adp} [mM]')

subplot(336) % katp
    plot(pHarray,katp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,katp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,katp,'.-','color','black')
title('k_{atp} [mM]')

subplot(337) % kg6p
    plot(pHarray,kg6p_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kg6p_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kg6p,'.-','color','black')
title('k_{g6p} [mM]')

subplot(338) % kglc
    plot(pHarray,kglc_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kglc_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kglc,'.-','color','black')
title('k_{glc} [mM]')

subplot(334)
plot(pHarray,xres_selected(5:end)','.-','color','black')
title('vm reference')

suptitle('HXK: parameter estimates vs pH')


%%
output_hxk.xres_selected = xres_selected;

output_hxk.pHarray = pHarray;
output_hxk.idxs2consider = idxs2consider;
output_hxk.DF = data.DF;

output_hxk.kadp = kadp;% = ones(numpHtested,1);
output_hxk.katp = katp;% = ones(numpHtested,1);
output_hxk.kg6p = kg6p;% = ones(numpHtested,1);
output_hxk.kglc = kglc;% = ones(numpHtested,1);
output_hxk.vm = vm;% = zeros(numpHtested,1);
output_hxk.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_hxk.kadp_up = kadp_up;% = ones(numpHtested,1);
output_hxk.katp_up = katp_up;% = ones(numpHtested,1);
output_hxk.kg6p_up = kg6p_up;% = ones(numpHtested,1);
output_hxk.kglc_up = kglc_up;% = ones(numpHtested,1);
output_hxk.vm_up = vm_up;% = zeros(numpHtested,1);
output_hxk.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_hxk.kadp_down = kadp_down;% = ones(numpHtested,1);
output_hxk.katp_down = katp_down;% = ones(numpHtested,1);
output_hxk.kg6p_down = kg6p_down;% = ones(numpHtested,1);
output_hxk.kglc_down = kglc_down;% = ones(numpHtested,1);
output_hxk.vm_down = vm_down;% = zeros(numpHtested,1);
output_hxk.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_hxk.keq_hxk = keq_hxk;% = setup.Keq_PGI;
output_hxk.keq_g6pdh = keq_g6pdh;% = setup.Keq_PGI;


%% saving output
output_hxk.Vmax_experimental = Vmax_experimental;
output_hxk.stDev_experimental = stDev_experimental;
saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_parEst.mat'];
save(saveName,'output_hxk');
set(105,'color','white'), savefig(105,['results/',setup.enzymeName,'/',setup.enzymeName, '_parEstimates.fig']);
% save('output_hxk.mat','output_hxk');

% % % % %% reverse to reload
% % % % load('output_hxk.mat','output_hxk');
% % % % 
% % % % xres_selected = output_hxk.xres_selected;
% % % % 
% % % % pHarray = output_hxk.pHarray;
% % % % 
% % % % kadp = output_hxk.kadp;%
% % % % katp = output_hxk.katp;%
% % % % kg6p = output_hxk.kg6p;%
% % % % kglc = output_hxk.kglc;%
% % % % vm = output_hxk.vm;%
% % % % 
% % % % kadp_up = output_hxk.kadp_up;%
% % % % katp_up = output_hxk.katp_up;%
% % % % kg6p_up = output_hxk.kg6p_up;%
% % % % kglc_up = output_hxk.kglc_up;%
% % % % vm_up = output_hxk.vm_up;%
% % % % 
% % % % kadp_down = output_hxk.kadp_down;%
% % % % katp_down = output_hxk.katp_down;%
% % % % kg6p_down = output_hxk.kg6p_down;%
% % % % kglc_down = output_hxk.kglc_down;%
% % % % vm_down = output_hxk.vm_down;%
% % % % 
% % % % keq_hxk = output_hxk.keq_hxk;%
% % % % keq_g6pdh = output_hxk.keq_g6pdh;%
% % % % x = 1; end

