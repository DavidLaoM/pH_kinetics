% % % % function x=dev_pEst_pgi
% % PEST_PGI.m
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
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 1;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    setup.caseStudyENO = 0;
    selectSetup_pH;
    % added
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_pgi = expData.pgi;
    
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
    
    pHarray = unique(import_pgi.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_pgi.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_pgi.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_pgi.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_pgi.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_pgi.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_pgi.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_pgi.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_pgi.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_pgi.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_pgi.rawData.absorbance_corrected{4,4};
    temp2 = import_pgi.rawData.absorbance_corrected{5,4};
    temp3 = import_pgi.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_pgi.rawData.time{1};
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Directly changing the concentration here, since the extinction
    % coefficient did not change.
    endPoint = zeros(numpHtested,DFs);
    for i = 1:DFs
        for j = 1:numpHtested
%             endPoint(j,i) = min(NADH{j,i});
            endPoint(j,i) = min([NADH{j,1}' NADH{j,2}' NADH{j,3}' NADH{j,4}']);
        end
    end
    for i = 1:DFs
        for j = 1:numpHtested
            dps = length(NADH{j,i});
            for k = 1:dps
                NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,i);
            end
        end
    end
    % making the last values increasing, zero
    for i = 1:DFs
        for j = 1:numpHtested
            [~,I] = min(NADH{j,i});
            % get the idx of the minimum, and tehen from the next idx to
            % end, make all zero
            for k = I:length(NADH{j,i})
                NADH{j,i}(k) = NADH{j,i}(I);
            end
        end
    end
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % % % % %     Addition PGI to delete the first part of the profile
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % NADH concentration 
        % pHs 1-3 | time 250 | dps 51
        % pHs 4 | time 200 | dps 41
        % pHs 5 | time 175 | dps 36
        % pHs 6 | time 125 | dps 26
        % pHs 7 | time 100 | dps 21
        % pHs 8-12 | time 75 | dps 16
    NADH2 = cell(size(NADH));
    for i = 1:DFs
        for j = 1:numpHtested
            if ((j == 1)||(j == 2)||(j == 3))
                NADH2{j,i} = NADH{j,i}(51:end);
            elseif j == 4
                NADH2{j,i} = NADH{j,i}(41:end);
            elseif j == 5
                NADH2{j,i} = NADH{j,i}(36:end);
            elseif j == 6
                NADH2{j,i} = NADH{j,i}(26:end);
            elseif j == 7
                NADH2{j,i} = NADH{j,i}(21:end);
            else
                NADH2{j,i} = NADH{j,i}(16:end);
            end
        end
    end
    % bring time to zero start
    for i = 1:DFs
        for j = 1:numpHtested
            dps = length(time{j,i});
            for k = 1:dps
                if ((j == 1)||(j == 2)||(j == 3))
                    time{j,i}(k) = time{j,i}(k) - 250;
                elseif j == 4
                    time{j,i}(k) = time{j,i}(k) - 200;
                elseif j == 5
                    time{j,i}(k) = time{j,i}(k) - 175;
                elseif j == 6
                    time{j,i}(k) = time{j,i}(k) - 125;
                elseif j == 7
                    time{j,i}(k) = time{j,i}(k) - 100;
                else
                    time{j,i}(k) = time{j,i}(k) - 75;
                end
            end
        end
    end
    % GPD reaction rate
    RRs2 = cell(size(RRs));
    for i = 1:DFs
        for j = 1:numpHtested
            if ((j == 1)||(j == 2)||(j == 3))
                RRs2{j,i} = RRs{j,i}(51:end);
            elseif j == 4
                RRs2{j,i} = RRs{j,i}(41:end);
            elseif j == 5
                RRs2{j,i} = RRs{j,i}(36:end);
            elseif j == 6
                RRs2{j,i} = RRs{j,i}(26:end);
            elseif j == 7
                RRs2{j,i} = RRs{j,i}(21:end);
            else
                RRs2{j,i} = RRs{j,i}(16:end);
            end
        end
    end
    data.RRs = RRs2;
    % time
    time2 = cell(size(time));
    for i = 1:DFs
        for j = 1:numpHtested
            if ((j == 1)||(j == 2)||(j == 3))
                time2{j,i} = time{j,i}(51:end);
            elseif j == 4
                time2{j,i} = time{j,i}(41:end);
            elseif j == 5
                time2{j,i} = time{j,i}(36:end);
            elseif j == 6
                time2{j,i} = time{j,i}(26:end);
            elseif j == 7
                time2{j,i} = time{j,i}(21:end);
            else
                time2{j,i} = time{j,i}(16:end);
            end
        end
    end
    clear NADH, NADH = NADH2; clear NADH2
    data.conc_mean = NADH;    
    clear time, time = time2; clear time2
    data.time = time;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    pHvals = unique(import_pgi.treatedData.pH_corrected);
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
        if setup.caseStudyPGI == 1
            ylim([0 0.15])
%             xlim([0 1000])
            xlim([0 800])
        end
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
    suptitle(suptitleName);
    
% % % %     figure
% % % %     plot(pHvals, Vmax(:,4),'.-')
% % % %     title('Starting estimate: Vmax [mM s-1] vs pH')    
end

        
%% (0.1) Experimental Vmax determination
setup.plotOutput = 0;
% %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 5; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
% select start point (this needs manual selection deppending on previous plots)
dp_start = 2 * ones(size(data.conc_mean));
dp_start(:,1) = 30 * ones(size(dp_start(:,1)));
dp_start(:,2) = 20 * ones(size(dp_start(:,2)));
dp_start(:,3) = 10 * ones(size(dp_start(:,3)));
dp_start(:,4) = 3 * ones(size(dp_start(:,4)));
% blank cell total length
total_len = zeros(size(dp_start));
% DFs considered
DFarray = [1/8 1/4 1/2 1/1];
% idxs2consider
% % % % idxs2consider = ones(size(DF));
idxs2consider = [0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1];

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
plength = 14; % Kms (2) + Vms (1) * numpH (12) (Keq is fixed to experimental data)
% plength = 15; % Kms (2) + Vms (1) * numpH (12) (Keq is fixed to experimental data) + link (1)
% plength = 26; % Kms (2) + Vms (1) * numpH (12) (Keq is fixed to experimental data) + link (12)
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
lb = -1*ones(1,plength);
%     lb(15) = -5;
%     lb(15:26) = -3;
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
% % parameter estimation
% tic
% [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t = toc;


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
selLambdaPos = 1; %13;%14; % 6, #1
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
% xres_selected = xres; %lambdalist based in 'ones', lam=0.3, loc=15.
% xres_selected = array_xres{15}; %lambdalist based in 'ones', lam=0.3, loc=15.
xres_selected = array_xres{selLambdaPos};

kg6p = ones(numpHtested,1);
kf6p = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    kg6p(i) = 1.0257 * 10.^xres_selected(1); %mM
    kf6p(i) = 0.307 * 10.^xres_selected(2); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+2);
end
keq_pgi = setup.Keq_PGI;
keq_pfk = setup.Keq_PFK;
keq_fba = setup.Keq_FBA;
keq_tpi = setup.Keq_TPI;
keq_gpd = setup.Keq_GPD;

% limits
kg6p_up = ones(numpHtested,1);
kf6p_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);

kg6p_down = ones(numpHtested,1);
kf6p_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);

for i = 1:numpHtested
    % up
    kg6p_up(i) = 1.0257 * 10.^(xres_selected(1)+stdp(1)); %mM
    kf6p_up(i) = 0.307 * 10.^(xres_selected(2)+stdp(2)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)+stdp(i+2)); %mM
    % down
    kg6p_down(i) = 1.0257 * 10.^(xres_selected(1)-stdp(1)); %mM
    kf6p_down(i) = 0.307 * 10.^(xres_selected(2)-stdp(2)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)-stdp(i+2)); %mM
end


%% (3.3) Study on paarameter values: plotting
vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;

% figure(105)
figure(106)

subplot(3,3,1) % vm
%     plot(pHarray,vm_up,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm_down,'.-','color',[0.5 0.5 0.5]), hold on, 
% plot(pHarray,vm,'.-','color','black')
% title('v_{m} [mM s^{1}]')
    plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_uChange,'.-','color','black')
title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})    
    
subplot(3,3,2) % kg6p
    plot(pHarray,kg6p_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kg6p_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kg6p,'.-','color','black')
    title('k_{g6p} [mM]')
subplot(3,3,3) % kf6p
    plot(pHarray,kf6p_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kf6p_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kf6p,'.-','color','black')
    title('k_{f6p} [mM]')

subplot(3,3,4) % vm(log_difference)
    plot(pHarray,xres_selected(3:14)','.-','color','black')
    title('v_{m} []')

subplot(3,3,5) % keq_pgi = setup.Keq_PGI;
    plot(pHarray,keq_pgi,'.-','color','black')
    title('k_{eq.PGI} [mM]')
subplot(3,3,6) % keq_pfk = setup.Keq_PFK;
    plot(pHarray,keq_pfk,'.-','color','black')
    title('k_{eq.PFK} [mM]')
subplot(3,3,7) % keq_fba = setup.Keq_FBA;
    plot(pHarray,keq_fba,'.-','color','black')
    title('k_{eq.FBA} [mM]')
subplot(3,3,8) % keq_tpi = setup.Keq_TPI;
    plot(pHarray,keq_tpi,'.-','color','black')
    title('k_{eq.TPI} [mM]')
% % subplot(3,3,9) % keq_gpd = setup.Keq_GPD;
% %     plot(pHarray,keq_gpd,'.-','color','black')
% %     title('k_{eq.GPD} [mM]')
    
subplot(339) % vm experimental
    errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
    title({'v_{m} experimental';'not normalized'})

suptitle('PGI: parameter estimates vs pH')


% %%
% set(105,'color','w');
% savefig(105,'PGI_pVals_pCis');
%%
output_pgi.xres_selected = xres_selected;

output_pgi.pHarray = pHarray;
output_pgi.idxs2consider = idxs2consider;
output_pgi.DF = data.DF;

output_pgi.kg6p = kg6p;% = ones(numpHtested,1);
output_pgi.kf6p = kf6p;% = ones(numpHtested,1);
output_pgi.vm = vm;% = zeros(numpHtested,1);
output_pgi.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_pgi.kg6p_up = kg6p_up;% = ones(numpHtested,1);
output_pgi.kf6p_up = kf6p_up;% = ones(numpHtested,1);
output_pgi.vm_up = vm_up;% = zeros(numpHtested,1);
output_pgi.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_pgi.kg6p_down = kg6p_down;% = ones(numpHtested,1);
output_pgi.kf6p_down = kf6p_down;% = ones(numpHtested,1);
output_pgi.vm_down = vm_down;% = zeros(numpHtested,1);
output_pgi.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_pgi.keq_pgi = keq_pgi;% = setup.Keq_PGI;
output_pgi.keq_pfk = keq_pfk;% = setup.Keq_PFK;
output_pgi.keq_fba = keq_fba;% = setup.Keq_FBA;
output_pgi.keq_tpi = keq_tpi;% = setup.Keq_TPI;
output_pgi.keq_gpd = keq_gpd;% = setup.Keq_GPD;

%% saving output
output_pgi.Vmax_experimental = Vmax_experimental;
output_pgi.stDev_experimental = stDev_experimental;
saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_parEst.mat'];
save(saveName,'output_pgi');
set(105,'color','white'), savefig(105,['results/',setup.enzymeName,'/',setup.enzymeName, '_parEstimates.fig']);


% % % % %% reverse to reload
% % % % load('output_pgi.mat','output_pgi');
% % % % 
% % % % % output_pgi.xres_selected = xres_selected;
% % % % 
% % % % xres_selected = output_pgi.xres_selected;
% % % % 
% % % % pHarray = output_pgi.pHarray;
% % % % 
% % % % kg6p = output_pgi.kg6p;% = output_pgi.ones(numpHtested,1);
% % % % kf6p = output_pgi.kf6p;% = ones(numpHtested,1);
% % % % vm = output_pgi.vm;% = zeros(numpHtested,1);
% % % % 
% % % % kg6p_up = output_pgi.kg6p_up;% = ones(numpHtested,1);
% % % % kf6p_up = output_pgi.kf6p_up;% = ones(numpHtested,1);
% % % % vm_up = output_pgi.vm_up;% = zeros(numpHtested,1);
% % % % 
% % % % kg6p_down = output_pgi.kg6p_down;% = ones(numpHtested,1);
% % % % kf6p_down = output_pgi.kf6p_down;% = ones(numpHtested,1);
% % % % vm_down = output_pgi.vm_down;% = zeros(numpHtested,1);
% % % % 
% % % % keq_pgi = output_pgi.keq_pgi;% = setup.Keq_PGI;
% % % % keq_pfk = output_pgi.keq_pfk;% = setup.Keq_PFK;
% % % % keq_fba = output_pgi.keq_fba;% = setup.Keq_FBA;
% % % % keq_tpi = output_pgi.keq_tpi;% = setup.Keq_TPI;
% % % % keq_gpd = output_pgi.keq_gpd;% = setup.Keq_GPD;
% % % % x = 1; end


