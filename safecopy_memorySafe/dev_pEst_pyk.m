% % % % function x=dev_pEst_pyk
% % PEST_PYK.m
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
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 1;
    setup.caseStudyTPI = 0;
    setup.caseStudyENO = 0;
    selectSetup_pH;
    % added
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_pyk = expData.pyk;
    
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
    
    pHarray = unique(import_pyk.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_pyk.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_pyk.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_pyk.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_pyk.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_pyk.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_pyk.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_pyk.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_pyk.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_pyk.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_pyk.rawData.absorbance_corrected{4,4};
    temp2 = import_pyk.rawData.absorbance_corrected{5,4};
    temp3 = import_pyk.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_pyk.rawData.time{1};
    
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
            NADH2{j,i} = NADH{j,i}(2:end);
        end
    end
    % bring time to zero start
    for i = 1:DFs
        for j = 1:numpHtested
            dps = length(time{j,i});
            for k = 1:dps
                time{j,i}(k) = time{j,i}(k) - 5;
            end
        end
    end
    % GPD reaction rate
    RRs2 = cell(size(RRs));
    for i = 1:DFs
        for j = 1:numpHtested
            RRs2{j,i} = RRs{j,i}(2:end);
        end
    end
    data.RRs = RRs2;
    % time
    time2 = cell(size(time));
    for i = 1:DFs
        for j = 1:numpHtested
            time2{j,i} = time{j,i}(2:end);
        end
    end
    clear NADH, NADH = NADH2; clear NADH2
    data.conc_mean = NADH;    
    clear time, time = time2; clear time2
    data.time = time;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    pHvals = unique(import_pyk.treatedData.pH_corrected);
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
    minwindow = 7; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
% select start point (this needs manual selection deppending on previous plots)
dp_start = ones(size(data.conc_mean));
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
plength = 18; % others (6) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
lb = -1*ones(1,plength);
options = optimset('Display','iter');
% %%
% % testing the costfunction
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
selLambdaPos = 1;%13;%17;%16;
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
% xres_selected = array_xres{14}; %lambdalist based in 'ones', lam=0.1, loc=5.
xres_selected = array_xres{selLambdaPos};

kadp = ones(numpHtested,1);
katp = ones(numpHtested,1);
kfbp = ones(numpHtested,1);
kpep = ones(numpHtested,1);
L = ones(numpHtested,1);
nHill = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    kadp(i) = 0.243 * 10.^xres_selected(1); %mM
    katp(i) = 9.3 * 10.^xres_selected(2); %mM
    kfbp(i) = 0.2 * 10.^xres_selected(3); %mM
    kpep(i) = 0.281 * 10.^xres_selected(4); %mM
    L(i) = 60000 * 10.^xres_selected(5);
    nHill(i) = 4 * 10.^xres_selected(6);
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+6);
end
keq_pyk = setup.Keq_PYK;
keq_ldh = setup.Keq_LDH;

% limits
kadp_up = ones(numpHtested,1);
katp_up = ones(numpHtested,1);
kfbp_up = ones(numpHtested,1);
kpep_up = ones(numpHtested,1);
L_up = ones(numpHtested,1);
nHill_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);

kadp_down = ones(numpHtested,1);
katp_down = ones(numpHtested,1);
kfbp_down = ones(numpHtested,1);
kpep_down = ones(numpHtested,1);
L_down = ones(numpHtested,1);
nHill_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    % up
    kadp_up(i) = 0.243 * 10.^(xres_selected(1)+stdp(1)); %mM
    katp_up(i) = 9.3 * 10.^(xres_selected(2)+stdp(2)); %mM
    kfbp_up(i) = 0.2 * 10.^(xres_selected(3)+stdp(3)); %mM
    kpep_up(i) = 0.281 * 10.^(xres_selected(4)+stdp(4)); %mM
    L_up(i) = 60000 * 10.^(xres_selected(5)+stdp(5)); %mM
    nHill_up(i) = 4 * 10.^(xres_selected(6)+stdp(6)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+6)+stdp(i+6)); %mM
    % down
    kadp_down(i) = 0.243 * 10.^(xres_selected(1)-stdp(1)); %mM
    katp_down(i) = 9.3 * 10.^(xres_selected(2)-stdp(2)); %mM
    kfbp_down(i) = 0.2 * 10.^(xres_selected(3)-stdp(3)); %mM
    kpep_down(i) = 0.281 * 10.^(xres_selected(4)-stdp(4)); %mM
    L_down(i) = 60000 * 10.^(xres_selected(5)-stdp(5)); %mM
    nHill_down(i) = 4 * 10.^(xres_selected(6)-stdp(6)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+6)-stdp(i+6)); %mM
end


%% (3.3) Study on paarameter values: plotting
vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;

figure(105)
% figure(106)

subplot(4,3,1) % vm
%     plot(pHarray,vm_up,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm_down,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm,'.-','color','black')
%     title('v_{m} [mM]')
    plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_uChange,'.-','color','black')
    title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})
subplot(4,3,4) % L
    plot(pHarray,L_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,L_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,L,'.-','color','black')
    title('k_{L} [mM]')
subplot(4,3,5) % nHill
    plot(pHarray,nHill_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,nHill_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,nHill,'.-','color','black')
    title('k_{nHill} [mM]')

subplot(4,3,9) % kadp
    plot(pHarray,kadp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kadp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kadp,'.-','color','black')
    title('k_{adp} [mM]')
subplot(4,3,6) % katp
    plot(pHarray,katp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,katp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,katp,'.-','color','black')
    title('k_{atp} [mM]')
subplot(4,3,7) % kfbp
    plot(pHarray,kfbp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kfbp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kfbp,'.-','color','black')
    title('k_{fbp} [mM]')
subplot(4,3,8) % kpep
    plot(pHarray,kpep_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kpep_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kpep,'.-','color','black')
    title('k_{pep} [mM]')

subplot(4,3,10) % keq_pyk
    plot(pHarray,keq_pyk,'.-','color','black')
    title('k_{eq.PYK} [mM]')
subplot(4,3,11) % keq_ldh
    plot(pHarray,keq_ldh,'.-','color','black')
    title('k_{eq.LDH} [mM]')

subplot(4,3,2)
plot(pHarray,xres_selected(7:18)','.-','color','black')
title('vm reference')

subplot(4,3,3) % vm experimental
    errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
    title({'v_{m} experimental';'not normalized'})

suptitle('parameter estimates vs pH')

% %%
% set(105,'color','w');
% savefig(105,'PYK_pVals_pCis');
%%
output_pyk.xres_selected = xres_selected;

output_pyk.pHarray = pHarray;
output_pyk.idxs2consider = idxs2consider;
output_pyk.DF = data.DF;

output_pyk.kadp = kadp;% = ones(numpHtested,1);
output_pyk.katp = katp;% = ones(numpHtested,1);
output_pyk.kfbp = kfbp;% = ones(numpHtested,1);
output_pyk.kpep = kpep;% = ones(numpHtested,1);
output_pyk.L = L;% = ones(numpHtested,1);
output_pyk.nHill = nHill;% = ones(numpHtested,1);
output_pyk.vm = vm;% = zeros(numpHtested,1);
output_pyk.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_pyk.kadp_up = kadp_up;% = ones(numpHtested,1);
output_pyk.katp_up = katp_up;% = ones(numpHtested,1);
output_pyk.kfbp_up = kfbp_up;% = ones(numpHtested,1);
output_pyk.kpep_up = kpep_up;% = ones(numpHtested,1);
output_pyk.L_up = L_up;% = ones(numpHtested,1);
output_pyk.nHill_up = nHill_up;% = ones(numpHtested,1);
output_pyk.vm_up = vm_up;% = zeros(numpHtested,1);
output_pyk.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_pyk.kadp_down = kadp_down;% = ones(numpHtested,1);
output_pyk.katp_down = katp_down;% = ones(numpHtested,1);
output_pyk.kfbp_down = kfbp_down;% = ones(numpHtested,1);
output_pyk.kpep_down = kpep_down;% = ones(numpHtested,1);
output_pyk.L_down = L_down;% = ones(numpHtested,1);
output_pyk.nHill_down = nHill_down;% = ones(numpHtested,1);
output_pyk.vm_down = vm_down;% = zeros(numpHtested,1);
output_pyk.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_pyk.keq_pyk = keq_pyk;% = setup.Keq_PYK;
output_pyk.keq_ldh = keq_ldh;% = setup.Keq_LDH;


%% saving output
output_pyk.Vmax_experimental = Vmax_experimental;
output_pyk.stDev_experimental = stDev_experimental;
saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_parEst.mat'];
save(saveName,'output_pyk');
set(105,'color','white'), savefig(105,['results/',setup.enzymeName,'/',setup.enzymeName, '_parEstimates.fig']);


% % % % %%
% % % % load('output_pyk.mat','output_pyk');
% % % % 
% % % % % output_pyk.xres_selected = xres_selected;
% % % % 
% % % % xres_selected = output_pyk.xres_selected;
% % % % 
% % % % pHarray = output_pyk.pHarray;
% % % % 
% % % % kadp = output_pyk.kadp;%
% % % % katp = output_pyk.katp;%
% % % % kfbp = output_pyk.kfbp;%
% % % % kpep = output_pyk.kpep;%
% % % % L = output_pyk.L;%
% % % % nHill = output_pyk.nHill;%
% % % % vm = output_pyk.vm;%
% % % % 
% % % % kadp_up = output_pyk.kadp_up;%
% % % % katp_up = output_pyk.katp_up;%
% % % % kfbp_up = output_pyk.kfbp_up;%
% % % % kpep_up = output_pyk.kpep_up;%
% % % % L_up = output_pyk.L_up;%
% % % % nHill_up = output_pyk.nHill_up;%
% % % % vm_up = output_pyk.vm_up;%
% % % % 
% % % % kadp_down = output_pyk.kadp_down;%
% % % % katp_down = output_pyk.katp_down;%
% % % % kfbp_down = output_pyk.kfbp_down;%
% % % % kpep_down = output_pyk.kpep_down;%
% % % % L_down = output_pyk.L_down;%
% % % % nHill_down = output_pyk.nHill_down;%
% % % % vm_down = output_pyk.vm_down;%
% % % % 
% % % % keq_pyk = output_pyk.keq_pyk;%
% % % % keq_ldh = output_pyk.keq_ldh;%
% % % % x = 1; end


