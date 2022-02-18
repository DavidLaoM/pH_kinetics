% % PEST_ENO.m
% Parameter estimation for the data in the Enolase assay.
% Vm are estimated changing with pH
% Kms are assumed constant
% Keq from the eQuilibrator

% % Contents
% AP 1  Estimation at pH #6
% AP 2  Estimation for all the other pH values
% AP 3  Regularization for the selected pH value

%% (0) Setup and data load
clear
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
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    setup.caseStudyENO = 1;
    selectSetup_pH;
    % added
    setup.saveOutput = 0;

    load('expData.mat','expData');
    import_eno = expData.eno;

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

    PEP = blankCell;
    Vmax = blank;
    for i = 1:(DFs*numpHtested)
        tempDiff = conc_mean{i} - conc_mean{i}(1); % all stoichiometries are 1-to-1.
        PEP{i} = conc_mean{i};
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
    data.chosenPEPini = 0.4;
    temp1 = import_eno.rawData.absorbance_corrected{4,4};
    temp2 = import_eno.rawData.absorbance_corrected{5,4};
    temp3 = import_eno.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_eno.rawData.time{1};

    pHvals = unique(import_eno.treatedData.pH_corrected);
    % visualize: check calculations made
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numpHtested
        subplot(3,4,i)
        for j = 1:DFs
            plot(time{i,j},PEP{i,j},'.-')
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
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': PEP concentration profile'];
    suptitle(suptitleName);
end


%% (1.1) Simple parameter fit. Parameter estimation. pH#6
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';

setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
setup.costfun = 1;

setup.weightData = 1;
setup.weightDataEsp = ones(1,12);
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% % Km fixed
% optfun = @costfun_Kmfixed;
% plength = 14; % Kms (2) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
% x_temp = zeros(1,plength);
% ub = 3*ones(1,plength);
% lb = -3*ones(1,plength);
% options = optimset('Display','iter');

% Vm fixed
optfun = @costfun_Vmfixed;
setup.startVm = 12;
plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');
% %%
% % testing the costfunction
% [error] = optfun(x_temp,data,setup);

%% parameter estimation
tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t = toc;

% %% (1.2) Simple parameter fit. Results Visualization
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
[error] = optfun(xres,data,setup);


%% (1.2) Simple parameter fit. Saving
set(101,'color','w');
set(102,'color','w');
savefig(101,'ENOfit_metabolites_Vmfixed.fig');
savefig(102,'ENOfit_rates_Vmfixed.fig');


%% (2.1) Simple parameter fit. Parameter estimation. All pH#
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
xres_cell = cell(numpHtested,1);
% loop to estimate each pH value
for i = 1:numpHtested
    % setup
    setup.startVm = i;
    setup.factMultiply = 1; %2;
    % parameter estimation
    tic
    [xres_temp,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;
    xres_cell{i} = xres_temp;
end
%% (2.2) plotting 1 by 1 + overall error
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
error_array = zeros(numpHtested,1);
for i = 1:numpHtested
    xres_temp = xres_cell{i};
    [error] = optfun(xres_temp,data,setup);
    error_array(i) = sum(abs(error));
end
%%
figure(1001)
semilogy(setup.fullpHarray, error_array, '.-')
hold on
line([6 8],[1E1 1E1],'color','blue','LineStyle','--')
xlabel('pH value')
ylabel('overall error simulation vs experimental data')
title({'Overall error obtained at each pH value';'\color{blue} error with km_{fixed} and vm_{variable} below 0.1.'})
ylim([1E0 1E4])
%%
set(1001,'color','w');
savefig(1001,'ENO_vmfixed_checkAllpHs.fig')
save('ENO_vmfixed_testAllpH.mat','xres_cell','error_array')


%% (3.1) Parameter estimation with regularization
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';

setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
setup.costfun = 1;

setup.weightData = 1;
setup.weightDataEsp = ones(1,12);
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0;

% % Km fixed
% optfun = @costfun_Kmfixed;
% plength = 14; % Kms (1) + Vms (1) * numpH (10) + (Keq is fixed to experimental data)
% x_temp = zeros(1,plength);
% ub = 3*ones(1,plength);
% lb = -3*ones(1,plength);
% options = optimset('Display','iter');

% Vm fixed
optfun = @costfun_Vmfixed;
setup.startVm = 12;
plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');

% % testing the costfunction
% [error] = optfun(x_temp,data,setup);
% 
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
%
% % lambdalist = [...
% %     1E-5,...
% %     1E-4,...
% %     1E-3,...
% %     1E-2,...
% %     1E-1,...
% %     1E0,...
% %     1E1,...
% %     1E2,...
% %     1E3,...
% %     1E4,...
% %     1E5,...
% %     ];
% lambdalist = [0 1E-5 1E0];
array_xres = cell(1,length(lambdalist));
array_eData = cell(1,length(lambdalist));
array_eParams = cell(1,length(lambdalist));
for i = 1:length(lambdalist)
    fprintf('pEst for lambda=%d\n',lambdalist(i));
    setup.selectedLambda = lambdalist(i);
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;
    setup.selectedLambda = 1;
    [error] = optfun(xres,data,setup);
    % seting in output arrays
    array_xres{i} = xres;
    array_eData{i} = error(1:end-14);
    array_eParams{i} = error(end-13:end);
end


% %% (2.2) Regularization. Results Visualization
eData = zeros(1,length(lambdalist));
eParameters = zeros(1,length(lambdalist));
for i = 1:length(lambdalist)
    eData(i) = sum(abs(array_eData{i}));
    eParameters(i) = sum(abs(array_eParams{i}));
end

figure(103)
yyaxis left
semilogx(lambdalist,eParameters,'o-')
hold on
yyaxis right
semilogx(lambdalist,eData,'o-')
legend('e_{Parameters}','e_{Data}')
suptitle('Regularization Enolase. Errors vs lambda')
%%
xres_selected = array_xres{17};
% xres_selected = xres;
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
[error] = optfun(xres_selected,data,setup);


%% (2.3) Simple parameter fit. Saving
set(112,'color','w');
set(212,'color','w');
set(103,'color','w');
savefig(112,'ENOfit_vmfixed_metabolites_lam5e_1.fig');
savefig(212,'ENOfit_vmfixed_rates_lam5e_1.fig');
savefig(103,'ENO_vmfixed_regularization.fig');


%% (3.1) Study on paarameter values: estimation with the lambda value
for case2 = 1 % getting the confidence intervals
    setup.ode = 'vanHeerden2014';
    setup.sourceVm = 'experimentalSlopes';
    setup.ode_pH = 'on';

    setup.plotResults = 0;
    setup.plotEachSimCF = 0;
    setup.simAllProfiles = 0;
    setup.plotEachSim = 1;

    setup.numpHtested = numpHtested;
    setup.DFstudy = 1:4;
    setup.costfun = 1;

    setup.weightData = 1;
    setup.weightDataEsp = ones(1,12);
    setup.weightHaldane = 0; % off for this case (Keq fixed)
    setup.selectedLambda = lambdalist(17); % by now just testing

%     % Km fixed
%     optfun = @costfun_Kmfixed;
%     plength = 14; % Kms (1) + Vms (1) * numpH (10) + (Keq is fixed to experimental data)
%     x_temp = zeros(1,plength);
%     ub = 3*ones(1,plength);
%     lb = -3*ones(1,plength);
%     options = optimset('Display','iter');
    
    % Vm fixed
    optfun = @costfun_Vmfixed;
    setup.startVm = 12;
    plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
    options = optimset('Display','iter');
    
    
    % % testing the costfunction
    % [error] = optfun(x_temp,data,setup);

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;
end
% %%
[error] = optfun(xres,data,setup);
for case2 = 1 % confidence intervals
    N = length(error);
    Jacobian = full(Jacobian);  
    varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
    stdp = sqrt(diag(varp));
end


%% (3.2) Study on paarameter values: recalculation
% parameter values
% xres_selected = array_xres{19}; % lambda value = 1
xres_selected = xres; % lambda value = 1
k2pg = ones(numpHtested,1);
kpep = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    k2pg(i) = 0.043 * 10.^xres_selected(1); %mM
    kpep(i) = 0.5 * 10.^xres_selected(2); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+2);
end
keq = setup.keq;
% limits
k2pg_up = ones(numpHtested,1);
kpep_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);
k2pg_down = ones(numpHtested,1);
kpep_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    k2pg_up(i) = 0.043 * 10.^(xres_selected(1)+stdp(1)); %mM
    kpep_up(i) = 0.5 * 10.^(xres_selected(2)+stdp(2)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)+stdp(i+2)); %mM s^[-1]
    k2pg_down(i) = 0.043 * 10.^(xres_selected(1)-stdp(1)); %mM
    kpep_down(i) = 0.5 * 10.^(xres_selected(2)-stdp(2)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)-stdp(i+2)); %mM s^[-1]
end


%% (3.3) Study on paarameter values: plotting
% figure(105)
vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;

figure(105)

subplot(231) % vm
%     plot(pHarray,vm_up,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm_down,'.-','color',[0.5 0.5 0.5]), hold on, 
% plot(pHarray,vm,'.-','color','black')
% title('v_{m} [mM s^{1}]')
    plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,vm_uChange,'.-','color','black')
title({'v_{m} [umol_{PEP} mg_{P}^{-1} min^{-1}]';'not normalized'})

subplot(232) % k2pg
    plot(pHarray,k2pg_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,k2pg_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,k2pg,'.-','color','black')
title('k_{p2g} [mM]')

subplot(233) % kpep
    plot(pHarray,kpep_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kpep_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kpep,'.-','color','black')
title('k_{pep} [mM]')

subplot(234) % vm(log_difference)
    plot(pHarray,xres_selected(3:14)','.-','color','black')
    title('v_{m} difference []')

subplot(236) % keq
plot(pHarray,keq,'.-','color','black')
title('k_{eq} [mM]')

suptitle('ENO: parameter estimates vs pH')


% %%
% set(105,'color','w');
% savefig(105,'ENO_pVals_pCis_vmfixed');
% %%
% output_eno.xres_selected = xres_selected;
% 
% output_eno.pHarray = pHarray;
% 
% output_eno.k2pg = k2pg;% = ones(numpHtested,1);
% output_eno.kpep = kpep;% = ones(numpHtested,1);
% output_eno.vm = vm;% = zeros(numpHtested,1);
% 
% output_eno.k2pg_up = k2pg_up;% = ones(numpHtested,1);
% output_eno.kpep_up = kpep_up;% = ones(numpHtested,1);
% output_eno.vm_up = vm_up;% = zeros(numpHtested,1);
% 
% output_eno.k2pg_down = k2pg_down;% = ones(numpHtested,1);
% output_eno.kpep_down = kpep_down;% = ones(numpHtested,1);
% output_eno.vm_down = vm_down;% = zeros(numpHtested,1);
% 
% output_eno.keq_eno = keq;% = setup.Keq_PGI;
% 
% % % % % save('output_eno.mat','output_eno');
