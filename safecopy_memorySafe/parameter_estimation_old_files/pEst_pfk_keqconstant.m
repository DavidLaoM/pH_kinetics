% % PEST_ALD.m
% Parameter estimation for the data in the Enolase assay.
% Vm are estimated changing with pH
% Kms are assumed constant
% Keq from the eQuilibrator


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
    
    figure
    plot(pHvals, Vmax(:,4),'.-')
    title('Starting estimate: Vmax [mM s-1] vs pH')    
end
% figure(11)
% plot(setup.fullpHarray, setup.Keq_HXK, '.-','MarkerSize',12,'LineWidth',1.2)
% xlabel('pH values','FontSize',16)
% ylabel('K_{eq.HXK}','FontSize',16)
% 
% figure(12)
% plot(setup.fullpHarray, setup.Keq_G6PDH, '.-','MarkerSize',12,'LineWidth',1.2)
% xlabel('pH values','FontSize',16)
% ylabel('K_{eq.G6PDH}','FontSize',16)


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
setup.costfun = 1;

setup.weightData = 1;
setup.weightDataEsp = ones(1,numpHtested);
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 1E-5; % by now just testing

% Km fixed
optfun = @costfun_Kmfixed;
% plength = 25; % vms (12) + other params (12)
% x_temp = zeros(1,plength);
plength = 12; % vms (12) + other params off
x_temp = -0.45*ones(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');
% %%
% % testing the costfunction
% setup.plotEachSimCF = 1;
% [error] = optfun(x_temp,data,setup);
% setup.plotEachSimCF = 0;


%% making Keq constant at the values of pH6.81
for i = 1:length(setup.Keq_FBA)
    setup.Keq_FBA(i) = setup.Keq_FBA(6);
    setup.Keq_TPI(i) = setup.Keq_TPI(6);
    setup.Keq_GPD(i) = setup.Keq_GPD(6);
    setup.Keq_PFK(i) = setup.Keq_PFK(6);
end


%%
% parameter estimation
tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t = toc;

% %% display result
disp(xres);
setup.plotEachSimCF = 1;
[error] = optfun(xres,data,setup);
% [error] = optfun(xres2,data,setup);
setup.plotEachSimCF = 0;


%%
% output system simulation
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
[error] = optfun(xres,data,setup);
setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

% %% calculation error
N = length(error);
Jacobian = full(Jacobian);  
varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
stdp = sqrt(diag(varp));


%% (3.2) Study on paarameter values: recalculation
% parameter values
xres_selected = xres; %lambdalist based in 'ones', lam=0.1, loc=5.
% xres_selected = array_xres{15}; %lambdalist based in 'ones', lam=0.1, loc=5.
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i);
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
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i)+stdp(i));
    %down
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i)-stdp(i));
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

subplot(332) % vm(log_difference)
    plot(pHarray,xres_selected','.-','color','black')
    title('v_{m} difference []')

subplot(337) % keq_pfk
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


%% final saving

% plots simulations (101 - metabolites, 102 - fluxes, 105 - parameters estimated)
% set(101,'color','w');
% set(102,'color','w');
% set(105,'color','w');
% savefig(101,'results/pfk/PFK_fit_metabolites.fig');
% savefig(102,'results/pfk/PFK_fit_rates.fig');
% savefig(105,'results/pfk/PFK_pVals_pCis');
output_pfk_keq_constant = output_pfk;
save('results/pfk/output_pfk_keq_constant.mat','output_pfk_keq_constant');

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

