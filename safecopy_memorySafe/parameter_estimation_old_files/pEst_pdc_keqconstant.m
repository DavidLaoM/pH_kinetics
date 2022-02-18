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
    setup.caseStudyPDC = 1;
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
    import_pdc = expData.pdc;
    
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
    
    pHarray = unique(import_pdc.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_pdc.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_pdc.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_pdc.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_pdc.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_pdc.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_pdc.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_pdc.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_pdc.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_pdc.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_pdc.rawData.absorbance_corrected{4,4};
    temp2 = import_pdc.rawData.absorbance_corrected{5,4};
    temp3 = import_pdc.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_pdc.rawData.time{1};
    
        % (1) Correct for minimum value
        % (2) Bring the minimum to zero (apply to all)
        % (3) In principle, use the 3 first dilution rates
        % (4) Watch out with the dilution factors (first 2 cases are
        % reversed)
        
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Adjusted to PDC
    
    % relocating well DF2,#7,8
    dps = length(NADH{1,1});
    for i = 3
        for j = 7:8
                for k = 1:dps
                    NADH{j,i}(k) = NADH{j,i}(k) - (NADH{j,3}(end)-NADH{j,4}(end));
                end
        end
    end
    
    % Directly changing the concentration here, since the extinction
    % coefficient did not change.
    endPoint = zeros(numpHtested,DFs);
    % locate the minimum
    for i = 1:DFs
        for j = 1:numpHtested
            endPoint(j,i) = min(NADH{j,i});
%             endPoint(j,i) = min([NADH{j,1}' NADH{j,2}' NADH{j,3}' NADH{j,4}']);
        end
    end
    % bringing the minimum to zero
    for i = 1:DFs
        for j = 1:numpHtested
            if((i==4)||(i==3)||((i==2)&&((j>=2)&&(j<=7))))
                for k = 1:dps
                    NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,i);
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
            for k = (tempidx+1):dps
                NADH{j,i}(k) = NADH{j,i}(tempidx);
            end
        end
    end
    
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    pHvals = unique(import_pdc.treatedData.pH_corrected);
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
% % % % setup.selectedLambda = 1E-1; % (TO PROPERLY SET)

% Km fixed
optfun = @costfun_Kmfixed;
plength = 14; % Kms (2) + Vms (1) * numpH (12)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');
% %%
% % testing the costfunction
% setup.plotEachSimCF = 1;
% [error] = optfun(x_temp,data,setup);
% setup.plotEachSimCF = 0;


%% making Keq constant at the values of pH6.81
for i = 1:length(setup.Keq_ADH)
    setup.Keq_ADH(i) = setup.Keq_ADH(6);
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
kpyr = ones(numpHtested,1);
hill = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    kpyr(i) = 8.5 * 10.^xres_selected(1); %mM
    hill(i) = 1.9 * 10.^xres_selected(2); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+2);
end
keq_adh = setup.Keq_ADH;

% limits
kpyr_up = ones(numpHtested,1);
hill_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);
kpyr_down = ones(numpHtested,1);
hill_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    %up
    kpyr_up(i) = 8.5 * 10.^(xres_selected(1)+stdp(1)); %mM
    hill_up(i) = 1.9 * 10.^(xres_selected(2)+stdp(2)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)+stdp(i+2));
    %down
    kpyr_down(i) = 8.5 * 10.^(xres_selected(1)-stdp(1)); %mM
    hill_down(i) = 1.9 * 10.^(xres_selected(2)-stdp(2)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)-stdp(i+2));
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
%     plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_uChange,'.-','color','black')
title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

subplot(332) % vm(log_difference)
    plot(pHarray,xres_selected(3:14)','.-','color','black')
    title('v_{m} difference []')

subplot(334) % kpyr
    plot(pHarray,kpyr_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kpyr_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,kpyr,'.-','color','black')
title('k_{pyr} [mM]')

subplot(335) % hill
    plot(pHarray,hill_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,hill_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pHarray,hill,'.-','color','black')
title('hill [mM]')

subplot(337) % keq_adh
plot(pHarray,keq_adh,'.-','color','black')
title('k_{eq.ADH} [mM]')

suptitle('PDC: parameter estimates vs pH')


%%
output_pdc.xres_selected = xres_selected;

output_pdc.pHarray = pHarray;

output_pdc.kpyr = kpyr;% = ones(numpHtested,1);
output_pdc.hill = hill;% = ones(numpHtested,1);
output_pdc.vm = vm;% = zeros(numpHtested,1);
output_pdc.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_pdc.kpyr_up = kpyr_up;% = ones(numpHtested,1);
output_pdc.hill_up = hill_up;% = ones(numpHtested,1);
output_pdc.vm_up = vm_up;% = zeros(numpHtested,1);
output_pdc.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_pdc.kpyr_down = hill_down;% = ones(numpHtested,1);
output_pdc.hill_down = hill_down;% = ones(numpHtested,1);
output_pdc.vm_down = vm_down;% = zeros(numpHtested,1);
output_pdc.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_pdc.keq_adh = keq_adh;% = setup.Keq_ADH;


%% final saving

% plots simulations (101 - metabolites, 102 - fluxes, 105 - parameters estimated)
% set(101,'color','w');
% set(102,'color','w');
% set(105,'color','w');
% savefig(101,'results/pdc/PDC_fit_metabolites.fig');
% savefig(102,'results/pdc/PDC_fit_rates.fig');
% savefig(105,'results/pdc/PDC_pVals_pCis');
% save('results/pdc/output_pdc.mat','output_pdc');
output_pdc_keq_constant = output_pdc;
save('results/pdc/output_pdc_keq_constant.mat','output_pdc_keq_constant');

% %% reverse to reload
% load('output_pdc.mat','output_pdc');
% 
% xres_selected = output_pdc.xres_selected;
% 
% pHarray = output_pdc.pHarray;
% 
% kfbp = output_pdc.kfbp;%
% hill = output_pdc.hill;%
% kdhap = output_pdc.kdhap;%
% vm = output_pdc.vm;%
% 
% kfbp_up = output_pdc.kfbp_up;%
% hill_up = output_pdc.hill_up;%
% kdhap_up = output_pdc.kdhap_up;%
% vm_up = output_pdc.vm_up;%
% 
% kfbp_down = output_pdc.kfbp_down;%
% hill_down = output_pdc.hill_down;%
% kdhap_down = output_pdc.kdhap_down;%
% vm_down = output_pdc.vm_down;%
% 
% keq_fba = output_pdc.keq_fba;%
% keq_tpi = output_pdc.keq_tpi;%
% keq_gpd = output_pdc.keq_gpd;%


