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

% Km fixed
optfun = @costfun_Kmfixed;
plength = 13; % Kms (3) + Vms (1) * numpH (8) + vm linking (2) + (Keq is fixed to experimental data)
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
for i = 1:length(setup.Keq_FBA)
    setup.Keq_FBA(i) = setup.Keq_FBA(6);
    setup.Keq_TPI(i) = setup.Keq_TPI(6);
    setup.Keq_GPD(i) = setup.Keq_GPD(6);
end


%%
% parameter estimation
tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t = toc;

% display result
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

subplot(331) % vm
%     plot(pHarray,vm_up,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm_down,'.-','color',[0.5 0.5 0.5]), hold on, 
%     plot(pHarray,vm,'.-','color','black')
% title('v_{m} [mM]')
    plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_uChange,'.-','color','black')
title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

subplot(332) % vm(log_difference)
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


% %%
% set(105,'color','w');
% savefig(105,'ALD_pVals_pCis');
% %%
output_ald.xres_selected = xres_selected;

output_ald.pHarray = pHarray;

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

% save('output_ald.mat','output_ald');


%% reverse to reload
% load('output_ald.mat','output_ald');

xres_selected = output_ald.xres_selected;

pHarray = output_ald.pHarray;

kfbp = output_ald.kfbp;%
kgap = output_ald.kgap;%
kdhap = output_ald.kdhap;%
vm = output_ald.vm;%

kfbp_up = output_ald.kfbp_up;%
kgap_up = output_ald.kgap_up;%
kdhap_up = output_ald.kdhap_up;%
vm_up = output_ald.vm_up;%

kfbp_down = output_ald.kfbp_down;%
kgap_down = output_ald.kgap_down;%
kdhap_down = output_ald.kdhap_down;%
vm_down = output_ald.vm_down;%

keq_fba = output_ald.keq_fba;%
keq_tpi = output_ald.keq_tpi;%
keq_gpd = output_ald.keq_gpd;%


%%
output_ald_keq_constant = output_ald;
save('results/ald/output_ald_keq_constant.mat','output_ald_keq_constant');

