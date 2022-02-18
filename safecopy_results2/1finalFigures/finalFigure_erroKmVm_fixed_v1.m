% % finalFigure_errorKmVm_fixed.m
% 1. Load data
% 2. Estimation vm variable
% 3. Estimation km variable
% 4. Plotting


%% (1) Setup and data load
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

% parameter estimation setup
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


%% (2) Estimation vm variable

% Km fixed
optfun = @costfun_Kmfixed;
plength = 14; % Kms (2) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');

% Estimation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t = toc;

% Error calculation
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
[error_varVm] = optfun(xres,data,setup);


%% (3) Estimation km variable

% Vm fixed
optfun = @costfun_Vmfixed;
setup.startVm = 12;
plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');

% Estimation
sourceVal = data.Vmax(1,4);
% numbersWanted = [2E-5 1E-4 2E-4 1E-3 2E-3 1E-2 2E-2];
% numbersWanted = data.Vmax(:,4)';
numbersWanted = [data.Vmax(1,4)*0.5, data.Vmax(1,4)*0.75, data.Vmax(:,4)', data.Vmax(12,4)*1.5, data.Vmax(12,4)*2, data.Vmax(12,4)*3, data.Vmax(12,4)*4, data.Vmax(12,4)*5,  data.Vmax(12,4)*6,  data.Vmax(12,4)*7,  data.Vmax(12,4)*8,  data.Vmax(12,4)*9, data.Vmax(12,4)*10];
factMultiply = numbersWanted ./ sourceVal;
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
xres_cell = cell(length(factMultiply),1);
% loop to estimate each pH value
for i = 1:length(factMultiply)
    % setup
%     setup.startVm = i;
    setup.startVm = 1;
%     setup.factMultiply = 1; %2;
    setup.factMultiply = factMultiply(i); %2;
    % parameter estimation
    tic
    [xres_temp,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;
    xres_cell{i} = xres_temp;
end

%% Error calculation
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
error_varKm = zeros(length(factMultiply),1);
for i = 1:length(factMultiply)
    xres_temp = xres_cell{i};
    [error] = optfun(xres_temp,data,setup);
    error_varKm(i) = sum(abs(error));
end


%%
save('tempData_finalFigure_errorKmVm_fixed.mat')


%% 
load('tempData_finalFigure_errorKmVm_fixed.mat')
%%
% selEst = 15;
selEst = 17; % selected estimate
% selEst = 18;
%%
simTime681 = cell(4,1);
simMet681 = cell(4,1);
simTime790 = cell(4,1);
simMet790 = cell(4,1);

DFstudy = setup.DFstudy;
obsMet = setup.observableMetabolite;

for sim6_81 = 6
    j = sim6_81;
    
    data.Vmaxs = data.Vmax(1,:).*factMultiply(selEst); % In case of Vmfixed, we have only one value for Vm here (took the 6th data point)
    data.PEP = data.conc_mean(j,:);
    data.Vprofs = data.RRs(j,:); 
    data.tempTime = data.time(j,:);
    % inputs to be selected
    data.chosenKeq = setup.keq(j);              
    data.i = j;
    % selecting the right parameters
    xassay = xres_cell{selEst}([7,19,1]);
%     xassay = xres_cell{selEst};

    % simulations
    for i = DFstudy
        % recall vmax for the specific value and simulate
        data.chosenDF = data.DF(j,i);
        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
        data.chosenPEPini = data.PEP{i}(1);
        % simulate metabolites
        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
        % cost function (NADHexp + vGAPDHr)
        simTime681{i} = simResult.t;
        simMet681{i} = simResult.y(:,obsMet);
    end
end

for sim7_90 = 12
    j = sim7_90;
    
    data.Vmaxs = data.Vmax(1,:).*factMultiply(selEst); % In case of Vmfixed, we have only one value for Vm here (took the 6th data point)
    data.PEP = data.conc_mean(j,:);
    data.Vprofs = data.RRs(j,:); 
    data.tempTime = data.time(j,:);
    % inputs to be selected
    data.chosenKeq = setup.keq(j);              
    data.i = j;
    % selecting the right parameters
    xassay = xres_cell{selEst}([13,25,1]);

    % simulations
    for i = DFstudy
        % recall vmax for the specific value and simulate
        data.chosenDF = data.DF(j,i);
        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
        data.chosenPEPini = data.PEP{i}(1);
        % simulate metabolites
        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
        % cost function (NADHexp + vGAPDHr)
        simTime790{i} = simResult.t;
        simMet790{i} = simResult.y(:,obsMet);
    end
end

%% get vm variable data
simTime681_vmvar = cell(4,1);
simMet681_vmvar = cell(4,1);
simTime790_vmvar = cell(4,1);
simMet790_vmvar = cell(4,1);

% selEst = 5; % selected estimate
DFstudy = setup.DFstudy;
obsMet = setup.observableMetabolite;

for sim6_81_vmvar = 6
    j = sim6_81_vmvar;
    
%     data.Vmaxs = data.Vmax(1,:).*factMultiply(selEst); % In case of Vmfixed, we have only one value for Vm here (took the 6th data point)
    data.Vmaxs = data.Vmax(j,:);
    data.PEP = data.conc_mean(j,:);
    data.Vprofs = data.RRs(j,:); 
    data.tempTime = data.time(j,:);
    % inputs to be selected
    data.chosenKeq = setup.keq(j);              
    data.i = j;
    % selecting the right parameters
%     xassay = xres_cell{selEst}([7,19,1]);
    xassay = xres([1,2,8]);
%     xassay = xres_cell{selEst};

    % simulations
    for i = DFstudy
        % recall vmax for the specific value and simulate
        data.chosenDF = data.DF(j,i);
        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
        data.chosenPEPini = data.PEP{i}(1);
        % simulate metabolites
        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
        % cost function (NADHexp + vGAPDHr)
        simTime681_vmvar{i} = simResult.t;
        simMet681_vmvar{i} = simResult.y(:,obsMet);
    end
end

for sim7_90_vmvar = 12
    j = sim7_90_vmvar;
    
%     data.Vmaxs = data.Vmax(1,:).*factMultiply(selEst); % In case of Vmfixed, we have only one value for Vm here (took the 6th data point)
    data.Vmaxs = data.Vmax(j,:);
    data.PEP = data.conc_mean(j,:);
    data.Vprofs = data.RRs(j,:); 
    data.tempTime = data.time(j,:);
    % inputs to be selected
    data.chosenKeq = setup.keq(j);              
    data.i = j;
    % selecting the right parameters
    xassay = xres([1,2,14]);

    % simulations
    for i = DFstudy
        % recall vmax for the specific value and simulate
        data.chosenDF = data.DF(j,i);
        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
        data.chosenPEPini = data.PEP{i}(1);
        % simulate metabolites
        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
        % cost function (NADHexp + vGAPDHr)
        simTime790_vmvar{i} = simResult.t;
        simMet790_vmvar{i} = simResult.y(:,obsMet);
    end
end

% % % % data.Vmaxs = data.Vmax(j,:);
% % % % 
% % % % % Error calculation
% % % % setup.plotEachSimCF = 1;
% % % % setup.plotEachSim = 0;
% % % % setup.simAllProfiles = 1;
% % % % [error_varVm] = optfun(xres,data,setup);

%%
figure(102)

% left-side error plot
subplot(2,3,[1 4])
semilogy(numbersWanted, error_varKm(1:length(factMultiply)), 'r.-')
hold on
line([numbersWanted(1) numbersWanted(end)],[sum(abs(error_varVm)) sum(abs(error_varVm))],'color','blue','LineStyle','--')
xlabel('starting vm value in optimization')
ylabel('error_{Data} = sim. - exp.')
legend('error_{Km variable}','error_{Vm variable}','Location','north','Orientation','horizontal')
title('error comparison')
xlim([numbersWanted(1) numbersWanted(end)])
% ylim([1E0 1E4])

% centre plot: vm variable
subplot(2,3,2)
data_PEP = data.conc_mean(6,:);
data_tempTime = data.time(6,:);
for i = 1:4
    plot(simTime681_vmvar{i}, simMet681_vmvar{i},'-','LineWidth',2)
    hold on
    plot(data_tempTime{i}, data_PEP{i},'k.','MarkerSize',4)
end
ylabel('pH6.81. PEP concentration [mM] vs assay time [s]')
xlabel('time [s]')
title('vm variable')
ylim([0.3 0.9])
xlim([0 300])

subplot(2,3,5)
data_PEP = data.conc_mean(12,:);
data_tempTime = data.time(12,:);
for i = 1:4
    plot(simTime790_vmvar{i}, simMet790_vmvar{i},'-','LineWidth',2)
    hold on
    plot(data_tempTime{i}, data_PEP{i},'k.','MarkerSize',4)
end
ylabel('pH7.90. PEP concentration [mM] vs assay time [s]')
xlabel('time [s]')
ylim([0.3 0.9])
xlim([0 300])


% right-side simulation plots: km variable
subplot(2,3,3)
data_PEP = data.conc_mean(6,:);
data_tempTime = data.time(6,:);
for i = 1:4
%     plot(data.tempTime{i}, simPEP{i,6},'-','LineWidth',2)
    plot(simTime681{i}, simMet681{i},'-','LineWidth',2)
    hold on
end
for i = 1:4
    plot(data_tempTime{i}, data_PEP{i},'k.','MarkerSize',4)
    hold on
end
ylabel('pH6.81. PEP concentration [mM] vs assay time [s]')
xlabel('time [s]')
title('km variable (minimal error found)')
legend('DF8','DF4','DF2','DF1')
ylim([0.3 0.9])
xlim([0 300])

subplot(2,3,6)
data_PEP = data.conc_mean(12,:);
data_tempTime = data.time(12,:);
for i = 1:4
    plot(simTime790{i}, simMet790{i},'-','LineWidth',2)
    hold on
    plot(data_tempTime{i}, data_PEP{i},'k.','MarkerSize',4)
end
ylabel('pH7.90. PEP concentration [mM] vs assay time [s]')
xlabel('time [s]')
ylim([0.3 0.9])
xlim([0 300])

set(gcf,'color','w')

% % % % %%
% % % % setup.factMultiply = factMultiply(selEst);
% % % % setup.startVm = 1;
% % % % 
% % % % 
% % % % setup.plotEachSimCF = 1;
% % % % setup.plotEachSim = 0;
% % % % setup.simAllProfiles = 1;
% % % %     xres_temp = xres_cell{17};
% % % %     [error] = optfun(xres_temp,data,setup);
% % % % setup.plotEachSimCF = 0;
% % % % setup.plotEachSim = 0;
% % % % setup.simAllProfiles = 0;
