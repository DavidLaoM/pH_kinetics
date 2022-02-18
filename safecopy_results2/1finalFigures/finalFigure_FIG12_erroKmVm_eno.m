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


%% (2a) Estimation vm variable

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
xres_varVm = xres;

% Error calculation
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
[error_varVm] = optfun(xres,data,setup);
sum_error_varVm = sum(abs(error_varVm));

%% (2b) Estimation km variable

% ntests = 12;
% error_varKm = cell(ntests,1);
% sum_error_varKm = zeros(ntests,1);
% xres_varKm = cell(ntests,1);
% for testVal = 1:ntests
%     % Km fixed
%     setup.fixedValues = data.Vmax(testVal,1:4);
%     optfun = @costfun_Vmfixed;
%     plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
%     x_temp = zeros(1,plength);
%     ub = 3*ones(1,plength);
%     lb = -3*ones(1,plength);
%     options = optimset('Display','iter');
% %     options = optimset('Display','iter','MaxIter',15);
% 
%     % Estimation
%     setup.plotEachSimCF = 0;
%     setup.plotEachSim = 0;
%     setup.simAllProfiles = 0;
%     tic
%     [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
%     t = toc;
% 
%     % Error calculation
% %     setup.startVm = 1;
% %     setup.plotEachSimCF = 1;
% %     setup.plotEachSim = 0;
% %     setup.simAllProfiles = 1;
%     [errorKm] = optfun(xres,data,setup);
%     error_varKm{testVal} = errorKm;
%     xres_varKm{testVal} = xres;
%     sum_error_varKm(testVal) = sum(abs(errorKm));
% end

ntests = 50;
testValues = linspace(-3,1,ntests);
data.Vmax2pre = data.Vmax(12,1:4);
data.Vmax2 = zeros(length(testValues),4);
for i = 1:ntests
    data.Vmax2(i,:) = data.Vmax2pre .* 10 .^ testValues(i);
end
% disp(data.Vmax2);

% first value
error_varKm = cell(ntests,1);
xres_varKm = cell(ntests,1);
sum_error_varKm = zeros(ntests,1);
for testVal = 1:ntests
    % Km fixed
    setup.fixedValues = data.Vmax2(testVal,:);
    optfun = @costfun_Vmfixed;
    plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
%     options = optimset('Display','iter');
    options = optimset('Display','iter','MaxIter',15);

    % Estimation
    setup.plotEachSimCF = 0;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0;
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;

    % Error calculation
    setup.startVm = 1;
    setup.plotEachSimCF = 1;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 1;
    [errorKm] = optfun(xres,data,setup);
    error_varKm{testVal} = errorKm;
    xres_varKm{testVal} = xres;
    sum_error_varKm(testVal) = sum(abs(errorKm));
end


%% first figure: error plot
x2 = sum_error_varKm;
x1 = sum_error_varVm * ones(size(sum_error_varKm));

figure

subplot(1,2,1) % boxplot
% boxplot([x1,x2],'Notch','on','Labels',{'vm variable, km fixed','vm fixed, km variable'})
boxplot([x1,x2],'Notch','on','Labels',{'vm variable km fixed','vm fixed, km variable'},...
    'PlotStyle','compact','LabelOrientation','horizontal')
xlabel('cases tested')
ylabel('error value')
ylim([0 250])

subplot(1,2,2)
loglog(data.Vmax2(:,end)',sum_error_varKm,'r.-')
hold on
line([data.Vmax2(1,end) data.Vmax2(end,end)],[sum_error_varVm sum_error_varVm],'Color','blue','LineStyle','--')
legend('Vmax variable','Km variable','location','south','orientation','horizontal')
xlabel('Sample starting vm value')
ylabel('error_{simulations - experimental data}')
ylim([1E0 1E4])
xlim([1E-6 1E-2])
% ylim([0 250])

% %%
% rng default  % For reproducibility
% x1 = normrnd(5,1,100,1);
% x2 = normrnd(6,1,100,1);
% figure
% boxplot([x1,x2],'Notch','on','Labels',{'mu = 5','mu = 6'})
% title('Compare Random Data from Different Distributions')


%% second figure: simulation vm variable, km constant
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
costfun_Kmfixed(xres_varVm,data,setup);


%% third figure: simulation vm constant, km variable (minimim error)

% detect the minimum
[~,idx_min] = min(sum_error_varKm);

% plot
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
setup.fixedValues = data.Vmax2(idx_min,:);
costfun_Vmfixed(xres_varKm{idx_min},data,setup);


%% third figure: simulation vm constant, km variable (average error)
% For more visibility

% plot
numInter = 15;
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
setup.fixedValues = data.Vmax2(numInter,:);
costfun_Vmfixed(xres_varKm{numInter},data,setup);


%% save resulting variables and plots
% to develop
% disp('Stop here, check and if good is reached, save.')

%%
% % save matrix
% save('20200810_finalFigure12_data2Run');

% % Error_Box_Plot
% set(1,'color','w');
% savefig(1,'20200810_finalFig12_errorBoxPlot.fig')

% % Vm_varaible
% set(101,'color','w');
% savefig(101,'20200810_finalFig12_PEPconcentration_Vm_varaible.fig')
% set(102,'color','w');
% savefig(102,'20200810_finalFig12_ENOrates_Vm_varaible.fig')

% % Vm_fixed_best_fit
% set(101,'color','w');
% savefig(101,'20200810_finalFig12_PEPconcentration_Vm_fixed_best_fit.fig')
% set(201,'color','w');
% savefig(201,'20200810_finalFig12_ENOrates_Vm_fixed_best_fit.fig')

% % Vm_fixed_intermediate_fit
% set(101,'color','w');
% savefig(101,'20200810_finalFig12_PEPconcentration_Vm_fixed_intermediate_fit.fig')
% set(201,'color','w');
% savefig(201,'20200810_finalFig12_ENOrates_Vm_fixed_intermediate_fit.fig')



