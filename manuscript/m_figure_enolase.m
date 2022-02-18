% % M_FIGURE_ENOLASE.M

if setup.fast_option == 0

    % % (1a) initial reload
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
    end

    % % (1b) units conversion
    minwindow = 60; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.6];  %Ylims plot conc
    dp_start = 6 * ones(size(data.conc_mean));
    total_len = zeros(size(dp_start));
    DFarray = [1/8 1/4 1/2 1/1];
    idxs2consider = ones(size(DF));
    
    setup.plotOutput = 0; % % 
    expRatesDetermination;

    % % (1c) pEst setup
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

    % % (2a) Estimation vm and km variable
    % Vm and Km variable
    optfun = @costfun_VmKmvariable;
    plength = 36; % (Kms (2) + Vms (1)) * numpH (12) + (Keq is fixed to experimental data)
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
    xres_varVmKm = xres;
    % Error calculation
    setup.plotEachSimCF = 1;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 1;
    [error_varVmKm] = optfun(xres,data,setup);
    sum_error_varVmKm = sum(abs(error_varVmKm));
    % 
    res_VmKm_pHdependent.xres_varVmKm = xres;%: [1�36 double]
    res_VmKm_pHdependent.resnorm = resnorm;%: 0.0717
    res_VmKm_pHdependent.residual = residual;%: [5811�1 double]
    res_VmKm_pHdependent.Jacobian = Jacobian;%: [5811�36 double]
    res_VmKm_pHdependent.error_varVmKm = error_varVmKm;%: [5811�1 double]
    res_VmKm_pHdependent.sum_error_varVmKm = sum_error_varVmKm;%: 14.8455        
    
    % % (2b) Estimation vm variable
    % Vmax variable
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
    % 
    res_Vm_pHdependent.xres_varVm = xres;%:
    res_Vm_pHdependent.resnorm = resnorm;%:
    res_Vm_pHdependent.residual = residual;%: 
    res_Vm_pHdependent.Jacobian = Jacobian;%:
    res_Vm_pHdependent.error_varVm = error_varVm;%: 
    res_Vm_pHdependent.sum_error_varVm = sum_error_varVm;%:   
    
    % % (2c) Estimation km variable
    % Km variable
    optfun = @costfun_Vmfixed;
    plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
    options = optimset('Display','iter','MaxIter',15);
    % create test arrays
    ntests = 50;
    testValues = linspace(-3,1,ntests);
    data.Vmax2pre = data.Vmax(12,1:4);
    data.Vmax2 = zeros(length(testValues),4);
    for i = 1:ntests
        data.Vmax2(i,:) = data.Vmax2pre .* 10 .^ testValues(i);
    end
    % first value
    error_varKm = cell(ntests,1);
    xres_varKm = cell(ntests,1);
    sum_error_varKm = zeros(ntests,1);
    resnorm_varKm = cell(ntests,1);
    residual_varKm = cell(ntests,1);
    Jacobian_varKm = cell(ntests,1);
    for testVal = 1:ntests
        % Vmax
        setup.fixedValues = data.Vmax2(testVal,:); % this could actually have been set inside the 'x_temp' input in the lsqnonlin call just below... ideally do in that one when time available... If so, thne remember to make the test as 'demanding' as the current one by increasing the boundaries for the Vmax parameter (x_temp(1)).
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
        % saving
        error_varKm{testVal} = errorKm;
        xres_varKm{testVal} = xres;
        sum_error_varKm(testVal) = sum(abs(errorKm));
        resnorm_varKm{testVal} = resnorm_varKm;
        residual_varKm{testVal} = residual_varKm;
        Jacobian_varKm{testVal} = Jacobian_varKm;
    end

    % saving in a bigger matrix
    res_Km_pHdependent.xres_varKm = xres_varKm;
    res_Km_pHdependent.resnorm = resnorm_varKm;
    res_Km_pHdependent.residual = residual_varKm;
    res_Km_pHdependent.Jacobian = Jacobian_varKm;
    res_Km_pHdependent.error_varKm = error_varKm;
    res_Km_pHdependent.sum_error_varKm = sum_error_varKm;

elseif setup.fast_option == 1 % directly recall the workspace from here
    load('workspace_enolase.mat');
end

%% PLOTTING
% %% get the simulations
tic
% %% first figure: simulation vm variable, km variable
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
[error_var_vmkm,simResult_var_vmkm] = simRes_costfun_VmKmvariable(xres_varVmKm,data,setup);

% %% second figure: simulation vm variable, km constant
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
[error_var_vm,simResult_var_vm] = simRes_costfun_Kmfixed(xres_varVm,data,setup);

% %% third figure: simulation vm constant, km variable (average error)
numInter = 15;
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
setup.fixedValues = data.Vmax2(numInter,:);
[error_var_km,simResult_var_km] = simRes_costfun_Vmfixed(xres_varKm{numInter},data,setup);


% %%
figure(102),
    
sp2 = subplot(122); 
    plot(res_Km_pHdependent.sum_error_varKm(1:40),'.-','Color',[1 0.5 0.5]),
    xlabel('sample number'), ylabel('error'),
    title({'focused view';'(first 40 sample)'}),
    hold on, line(sp2.XLim,ones(1,2)*mean(res_VmKm_pHdependent.sum_error_varVmKm),'LineStyle','--','color',[0.9 0.9 0.9])
    hold on, line(sp2.XLim,ones(1,2)*mean(res_Vm_pHdependent.sum_error_varVm),'LineStyle','--','color',[0.5 0.5 1])
    
subplot(121), 
    area(sp2.XLim,ones(1,2)*sp2.YLim(2),'FaceColor',[0.95 0.95 0.95])
    hold on
    plot(res_Km_pHdependent.sum_error_varKm,'.-','Color',[1 0.5 0.5]),
    xlabel('sample number'), ylabel('error'),
    title('general view'),
    
set(gcf,'color','w');

%% create the figure

% specific color palette
c_royalBlue = [65	105	225]/255; % royalblue
c_midnightblue = [25	25	112]/255; % midnightblue
c_CCCCCC = [204	204	204]/255; % #CCCCCC
c_E5E5E5 = [229 229 229]/255; % #E5E5E5
c_0f1076 = [15	16	118]/255; % #0f1076
c_chocolate = [210	105	30]/255; % (#e59400 temp orange)
% 
resizeVal = 2/3;

% 
x1 = res_VmKm_pHdependent.sum_error_varVmKm;
x2 = res_Vm_pHdependent.sum_error_varVm;
x3 = res_Km_pHdependent.sum_error_varKm;

selValues = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
[idxs,~] = find(simResult_var_vm.data.tempTime{1} == selValues);
idxs2 = [idxs; 21; 51; 101];
idxs3 = sort(idxs2);
changed_j = 9:20;

% clf(101)
if exist('h101','var')
    clf(101)
end
h101 = figure(101);
dim = [.04 .02 .92 .59];
str = ' ';
an_h = annotation('textbox',dim,'String',str,...
    'LineWidth',1);

% % subplot(3,6,[1,2,7,8]) % boxplot
sp1_h = subplot(5,4,[2,3,6,7]); % boxplot
hb = bar([1 2],[mean(x2) mean(x3(1:40))]);
hb.FaceColor = 'flat';
hb.CData = [c_midnightblue; c_chocolate];
hb.EdgeColor = 'none';
set(gca, 'FontSize',14*resizeVal,'XTick',[1 2],'XTickLabel',{'V_{max}','K_{m}'});
fix_xticklabels(gca,0.1,{'FontSize',14*resizeVal});
ylabel({'Error (mM)';'simulation - experimental'}, 'FontSize',18*resizeVal)%12)
errbar = [std(x2); std(x3(1:40))];   % CREATE �errbar� MATRIX
yd = [mean(x2) mean(x3(1:40))]';
hold on
for k1 = 1:2
    errorbar(k1,  yd(k1,:),  errbar(k1,:), '.k', 'LineWidth',2)
end
ylim([0 200])
hold on
text(0, sp1_h.YLim(2)*0.9, 'A', 'FontSize', 20)
text(-2.7, sp1_h.YLim(2)*-0.3, 'B', 'FontSize', 20)

% 
set(0,'CurrentFigure', h101)
for j = 1:12
    sp2_h(j) = subplot(5,4,changed_j(j));   
    for i = simResult_var_vm.DFstudy
 
        x3 = [simResult_var_vm.data.tempTime{i}(idxs3)', fliplr(simResult_var_vm.data.tempTime{i}(idxs3)')];
            inBetween = [simResult_var_vm.simPEP{i,j}(idxs3)', fliplr(simResult_var_km.simPEP{i,j}(idxs3)')];
            f1 = fill(x3, inBetween, 'g');
            f1.FaceColor = [228	165	120]/255;
            f1.EdgeColor = 'none';
            hold on
        plot(simResult_var_vm.data.tempTime{i}(idxs3), simResult_var_vm.simPEP{i,j}(idxs3),'-','LineWidth',1.5,'color', hb.CData(1,:))
        hold on
        plot(simResult_var_vm.data.tempTime{i}(idxs3), simResult_var_km.simPEP{i,j}(idxs3),'-','LineWidth',1.5,'color', hb.CData(2,:))
        hold on
        csp = plot(simResult_var_vm.data.tempTime{i}(idxs3), simResult_var_vm.expPEP{i,j}(idxs3),'k.','MarkerSize',10);
        hold on
        ylim([0 1.2])
        tempText = erase(sprintf('pH %d', setup.fullpHarray(j)),"0000e+00");
        text(30, 1.2*0.9, tempText,'FontSize',14 * resizeVal)
    end
    if j == 9
        xlab = xlabel('assay time (s)','FontSize',18 * resizeVal);
        xlab.Position = xlab.Position + [0 0 0];
        ylab = ylabel('PEP concentration (mM)','FontSize',18 * resizeVal);
        ylab.Position = ylab.Position + [0 1.5 0];
    end
    ax = gca;
    ax.FontSize = 14 * resizeVal;
    xlim([0 600])
    ax.XTick = [0 300 600];
end
hold off
set(gcf,'color','w');
set(101, 'Position', [100 100 900 900])

% edits for the boxes and labels A and B
sp1_h.Position = sp1_h.Position + [0 0.025 0 0]; % bit up with A
for j = 1:12
    sp2_h(j).Position = sp2_h(j).Position + [0 0 0 0];
end
% box
hold off

%% latest additions
sp_temp = h101.Children(2);
hL1 = legend(sp_temp.Children([4 3 2]),'pH-dependent V_{max}',...
    'pH-dependent K_{m}','experimental data');
hL1.Orientation = 'horizontal';
hL1.Box = 'off';
hL1.FontSize = 10;
hL1.Position = [0.3    0.02    0.3960    0.0340];

%% comments 2021 11 15
% delete(temp)
dim = [.365 .338 .3 .3];
str = 'pH dependency';
temp = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'EdgeColor','none','HorizontalAlignment','center');

%% needed stop not to get truncated... #matlabUselessSecrets
% save
if setup.saveOutput == 1
    savefig(101,'pHmanus_enolase')
    % specs printing (method 3)
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    % print -dpdf -painters panel_11
    print -dpdf -painters enolase
end

