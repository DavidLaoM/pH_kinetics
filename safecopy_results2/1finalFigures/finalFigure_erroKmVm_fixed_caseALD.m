% % finalFigure_errorKmVm_fixed_caseALD.m
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
setup.weightDataEsp = ones(1,numpHtested);
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 1E-5; % by now just testing


%% (2a) Estimation vm variable

% Km fixed
optfun = @costfun_Kmfixed;
plength = 13; % Kms (4) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
lb = -1*ones(1,plength);
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


%% (2b) Estimation km variable. Quick comparison
% first value
% ntests = 12;
ntests = 1;
error_varKm = cell(ntests,1);
for testVal = 1:ntests
% % % % numtests = [1,12];
% % % % error_varKm = cell(length(numtests),1);
% % % % for testVal = numtests
    % Km fixed
    setup.fixedValues = data.Vmax(testVal,1:4);
    optfun = @costfun_Vmfixed;
    plength = 37; % Kms (4) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 1*ones(1,plength);
    lb = -1*ones(1,plength);
%     options = optimset('Display','iter');
    options = optimset('Display','iter','MaxIter',20);

    % Estimation
    setup.plotEachSimCF = 0;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0;
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;

    % Error calculation
%     setup.startVm = 1;
    setup.plotEachSimCF = 1;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 1;
    [errorKm] = optfun(xres,data,setup);
    error_varKm{testVal} = errorKm;
end
% %%
errorvm4plot = sum(abs(error_varVm));
errorkm4plot = zeros(1,ntests);
for i = 1:numtests
% errorkm4plot = zeros(1,length(numtests));
% for i = 1:length(numtests)
    errorkm4plot(i) = sum(abs(error_varKm{i}));
end
xlims = [6 8];

figure(301)
line(xlims,[errorvm4plot errorvm4plot],'Color','blue','LineStyle','--')
hold on
% plot(setup.fullpHarray,errorkm4plot,'.-')
plot([setup.fullpHarray(1) setup.fullpHarray(end)],errorkm4plot,'.-')
legend('Vmax variable','Km variable')


% % % % %% (2c) Estimation km variable. Quick comparison
% % % % ntests = 50;
% % % % testValues = linspace(-3,1,ntests);
% % % % data.Vmax2pre = data.Vmax(12,1:4);
% % % % data.Vmax2 = zeros(length(testValues),4);
% % % % for i = 1:ntests
% % % %     data.Vmax2(i,:) = data.Vmax2pre .* 10 .^ testValues(i);
% % % % end
% % % % disp(data.Vmax2);
% % % % 
% % % % % first value
% % % % error_varKm = cell(ntests,1);
% % % % xres_varKm = cell(ntests,1);
% % % % for testVal = 1:ntests
% % % %     % Km fixed
% % % %     setup.fixedValues = data.Vmax2(testVal,:);
% % % %     optfun = @costfun_Vmfixed;
% % % %     plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
% % % %     x_temp = zeros(1,plength);
% % % %     ub = 3*ones(1,plength);
% % % %     lb = -3*ones(1,plength);
% % % % %     options = optimset('Display','iter');
% % % %     options = optimset('Display','iter','MaxIter',15);
% % % % 
% % % %     % Estimation
% % % %     setup.plotEachSimCF = 0;
% % % %     setup.plotEachSim = 0;
% % % %     setup.simAllProfiles = 0;
% % % %     tic
% % % %     [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% % % %     t = toc;
% % % %     xres_varKm{testVal} = xres;
% % % % 
% % % %     % Error calculation
% % % %     setup.startVm = 1;
% % % %     setup.plotEachSimCF = 1;
% % % %     setup.plotEachSim = 0;
% % % %     setup.simAllProfiles = 1;
% % % %     [errorKm] = optfun(xres,data,setup);
% % % %     error_varKm{testVal} = errorKm;
% % % % end
% % % % %%
% % % % errorvm4plot = sum(abs(error_varVm));
% % % % errorkm4plot = zeros(1,ntests);
% % % % for i = 1:ntests
% % % %     errorkm4plot(i) = sum(abs(error_varKm{i}));
% % % % end
% % % % % xlims = [6 8];
% % % % % %% test vales 6 and 7
% % % % % first value
% % % % testsNums = 6:7;
% % % % 
% % % % error_varKm_test67 = cell(length(testsNums),1);
% % % % xres_varKm_test67 = cell(length(testsNums),1);
% % % % simPEP_full = cell(length(testsNums),1);
% % % % simENO_full = cell(length(testsNums),1);
% % % % expPEP_full = cell(length(testsNums),1);
% % % % expENO_full = cell(length(testsNums),1);
% % % % 
% % % % for testVal = testsNums
% % % %     % Km fixed
% % % %     setup.fixedValues = data.Vmax2(testVal,:);
% % % %     optfun = @costfun_Vmfixed;
% % % %     plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
% % % %     x_temp = zeros(1,plength);
% % % %     ub = 3*ones(1,plength);
% % % %     lb = -3*ones(1,plength);
% % % % %     options = optimset('Display','iter');
% % % %     options = optimset('Display','iter','MaxIter',15);
% % % % 
% % % %     % Estimation
% % % %     setup.plotEachSimCF = 0;
% % % %     setup.plotEachSim = 0;
% % % %     setup.simAllProfiles = 0;
% % % %     tic
% % % %     [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% % % %     t = toc;
% % % %     xres_varKm_test67{testVal-5} = xres;
% % % % 
% % % %     % Error calculation
% % % %     setup.startVm = 1;
% % % %     setup.plotEachSimCF = 1;
% % % %     setup.plotEachSim = 0;
% % % %     setup.simAllProfiles = 1;
% % % % %     [errorKm] = optfun(xres,data,setup);
% % % %     [errorKm,simPEP,simENO,expPEP,expENO] = costfun_Vmfixed_simsENO(xres,data,setup);
% % % %     error_varKm{testVal-5} = errorKm;
% % % %     simPEP_full{testVal-5} = simPEP;
% % % %     simENO_full{testVal-5} = simENO;
% % % %     expPEP_full{testVal-5} = expPEP;
% % % %     expENO_full{testVal-5} = expENO;
% % % % end
% % % % 
% % % % 
% % % % %%
% % % % figure(402)
% % % % % line([0 0.1],[errorvm4plot errorvm4plot],'Color','blue','LineStyle','--')
% % % % % hold on
% % % % subplot(1,3,[1 2])
% % % % loglog(data.Vmax2(:,end)',errorkm4plot,'r.-')
% % % % hold on
% % % % line([data.Vmax2(1,end) data.Vmax2(end,end)],[errorvm4plot errorvm4plot],'Color','blue','LineStyle','--')
% % % % legend('Vmax variable','Km variable','location','south','orientation','horizontal')
% % % % xlabel('Sample starting vm value')
% % % % ylabel('error_{simulations - experimental data}')
% % % % ylim([1E0 1E4])
% % % % xlim([1E-6 1E-2])
% % % % 
% % % % % case low pH values (1)
% % % % subplot(2,3,3)
% % % % for j = 1
% % % %     for i = setup.DFstudy
% % % %         plot(data.time{j,i}, simPEP_full{1}{i,j},'-','LineWidth',2)
% % % %         hold on
% % % %         plot(data.time{j,i}, expPEP_full{1}{i,j},'k.','MarkerSize',4)
% % % %         ylim([0 1.2])
% % % %     end
% % % % end
% % % % yaxisname = erase(sprintf('PEP [mM] concentration fit @pH%d',setup.fullpHarray(j)),"0000e+00");
% % % % ylabel(yaxisname)
% % % % xlabel('time [s]')
% % % % 
% % % % % case high pH values
% % % % subplot(2,3,6)
% % % % for j = 12
% % % %     for i = setup.DFstudy
% % % %         plot(data.time{j,i}, simPEP_full{1}{i,j},'-','LineWidth',2)
% % % %         hold on
% % % %         plot(data.time{j,i}, expPEP_full{1}{i,j},'k.','MarkerSize',4)
% % % %         ylim([0 1.2])
% % % %     end
% % % % end
% % % % yaxisname = erase(sprintf('PEP [mM] concentration fit @pH%d',setup.fullpHarray(j)),"0000e+00");
% % % % ylabel(yaxisname)
% % % % xlabel('time [s]')



