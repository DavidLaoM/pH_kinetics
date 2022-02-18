% % DoE_ENO.m
% design of experiments for the enolase enzyme.
% The idea of this code is to find out the conditions in which parameters
% Km may also become identifiable.

% Structure
% 1. Recall Eno estimation results
% 2. Test PSA for changing Km: (+/-) 2 orders of magnitude (up/down)
% 3. Test PSA for change in concentrations: (s/p) substrate/prooduct
% 4. Create mock dataset
% 5. Estimate back
% 6. Other enzymes
    % ENO highly recommended. Cannot be missed. Because:
        % a. Little parameters
        % b. Directly measured
        % c. Keq does not change.
    % GAPDH enzyme of interest
        % a. One reaction is directly measures
        % -. Complicated topology
        % -. Keq change also to take into account
        % -. Apparently difficulty in Vmax measurement
        % +. It could actually help visualize Vmax
        
        
%% 1. Recall Eno estimation results
clear, close all
set_paths_pHstudy;
dbstop if error
for section1 = 1
    
    % (0) Setup and data load
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
    
    % (0.1) Experimental Vmax determination
    % intial things that could be in the setup
    minwindow = 60; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.6];  %Ylims plot conc
    % select start point (this needs manual selection deppending on previous plots)
    dp_start = 6 * ones(size(data.conc_mean));
    % % % % dp_start = 10 * ones(size(data.conc_mean));
    % blank cell total length
    total_len = zeros(size(dp_start));
    % DFs considered
    DFarray = [1/8 1/4 1/2 1/1];
    % idxs2consider
    idxs2consider = ones(size(DF));
    % Experimental rates determination and plotting
    expRatesDetermination;
    
    % (1.1) Simple parameter fit. Parameter estimation
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
    plength = 14; % Kms (2) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
    options = optimset('Display','iter');
    
    % (1.2) Check that we get the right result parameter set
    close all
    load('eno_parEst.mat')
    xres = output_eno.xres_selected;
    % testing the costfunction
    setup.plotEachSimCF = 1;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 1;
    [error] = optfun(xres,data,setup);
    setup.plotEachSimCF = 0;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0;
end


%% 1.2. Adapted simulation setup
% simulation time extended
% output to select
% plotting right afterwards

% % simulation
% setup.plotEachSimCF = 0; setup.plotEachSim = 0; setup.simAllProfiles = 0;
% nPSAvals = 1;
% simResults = cell(1,nPSAvals);
% for i = 1:nPSAvals
%     [error,simResult] = simRes_costfun_Kmfixed(xres,data,setup);
%     simResults{i} = simResult;
% end
% % visualization
% numpH = setup.numpHtested;
% for j = 1:numpH   
%     simulationVisualization_mainPSA;
% end    


%% 2. Test PSA for changing Km: (+/-) 2 orders of magnitude (up/down)
for tempTest = 1
%     % % first test on a pretty big range
%     % prior results: going too high or too down qill block.
%     x_var = linspace(-5,5,21);
%     tempTickLabels = cell(6,1);
%         tempTickLabels{1} = '-5';
%         tempTickLabels{2} = '-3';
%         tempTickLabels{3} = '-1';
%         tempTickLabels{4} = '1';
%         tempTickLabels{5} = '3';
%         tempTickLabels{6} = '5';
%     n = length(x_var);
%     m = length(xres);
% 
%     x_K2pg = zeros(n,m);
%     x_Kpep = zeros(n,m);
%     for i = 1:n
%         x_K2pg(i,:) = xres;
%         x_K2pg(i,1) = x_K2pg(i,1) + x_var(i);
% 
%         x_Kpep(i,:) = xres;
%         x_Kpep(i,2) = x_K2pg(i,2) + x_var(i);
%     end
% 
%     % % K2pg
%     % simulation
%     setup.plotEachSimCF = 0; setup.plotEachSim = 0; setup.simAllProfiles = 0;
%     nPSAvals = n;
%     simResults = cell(1,nPSAvals);
%     for i = 1:nPSAvals
%         xres = x_K2pg(i,:);
%         [error,simResult] = simRes_costfun_Kmfixed(xres,data,setup);
%         simResults{i} = simResult;
%     end
% 
%     % visualization
%     numpH = setup.numpHtested;
%     for j = 1:numpH   
%         simulationVisualization_mainPSA;
%     end   
end


%% % second test on a focused range
% % % kms to zero
% % xres(1:2) = zeros(1,2);
% % % prior results: going too high or too down qill block.
% % setup.PSAstudy = 1; % if isfield
% % setup.PSAstudy_ENO_k2pg = 1; % if isfield + if ==
% % setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
% % 
% % % % PSA ENO_ini
% % % simulation
% % setup.plotEachSimCF = 0;
% % setup.plotEachSim = 0;
% % setup.simAllProfiles = 0;
% % 
% % nPSAvals = length(setup.PSAvals.ENO_k2pg_ini);
% % simResults = cell(1,nPSAvals);
% % for i = 1:nPSAvals
% %     setup.idx = i;
% %     [error,simResult] = simRes_costfun_Kmfixed_2(xres2,data,setup);
% %     simResults{i} = simResult;
% % end
% % %%
% % % visualization
% % numpH = setup.numpHtested;
% % for j = 1:numpH   
% %     simulationVisualization_mainPSA;
% % end

%% temp
% xres3 = xres; xres3(1:2) = zeros(1,2);
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error,simResult] = simRes_costfun_Kmfixed_2(xres3,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;


%% (2A) initial testing of Laura's protocol concentrations + longer timespan
xres3 = xres; xres3(1:2) = zeros;
% kms to zero

% just the DF 1 is simulated
setup.DFstudy = 4;

% prior results: going too high or too down qill block.
setup.PSAstudy = 1; % if isfield
setup.PSAstudy_ENO_k2pg = 1; % if isfield + if ==
% setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
setup.PSAvals.ENO_k2pg_ini = [6E0 6E-1 5E-2];
setup.supTitleText = 'PSA.ENO.Km2PG';
setup.legendLocation =[0.35 0.925 0.1 0.05];
    A = string(setup.PSAvals.ENO_k2pg_ini);
setup.legNames = cellstr(A);
setup.addedTime = linspace(605,1000,60)';

% % PSA ENO_ini
% simulation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.PSAvals.ENO_k2pg_ini);
simResults = cell(1,nPSAvals);
for i = 1:nPSAvals
    setup.idx = i;
    [error,simResult] = simRes_costfun_Kmfixed_2(xres3,data,setup);
    simResults{i} = simResult;
end

% visualization
numpH = setup.numpHtested;
for j = 1:numpH   
    simulationVisualization_mainPSA;
end


%% (2B) testing the other region in Km (in case the K2pg value changes much
xres4 = xres; xres4(1:2) = [1.1990 0];
% kms to zero

% just the DF 1 is simulated
setup.DFstudy = 4;

% prior results: going too high or too down qill block.
setup.PSAstudy = 1; % if isfield
setup.PSAstudy_ENO_k2pg = 1; % if isfield + if ==
% setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
setup.PSAvals.ENO_k2pg_ini = [6E0 6E-1 5E-2];
setup.supTitleText = 'PSA.ENO.Km2PG';
setup.legendLocation =[0.35 0.925 0.1 0.05];
    A = string(setup.PSAvals.ENO_k2pg_ini);
setup.legNames = cellstr(A);
setup.addedTime = linspace(605,1000,60)';

% % PSA ENO_ini
% simulation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.PSAvals.ENO_k2pg_ini);
simResults = cell(1,nPSAvals);
for i = 1:nPSAvals
    setup.idx = i;
    [error,simResult] = simRes_costfun_Kmfixed_2(xres4,data,setup);
    simResults{i} = simResult;
end

% visualization
numpH = setup.numpHtested;
for j = 1:numpH   
    simulationVisualization_mainPSA;
end


%% (2C) which values
xres3 = xres; xres3(1:2) = zeros; xres3(1:2) = [1.1990 0];
% kms to zero

% just the DF 1 is simulated
setup.DFstudy = 4;

% prior results: going too high or too down qill block.
setup.PSAstudy = 1; % if isfield
setup.PSAstudy_ENO_k2pg = 1; % if isfield + if ==
% setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
setup.PSAvals.ENO_k2pg_ini = [6E0 6.*10.^-0.5 6E-1 6.*10.^-1.5 5E-2]; %2,0.2
setup.supTitleText = 'PSA.ENO.Km2PG';
setup.legendLocation =[0.5 0.925 0.1 0.05];
    A = string(setup.PSAvals.ENO_k2pg_ini);
setup.legNames = cellstr(A);
setup.addedTime = linspace(605,1000,60)';

% % PSA ENO_ini
% simulation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.PSAvals.ENO_k2pg_ini);
simResults = cell(1,nPSAvals);
for i = 1:nPSAvals
    setup.idx = i;
    [error,simResult] = simRes_costfun_Kmfixed_2(xres3,data,setup);
    simResults{i} = simResult;
end

% visualization
numpH = setup.numpHtested;
for j = 1:numpH   
    simulationVisualization_mainPSA;
end


%% (2D) what if Kmp2g even fall further away
xres3 = xres; xres3(1:2) = zeros; xres3(1:2) = [2.1990 0];
% xres3 = xres; xres3(1:2) = zeros; xres3(1:2) = [-1 0];
% kms to zero

% just the DF 1 is simulated
setup.DFstudy = 4;

% prior results: going too high or too down qill block.
setup.PSAstudy = 1; % if isfield
setup.PSAstudy_ENO_k2pg = 1; % if isfield + if ==
% setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
setup.PSAvals.ENO_k2pg_ini = [6E0 6E-1 5E-2];
setup.supTitleText = 'PSA.ENO.Km2PG';
setup.legendLocation =[0.35 0.925 0.1 0.05];
    A = string(setup.PSAvals.ENO_k2pg_ini);
setup.legNames = cellstr(A);
setup.addedTime = linspace(605,1000,60)';

% % PSA ENO_ini
% simulation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.PSAvals.ENO_k2pg_ini);
simResults = cell(1,nPSAvals);
for i = 1:nPSAvals
    setup.idx = i;
    [error,simResult] = simRes_costfun_Kmfixed_2(xres3,data,setup);
    simResults{i} = simResult;
end

% visualization
numpH = setup.numpHtested;
for j = 1:numpH   
    simulationVisualization_mainPSA;
end


%% (2E) what about the vaues for Kpep
% Keep PEP as low as possible to avoid it's effect.
xres3 = xres; xres3(1:2) = zeros;
% xres3 = xres; xres3(1:2) = zeros; xres3(1:2) = [-1 0];
% kms to zero

% just the DF 1 is simulated
setup.DFstudy = 4;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% prior results: going too high or too down qill block.
setup.PSAstudy = 1; % if isfield

% setup.PSAstudy_ENO_k2pg = 1; % if isfield + if ==
% % setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
% setup.PSAvals.ENO_k2pg_ini = [6E0 6E-1 5E-2];
% setup.supTitleText = 'PSA.ENO.Km2PG';

setup.PSAstudy_ENO_k2pg = 0; % if isfield + if ==
setup.PSAstudy_ENO_kpep = 1; % if isfield + if ==
% setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
% setup.PSAvals.ENO_k2pg_ini = [6E0 6E-1 5E-2];
setup.PSAvals.ENO_kpep = [-1 -0.5 0 0.5 1];
setup.supTitleText = 'PSA.ENO.KmPEP';

setup.legendLocation =[0.5 0.925 0.1 0.05];
    A = string(setup.PSAvals.ENO_kpep);
setup.legNames = cellstr(A);
setup.addedTime = linspace(605,1000,60)';
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % PSA ENO_ini
% simulation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.PSAvals.ENO_kpep);
simResults = cell(1,nPSAvals);
for i = 1:nPSAvals
    if setup.PSAstudy_ENO_kpep == 1
        xres4 = xres3;
        xres4(2) = setup.PSAvals.ENO_kpep(i);
    end
    [error,simResult] = simRes_costfun_Kmfixed_2(xres4,data,setup);
    simResults{i} = simResult;
end

% visualization
numpH = setup.numpHtested;
for j = 1:numpH   
    simulationVisualization_mainPSA;
end


%% 3. Test PSA for change in concentrations: (s/p) substrate/prooduct
%% 4. Create mock dataset
%% 5. Estimate back
%% 6. Other enzymes

%% memoryDump


%%
% for i = 1:length(simResults)
%     % recall
%     simResult2 = simResults{i};
%     nfig = 1000 + i;
%     % plot
%     figure(nfig)
%     subplot(131), plot(simResult2.t, simResult2.y(:,1), '.-'), title('p2g')
%     subplot(132), plot(simResult2.t, simResult2.y(:,2), '.-'), title('pep') 
%     subplot(133), plot(simResult2.t, simResult2.v, '.-')
% end
    
    
% figure(1001)
% subplot(131), plot(simResult2.t, simResult2.y(:,1), '.-'), title('p2g')
% subplot(132), plot(simResult2.t, simResult2.y(:,2), '.-'), title('pep') 
% subplot(133), plot(simResult2.t, simResult2.v, '.-'),  title('vENO')
% figure(1002)
% subplot(131), plot(simResult2.t, simResult2.y(:,1), '.-'), title('p2g')
% subplot(132), plot(simResult2.t, simResult2.y(:,2), '.-'), title('pep') 
% subplot(133), plot(simResult2.t, simResult2.v, '.-'),  title('vENO')
% figure(1003)
% subplot(131), plot(simResult2.t, simResult2.y(:,1), '.-'), title('p2g')
% subplot(132), plot(simResult2.t, simResult2.y(:,2), '.-'), title('pep') 
% subplot(133), plot(simResult2.t, simResult2.v, '.-'),  title('vENO')
