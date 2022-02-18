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

    %     figure
    %     plot(pHvals, Vmax(:,4),'.-')
    %     title('Starting estimate: Vmax [mM s-1] vs pH')    
    end

    % % (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 70; % minimum size of the window
    limRates = [0 1.5E-4]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
    % select start point (this needs manual selection deppending on previous plots)
    dp_start = 6 * ones(size(data.conc_mean));
    % blank cell total length
    total_len = zeros(size(dp_start));
    % DFs considered
    DFarray_pre = 1./data.DF;
    DFarray = DFarray_pre(1,:);
    % idxs2consider
    idxs2consider = [0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1;
                    0 1 1 1];
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
    setup.problemStudy = 'fullKinetics_paramsPartFixed';
    problemStudy = setup.problemStudy;
    optfun = @costfun_Kmfixed;
    switch problemStudy
        case 'onlyVmax'
            plength = 12; % only the vmax change with pH. Other params OFF
        case 'fullKinetics_paramsPartFixed'
            plength = 25; % only the vmax change with pH. Other params ON
        case 'fullKinetics_paramsAllFlexible'
            plength = 168; % All parameters change with pH. Other params ON
        otherwise
            disp('Warning: problem study (reaction kinetics) have not been selected');
    end
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
    options = optimset('Display','iter');
    
    % (1.2) Check that we get the right result parameter set
    close all
    load('pfk_parEst.mat')
    xres = output_pfk.xres_selected;
    % testing the costfunction
    setup.plotEachSimCF = 1;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 1;
    [error] = optfun(xres,data,setup);
    setup.plotEachSimCF = 0;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0;
end


%% temp
% xres2 = xres; %xres2(1:4) = zeros(1,4);
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
[error,simResult] = simRes_costfun_Kmfixed_2(xres,data,setup);
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;


%% (2A) initial testing of Laura's protocol concentrations + longer timespan
% little change
xres3 = xres; %xres3(1:2) = zeros;
% kms to zero

% just the DF 1 is simulated
setup.DFstudy = 4;
% prior results: going too high or too down qill block.
setup.PSAstudy = 1; % if isfield
setup.PSAstudy_PFK_kf6p = 1; % if isfield + if ==
% setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
setup.PSAvals.PFK_kf6p_ini = [0 1E1 1E0 1E-1];
setup.supTitleText = 'PSA.PFK.KmF6P';
setup.legendLocation =[0.35 0.925 0.1 0.05];
    A = string(setup.PSAvals.PFK_kf6p_ini);
setup.legNames = cellstr(A);
setup.addedTime = [];%linspace(605,1000,60)';

% % PSA ENO_ini
% simulation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.PSAvals.PFK_kf6p_ini);
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


%% (2B) Maybe more values tested
% decreases the rate slightly and more visual at the rate level.
% But seems to be limited by other metabolites

% little change
xres3 = xres; %xres3(1:2) = zeros;
% kms to zero

% just the DF 1 is simulated
setup.DFstudy = 4;
% prior results: going too high or too down qill block.
setup.PSAstudy = 1; % if isfield
setup.PSAstudy_PFK_kf6p = 1; % if isfield + if ==
% setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
% setup.PSAvals.PFK_kf6p_ini = [1E1 1E0 1E-1 1E-2 1E-3];
setup.PSAvals.PFK_kf6p_ini = [1E3 1E2 1E1 1E0 1E-1];
setup.supTitleText = 'PSA.PFK.KmF6P';
setup.legendLocation = [0.35 0.925 0.1 0.05];
    A = string(setup.PSAvals.PFK_kf6p_ini);
setup.legNames = cellstr(A);
setup.addedTime = [];%linspace(605,1000,60)';

% % PSA ENO_ini
% simulation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.PSAvals.PFK_kf6p_ini);
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


%% (2C) What if reference parameters?
% little change
xres3 = zeros(size(xres));
% kms to zero

% just the DF 1 is simulated
setup.DFstudy = 4;
% prior results: going too high or too down qill block.
setup.PSAstudy = 1; % if isfield
setup.PSAstudy_PFK_kf6p = 1; % if isfield + if ==
% setup.PSAvals.ENO_k2pg_ini = [6E1 6E0 6E-1 5E-2 5E-3];
setup.PSAvals.PFK_kf6p_ini = [1E1 1E0 1E-1];
setup.supTitleText = 'PSA.PFK.KmF6P';
setup.legendLocation =[0.35 0.925 0.1 0.05];
    A = string(setup.PSAvals.PFK_kf6p_ini);
setup.legNames = cellstr(A);
setup.addedTime = [];%linspace(605,1000,60)';

% % PSA ENO_ini
% simulation
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

nPSAvals = length(setup.PSAvals.PFK_kf6p_ini);
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




