% % PEST_GAPDH_FWD_REV
% In this code, we test again if we can combine gapdh_fwd and gapdh_rev to
% estimate reaction parameters
clear, close all
set_paths_pHstudy;
dbstop if error


%% Load and relocate data

% forward direction
for loading_gapdhFWD = 1
    % select specific case and recall data
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 1;
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
    setup.plotOutput = 0;
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_gapdh = expData.gapdh;
    % changes in the dataset for gapdh_fwd
        %         excelSheet_corrected: [2×23 double]
        import_gapdh.treatedData.excelSheet_corrected(1,24) = import_gapdh.treatedData.excelSheet_corrected(1,23);
        import_gapdh.treatedData.excelSheet_corrected(2,24) = import_gapdh.treatedData.excelSheet_corrected(2,23);
        %             pH_corrected: [2×23 double]
        import_gapdh.treatedData.pH_corrected(1,24) = import_gapdh.treatedData.pH_corrected(1,23);
        import_gapdh.treatedData.pH_corrected(2,24) = import_gapdh.treatedData.pH_corrected(2,23);
        %       dilution_corrected: [2×23 double]
        import_gapdh.treatedData.dilution_corrected(1,24) = import_gapdh.treatedData.dilution_corrected(1,23);
        import_gapdh.treatedData.dilution_corrected(2,24) = import_gapdh.treatedData.dilution_corrected(2,23);
        %          absorbance_mean: {2×23 cell}
        import_gapdh.treatedData.absorbance_mean{1,24} = import_gapdh.treatedData.absorbance_mean{1,23};
        import_gapdh.treatedData.absorbance_mean{2,24} = import_gapdh.treatedData.absorbance_mean{2,23};
        %       absorbance_samples: {2×23 cell}
        import_gapdh.treatedData.absorbance_samples{1,24} = import_gapdh.treatedData.absorbance_samples{1,23};
        import_gapdh.treatedData.absorbance_samples{2,24} = import_gapdh.treatedData.absorbance_samples{2,23};
        %           absorbance_std: {2×23 cell}
        import_gapdh.treatedData.absorbance_std{1,24} = import_gapdh.treatedData.absorbance_std{1,23};
        import_gapdh.treatedData.absorbance_std{2,24} = import_gapdh.treatedData.absorbance_std{2,23};
        %       concentration_mean: {2×23 cell}
        import_gapdh.treatedData.concentration_mean{1,24} = import_gapdh.treatedData.concentration_mean{1,23};
        import_gapdh.treatedData.concentration_mean{2,24} = import_gapdh.treatedData.concentration_mean{2,23};
        %        concentration_std: {2×23 cell}
        import_gapdh.treatedData.concentration_std{1,24} = import_gapdh.treatedData.concentration_std{1,23};
        import_gapdh.treatedData.concentration_std{2,24} = import_gapdh.treatedData.concentration_std{2,23};
        %                     time: {2×23 cell}
        import_gapdh.treatedData.time{1,24} = import_gapdh.treatedData.time{1,23};
        import_gapdh.treatedData.time{2,24} = import_gapdh.treatedData.time{2,23};
        %            reaction_rate: {2×23 cell}
        import_gapdh.treatedData.reaction_rate{1,24} = import_gapdh.treatedData.reaction_rate{1,23};
        import_gapdh.treatedData.reaction_rate{2,24} = import_gapdh.treatedData.reaction_rate{2,23};
    
    
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
    
    pHarray = unique(import_gapdh.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_gapdh.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_gapdh.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_gapdh.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_gapdh.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_gapdh.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_gapdh.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_gapdh.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_gapdh.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_gapdh.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_gapdh.rawData.absorbance_corrected{4,4};
    temp2 = import_gapdh.rawData.absorbance_corrected{5,4};
    temp3 = import_gapdh.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_gapdh.rawData.time{1};
    
        % (1) Correct for minimum value
        % (2) Bring the minimum to zero (apply to all)
        % (3) In principle, use the 3 first dilution rates
        % (4) Watch out with the dilution factors (first 2 cases are
        % reversed)
        
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Directly changing the concentration here, sicne the extinction
    % coefficient did not change.
    dps = length(NADH{1,1});
    for i = 1:DFs
        for j = 1:numpHtested
            dps = length(NADH{j,i});
            k2 = 1:dps;
            for k = 1:dps
                k3 = k2(end-k+1);
%                 switch j
%                     case {1,2,3,4,5,6,7}
%                         NADH{j,i}(k) = NADH{j,i}(k) - NADH{j,DFs}(dps);
%                         NADH{j,i}(k) = NADH{j,i}(k) - NADH{j,DFs}(1);
                        NADH{j,i}(k3) = NADH{j,i}(k3) - NADH{j,i}(1);
%                     case {8,9,10,11,12}
%                         NADH{j,i}(k) = NADH{j,i}(k) - 0.0453;
%     %                     NADH{j,i}(k) = NADH{j,i}(k) - NADH{6,DFs}(dps);
%                 end
            end
        end
    end
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    pHvals = unique(import_gapdh.treatedData.pH_corrected);
    % visualize: check calculations made
    if setup.plotOutput == 1
        figure('units','normalized','outerposition',[0 0 1 1])
        for i = 1:numpHtested
            subplot(3,4,i)
            for j = 1:DFs
                plot(time{i,j},NADH{i,j},'.-')
                hold on
            end
%             title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
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
            if setup.caseStudyGAPDHr == 1
                ylim([0 0.15])
            end
        end
        suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
        suptitle(suptitleName);
    end
%     figure
%     plot(pHvals, Vmax(:,4),'.-')
%     title('Starting estimate: Vmax [mM s-1] vs pH')  

    % experimental Vmax determination
    % intial things that could be in the setup
    minwindow = 6; % minimum size of the window
    limRates = [0 1E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.05];  %Ylims plot conc
    % select start point (this needs manual selection deppending on previous plots)
    dp_start = ones(size(data.conc_mean));
    % blank cell total length
    total_len = zeros(size(dp_start));
    % DFs considered
    DFarray = [1/8 1/4 1/2 1/1];
    % idxs2consider
    idxs2consider = [1 1 1 1;
                    1 1 1 1;
                    1 1 1 0;
                    1 1 1 1;
                    1 1 1 1;
                    1 1 1 1;
                    1 1 1 1;
                    1 1 1 1;
                    1 1 1 1;
                    1 1 1 1;
                    1 1 1 1;
                    1 1 1 0];
    setup.idxs2consider = idxs2consider;

    % Experimental rates determination and plotting
    expRatesDetermination;
    
end
dataFWD = data;
setupFWD = setup;
clear data setup idxs2consider

% reverse direction
for loading_gapdhREV = 1
    % select specific case and recall data
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 1;
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
    setup.plotOutput = 0;
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_gapdhr = expData.gapdhr;
    
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
    
    pHarray = unique(import_gapdhr.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_gapdhr.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_gapdhr.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_gapdhr.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_gapdhr.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_gapdhr.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_gapdhr.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_gapdhr.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_gapdhr.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_gapdhr.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_gapdhr.rawData.absorbance_corrected{4,4};
    temp2 = import_gapdhr.rawData.absorbance_corrected{5,4};
    temp3 = import_gapdhr.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_gapdhr.rawData.time{1};
    
        % (1) Correct for minimum value
        % (2) Bring the minimum to zero (apply to all)
        % (3) In principle, use the 3 first dilution rates
        % (4) Watch out with the dilution factors (first 2 cases are
        % reversed)
        
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Directly changing the concentration here, sicne the extinction
    % coefficient did not change.
    dps = length(NADH{1,1});
    for i = 1:DFs
        for j = 1:numpHtested
            for k = 1:dps
                switch j
                    case {1,2,3,4,5,6,7}
                        NADH{j,i}(k) = NADH{j,i}(k) - NADH{j,DFs}(dps);
                    case {8,9,10,11,12}
                        NADH{j,i}(k) = NADH{j,i}(k) - 0.0453;
    %                     NADH{j,i}(k) = NADH{j,i}(k) - NADH{6,DFs}(dps);
                end
            end
        end
    end
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    pHvals = unique(import_gapdhr.treatedData.pH_corrected);
    % visualize: check calculations made
    if setup.plotOutput == 1
        figure('units','normalized','outerposition',[0 0 1 1])
        for i = 1:numpHtested
            subplot(3,4,i)
            for j = 1:DFs
                plot(time{i,j},NADH{i,j},'.-')
                hold on
            end
%             title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
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
            if setup.caseStudyGAPDHr == 1
                ylim([0 0.15])
            end
        end
        suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
        suptitle(suptitleName);
    end
    
%     figure
%     plot(pHvals, Vmax(:,4),'.-')
%     title('Starting estimate: Vmax [mM s-1] vs pH') 

    % %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 6; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
    % select start point (this needs manual selection deppending on previous plots)
    dp_start = ones(size(data.conc_mean));
    % blank cell total length
    total_len = zeros(size(dp_start));
    % DFs considered
    DFarray = [1/8 1/4 1/2 1/1];
    % idxs2consider
    idxs2consider = ones(size(DF));
    setup.idxs2consider = idxs2consider;

    % Experimental rates determination and plotting
    expRatesDetermination;

end
dataREV = data;
setupREV = setup;
clear data setup idxs2consider

% relocate
data.FWD = dataFWD;
data.REV = dataREV;
setup.FWD = setupFWD;
setup.REV = setupREV;
clear dataFWD dataREV setupFWD setupREV


%% setup cost function and test run

% specific setup
setup.FWD.ode = 'gapdh_fwd_simplified';
setup.FWD.sourceVm = 'experimentalSlopes';
setup.FWD.ode_pH = 'on';
setup.FWD.weightDataEsp = setup.FWD.idxs2consider;
setup.REV.ode = 'gapdh_rev_simplified';
setup.REV.sourceVm = 'experimentalSlopes';
setup.REV.ode_pH = 'on';
setup.REV.weightDataEsp = setup.REV.idxs2consider;

% common setup
setup.caseKm = 'pH_dependent'; % setup.caseKm = 'pH_independent';
setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;
setup.numpHtested = length(setup.FWD.pHtested);
setup.DFstudy = 1:4;
setup.costfun = 3;
setup.weightData = 1;
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% cost function
caseKm = setup.caseKm;
switch caseKm
    case 'pH_independent'
        plength = 28; % Kms (4) + Vms (2) * numpH (12)
    case 'pH_dependent'
        plength = 72; % Kms (4) * numpH (12) + Vms (1) * numpH (12)
    otherwise
        disp('Warning: no specification has been made on Km being pH dependent or independent');
end
optfun = @costfun_gapdhFWDREV;
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');

% % testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(x_temp,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;


%% Parameter Estimation
setup.enzymeName = 'gapdh_FWD_REV';
% parameter estimation
lambdalist = 1E-5;
% lambdalist = [...
%     1E-5, 2E-5, 5E-5,...
%     1E-4, 2E-4, 5E-4,...
%     1E-3, 2E-3, 5E-3,...
%     1E-2, 2E-2, 5E-2,...
%     1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 7E-1,... %area of change
%     1E0, 2E0, 3E0, 4E0, 5E0, 7E0,...%area of change
%     1E1, 2E1, 5E1,...
%     1E2, 2E2, 5E2,...
%     1E3, 2E3, 5E3,...
%     1E4, 2E4, 5E4,...
%     1E5, 2E5, 5E5,...
%     ];
haldaneList = [1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1E0, ...
    1E1, 1E2, 1E3, 1E4, 1E5,...
    1E6, 1E7, 1E8, 1E9, 1E10];
num = 7;
haldaneList = haldaneList(num);

% parameter estimation
array_xres = cell(length(haldaneList),length(lambdalist));

array_eData = cell(length(haldaneList),length(lambdalist));
array_eHaldane = cell(length(haldaneList),length(lambdalist));
array_eParams = cell(length(haldaneList),length(lambdalist));

array_resnorm = cell(length(haldaneList),length(lambdalist));
array_residual = cell(length(haldaneList),length(lambdalist));
array_Jacobian = cell(length(haldaneList),length(lambdalist));
array_raw_error = cell(length(haldaneList),length(lambdalist));

for j = 1:length(haldaneList)
    setup.weightHaldane = haldaneList(j);
    for i = 1:length(lambdalist)
        fprintf('pEst for lambda=%d, haldane=%d\n',lambdalist(i),haldaneList(j));
        setup.selectedLambda = lambdalist(i);
        % parameter estimation
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
        t = toc;
        [raw_error] = optfun(xres,data,setup);
        setup.selectedLambda = 1;
        setup.weightHaldane = 1;
        [error] = optfun(xres,data,setup);
        % seting in output arrays
        array_xres{j,i} = xres;
        array_eData{j,i} = error(1:end-61);
        array_eHaldane{j,i} =  error(end-60:end-48);
        array_eParams{j,i} = error(end-47:end);
        array_resnorm{j,i} = resnorm;
        array_residual{j,i} = residual;
        array_Jacobian{j,i} = Jacobian;
        array_raw_error{j,i} = raw_error;
    end

end

% save
saveName = sprintf('tempRes_A%d.mat',num);
save(saveName);

