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
% options = optimset('Display','iter','MaxFunEvals',1);

% %% testing the costfunction
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
% lambdalist = 1E-5;
lambdalist = [...
    1E-5,...
    1E-4,...
    1E-3,...
    1E-2,...
    1E-1,...
    1E0,...
    1E1,...
    1E2,...
    1E3,...
    1E4,...
    1E5,...
    ];
haldaneList = [1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1E0, ...
    1E1, 1E2, 1E3, 1E4, 1E5,...
    1E6, 1E7, 1E8, 1E9, 1E10];
haldaneList = haldaneList(8);

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


% %% testing haldane relationship effect and simulations
% selHaldanePos = 8;
% i = 1;
% 
% % Haldane regularization
% eData = zeros(length(haldaneList),1);
% eHaldane = zeros(length(haldaneList),1);
% for j = 1:length(haldaneList)
%     eData(j) = sum(abs(array_eData{j,i}));
%     eHaldane(j) = sum(abs(array_eHaldane{j,i}));
% end
% % plotting
% f1 = figure(103);
% yyaxis left
% s1 = semilogx(haldaneList,eHaldane,'o-','MarkerSize',6);
% ylabel('error_{Haldane} []')
% xlabel(['Haldane reg. factor []: ',erase(sprintf('selected reg. = %d',haldaneList(selHaldanePos)),".000000")])
% hold on
% yyaxis right
% s2 = semilogx(haldaneList,eData,'o-','MarkerSize',6);
% ylabel('error_{Data} []')
% textHere = {['Regulatization ',setup.enzymeName,': '];'Errors in parameter values and data fit in blue and red, ';'respectively, in the y-axes. Regularization factor in the x-axis'};
% suptitle(textHere)
% l1 = line([haldaneList(selHaldanePos) haldaneList(selHaldanePos)],[s1.Parent.YLim(1) s1.Parent.YLim(2)]);
%     l1.Color = 'black';
%     l1.LineStyle = '--';
% set(103,'color','white')
% hold off
% 
% % simulation
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(array_xres{selHaldanePos},data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;


%% in the case of regularization at a given haldane = 100 (#8)
selHaldanePos = 8; %8
j = 1;
selLambdaPos = 5;

% Haldane regularization
eData = zeros(length(lambdalist),1);
eParameters = zeros(length(lambdalist),1);
for i = 1:length(lambdalist)
    eData(i) = sum(abs(array_eData{j,i}));
    eParameters(i) = sum(abs(array_eParams{j,i}));
end
% plotting
f1 = figure(103);
yyaxis left
s1 = semilogx(lambdalist,eParameters,'o-','MarkerSize',6);
ylabel('error_{Parameters} []')
xlabel(['Regularization factor []: ',erase(sprintf('selected reg. = %d',lambdalist(selLambdaPos)),".000000")])
hold on
yyaxis right
s2 = semilogx(lambdalist,eData,'o-','MarkerSize',6);
ylabel('error_{Data} []')
% textHere = {['Regulatization: '];'Errors in parameter values and data fit in blue and red, ';'respectively, in the y-axes. Regularization factor in the x-axis'};
% suptitle(textHere)
l1 = line([lambdalist(selLambdaPos) lambdalist(selLambdaPos)],[s1.Parent.YLim(1) s1.Parent.YLim(2)]);
    l1.Color = 'black';
    l1.LineStyle = '--';
set(103,'color','white')
hold off

% simulation
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
[error] = optfun(array_xres{selLambdaPos},data,setup);
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;


% % %% (2.3) Saving data after regularization
% if setup.plotOutput == 1
%     saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_regularizationResults.mat'];
% %     save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected')
%         save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected','array_resnorm','array_residual','array_Jacobian','array_raw_error')
%     set(101,'color','white'), savefig(101,['results/',setup.enzymeName,'/',setup.enzymeName, '_trainData_fit_metabolites_reg.fig']);
%     set(102,'color','white'), savefig(102,['results/',setup.enzymeName,'/',setup.enzymeName, '_testData_fit_fluxes_reg.fig']);
%     set(103,'color','white'), savefig(103,['results/',setup.enzymeName,'/',setup.enzymeName, '_regularization.fig']);
% end

% %% (2.4) Studying confidence intervals vs regularization
% plot_ParsCIs_lambdalist;
% if setup.plotOutput == 1
%     savefig(104,['results/',setup.enzymeName,'/',setup.enzymeName, '_ParsCIs_lambdalist.fig']);
% end


%% (3.1) Study on paarameter values: estimation with the lambda value
% % % % % recall specific values
% % % % setup.selectedLambda = lambdalist(selLambdaPos);
% % % % setup.weightHaldane = haldaneList(selHaldanePos);
% % % % error = array_raw_error{selHaldanePos};
% % % % % N = cell_N{selLambdaPos};
% % % % N = length(array_raw_error{selHaldanePos});
% % % % Jacobian = array_Jacobian{selHaldanePos};
% % % % resnorm = array_resnorm{selHaldanePos}; 
% % % % varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
% % % % stdp = sqrt(diag(varp));
% % % % % varp = cell_varp{selHaldanePos};
% % % % % stdp = cell_stdp{selHaldanePos};

% recall specific values
setup.selectedLambda = lambdalist(selLambdaPos);
setup.weightHaldane = haldaneList(1);
error = array_raw_error{selLambdaPos};
% N = cell_N{selLambdaPos};
N = length(array_raw_error{selLambdaPos});
Jacobian = array_Jacobian{selLambdaPos};
resnorm = array_resnorm{selLambdaPos}; 
varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
stdp = sqrt(diag(varp));
% varp = cell_varp{selHaldanePos};
% stdp = cell_stdp{selHaldanePos};


%% ENZYME-SPECIFIC
%% (3.2) Study on paarameter values: recalculation

% parameter values (pre check)
xres_zero = array_xres{1};
% % % % xres_selected = array_xres{selHaldanePos};
xres_selected = array_xres{selLambdaPos};

figure(201),

subplot(231)
plot(xres_zero(1:12),'.-')
hold on
plot(xres_selected(1:12),'.-')
title('k_{gap}')
xlim([0 12])
hold off

subplot(232)
plot(xres_zero(13:24),'.-')
hold on
plot(xres_selected(13:24),'.-')
title('k_{bpg}')
xlim([0 12])
hold off

subplot(233)
plot(xres_zero(25:36),'.-')
hold on
plot(xres_selected(25:36),'.-')
title('k_{nad}')
xlim([0 12])
hold off

subplot(234)
plot(xres_zero(37:48),'.-')
hold on
plot(xres_selected(37:48),'.-')
title('k_{nadh}')
xlim([0 12])
legend('xres_{zero}','xres_{selected}','location','south')
hold off

subplot(235)
plot(xres_zero(49:60),'.-')
hold on
plot(xres_selected(49:60),'.-')
title('v_{max.fwd}')
xlim([0 12])
hold off

subplot(236)
plot(xres_zero(61:72),'.-')
hold on
plot(xres_selected(61:72),'.-')
title('v_{max.rev}')
xlim([0 12])
hold off


%% parameter especification and c.i.
numpHtested = 12;
    pH_fwd = data.FWD.pH(:,4);
    pH_rev = data.REV.pH(:,4);
    vmax_rev = data.REV.Vmax(:,4);
    vmax_rev_interp = interp1(pH_rev,vmax_rev,pH_fwd,'pchip');

kgap = ones(numpHtested,1);
kbpg = ones(numpHtested,1);
knad = ones(numpHtested,1);
knadh = ones(numpHtested,1);
vm_fwd = zeros(numpHtested,1);
vm_rev = zeros(numpHtested,1);
for i = 1:numpHtested
    kgap(i) = 2.48 * 10.^xres_selected(i+0); %mM
    kbpg(i) = 1.18 * 10.^xres_selected(i+12); %mM
    knad(i) = 2.92 * 10.^xres_selected(i+24); %mM
    knadh(i) = 0.022 * 10.^xres_selected(i+36); %mM
    vm_fwd(i) = data.FWD.Vmax(i) * 10 .^ xres_selected(i+48);
    vm_rev(i) = vmax_rev_interp(i) * 10 .^ xres_selected(i+60);
end
keq_gapdh = setup.FWD.pH_Keq_gapdh_eQ_fwd;
keq_pgk = setup.FWD.pH_Keq_pgk_fwd;

% limits
kgap_up = ones(numpHtested,1);
kbpg_up = ones(numpHtested,1);
knad_up = ones(numpHtested,1);
knadh_up = ones(numpHtested,1);
vm_fwd_up = zeros(numpHtested,1);
vm_rev_up = zeros(numpHtested,1);
kgap_down = ones(numpHtested,1);
kbpg_down = ones(numpHtested,1);
knad_down = ones(numpHtested,1);
knadh_down = ones(numpHtested,1);
vm_fwd_down = zeros(numpHtested,1);
vm_rev_down = zeros(numpHtested,1);
for i = 1:numpHtested
    %up
    kgap_up(i) = 2.48 * 10.^(xres_selected(i+0) + stdp(i+0)); %mM
    kbpg_up(i) = 1.18 * 10.^(xres_selected(i+12) + stdp(i+12)); %mM
    knad_up(i) = 2.92 * 10.^(xres_selected(i+24) + stdp(i+24)); %mM
    knadh_up(i) = 0.022 * 10.^(xres_selected(i+36) + stdp(i+36)); %mM
    vm_fwd_up(i) = data.FWD.Vmax(i) * 10.^(xres_selected(i+48) + stdp(i+48));
    vm_rev_up(i) = vmax_rev_interp(i) * 10.^(xres_selected(i+60) + stdp(i+60));
    %down
    kgap_down(i) = 2.48 * 10.^(xres_selected(i+0) - stdp(i+0)); %mM
    kbpg_down(i) = 1.18 * 10.^(xres_selected(i+12) - stdp(i+12)); %mM
    knad_down(i) = 2.92 * 10.^(xres_selected(i+24) - stdp(i+24)); %mM
    knadh_down(i) = 0.022 * 10.^(xres_selected(i+36) - stdp(i+36)); %mM
    vm_fwd_down(i) = data.FWD.Vmax(i) * 10.^(xres_selected(i+48) - stdp(i+48));
    vm_rev_down(i) = vmax_rev_interp(i) * 10.^(xres_selected(i+60) - stdp(i+60));
end

vm_fwd_up_uChange = vm_fwd_up .* 60 .* 60 ./ setup.FWD.concProtein;
vm_fwd_down_uChange = vm_fwd_down .* 60 .* 60 ./ setup.FWD.concProtein;
vm_fwd_uChange = vm_fwd .* 60 .* 60 ./ setup.FWD.concProtein;

vm_rev_up_uChange = vm_rev_up .* 60 .* 60 ./ setup.REV.concProtein;
vm_rev_down_uChange = vm_rev_down .* 60 .* 60 ./ setup.REV.concProtein;
vm_rev_uChange = vm_rev .* 60 .* 60 ./ setup.REV.concProtein;


%% (3.3) Study on paarameter values: plotting
Vmax_fwd_experimental = setup.FWD.exp_vmax_gapdhf*60*60/setup.FWD.concProtein;
Vmax_rev_experimental = setup.FWD.exp_vmax_gapdhr*60*60/setup.REV.concProtein;
pH_fwd = setup.FWD.pH_vals;
pH_rev = setup.REV.pH_vals;
keq_gapdh = setup.FWD.pH_Keq_gapdh_eQ_fwd;
keq_pgk = setup.FWD.pH_Keq_pgk_fwd;

figure(105)

subplot(341) % vm_fwd
    plot(pH_fwd,vm_fwd_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pH_fwd,vm_fwd_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pH_fwd,vm_fwd_uChange,'.-','color','black')
    title({'v_{m.fwd} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

subplot(342) % vm_rev
    plot(pH_fwd,vm_rev_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pH_fwd,vm_rev_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pH_fwd,vm_rev_uChange,'.-','color','black')
	title({'v_{m.rev} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

subplot(343) % vm_fwd_exp
    plot(pH_rev,Vmax_fwd_experimental,'.-')
    title({'v_{m.fwd} experimental';'not normalized'})

subplot(344) % vm_rev_exp
    plot(pH_rev,Vmax_rev_experimental,'.-')
    title({'v_{m.rev} experimental';'not normalized'})
    
subplot(345) % vm_fwd(log_difference)
    plot(pH_fwd,xres_selected(49:60)','.-','color','black')
    title('v_{m.fwd} difference []')
    
subplot(346) % vm_rev(log_difference)
    plot(pH_fwd,xres_selected(61:72)','.-','color','black')
    title('v_{m.rev} difference []')

subplot(347) % keq_gapdh
    plot(pH_fwd,keq_gapdh,'.-','color','black')
    title('k_{eq.GAPDH} [mM]')

subplot(348) % keq_pgk
    plot(pH_fwd,keq_pgk,'.-','color','black')
    title('k_{eq.PGK} [mM]')
    
subplot(349) % kpyr
    plot(pH_fwd,kgap_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pH_fwd,kgap_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pH_fwd,kgap,'.-','color','black')
title('k_{gap} [mM]')

subplot(3,4,10) % kbpg
    plot(pH_fwd,kbpg_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pH_fwd,kbpg_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pH_fwd,kbpg,'.-','color','black')
title('k_{bpg} [mM]')

subplot(3,4,11) % knad
    plot(pH_fwd,knad_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pH_fwd,knad_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pH_fwd,knad,'.-','color','black')
title('k_{nad} [mM]')

subplot(3,4,12) % knadh
    plot(pH_fwd,knadh_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pH_fwd,knadh_down,'.-','color',[0.5 0.5 0.5]), hold on, 
plot(pH_fwd,knadh,'.-','color','black')
title('k_{nadh} [mM]')

suptitle('GAPDH_{fwd+rev}: parameter estimates vs pH')


% %%
% output_gapdh.xres_selected = xres_selected;
% 
% output_gapdh.pHarray = pHarray;
% 
% output_gapdh.kgap = kgap;% = ones(numpHtested,1);
% output_gapdh.kbpg = kbpg;% = ones(numpHtested,1);
% output_gapdh.knad = knad;% = ones(numpHtested,1);
% output_gapdh.knadh = knadh;% = ones(numpHtested,1);
% output_gapdh.vm = vm;% = zeros(numpHtested,1);
% output_gapdh.vm_uChange = vm_uChange;% = zeros(numpHtested,1);
% 
% output_gapdh.kgap_up = kgap_up;% = ones(numpHtested,1);
% output_gapdh.kbpg_up = kbpg_up;% = ones(numpHtested,1);
% output_gapdh.knad_up = knad_up;% = ones(numpHtested,1);
% output_gapdh.knadh_up = knadh_up;% = ones(numpHtested,1);
% output_gapdh.vm_up = vm_up;% = zeros(numpHtested,1);
% output_gapdh.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);
% 
% output_gapdh.kgap_down = kgap_down;% = ones(numpHtested,1);
% output_gapdh.kbpg_down = kbpg_down;% = ones(numpHtested,1);
% output_gapdh.knad_down = knad_down;% = ones(numpHtested,1);
% output_gapdh.knadh_down = knadh_down;% = ones(numpHtested,1);
% output_gapdh.vm_down = vm_down;% = zeros(numpHtested,1);
% output_gapdh.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);
% 
% output_gapdh.keq_gapdh = keq_gapdh;% = setup.Keq_ADH;
% output_gapdh.keq_pgk = keq_pgk;% = setup.Keq_ADH;
% 
% 
% %% saving output
% output_gapdh.Vmax_experimental = Vmax_experimental;
% output_gapdh.stDev_experimental = stDev_experimental;
% saveName = ['results/',setup.enzymeName,'/',setup.enzymeName, '_parEst.mat'];
% save(saveName,'output_gapdh');
% set(105,'color','white'), savefig(105,['results/',setup.enzymeName,'/',setup.enzymeName, '_parEstimates.fig']);
% 
% 
% % % % % x = 1; end
% 
% 
% %% memoryDump
% % %% copying files
% % % copyfile dev_pEst_gapdh_fwd_rev.m A1.m
% % %          dev_pEst_gapdh_fwd_rev
% % % copyfile('dev_pEst_gapdh_fwd_rev.m','test1.m')
% % for i = 2:16
% %     fid = fopen('A1.m','rt');
% %     X = fread(fid);
% %     fclose(fid);
% %     X = char(X.');
% %     str2rep = sprintf('num = %d;',i);
% %     Y = strrep(X,'num = 1;',str2rep);
% %     tempName = sprintf('A%d.m',i);
% %     fid2 = fopen(tempName,'wt') ;
% %     fwrite(fid2,Y) ;
% %     fclose(fid2);
% % end
% % %% putting arrays together
% % i = 1;
% % for j = 1:length(haldaneList)
% %     % loading data
% %     tempName = sprintf('tempRes_A%d.mat',j);
% %     load(tempName,'xres','error','resnorm','residual','Jacobian','raw_error');
% %     % locating up data again
% %     array_xres{j,i} = xres;
% %     array_eData{j,i} = error(1:end-61);
% %     array_eHaldane{j,i} =  error(end-60:end-48);
% %     array_eParams{j,i} = error(end-47:end);
% %     array_resnorm{j,i} = resnorm;
% %     array_residual{j,i} = residual;
% %     array_Jacobian{j,i} = Jacobian;
% %     array_raw_error{j,i} = raw_error;
% % end

