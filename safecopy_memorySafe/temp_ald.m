%% (0) Setup and data load

%% (0) Setup and data load
clear, close all
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
    
% % % %     figure
% % % %     plot(pHvals, Vmax(:,4),'.-')
% % % %     title('Starting estimate: Vmax [mM s-1] vs pH')    
end


% %% (0.1) Experimental Vmax determination
setup.plotOutput = 0;
% correcting DF, specific for pdc case
data.DF(7:8,3) = [4;4]; 
DF(7:8,3) = [4;4];
% %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 6; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
% select start point (this needs manual selection deppending on previous plots)
dp_start = ones(size(data.conc_mean));
dp_start(:,1) = 10 * ones(size(dp_start(:,1)));
dp_start(:,2) = 10 * ones(size(dp_start(:,2)));
dp_start(:,3) = 5 * ones(size(dp_start(:,3)));
dp_start(:,4) = 2 * ones(size(dp_start(:,4)));
% blank cell total length
total_len = zeros(size(dp_start));
% DFs considered
DFarray = [1/8 1/4 1/2 1/1];
% idxs2consider
idxs2consider = [0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 1 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1];

% Experimental rates determination and plotting
expRatesDetermination;

% %% (1.1) Simple parameter fit. Parameter estimation
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';

setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
setup.costfun = 3;

setup.weightData = 1;
setup.weightDataEsp = idxs2consider;
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% Km fixed
optfun = @costfun_Kmfixed;
plength = 14; % Kms (2) + Vms (1) * numpH (12)
x_temp = zeros(1,plength);
% % % % ub = 3*ones(1,plength);
% % % % lb = -3*ones(1,plength);
ub = 1*ones(1,plength);
lb = -1*ones(1,plength);
options = optimset('Display','iter');
% %%

% %% (2.1) Parameter estimation with regularization
% parameter estimation
lambdalist = 1;
parameterEstimation_lambdalist;

% %% (2.2) Regularization. Results Visualization
selLambdaPos = 1;%15;%1;%17;%14;
regularizationVisualization;

%
xres_Keq_pH_dependent = array_xres{selLambdaPos};

%
for i = 1:length(setup.Keq_ADH)
    setup.Keq_ADH(i) = setup.Keq_ADH(6);
end

% %% parameter estimation pH independent
tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t = toc;
% 
xres_Keq_pH_independent = xres;

% parameter values
vm_keq_pH_dependent = zeros(numpHtested,1);
vm_keq_pH_independent = zeros(numpHtested,1);
for i = 1:numpHtested
    vm_keq_pH_dependent(i) = data.Vmax(i,4) * 10.^xres_Keq_pH_dependent(i+2);
    vm_keq_pH_independent(i) = data.Vmax(i,4) * 10.^xres_Keq_pH_independent(i+2);
end
vm_keq_pH_dependent_uChange = vm_keq_pH_dependent .* 60 .* 60 ./ setup.concProtein;
vm_keq_pH_independent_uChange = vm_keq_pH_independent .* 60 .* 60 ./ setup.concProtein;
% 
updated_keqconst_pdc.vm_keq_pH_dependent_uChange = vm_keq_pH_dependent_uChange;
updated_keqconst_pdc.vm_keq_pH_independent_uChange = vm_keq_pH_independent_uChange;
updated_keqconst_pdc.pH = pH;




