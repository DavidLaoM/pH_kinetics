% % PEST_TPI.m

%% (0) Setup and data load
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
    setup.caseStudyTPI = 1;
    setup.caseStudyENO = 0;
    selectSetup_pH;
    % added
%     setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_tpi = expData.tpi;
    
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
    
    pHarray = unique(import_tpi.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_tpi.treatedData.pH_corrected==pHval);
        if pHval == 6.19
            tempID = tempID(1:4);
        else
        end
        pHTemp(:,i) = import_tpi.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_tpi.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_tpi.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_tpi.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_tpi.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_tpi.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_tpi.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_tpi.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_tpi.rawData.absorbance_corrected{4,4};
    temp2 = import_tpi.rawData.absorbance_corrected{5,4};
    temp3 = import_tpi.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_tpi.rawData.time{1};
    
        % (1) Correct for minimum value
        % (2) Bring the minimum to zero (apply to all)
        % (3) In principle, use the 3 first dilution rates
        % (4) Watch out with the dilution factors (first 2 cases are
        % reversed)
        
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    dps = length(NADH{1,1});
    stopVals_dps = [ 40, 23, 14, 9;...
                     40, 24, 15, 10;...
                     35, 20, 11, 8;...
                     34, 19, 11, 8;...
                     34, 19, 12, 10;...
                     30, 16, 10, 8;...
                     23, 13, 9, 8];

    % Adjusted to TPI
    % Directly changing the concentration here, since the extinction
    % coefficient did not change.
    endPoint = zeros(numpHtested,DFs);
    % locate the minimum
    for i = 1:DFs
        for j = 1:numpHtested
            endPoint(j,i) = NADH{j,i}(stopVals_dps(j,i));
        end
    end
    % from that value onwards, copy the value
     for i = 1:DFs
        for j = 1:numpHtested
            inival = stopVals_dps(j,i)+1;
            for k = inival:dps
                NADH{j,i}(k) = NADH{j,i}(inival-1);
            end
        end
    end
    
    % bringing the minimum to zero
    for i = 1:DFs
        for j = 1:numpHtested
                for k = 1:dps
                    NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,i);
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

    pHvals = unique(import_tpi.treatedData.pH_corrected);
    % visualize: check calculations made
    if(setup.plotOutput == 1)
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
    end  
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % ADDED 2020/11/17  % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % setup.plotOutput = 0;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%% (0.1) Experimental Vmax determination
% setup.plotOutput = 0;
% %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 5; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
% select start point (this needs manual selection deppending on previous plots)
dp_start = 2 * ones(size(data.conc_mean));
% blank cell total length
total_len = zeros(size(dp_start));
% DFs considered
DFarray = [1/8 1/4 1/2 1/1];
% idxs2consider
idxs2consider = [0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1;
                0 0 1 1];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % ADDED 2020/11/17  % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % idxs2consider = [0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0;
% % % %                 0 0 1 0];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Experimental rates determination and plotting
expRatesDetermination;

% %% (0.3) saving
save_initial_recall

%% (1.1) Simple parameter fit. Parameter estimation
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';

setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 0;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
% % % % setup.costfun = 1;
setup.costfun = 3;

setup.weightData = 1;
setup.weightDataEsp = idxs2consider;
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% Km fixed
optfun = @costfun_Kmfixed;
plength = 9; % Kms (2) + Vms (1) * numpH (7)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');
% %%
% % testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% % [error] = optfun(x_temp,data,setup);
% [error] = optfun(xres,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;
% %%
% % parameter estimation
% tic
% [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t = toc;


%% (2.1) Parameter estimation with regularization
% parameter estimation
lambdalist = [...
    1E-5, 2E-5, 5E-5,...
    1E-4, 2E-4, 5E-4,...
    1E-3, 2E-3, 5E-3,...
    1E-2, 2E-2, 5E-2,...
    1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 7E-1,... %area of change
    1E0, 2E0, 3E0, 4E0, 5E0, 7E0,...%area of change
    1E1, 2E1, 5E1,...
    1E2, 2E2, 5E2,...
    1E3, 2E3, 5E3,...
    1E4, 2E4, 5E4,...
    1E5, 2E5, 5E5,...
    ];
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % ADDED 2020/11/17  % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % lambdalist = 0;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
parameterEstimation_lambdalist;


%% (2.2) Regularization. Results Visualization
selLambdaPos = 1;%28;%18;%1;%16;%14;
regularizationVisualization;
% %% (2.3) Saving data after regularization
save_regularization


%% (2.3) Studying confidence intervals vs regularization
plot_ParsCIs_lambdalist;
save_pars_ci


%% (3.1) Study on paarameter values: estimation with the lambda value
% % select lambda
% setup.selectedLambda = lambdalist(selLambdaPos);
% % estimate parameters
% tic
% [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t = toc;
% % simulate for error
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;
% [error] = optfun(xres,data,setup);
%  % confidence intervals
% N = length(error);
% Jacobian = full(Jacobian);  
% varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
% stdp = sqrt(diag(varp));

% recall specific values
setup.selectedLambda = lambdalist(selLambdaPos);
error = array_raw_error{selLambdaPos};
N = cell_N{selLambdaPos};
Jacobian = cell_Jacobian{selLambdaPos};
varp = cell_varp{selLambdaPos};
stdp = cell_stdp{selLambdaPos};



%% ENZYME-SPECIFIC
%% (3.2) Study on paarameter values: recalculation
% parameter values
% xres_selected = xres; %lambdalist based in 'ones', lam=0.1, loc=5.
% xres_selected = array_xres{15}; %lambdalist based in 'ones', lam=0.1, loc=5.
xres_selected = array_xres{selLambdaPos};

kdhap = ones(numpHtested,1);
kgap = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    kdhap(i) = 6.45 * 10.^xres_selected(1); %mM
    kgap(i) = 5.25 * 10.^xres_selected(2); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+2);
end
keq_tpi = setup.Keq_TPI;
keq_gpd = setup.Keq_GPD;

% limits
kdhap_up = ones(numpHtested,1);
kgap_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);
kdhap_down = ones(numpHtested,1);
kgap_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    %up
    kdhap_up(i) = 0.08 * 10.^(xres_selected(1)+stdp(1)); %mM
    kgap_up(i) = 1.2 * 10.^(xres_selected(2)+stdp(2)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)+stdp(i+2));
    %down
    kdhap_down(i) = 0.08 * 10.^(xres_selected(1)-stdp(1)); %mM
    kgap_down(i) = 1.2 * 10.^(xres_selected(2)-stdp(2)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)-stdp(i+2));
end

vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;
% % % % % Test resulting value
% % % % vm = data.Vmax(:,4) .* (10 .^ xres(4:11)') .* 60 .* 60 ./ 1.7781;
% % % % vm = mean([data.Vmax(:,1)*8 data.Vmax(:,2)*4],2) .* (10 .^ zeros(8,1)) .* 60 .* 60 ./ 1.7781;


%% (3.3) Study on paarameter values: plotting
if(setup.plotOutput == 1)
    figure(105)

    subplot(331) % vm
        plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_uChange,'.-','color','black')
    title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

    subplot(332) % vm experimental
        errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
        title({'v_{m} experimental';'not normalized'})

    subplot(333) % vm(log_difference)
        plot(pHarray,xres_selected(3:9)','.-','color','black')
        title('v_{m} difference []')

    subplot(334) % kdhap
        plot(pHarray,kdhap_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kdhap_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kdhap,'.-','color','black')
    title('k_{dhap} [mM]')

    subplot(335) % kgap
        plot(pHarray,kgap_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kgap_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kgap,'.-','color','black')
    title('k_{gap} [mM]')

    subplot(337) % keq_tpi
    plot(setup.fullpHarray,keq_tpi,'.-','color','black')
    title('k_{eq.TPI} [mM]')

    subplot(338) % keq_gpd
    plot(setup.fullpHarray,keq_gpd,'.-','color','black')
    title('k_{eq.ENO} [mM]')

    suptitle('TPI: parameter estimates vs pH')
end

%%
output_tpi.xres_selected = xres_selected;

output_tpi.pHarray = pHarray;
output_tpi.idxs2consider = idxs2consider;
output_tpi.DF = data.DF;

output_tpi.kdhap = kdhap;% = ones(numpHtested,1);
output_tpi.kgap = kgap;% = ones(numpHtested,1);
output_tpi.vm = vm;% = zeros(numpHtested,1);
output_tpi.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_tpi.kdhap_up = kdhap_up;% = ones(numpHtested,1);
output_tpi.kgap_up = kgap_up;% = ones(numpHtested,1);
output_tpi.vm_up = vm_up;% = zeros(numpHtested,1);
output_tpi.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_tpi.kdhap_down = kgap_down;% = ones(numpHtested,1);
output_tpi.kgap_down = kgap_down;% = ones(numpHtested,1);
output_tpi.vm_down = vm_down;% = zeros(numpHtested,1);
output_tpi.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_tpi.keq_tpi = keq_tpi;% = setup.Keq_TPI;
output_tpi.keq_gpd = keq_gpd;% = setup.Keq_GPD;


%% saving output
output_tpi.Vmax_experimental = Vmax_experimental;
output_tpi.stDev_experimental = stDev_experimental;
save_output

% %% reverse to reload
% load('output_tpi.mat','output_tpi');
% 
% xres_selected = output_tpi.xres_selected;
% 
% pHarray = output_tpi.pHarray;
% 
% kfbp = output_tpi.kfbp;%
% kgap = output_tpi.kgap;%
% kdhap = output_tpi.kdhap;%
% vm = output_tpi.vm;%
% 
% kfbp_up = output_tpi.kfbp_up;%
% kgap_up = output_tpi.kgap_up;%
% kdhap_up = output_tpi.kdhap_up;%
% vm_up = output_tpi.vm_up;%
% 
% kfbp_down = output_tpi.kfbp_down;%
% kgap_down = output_tpi.kgap_down;%
% kdhap_down = output_tpi.kdhap_down;%
% vm_down = output_tpi.vm_down;%
% 
% keq_fba = output_tpi.keq_fba;%
% keq_tpi = output_tpi.keq_tpi;%
% keq_gpd = output_tpi.keq_gpd;%
% % % % x = 1; end