% % PEST_PDC.m

%% (0) Setup and data load
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
%     setup.saveOutput = 0;
    
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


%% (0.1) Experimental Vmax determination
% setup.plotOutput = 0;
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

% %% (0.3) saving
save_initial_recall


%% (1.1) Simple parameter fit. Parameter estimation
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
% testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(x_temp,data,setup);
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
parameterEstimation_lambdalist;


%% (2.2) Regularization. Results Visualization
selLambdaPos = 1;%15;%1;%17;%14;
regularizationVisualization;
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

kpyr = ones(numpHtested,1);
hill = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    kpyr(i) = 8.5 * 10.^xres_selected(1); %mM
    hill(i) = 1.9 * 10.^xres_selected(2); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+2);
end
keq_adh = setup.Keq_ADH;

% limits
kpyr_up = ones(numpHtested,1);
hill_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);
kpyr_down = ones(numpHtested,1);
hill_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    %up
    kpyr_up(i) = 8.5 * 10.^(xres_selected(1)+stdp(1)); %mM
    hill_up(i) = 1.9 * 10.^(xres_selected(2)+stdp(2)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)+stdp(i+2));
    %down
    kpyr_down(i) = 8.5 * 10.^(xres_selected(1)-stdp(1)); %mM
    hill_down(i) = 1.9 * 10.^(xres_selected(2)-stdp(2)); %mM
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
    % figure(106)

    subplot(331) % vm
        plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_uChange,'.-','color','black')
    title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

    subplot(332) % vm experimental
        errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
        title({'v_{m} experimental';'not normalized'})

    subplot(333) % vm(log_difference)
        plot(pHarray,xres_selected(3:14)','.-','color','black')
        title('v_{m} difference []')

    subplot(334) % kpyr
        plot(pHarray,kpyr_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kpyr_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kpyr,'.-','color','black')
    title('k_{pyr} [mM]')

    subplot(335) % hill
        plot(pHarray,hill_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,hill_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,hill,'.-','color','black')
    title('hill [mM]')

    subplot(337) % keq_adh
    plot(pHarray,keq_adh,'.-','color','black')
    title('k_{eq.ADH} [mM]')

    suptitle('PDC: parameter estimates vs pH')
end

%%
output_pdc.xres_selected = xres_selected;

output_pdc.pHarray = pHarray;
output_pdc.idxs2consider = idxs2consider;
output_pdc.DF = data.DF;

output_pdc.kpyr = kpyr;% = ones(numpHtested,1);
output_pdc.hill = hill;% = ones(numpHtested,1);
output_pdc.vm = vm;% = zeros(numpHtested,1);
output_pdc.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_pdc.kpyr_up = kpyr_up;% = ones(numpHtested,1);
output_pdc.hill_up = hill_up;% = ones(numpHtested,1);
output_pdc.vm_up = vm_up;% = zeros(numpHtested,1);
output_pdc.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_pdc.kpyr_down = hill_down;% = ones(numpHtested,1);
output_pdc.hill_down = hill_down;% = ones(numpHtested,1);
output_pdc.vm_down = vm_down;% = zeros(numpHtested,1);
output_pdc.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_pdc.keq_adh = keq_adh;% = setup.Keq_ADH;


%% saving output
output_pdc.Vmax_experimental = Vmax_experimental;
output_pdc.stDev_experimental = stDev_experimental;
save_output

% %% reverse to reload
% load('output_pdc.mat','output_pdc');
% 
% xres_selected = output_pdc.xres_selected;
% 
% pHarray = output_pdc.pHarray;
% 
% kfbp = output_pdc.kfbp;%
% hill = output_pdc.hill;%
% kdhap = output_pdc.kdhap;%
% vm = output_pdc.vm;%
% 
% kfbp_up = output_pdc.kfbp_up;%
% hill_up = output_pdc.hill_up;%
% kdhap_up = output_pdc.kdhap_up;%
% vm_up = output_pdc.vm_up;%
% 
% kfbp_down = output_pdc.kfbp_down;%
% hill_down = output_pdc.hill_down;%
% kdhap_down = output_pdc.kdhap_down;%
% vm_down = output_pdc.vm_down;%
% 
% keq_fba = output_pdc.keq_fba;%
% keq_tpi = output_pdc.keq_tpi;%
% keq_gpd = output_pdc.keq_gpd;%

% % % % x = 1; end


