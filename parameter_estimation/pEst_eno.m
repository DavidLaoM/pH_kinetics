% % PEST_ENO.m

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
    setup.caseStudyTPI = 0;
    setup.caseStudyENO = 1;
    selectSetup_pH;
%     % added
%     setup.saveOutput = 0;

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
    if(setup.plotOutput == 1)
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
end


%% (0.1) Experimental Vmax determination
% setup.plotOutput = 0;
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
% % % % setup.weightDataEsp = ones(1,12);
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
% %%
% % testing the costfunction
% setup.plotEachSimCF = 1;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 1;
% [error] = optfun(x_temp,data,setup);
% setup.plotEachSimCF = 0;
% setup.plotEachSim = 0;
% setup.simAllProfiles = 0;
% %%
% parameter estimation
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
selLambdaPos = 1; %12;
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
% xres_selected = array_xres{19}; % lambda value = 1
% xres_selected = xres; % lambda value = 1
xres_selected = array_xres{selLambdaPos};

k2pg = ones(numpHtested,1);
kpep = ones(numpHtested,1);
vm = zeros(numpHtested,1);
for i = 1:numpHtested
    k2pg(i) = 0.043 * 10.^xres_selected(1); %mM
    kpep(i) = 0.5 * 10.^xres_selected(2); %mM
    vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+2);
end
keq = setup.keq;
% limits
k2pg_up = ones(numpHtested,1);
kpep_up = ones(numpHtested,1);
vm_up = zeros(numpHtested,1);
k2pg_down = ones(numpHtested,1);
kpep_down = ones(numpHtested,1);
vm_down = zeros(numpHtested,1);
for i = 1:numpHtested
    k2pg_up(i) = 0.043 * 10.^(xres_selected(1)+stdp(1)); %mM
    kpep_up(i) = 0.5 * 10.^(xres_selected(2)+stdp(2)); %mM
    vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)+stdp(i+2)); %mM s^[-1]
    k2pg_down(i) = 0.043 * 10.^(xres_selected(1)-stdp(1)); %mM
    kpep_down(i) = 0.5 * 10.^(xres_selected(2)-stdp(2)); %mM
    vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+2)-stdp(i+2)); %mM s^[-1]
end


%% (3.3) Study on paarameter values: plotting
vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;

if(setup.plotOutput == 1)
    figure(105)

    subplot(231) % vm
    %     plot(pHarray,vm_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    %     plot(pHarray,vm_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    % plot(pHarray,vm,'.-','color','black')
    % title('v_{m} [mM s^{1}]')
        plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,vm_uChange,'.-','color','black')
    title({'v_{m} [umol_{PEP} mg_{P}^{-1} min^{-1}]';'not normalized'})

    subplot(232) % vm experimental
        errorbar(pHarray,Vmax_experimental,stDev_experimental,'o-')
        title({'v_{m} experimental';'not normalized'})

    subplot(233) % vm(log_difference)
        plot(pHarray,xres_selected(3:14)','.-','color','black')
        title('v_{m} difference []')

    subplot(234) % k2pg
        plot(pHarray,k2pg_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,k2pg_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,k2pg,'.-','color','black')
    title('k_{p2g} [mM]')

    subplot(235) % kpep
        plot(pHarray,kpep_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kpep_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(pHarray,kpep,'.-','color','black')
    title('k_{pep} [mM]')

    subplot(236) % keq
    plot(pHarray,keq,'.-','color','black')
    title('k_{eq} [mM]')

    suptitle('ENO: parameter estimates vs pH')
end
% %%
% set(105,'color','w');
% savefig(105,'ENO_pVals_pCis');
%%
output_eno.xres_selected = xres_selected;

output_eno.pHarray = pHarray;
output_eno.idxs2consider = idxs2consider;
output_eno.DF = data.DF;

output_eno.k2pg = k2pg;% = ones(numpHtested,1);
output_eno.kpep = kpep;% = ones(numpHtested,1);
output_eno.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

output_eno.k2pg_up = k2pg_up;% = ones(numpHtested,1);
output_eno.kpep_up = kpep_up;% = ones(numpHtested,1);
output_eno.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

output_eno.k2pg_down = k2pg_down;% = ones(numpHtested,1);
output_eno.kpep_down = kpep_down;% = ones(numpHtested,1);
output_eno.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

output_eno.keq_eno = keq;% = setup.Keq_PGI;


%% saving output
output_eno.Vmax_experimental = Vmax_experimental;
output_eno.stDev_experimental = stDev_experimental;
save_output

