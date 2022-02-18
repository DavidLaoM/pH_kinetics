% % % % function x=dev_pEst_gapdh_rev
% % PEST_ALD.m
% Parameter estimation for the data in the Enolase assay.
% Vm are estimated changing with pH
% Kms are assumed constant
% Keq from the eQuilibrator


%% (0) Setup and data load
clear, close all
set_paths_pHstudy;
dbstop if error
for step0 = 1
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
    
%     figure
%     plot(pHvals, Vmax(:,4),'.-')
%     title('Starting estimate: Vmax [mM s-1] vs pH')    
end


%% (0.1) Experimental Vmax determination

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

% Experimental rates determination and plotting
expRatesDetermination;


%% (1.1) Simple parameter fit. Parameter estimation
% % % % setup.ode = 'vanHeerden2014';
setup.ode = 'gapdh_rev_simplified';
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
% setup.weightDataEsp = ones(1,numpHtested);
setup.weightDataEsp = idxs2consider;
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% Km fixed
optfun = @costfun_Kmfixed;
% plength = 10; % Kms (2) + Vms (1) * numpH (12) % overly simplified 2020-09-30
plength = 14; % Kms (2) + Vms (1) * numpH (12) % still keeping kms
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');


% % % % %% PSA on the effect of Keq for gapdh_reverse (run all of them)
% % % % setup.KeqPSA = 1;
% % % % x_temp = [-0.0284   -1.4314   -0.0060   -0.1678    0.3818    0.3831...
% % % %            0.3885    0.4256    0.3735    0.3420    0.3257    0.3188...
% % % %            0.2959    0.3072];
% % % % setup.plotEachSimCF = 0;
% % % % setup.plotEachSim = 0;
% % % % setup.simAllProfiles = 0;
% % % % 
% % % % nPSAvals = length(setup.pH_Keq_gapdh);
% % % % simResults = cell(1,nPSAvals);
% % % % for i = 1:nPSAvals
% % % %     setup.nPSA = i;
% % % %     [error,simResult] = costfun_Kmfixed_simOutput(x_temp,data,setup);
% % % %     simResults{i} = simResult;
% % % % end
% % % % 
% % % % numpH = setup.numpHtested;
% % % % for j = 1:numpH   
% % % %     simulationVisualization_PSAkeq;
% % % % end
% % % % 
% % % % % if setup.plotOutput == 1
% % % % %     set(1,'color','white'), savefig(1,['results/',setup.enzymeName,'/',setup.enzymeName, '_concentrations_basezero.fig']);
% % % % %     set(2,'color','white'), savefig(2,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_df.fig']);
% % % % % end


%% PSA on the effect of Keq for gapdh_reverse (simple run)
setup.KeqPSA = 1;

setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

for simulation = 1
    nPSAvals = length(setup.pH_Keq_gapdh);
    simResults = cell(1,nPSAvals);
    for ix = 1:nPSAvals
        setup.nPSA = ix;

        DFs = 4; numpH = 10; obsMet = 8;
        simNADH = cell(DFs,numpH);
        expNADH = cell(DFs,numpH);
        simGAPDHr = cell(DFs,numpH);
        expGAPDHr = cell(DFs,numpH);

        % simulations loop for each pH value
        for j = 5
            % select required data
            data.NADH = data.conc_mean(j,:);
            data.tempTime = data.time(j,:);

            % added in case of the Keq_PSA
            eQvals = [1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1];
            setup.pH_Keq_gapdh_eQ = 0.05 * (eQvals);
            % Keq_pgk obtained with eQuilibrator
            setup.pH_Keq_pgk = [1/(7.4E-4),    1/(7.2E-4),    1/(6.9E-4),    1/(6.5E-4),    1/(6.1E-4),    1/(5.7E-4),    1/(5.5E-4),    1/(5.3E-4),    1/(5.2E-4),   1/(5.2E-4)];

            if isfield(setup,'KeqPSA')
                if setup.KeqPSA == 1
                    nPSA = setup.nPSA;
                    % inputs to be selected
                    data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(nPSA);
                    data.chosenKeqPGK = setup.pH_Keq_pgk(nPSA);
                end
            end
            
            % simulations
            for i = 4

                % still keeping kms
                odefun = @odeGAPDHrev_simplified;

                p = struct;
                p.TDH1_Kgap = 2.3230; %10.^-0.0284.*2.48; % mM
                p.TDH1_Kbpg = 0.0437; %10.^-1.4314.*1.18; % mM
                p.TDH1_Knad = 2.8799; %10.^-0.0060.*2.92; %mM
                p.TDH1_Knadh = 0.0149; %10.^-0.1678.*0.022; % mM
                p.TDH1_Vmr = 0.0038;%0.0016 * 10 .^ 0.3735; %(UNIT!) 
                p.GAPDH_Keq = data.chosenKeqGAPDH;
                p.PGK_Vm = 1306.45 / 60;
                p.PGK_Keq = data.chosenKeqPGK;

                
                % set initial conditions and timespan
                P3Go = 5;
                ATPo = 1;
                BPGo = 0;
                ADPo = 0;
                NADo = 0;
                GAPo = 0;
                PHOSo = 500;
                NADHo = 0.0924;
                tspan   = [0 300];  % time [300 s (6 min)]
                y0 = [P3Go ATPo BPGo ADPo NADo GAPo PHOSo NADHo];
                options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));

                % run ode15s
                f = 1; % not being used for PSA
                [t,y] = ode15s(odefun,tspan,y0,options,p,f,data,setup);

                simResult.t = t;
                simResult.y = y;

                % cost function (NADHexp + vGAPDHr)
                simTime = simResult.t;
                simMet = simResult.y(:,obsMet);

                simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                expNADH{i,j} = data.NADH{i};
            end

       end
        % save simResult
        simResult.simNADH = simNADH;
        simResult.expNADH = expNADH;
        % put them all together
        simResults{ix} = simResult;
        
    end
end

% done
for plotting = 1
    numpH = setup.numpHtested;
    for j = 5
        data.tempTime = data.time(j,:);
        
        h101 = figure(101);

        set(0,'CurrentFigure', h101);
        p101 = subplot(3,4,j);
%         p101 = subplot(1,1,1);
        % % % % % Added for the Keq
        colArray = jet(length(simResults));
        for k = 1:length(simResults)
            % load back simResult
        %     simNADH = simResults{k}.simNADH';
        %     expNADH = simResults{k}.expNADH';
            simNADH = simResults{k}.simNADH;
            expNADH = simResults{k}.expNADH;
            % 
            DFstudy = [1 2 3 4];
            data.tempTime = data.time(j,:);

            for i = 4
                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',colArray(k,:))
                    hold on
                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                    ylim([0 0.15])
            end
            
            % ylabel left
            ylabel('NADH concentration [mM]')        
            % xlabel
            xlabel('assay time [s]')
            % pH text box
            tempText = erase(sprintf('pH %d', data.pH(j,1)),"0000e+00");
            text(30, p101.YLim(2)*0.9, tempText);
        end
    end
end

