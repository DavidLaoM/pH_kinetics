%% (0) Setup
clear
set_paths_pHstudy;
dbstop if error

%% (0.1) Setup and data load: gapdh_forward
for reducingVisualLines = 1
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
    selectSetup_pH;
    setup.saveOutput = 0;

    load('expData.mat','expData');
    import_gapdhF_temp = expData.gapdh;
    % avoid extra 2 pH values
    import_gapdhF.rawData = import_gapdhF_temp.rawData;
    
    import_gapdhF.treatedData.excelSheet_corrected = import_gapdhF_temp.treatedData.excelSheet_corrected(:,1:20);
    import_gapdhF.treatedData.pH_corrected = import_gapdhF_temp.treatedData.pH_corrected(:,1:20);
    import_gapdhF.treatedData.dilution_corrected = import_gapdhF_temp.treatedData.dilution_corrected(:,1:20);
    import_gapdhF.treatedData.absorbance_mean = import_gapdhF_temp.treatedData.absorbance_mean(:,1:20);
    import_gapdhF.treatedData.absorbance_samples = import_gapdhF_temp.treatedData.absorbance_samples(:,1:20);
    import_gapdhF.treatedData.absorbance_std = import_gapdhF_temp.treatedData.absorbance_std(:,1:20);
    import_gapdhF.treatedData.concentration_mean = import_gapdhF_temp.treatedData.concentration_mean(:,1:20);
    import_gapdhF.treatedData.concentration_std = import_gapdhF_temp.treatedData.concentration_std(:,1:20);
    import_gapdhF.treatedData.time = import_gapdhF_temp.treatedData.time(:,1:20);
    import_gapdhF.treatedData.reaction_rate = import_gapdhF_temp.treatedData.reaction_rate(:,1:20);
    import_gapdhF.treatedData.unitsTime = import_gapdhF_temp.treatedData.unitsTime;
    import_gapdhF.treatedData.unitsAbsorbance = import_gapdhF_temp.treatedData.unitsAbsorbance;
    import_gapdhF.treatedData.unitsConcentration = import_gapdhF_temp.treatedData.unitsConcentration;
    import_gapdhF.treatedData.unitsRates = import_gapdhF_temp.treatedData.unitsRates;
    import_gapdhF.treatedData.concProtein = import_gapdhF_temp.treatedData.concProtein;
    import_gapdhF.treatedData.unitsProtein = import_gapdhF_temp.treatedData.unitsProtein;
    import_gapdhF.treatedData.protDF = import_gapdhF_temp.treatedData.protDF;
    
%     DFs = setup.DFactorsTotal;
    DFs = 4;
%     pHtested = setup.pHtested;
    pHtested = [1 1 0 1 1 1 1 1 1 1 1 0];
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

    pHarray = unique(import_gapdhF.treatedData.pH_corrected);
    for i = 1:numpHtested
%         disp(i);
        pHval = pHarray(i);
        tempID = find(import_gapdhF.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_gapdhF.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_gapdhF.treatedData.dilution_corrected(tempID);
        for j = 1:4
%             disp(j);
            abs_meanTemp{j,i} = import_gapdhF.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_gapdhF.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_gapdhF.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_gapdhF.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_gapdhF.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_gapdhF.treatedData.reaction_rate{tempID(j)};
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
        Vmax(i) = max(abs(RRs{i}));
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
    data_gapdhF.pH = pH;
    data_gapdhF.DF = DF;
    data_gapdhF.abs_mean = abs_mean;
    data_gapdhF.abs_std = abs_std;
        conc_mean{1,1} = conc_mean{1,1} * 2;
        conc_mean{1,2} = conc_mean{1,2} * 2;
        conc_mean{1,3} = conc_mean{1,3} * 2;
        conc_mean{1,4} = conc_mean{1,4} * 2;
        conc_mean{2,1} = conc_mean{2,1} * 2;
        conc_mean{2,2} = conc_mean{2,2} * 2;
        conc_mean{2,3} = conc_mean{2,3} * 2;
        conc_mean{2,4} = conc_mean{2,4} * 2;
    data_gapdhF.conc_mean = conc_mean;              % 2
        conc_std{1,1} = conc_std{1,1} * 2;
        conc_std{1,2} = conc_std{1,2} * 2;
        conc_std{1,3} = conc_std{1,3} * 2;
        conc_std{1,4} = conc_std{1,4} * 2;
        conc_std{2,1} = conc_std{2,1} * 2;
        conc_std{2,2} = conc_std{2,2} * 2;
        conc_std{2,3} = conc_std{2,3} * 2;
        conc_std{2,4} = conc_std{2,4} * 2;
    data_gapdhF.conc_std = conc_std;                % 2
    data_gapdhF.time = time;
        RRs{1,1} = RRs{1,1} * 2;
        RRs{1,2} = RRs{1,2} * 2;
        RRs{1,3} = RRs{1,3} * 2;
        RRs{1,4} = RRs{1,4} * 2;
        RRs{2,1} = RRs{2,1} * 2;
        RRs{2,2} = RRs{2,2} * 2;
        RRs{2,3} = RRs{2,3} * 2;
        RRs{2,4} = RRs{2,4} * 2;
    data_gapdhF.RRs = RRs;                          % 2
        Vmax(1,1) = Vmax(1,1) * 2;
        Vmax(1,2) = Vmax(1,2) * 2;
        Vmax(1,3) = Vmax(1,3) * 2;
        Vmax(1,4) = Vmax(1,4) * 2;
        Vmax(2,1) = Vmax(2,1) * 2;
        Vmax(2,2) = Vmax(2,2) * 2;
        Vmax(2,3) = Vmax(2,3) * 2;
        Vmax(2,4) = Vmax(2,4) * 2;
    data_gapdhF.Vmax = Vmax;                        % 2
    data_gapdhF.chosenVmax = max(max(Vmax));
    data_gapdhF.chosenNADini = 0.05;
    temp1 = import_gapdhF.rawData.absorbance_corrected{4,4};
    temp2 = import_gapdhF.rawData.absorbance_corrected{5,4};
    temp3 = import_gapdhF.rawData.absorbance_corrected{6,4};
    data_gapdhF.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data_gapdhF.raw.time = import_gapdhF.rawData.time{1};

    pHvals = unique(import_gapdhF.treatedData.pH_corrected);
    % visualize: check calculations made
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numpHtested
        subplot(3,4,i)
        for j = 1:DFs
            plot(time{i,j},NADH{i,j},'.-')
            hold on
        end
        title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
        if((i == 1)||(i == 2))
            text(200,0.039,'DF.2,4,8,16');
        end
        xlim([0 600])
        ylim([0 0.2])
        if i == numpHtested
            if setup.caseStudyGAPDHr == 1
                legend('DF 8','DF 4','DF 2','DF 1')
            end
        end
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
    suptitle(suptitleName);
end

%% (0.2) Setup and data load: gapdh_reverse
for reducingVisualLines = 1
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
    selectSetup_pH;
    % added
    setup.saveOutput = 0;

    load('expData.mat','expData');
    import_gapdhR = expData.gapdhr;

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

    pHarray = unique(import_gapdhR.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_gapdhR.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_gapdhR.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_gapdhR.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_gapdhR.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_gapdhR.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_gapdhR.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_gapdhR.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_gapdhR.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_gapdhR.treatedData.reaction_rate{tempID(j)};
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
        Vmax(i) = max(abs(RRs{i}));
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
    data_gapdhR.pH = pH;
    data_gapdhR.DF = DF;
    data_gapdhR.abs_mean = abs_mean;
    data_gapdhR.abs_std = abs_std;
    data_gapdhR.conc_mean = conc_mean;
    data_gapdhR.conc_std = conc_std;
    data_gapdhR.time = time;
    data_gapdhR.RRs = RRs;
    data_gapdhR.Vmax = Vmax;
    data_gapdhR.chosenVmax = max(max(Vmax));
    data_gapdhR.chosenNADini = 0.15;
    temp1 = import_gapdhR.rawData.absorbance_corrected{4,4};
    temp2 = import_gapdhR.rawData.absorbance_corrected{5,4};
    temp3 = import_gapdhR.rawData.absorbance_corrected{6,4};
    data_gapdhR.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data_gapdhR.raw.time = import_gapdhR.rawData.time{1};

    pHvals = unique(import_gapdhR.treatedData.pH_corrected);
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
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
    suptitle(suptitleName);
end

%% simple simulations
setup.params =  {'v_{maxFWD}^{REVdata} [mM s^{-1}]'; 'K_{gap} [mM]'; 'K_{bpg} [mM]'; 'K_{nad} [mM]'; 'K_{nadh} [mM]'; 'v_{maxREV}^{REVdata} [mM s^{-1}]'; 'v_{maxFWD}^{FWDdata} [mM s^{-1}]'; 'v_{maxREV}^{FWDdata} [mM s^{-1}]'};
xtemp = zeros(8,1);

% Keq    
pH_Keq_gapdh = 0.05 *[1.7E-3, 2.3E-3, 4.0E-3, 8.0E-3, 1.62E-2, 3.46E-2, 6.56E-2, 1.16E-1, 1.78E-1, 2.45E-1];
pH_Keq_pgk = [1/(7.4E-4),    1/(7.2E-4),    1/(6.9E-4),    1/(6.5E-4),    1/(6.1E-4),    1/(5.7E-4),    1/(5.5E-4),    1/(5.3E-4),    1/(5.2E-4),   1/(5.2E-4)];
simRes = cell(1,10);
iter = 1:10;
iter = 6;
for i = iter
    setup.chosen_keq_gapdh = pH_Keq_gapdh(i);
    setup.chosen_keq_pgk = pH_Keq_pgk(i);
    [testSim] = simSys_temp(xtemp,data,setup);
    simRes{i} = testSim;
end

figure(1)
for i = iter
    t = simRes{i}.t;
    y = simRes{i}.y;
    plot(t,y(:,8),'b','Linewidth',1.2)
    hold on
end
ylabel('NADH [mM]')
xlabel('time [s]')
suptitle('Progression curve: NADH vs time')

function [simResult] = simSys_temp(xtemp,data,setup)
    % set parameters
    % GAPDH   
    p.TDH1_Vm = 10.^xtemp(1).*30.98/1000; % mM/s
    p.TDH1_Kgap = 10.^xtemp(2).*2.48; % mM
    p.TDH1_Kbpg = 10.^xtemp(3).*1.18; % mM
    p.TDH1_Knad = 10.^xtemp(4).*2.92; %mM
    p.TDH1_Knadh = 10.^xtemp(5).*0.022; % mM
%     p.TDH1_Keq = 10.^xtemp(6).*0.00173; % []
% % % %     p.TDH1_Keq = 10.^xtemp(6).*0.00149; %/0.05; % []
    % p.TDH1_Keq = data.chosenKeqGAPDH; % []
    p.TDH1_Keq = setup.chosen_keq_gapdh;
    % PGK
    p.PGK_Katp = 0.3; % mM
    p.PGK_Kp3g = 0.53; % mM
    p.PGK_Kbpg = 0.003; % mM
    p.PGK_Kadp = 0.2; % mM
%     p.PGK_Keq = 1754.386; % [] %/10
% % % %     p.PGK_Keq = 1724; %3200; %1724; % [] %/10
    p.PGK_Keq = setup.chosen_keq_pgk;
%     p.PGK_Keq = 1000;
    p.PGK_Vm = 44.5; %1306.45; % mM/s
    % p.PGK_Vm = 1306.45 / 60 * setup.excessPGK / data.chosenDF; % mM s^{-1} % corrected to make it appear in excess
        
    % set initial conditions, concentrations and timespan 
    P3Go = 5; %*1e-3;% [mM]
    ATPo = 1; %*1e-3;
    BPGo = 0; %*1e-3;
    ADPo = 0; %*1e-3;
    NADo = 0; %*1e-3;
    GAPo = 0; %*1e-3;
    PHOSo = 50; %*1e-3;
    NADHo = 0.15; %*1e-3;
%     % test PGK
%     P3Go = 0; %*1e-3;% [mM]
%     ATPo = 0; %*1e-3;
%     BPGo = 1; %*1e-3;
%     ADPo = 1; %*1e-3;
%     NADo = 0; %*1e-3;
%     GAPo = 0; %*1e-3;
%     PHOSo = 0; %*1e-3;
%     NADHo = 0; %*1e-3;    
    
    
    y0 = [P3Go ATPo BPGo ADPo NADo GAPo PHOSo NADHo];  
    tspan = [0 300];
    options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(y0)));
    % run ode15s
    f = 1; % not being used for PSA
    
    % simulations
    [t,y] = ode15s(@ode_temp,tspan,y0,options,p,f,data,setup);
    simResult.t = t;
    simResult.y = y;

    % plotting
    figure(1)
    plot(t,y(:,8),'b','Linewidth',1.2)
    ylabel('NADH [mM]')
    xlabel('time [s]')
    suptitle('Progression curve: NADH vs time')
    
%     figure
%     for i = 1:n12
%         subplot(3,3,i)
%         plot(t,y(:,i),'b')
%         title(setup.PSAmets{i})
%         hold on
%     end
%     suptitle('System simulation: Reverse direction')
    
%     % test PGK
%     titleName = {'P3G','ATP','BPG','ADP'};
%     figure
%     for i = 1:4
%         subplot(2,2,i)
%         plot(t,y(:,i),'b','Linewidth',1.2)
%         title(titleName{i})
%     end

end 

function [v] = ode_temp(tspan,y0,p,f,data,setup)
ode_pH = setup.ode_pH;
typeVm = setup.typeVm;

% select initial points
% reverse
P3G = y0(1);
ATP = y0(2);
BPG = y0(3);
ADP = y0(4);
NAD = y0(5);
GAP = y0(6);
PHOS = y0(7);
NADH = y0(8);

% rateEquations
% v_GAPDH = p.TDH1_Vm.*(GAP.*NAD.* PHOS - BPG .* NADH ./ p.TDH1_Keq);
v_GAPDH = p.TDH1_Vm.*(GAP .* NAD - BPG .* NADH ./ p.TDH1_Keq);
v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq);
% v_PGK = p.PGK_Vm .* (BPG .* ADP .* p.PGK_Keq - P3G .* ATP);
% count Pi effect?
% v_GAPDH = ((p.TDH1_Vm./(p.TDH1_Kgap.*p.TDH1_Knad)).*(GAP.*NAD.*PHOS-(BPG.*NADH)./p.TDH1_Keq))./...
%     ((1+GAP./p.TDH1_Kgap).*(1+NAD./p.TDH1_Knad)+(1+BPG./p.TDH1_Kbpg).*(1+NADH./p.TDH1_Knadh)-1);
% v_PGK = p.PGK_VmPGK.*((p.PGK_Keq.*BPG.*ADP)-ATP.*P3G)./...
%       (p.PGK_Katp*p.PGK_Kp3g.*(1 + ADP./p.PGK_Kadp + ATP./p.PGK_Katp).*(1 + BPG/p.PGK_Kbpg + P3G/p.PGK_Kp3g));
% v_GAPDH = p.TDH1_Vmf_R .* ( rNAD .* rGAP - p.TDH1_Vmr_R .* rBPG .* rNADH)./((p.TDH1_Kgap .* p.TDH1_Knad).*((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap)));
% v_PGK = p.PGK_Vm .* (rBPG .* rADP - rP3G .* rATP ./ p.PGK_Keq);
% v_GAPDH = 0;

% mass balances (assumed ideally mixed batch system configuration)
% reverse
v(1) = + v_PGK; %P3G
v(2) = + v_PGK; %ATP
v(3) = - v_PGK + v_GAPDH; %BPG
v(4) = - v_PGK; %ADP
v(5) = - v_GAPDH; %NAD
v(6) = - v_GAPDH; %GAP
v(7) = - v_GAPDH; %PHOS
v(8) = + v_GAPDH; % NADH
v=v';

end


%%




