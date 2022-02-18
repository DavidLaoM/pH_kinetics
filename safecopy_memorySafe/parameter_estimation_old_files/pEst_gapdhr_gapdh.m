% % PEST_GAPDHR_GAPDH.M 
% Parameter estimation for GAPDH was tested with the gapdh_r dataset in the
% code 'pEst_gapdhr_updatedTopology'. In this code, the exercise is 
% supplemented with data for the GAPDH (forward sense) reaction.

% The idea is to see if the parameters can be well-quantified with both
% data sets at the same time and see if this adds information to the
% problem. In principle, the Haldane constraint will not be considered (the
% weight will be set to zero) since the Keq approach was not solved in the
% previous assay, and needs to be further discussed.

% Sections
% (0) Setup: The data sets to be used are loaded.
% (1) Test for simulation file
% (2) Test for parameter the cost function
% (3) Parameter estimation: no regularization
% (4) Parameter estimation: added regularization


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


%% (1.1) Simple system simulation @pH7
data.gapdhF = data_gapdhF;
data.gapdhR = data_gapdhR;
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode_pH = 'on';
setup.typeVm = 'common';

xtemp = zeros(6,1);
data.chosenDF = 1;
data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(6);
data.chosenKeqPGK = setup.pH_Keq_pgk(6);
data.i = 6;
setup.excessPGK = 1;
setup.plotEachSim = 1;
[testSim] = simSys_FandR(xtemp,data,setup);


%% (1.2) PSA: NADH, vappPGK and vappGAPDH plotted
data.gapdhF = data_gapdhF;
data.gapdhR = data_gapdhR;
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode_pH = 'on';
setup.typeVm = 'common';

n = 21;
% xval = linspace(-1,1,n);
xval = linspace(-3,3,n);
c = cool(n);
xobs = 8;
np = 6;

allSims2 = cell(n,np);
for k = 1:np
    for i = 1:n
        % inputs before simSys required before running
        xtemp = zeros(6,1);
        xtemp(k) = xval(i);
        data.chosenDF = 1; % reference: DF = 1
        data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(6); %setup.pH_Keq_gapdh(6); % reference from pH = 7.06
        data.chosenKeqPGK = setup.pH_Keq_pgk(6); % reference from pH = 7.06
        data.i = 6;
        setup.excessPGK = 1; % excess of PGK protein for faster linking reaction
        setup.plotEachSim = 0;

        % simulation
        [tempSim] = simSys_FandR(xtemp,data,setup);
        
        % calculate reaction rates in the system
        % GAPDH
        p.TDH1_Vmf = 10.^xtemp(1).*1184.52/60 / data.chosenDF;% mM s^{-1}        
        p.TDH1_Kgap = 10.^xtemp(2).*2.48; % mM
        p.TDH1_Kbpg = 10.^xtemp(3).*1.18; % mM
        p.TDH1_Knad = 10.^xtemp(4).*2.92; %mM
        p.TDH1_Knadh = 10.^xtemp(5).*0.022; % mM
        p.TDH1_Vmr = 10.^xtemp(6).*6549.8/60 / data.chosenDF; % mM s^{-1}
        p.TDH1_Keq = data.chosenKeqGAPDH; % []
        % PGK
        p.PGK_Keq = data.chosenKeqPGK; % [] %/10
        p.PGK_Vm = 1306.45 / 60 * setup.excessPGK / data.chosenDF; % mM s^{-1} % corrected to make it appear in excess
        p.PGK_Katp = 0.3; % mM
        p.PGK_Kp3g = 0.53; % mM
        p.PGK_Kbpg = 0.003; % mM
        p.PGK_Kadp = 0.2; % mM
        % select initial points
        % reverse
        rP3G = tempSim.y(:,1);
        rATP = tempSim.y(:,2);
        rBPG = tempSim.y(:,3);
        rADP = tempSim.y(:,4);
        rNAD = tempSim.y(:,5);
        rGAP = tempSim.y(:,6);
        rPHOS = tempSim.y(:,7);
        rNADH = tempSim.y(:,8);
        % forward
        fP3G = tempSim.y(:,9);
        fATP = tempSim.y(:,10);
        fBPG = tempSim.y(:,11);
        fADP = tempSim.y(:,12);
        fNAD = tempSim.y(:,13);
        fGAP = tempSim.y(:,14);
        fPHOS = tempSim.y(:,15);
        fNADH = tempSim.y(:,16);
        % calulate v (rateEquations)
        H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
        v_GAPDHr = (-(p.TDH1_Vmr .* rBPG .* rNADH .* H ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* rNAD .* rGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + rNAD ./ p.TDH1_Knad + rNADH ./ p.TDH1_Knadh) .* (1 + rBPG ./ p.TDH1_Kbpg + rGAP ./ p.TDH1_Kgap));
        v_GAPDHf = (-(p.TDH1_Vmr .* fBPG .* fNADH .* H ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* fNAD .* fGAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + fNAD ./ p.TDH1_Knad + fNADH ./ p.TDH1_Knadh) .* (1 + fBPG ./ p.TDH1_Kbpg + fGAP ./ p.TDH1_Kgap));
        v_PGKr = p.PGK_Vm .* (rBPG .* rADP - rP3G .* rATP ./ p.PGK_Keq);
        v_PGKf = p.PGK_Vm .* (fBPG .* fADP - fP3G .* fATP ./ p.PGK_Keq);
        tempSim.v = [v_GAPDHr, v_PGKr, v_GAPDHf, v_PGKf];
%         v_GAPDH = (-(p.TDH1_Vmr .* BPG .* NADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* NAD .* GAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./ p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
%         v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq);
%         tempSim.v = [v_GAPDH, v_PGK];
        % write down in output matrix
        allSims2{i,k} = tempSim;
    end
end

% Plotting results: reverse
figure
for k = 1:np
    % plots NADH (3,6,[1:6])
    subplot(3,6,k)
    for i = 1:n
        % select specific data
        Tdata = allSims2{i,k}.t;
        Ydata = allSims2{i,k}.y;
        plot(Tdata,Ydata(:,8),'Color',[47/255 126/255 178/255])
%         ylim([0 0.15])
        if k == 1
            ylabel('NADH concentration [mM]','FontWeight','bold')
        else
            set(gca,'ytick',[])
        end
        set(gca,'xtick',[])
        xlabel(setup.params{k});
        set(gca,'xaxisLocation','top','FontWeight','bold')
        hold on
    end
    % plots v_{gapdh} (3,6,[7:12])
    subplot(3,6,k+6)
    for i = 1:n
        % select specific data
        Tdata = allSims2{i,k}.t;
        Vdata = allSims2{i,k}.v;
        plot(Tdata,Vdata(:,1),'Color',[47/255 126/255 178/255])
%         ylim([-0.08 0.02])
        if k == 1
            ylabel('v_{GAPDH} [mM s^{-1}]','FontWeight','bold')
        else
            set(gca,'ytick',[])
        end
        set(gca,'xtick',[])
        hold on
    end
    % plots v_{pgk} (3,6,[13:18])
    subplot(3,6,k+12)
    for i = 1:n
        % select specific data
        Tdata = allSims2{i,k}.t;
        Vdata = allSims2{i,k}.v;
        plot(Tdata,Vdata(:,2),'Color',[47/255 126/255 178/255])
%         ylim([-0.08 0.02])
        if k == 1
            ylabel('v_{PGK} [mM s^{-1}]','FontWeight','bold')
        else
            set(gca,'ytick',[])
            set(gca,'xtick',[])
        end
        hold on
    end
end
suptitle('PSA reverse: Each parameter is changed 3 orders of magnitude up and down');

% Plotting results: forward
figure
for k = 1:np
    % plots NADH (3,6,[1:6])
    subplot(3,6,k)
    for i = 1:n
        % select specific data
        Tdata = allSims2{i,k}.t;
        Ydata = allSims2{i,k}.y;
        plot(Tdata,Ydata(:,16),'Color',[47/255 126/255 178/255])
%         ylim([0 0.15])
        if k == 1
            ylabel('NADH concentration [mM]','FontWeight','bold')
        else
            set(gca,'ytick',[])
        end
        set(gca,'xtick',[])
        xlabel(setup.params{k});
        set(gca,'xaxisLocation','top','FontWeight','bold')
        hold on
    end
    % plots v_{gapdh} (3,6,[7:12])
    subplot(3,6,k+6)
    for i = 1:n
        % select specific data
        Tdata = allSims2{i,k}.t;
        Vdata = allSims2{i,k}.v;
        plot(Tdata,Vdata(:,3),'Color',[47/255 126/255 178/255])
%         ylim([-0.08 0.02])
        if k == 1
            ylabel('v_{GAPDH} [mM s^{-1}]','FontWeight','bold')
        else
            set(gca,'ytick',[])
        end
        set(gca,'xtick',[])
        hold on
    end
    % plots v_{pgk} (3,6,[13:18])
    subplot(3,6,k+12)
    for i = 1:n
        % select specific data
        Tdata = allSims2{i,k}.t;
        Vdata = allSims2{i,k}.v;
        plot(Tdata,Vdata(:,4),'Color',[47/255 126/255 178/255])
%         ylim([-0.08 0.02])
        if k == 1
            ylabel('v_{PGK} [mM s^{-1}]','FontWeight','bold')
        else
            set(gca,'ytick',[])
            set(gca,'xtick',[])
        end
        hold on
    end
end
suptitle('PSA forward: Each parameter is changed 3 orders of magnitude up and down');


%% (2.1a) Parameter estimation: common Vm. Data-based optimization (wD = 1, wH = 0, wR = 0)
data.gapdhF = data_gapdhF;
data.gapdhR = data_gapdhR;
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode_pH = 'on';
setup.typeVm = 'common';

setup.DFstudy = 4;
setup.costfun = 1;

setup.plotResults = 0;
setup.plotEachSim = 0;
setup.plotEachSimCF = 0;
setup.weightDataR = 1;
setup.weightDataF = 0;%0.01; %0.5;
setup.weightHaldane = 0;
setup.selectedLambda = 0;
    temp_wDR = setup.weightDataR;
    temp_wDF = setup.weightDataF;
    temp_wH = setup.weightHaldane;
    temp_wL = setup.selectedLambda;

ode_pH = setup.ode_pH;
numpH = numpHtested;
plength = length(setup.params);
pvals = zeros(numpH,plength);
pcis = zeros(numpH,plength);
x_temp = zeros(1,plength);
ub = 3*ones(1,plength); ub([1,6]) = 6;
lb = -3*ones(1,plength); lb([1,6]) = -6;
options = optimset('Display','iter');

errorData = zeros(numpH,1);
errorHaldane = zeros(numpH,1);
errorRegpars = zeros(numpH,1);
for i = 1:numpH 
% % % % for i = 1
    % inputs to be selected
    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(i);
    setup.excessPGK = 1;
    data.gapdhR.NADH = data.gapdhR.conc_mean(i,:);
    data.gapdhR.Vprofs = data.gapdhR.RRs(i,:);
    data.gapdhR.tempTime = data.gapdhR.time(i,:);
    data.gapdhF.NADH = data.gapdhF.conc_mean(i,:);
    data.gapdhF.Vprofs = data.gapdhF.RRs(i,:);
    data.gapdhF.tempTime = data.gapdhF.time(i,:);
    data.i = i;
    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_FandR,x_temp,lb,ub,options,data,setup);
    t = toc;
%     xres = pvals_VMconst(i,1:6);
%     xres = pvals_KMconst(i,1:6);
%     xres = pvals_Allconst(1,1:6);
    pvals(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    % confidence intervals estimated from covar./FIM. Only experimental
    % datapoins are considered for total N, and not regularization.
    lN = length(setup.DFstudy);
    switch lN
        case 1
            N = length(data.gapdhR.NADH{4});
        case 2
            N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3});
        case 4
            N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3}) + length(data.gapdhR.NADH{2}) + length(data.gapdhR.NADH{1});
        otherwise
            disp('No N has been selected');
    end
    Jacobian = full(Jacobian);  
    varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
    stdp = sqrt(diag(varp));
    pcis(i,:) = stdp; % confidence intervals

    % simulate and plot results
%     setup.weightData = 1;
%     setup.weightHaldane = 1; %0.1; %1;
%     setup.selectedLambda = 0;
    setup.plotEachSimCF = 1; %1;
    setup.simAllProfiles = 0; %1;
    setup.weightDataR = 1;
    setup.weightDataF = 1;
    setup.weightHaldane = 1;
    [error] = costfun_pH_FandR(xres,data,setup);
    setup.weightDataR = temp_wDR;
    setup.weightDataF = temp_wDF;
    setup.weightHaldane = temp_wH;
    setup.simAllProfiles = 0;
    setup.plotEachSimCF = 0;
    % calculating errors
    errorData(i) = sum(abs(error(1:end-7)));
    errorHaldane(i) = sum(abs(error(end-6)));
    errorRegpars(i) = sum(abs(error(end-5,end)));
end


%% (2.1b) parameter visualization + add confidence intervals
sourceVm = setup.sourceVm;
figure
for i = 1:plength
    % plot parameter values
    subplot(3,3,i)
%     plot([pHvals(1:2);pHvals(4:end)],[pvals(1:2,i);pvals(4:end,i)],'.-')
    plot(pHvals,pvals(:,i),'.-')
%     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    titleName = setup.params{i};
    title(titleName);
    % plot errors
    if i == plength
        subplot(3,3,i+1)
        plot(pHvals,errorData,'.-')
        hold on
        plot(pHvals,errorHaldane,'.-')
        hold on
        plot(pHvals,errorRegpars,'.-')
        legend('error_{Data}','error_{Haldane}','error_{Regpars}','location','southoutside','orientation','horizontal')
        ylim([0 0.1])
    end
    % plot haldaner relationship
    if i == plength
        Keq_haldane_estimated = zeros(1,numpH);
        for j = 1:numpH
            data.i = j;
            switch sourceVm
                case 'literature'
                    vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                    vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                case 'experimentalSlopes'
                    vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                    vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                case 'experimentalSlopesFixed'
                    vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                    vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                otherwise
                    disp('No source for vmax has been selected');
            end
            ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
            ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
            kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
            kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
            switch ode_pH
                case 'on'
                    H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
%                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                otherwise
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
            end
        end
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        subplot(3,3,i+3)
        semilogy(pHvals,Keq_haldane_estimated)
        hold on
        semilogy(pHvals,Keq_haldane_theory,'k+')
        legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')


%% (2.1c) Pareto front: common Vm. Data-based optimization (wD = 1, wH = 0, wR = 0)
data.gapdhF = data_gapdhF;
data.gapdhR = data_gapdhR;
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode_pH = 'on';
setup.typeVm = 'common';
typeVm = setup.typeVm;
switch typeVm
    case 'common'
    setup.params =  {'v_{maxFWD} [mM s^{-1}]';
                    'K_{gap} [mM]';
                    'K_{bpg} [mM]';
                    'K_{nad} [mM]';
                    'K_{nadh} [mM]';
                    'v_{maxREV} [mM s^{-1}]'};
    case 'specific'
    setup.params =  {'v_{maxFWD}^{REVdata} [mM s^{-1}]';
                    'K_{gap} [mM]';
                    'K_{bpg} [mM]';
                    'K_{nad} [mM]';
                    'K_{nadh} [mM]';
                    'v_{maxREV}^{REVdata} [mM s^{-1}]';
                    'v_{maxFWD}^{FWDdata} [mM s^{-1}]';
                    'v_{maxREV}^{FWDdata} [mM s^{-1}]'};
    otherwise
        disp('No parameters have been selected.');
end

setup.DFstudy = 4;
setup.costfun = 1;

pvals_cell = cell(10,1);
pcis_cell = cell(10,1);
errorData_cell = cell(10,1);
errorDataReverse_cell = cell(10,1);
errorDataForward_cell = cell(10,1);
errorHaldane_cell = cell(10,1);
errorRegpars_cell = cell(10,1);
% the weights tested now start at a low value for the experiment forward
% and then increase.
weightTest = [0 1E-4 2E-4 1E-3 2E-3 1E-2 2E-2 1E-1 2E-1 1E0];
for o = 1:10   
    setup.plotResults = 0;
    setup.plotEachSim = 0;
    setup.plotEachSimCF = 0;
    setup.weightDataR = 1;
    setup.weightDataF = weightTest(o);
    setup.weightHaldane = 0;
    setup.selectedLambda = 0;
        temp_wDR = setup.weightDataR;
        temp_wDF = setup.weightDataF;
        temp_wH = setup.weightHaldane;
        temp_wL = setup.selectedLambda;

    ode_pH = setup.ode_pH;
    numpH = numpHtested;
    plength = length(setup.params);
    pvals = zeros(numpH,plength);
    pcis = zeros(numpH,plength);
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength); ub([1,6]) = 6;
    lb = -3*ones(1,plength); lb([1,6]) = -6;
    options = optimset('Display','iter');
    
    errorDataReverse = zeros(numpH,1);
    errorDataForward = zeros(numpH,1);
    errorHaldane = zeros(numpH,1);
    errorRegpars = zeros(numpH,1);
    for i = 1:numpH 
    % for i = 1
        % inputs to be selected
        data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
        data.KeqPGK = setup.pH_Keq_pgk(i);
        setup.excessPGK = 1;
        data.gapdhR.NADH = data.gapdhR.conc_mean(i,:);
        data.gapdhR.Vprofs = data.gapdhR.RRs(i,:);
        data.gapdhR.tempTime = data.gapdhR.time(i,:);
        data.gapdhF.NADH = data.gapdhF.conc_mean(i,:);
        data.gapdhF.Vprofs = data.gapdhF.RRs(i,:);
        data.gapdhF.tempTime = data.gapdhF.time(i,:);
        data.i = i;
        % parameter estimation
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_FandR,x_temp,lb,ub,options,data,setup);
        t = toc;
        pvals(i,:) = xres;
        fprintf('Pest finished for pH #%d, time %d s\n',i,t);
        % confidence intervals estimated from covar./FIM. Only experimental
        % datapoins are considered for total N, and not regularization.
        lN = length(setup.DFstudy);
        switch lN
            case 1
                N = length(data.gapdhR.NADH{4});
            case 2
                N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3});
            case 4
                N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3}) + length(data.gapdhR.NADH{2}) + length(data.gapdhR.NADH{1});
            otherwise
                disp('No N has been selected');
        end
        Jacobian = full(Jacobian);  
        varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
        stdp = sqrt(diag(varp));
        pcis(i,:) = stdp; % confidence intervals

        % simulate and plot results
        setup.plotEachSimCF = 0; %1;
        setup.simAllProfiles = 0; %1;
        setup.weightDataR = 1;
        setup.weightDataF = 1;
        setup.weightHaldane = 1;
        [error] = costfun_pH_FandR(xres,data,setup);
        setup.weightDataR = temp_wDR;
        setup.weightDataF = temp_wDF;
        setup.weightHaldane = temp_wH;
        setup.simAllProfiles = 0;
        setup.plotEachSimCF = 0;
        % calculating errors
        errorData(i) = sum(abs(error(1:end-7)));
        errorDataReverse(i) = sum(abs(error(1:26)));
        errorDataForward(i) = sum(abs(error(27:end-7)));
        errorHaldane(i) = sum(abs(error(end-6)));
        errorRegpars(i) = sum(abs(error(end-5,end)));
    end
    pvals_cell{o} = pvals;
    pcis_cell{o} = pcis;
    errorData_cell{o} = errorData;
    errorDataReverse_cell{o} = errorDataReverse;
    errorDataForward_cell{o} = errorDataForward;
    errorHaldane_cell{o} = errorHaldane;
    errorRegpars_cell{o} = errorRegpars;

    % %% parameter visualization + add confidence intervals
    sourceVm = setup.sourceVm;
    figure
    for i = 1:plength
        % plot parameter values
        subplot(3,3,i)
        plot(pHvals,pvals(:,i),'.-')
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
        titleName = setup.params{i};
        title(titleName);
        % plot errors
        if i == plength
            subplot(3,3,i+1)
            plot(pHvals,errorData,'.-')
            hold on
            plot(pHvals,errorHaldane,'.-')
            hold on
            plot(pHvals,errorRegpars,'.-')
            legend('error_{Data}','error_{Haldane}','error_{Regpars}','location','southoutside','orientation','horizontal')
            ylim([0 0.1])
        end
        % plot haldaner relationship
        if i == plength
            Keq_haldane = zeros(1,numpH);
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
                if setup.ode_pH == 'on'
                    H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
                else
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated)
            hold on
            semilogy(pHvals,Keq_haldane_theory,'k+')
            legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
        end

    end
    suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')

end


% %% (2.1d) [Visualization all p-estimates] Pareto front
% parameter visualization + add confidence intervals
sourceVm = setup.sourceVm;
c = cool(10);
figure
for i = 1:plength
    
    % plot parameter values
    subplot(3,3,i)
    for o = 1:10
%     for o = 1
        pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
%         errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
%         hold on
    end
    titleName = setup.params{i};
    title(titleName);
    
    % plot error Data
    if i == plength
        subplot(3,3,i+1)
        for o = 1:10
            errorData = errorData_cell{o};
            plot(pHvals,errorData,'.-','color',c(o,:))
            hold on
        end
        title('error_{Data}')
    end
    
    % plot haldane relationship
    if i == plength
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        Keq_haldane_estimated = zeros(1,numpH);
        for o = 1:10
            pvals = pvals_cell{o};
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
                if setup.ode_pH == 'on'
                    H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
                else
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated,'.-','color',c(o,:))
            hold on
            if o == 10
                semilogy(pHvals,Keq_haldane_theory,'k.','MarkerSize',10)
            end
        end
%         legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('Fitting the experimental data or the haldane relationship (k_{eq})')


% %% (2.1e) [visualization] Pareto front
sumErrorDataReverse = zeros(10,1);
sumErrorDataForward = zeros(10,1);
for o = 1:10
    sumErrorDataReverse(o) = sum(abs(errorDataReverse_cell{o}));
    sumErrorDataForward(o) = sum(abs(errorDataForward_cell{o}));
end

figure
plot(sumErrorDataForward, sumErrorDataReverse,'k.-')
xlabel('Error Data Forward (NADH)')
ylabel('Error Data Reverse (NADH)')
% adding labels (locations 3 to 6)
for i = 1:10 %3:6
    hold on
    txt = erase(sprintf('w_{DF}=%d, w_{DR}=%d',weightTest(i),1),".000000");
%     if i == 3, sumErrorHaldane(i) = sumErrorHaldane(i)-0.4; end
    text(sumErrorDataForward(i),sumErrorDataReverse(i),txt,'color','blue')
end


%% (2.2a) specific Vm. Data-based optimization (wD = 1, wH = 0, wR = 0)
data.gapdhF = data_gapdhF;
data.gapdhR = data_gapdhR;
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode_pH = 'on';
setup.typeVm = 'specific';
typeVm = setup.typeVm;
switch typeVm
    case 'common'
    setup.params =  {'v_{maxFWD} [mM s^{-1}]';
                    'K_{gap} [mM]';
                    'K_{bpg} [mM]';
                    'K_{nad} [mM]';
                    'K_{nadh} [mM]';
                    'v_{maxREV} [mM s^{-1}]'};
    case 'specific'
    setup.params =  {'v_{maxFWD}^{REVdata} [mM s^{-1}]';
                    'K_{gap} [mM]';
                    'K_{bpg} [mM]';
                    'K_{nad} [mM]';
                    'K_{nadh} [mM]';
                    'v_{maxREV}^{REVdata} [mM s^{-1}]';
                    'v_{maxFWD}^{FWDdata} [mM s^{-1}]';
                    'v_{maxREV}^{FWDdata} [mM s^{-1}]'};
    otherwise
        disp('No parameters have been selected.');
end


setup.DFstudy = 4;
setup.costfun = 1;

pvals_cell = cell(10,1);
pcis_cell = cell(10,1);
errorData_cell = cell(10,1);
errorDataReverse_cell = cell(10,1);
errorDataForward_cell = cell(10,1);
errorHaldane_cell = cell(10,1);
errorRegpars_cell = cell(10,1);
% the weights tested now start at a low value for the experiment forward
% and then increase.
weightTest = [0 1E-4 2E-4 1E-3 2E-3 1E-2 2E-2 1E-1 2E-1 1E0];
% for o = 1:10 
for o = 10
    setup.plotResults = 0;
    setup.plotEachSim = 0;
    setup.plotEachSimCF = 0;
    setup.weightDataR = 1;
% % % %     setup.weightDataF = weightTest(o);
    setup.weightDataF = 0;
    setup.weightHaldane = 0;
    setup.selectedLambda = 0;
% % % %     setup.selectedLambda = 0.005; %0.01;
    
        temp_wDR = setup.weightDataR;
        temp_wDF = setup.weightDataF;
        temp_wH = setup.weightHaldane;
        temp_wL = setup.selectedLambda;

    ode_pH = setup.ode_pH;
    numpH = numpHtested;
    plength = length(setup.params);
    pvals = zeros(numpH,plength);
    pcis = zeros(numpH,plength);
    x_temp = zeros(1,plength);
% %     ub = 3*ones(1,plength); ub([1,6,7,8]) = 6;
% %     lb = -3*ones(1,plength); lb([1,6,7,8]) = -6;
    ub = 1*ones(1,plength); ub([1,6,7,8]) = 3;
    lb = -1*ones(1,plength); lb([1,6,7,8]) = -3;
    options = optimoptions('lsqnonlin','Display','iter');
%     options = optimoptions('lsqnonlin','Display','iter','FunctionTolerance',1e-10);
%     options = optimset('Display','iter');
    
    errorDataReverse = zeros(numpH,1);
    errorDataForward = zeros(numpH,1);
    errorHaldane = zeros(numpH,1);
    errorRegpars = zeros(numpH,1);
    for i = 1:numpH 
%     for i = numpH:1
    % for i = 1
% % % %     for i = numpH
        % inputs to be selected
        data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
        data.KeqPGK = setup.pH_Keq_pgk(i);
        setup.excessPGK = 1;
        data.gapdhR.NADH = data.gapdhR.conc_mean(i,:);
        data.gapdhR.Vprofs = data.gapdhR.RRs(i,:);
        data.gapdhR.tempTime = data.gapdhR.time(i,:);
        data.gapdhF.NADH = data.gapdhF.conc_mean(i,:);
        data.gapdhF.Vprofs = data.gapdhF.RRs(i,:);
        data.gapdhF.tempTime = data.gapdhF.time(i,:);
        data.i = i;
        % parameter estimation
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_FandR,x_temp,lb,ub,options,data,setup);
        t = toc;
        pvals(i,:) = xres;
        fprintf('Pest finished for pH #%d, time %d s\n',i,t);
        % confidence intervals estimated from covar./FIM. Only experimental
        % datapoins are considered for total N, and not regularization.
        lN = length(setup.DFstudy);
        switch lN
            case 1
                N = length(data.gapdhR.NADH{4});
            case 2
                N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3});
            case 4
                N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3}) + length(data.gapdhR.NADH{2}) + length(data.gapdhR.NADH{1});
            otherwise
                disp('No N has been selected');
        end
        Jacobian = full(Jacobian);  
        varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
        stdp = sqrt(diag(varp));
        pcis(i,:) = stdp; % confidence intervals

        % simulate and plot results
        setup.plotEachSimCF = 1; %1;
        setup.simAllProfiles = 0; %1;
        setup.weightDataR = 1;
        setup.weightDataF = 1;
        setup.weightHaldane = 1;
        [error] = costfun_pH_FandR(xres,data,setup);
        setup.weightDataR = temp_wDR;
        setup.weightDataF = temp_wDF;
        setup.weightHaldane = temp_wH;
        setup.simAllProfiles = 0;
        setup.plotEachSimCF = 0;
        % calculating errors
        errorData(i) = sum(abs(error(1:end-7)));
        errorDataReverse(i) = sum(abs(error(1:26)));
        errorDataForward(i) = sum(abs(error(27:end-7)));
        errorHaldane(i) = sum(abs(error(end-6)));
        errorRegpars(i) = sum(abs(error(end-5,end)));
    end
    pvals_cell{o} = pvals;
    pcis_cell{o} = pcis;
    errorData_cell{o} = errorData;
    errorDataReverse_cell{o} = errorDataReverse;
    errorDataForward_cell{o} = errorDataForward;
    errorHaldane_cell{o} = errorHaldane;
    errorRegpars_cell{o} = errorRegpars;

    % %% parameter visualization + add confidence intervals
    sourceVm = setup.sourceVm;
    figure
    for i = 1:plength
        % plot parameter values
        subplot(4,3,i)
        plot(pHvals,pvals(:,i),'.-')
%         errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
        titleName = setup.params{i};
        title(titleName);
        % plot errors
        if i == plength
            subplot(4,3,i+1)
            plot(pHvals,errorData,'.-')
            hold on
            plot(pHvals,errorHaldane,'.-')
            hold on
            plot(pHvals,errorRegpars,'.-')
            legend('error_{Data}','error_{Haldane}','error_{Regpars}','location','southoutside','orientation','horizontal')
%             ylim([0 0.1])
        end
        % plot haldaner relationship
        if i == plength
            Keq_haldane = zeros(1,numpH);
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
                if setup.ode_pH == 'on'
                    H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
                else
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
            subplot(4,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated)
            hold on
            semilogy(pHvals,Keq_haldane_theory,'k+')
            legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
        end

    end
    suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')

    figure
    for i = 1:plength
        % plot parameter values
        subplot(4,3,i)
        errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
        titleName = setup.params{i};
        title(titleName);
        ylim([-3 3])
        line([6 8],[1 1],'Color','black','LineStyle','--')
        line([6 8],[-1 -1],'Color','black','LineStyle','--')
    end
    
end


%% (2.2b) Last resource/test: MultiStart. specific Vm. Data-based optimization (wD = 1, wH = 0, wR = 
data.gapdhF = data_gapdhF;
data.gapdhR = data_gapdhR;
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode_pH = 'on';
setup.typeVm = 'specific';
typeVm = setup.typeVm;
switch typeVm
    case 'common'
    setup.params =  {'v_{maxFWD} [mM s^{-1}]';
                    'K_{gap} [mM]';
                    'K_{bpg} [mM]';
                    'K_{nad} [mM]';
                    'K_{nadh} [mM]';
                    'v_{maxREV} [mM s^{-1}]'};
    case 'specific'
    setup.params =  {'v_{maxFWD}^{REVdata} [mM s^{-1}]';
                    'K_{gap} [mM]';
                    'K_{bpg} [mM]';
                    'K_{nad} [mM]';
                    'K_{nadh} [mM]';
                    'v_{maxREV}^{REVdata} [mM s^{-1}]';
                    'v_{maxFWD}^{FWDdata} [mM s^{-1}]';
                    'v_{maxREV}^{FWDdata} [mM s^{-1}]'};
    otherwise
        disp('No parameters have been selected.');
end


setup.DFstudy = 4;
setup.costfun = 1;

pvals_cell = cell(10,1);
pcis_cell = cell(10,1);
errorData_cell = cell(10,1);
errorDataReverse_cell = cell(10,1);
errorDataForward_cell = cell(10,1);
errorHaldane_cell = cell(10,1);
errorRegpars_cell = cell(10,1);
% the weights tested now start at a low value for the experiment forward
% and then increase.
weightTest = [0 1E-4 2E-4 1E-3 2E-3 1E-2 2E-2 1E-1 2E-1 1E0];
% for o = 1:10 
for o = 10
    setup.plotResults = 0;
    setup.plotEachSim = 0;
    setup.plotEachSimCF = 0;
    setup.weightDataR = 1;
    setup.weightDataF = weightTest(o);
    setup.weightHaldane = 0.1;
    setup.selectedLambda = 0;
        temp_wDR = setup.weightDataR;
        temp_wDF = setup.weightDataF;
        temp_wH = setup.weightHaldane;
        temp_wL = setup.selectedLambda;

    ode_pH = setup.ode_pH;
    numpH = numpHtested;
    plength = length(setup.params);
    pvals = zeros(numpH,plength);
    pcis = zeros(numpH,plength);
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength); ub([1,6,7,8]) = 6;
    lb = -3*ones(1,plength); lb([1,6,7,8]) = -6;
    options = optimoptions('lsqnonlin','Display','iter');
    options = optimoptions('lsqnonlin','Display','iter','FunctionTolerance',1e-10);
%     options = optimset('Display','iter');
    
    errorDataReverse = zeros(numpH,1);
    errorDataForward = zeros(numpH,1);
    errorHaldane = zeros(numpH,1);
    errorRegpars = zeros(numpH,1);
    for i = 1:numpH 
    % for i = 1
% % % %     for i = numpH
        % inputs to be selected
        data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
        data.KeqPGK = setup.pH_Keq_pgk(i);
        setup.excessPGK = 1;
        data.gapdhR.NADH = data.gapdhR.conc_mean(i,:);
        data.gapdhR.Vprofs = data.gapdhR.RRs(i,:);
        data.gapdhR.tempTime = data.gapdhR.time(i,:);
        data.gapdhF.NADH = data.gapdhF.conc_mean(i,:);
        data.gapdhF.Vprofs = data.gapdhF.RRs(i,:);
        data.gapdhF.tempTime = data.gapdhF.time(i,:);
        data.i = i;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        nMS = 10;
        model = @(x_temp)costfun_pH_FandR(x_temp,data,setup);
        problem = createOptimProblem('lsqnonlin', 'objective', model, ...
            'xdata', data, 'ydata', data, 'x0', x_temp, 'lb', lb, ...
            'ub', ub);
        ms = MultiStart('Display','iter');
        tic
        [xres,fval,exitflag,output,solutions] = run(ms, problem, nMS);
        t = toc;
        pvals(i,:) = xres;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % %         % parameter estimation
% % % %         tic
% % % %         [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_FandR,x_temp,lb,ub,options,data,setup);
% % % %         t = toc;
% % % %         pvals(i,:) = xres;
% % % %         fprintf('Pest finished for pH #%d, time %d s\n',i,t);
% % % %         % confidence intervals estimated from covar./FIM. Only experimental
% % % %         % datapoins are considered for total N, and not regularization.
% % % %         lN = length(setup.DFstudy);
% % % %         switch lN
% % % %             case 1
% % % %                 N = length(data.gapdhR.NADH{4});
% % % %             case 2
% % % %                 N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3});
% % % %             case 4
% % % %                 N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3}) + length(data.gapdhR.NADH{2}) + length(data.gapdhR.NADH{1});
% % % %             otherwise
% % % %                 disp('No N has been selected');
% % % %         end
% % % %         Jacobian = full(Jacobian);  
% % % %         varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
% % % %         stdp = sqrt(diag(varp));
% % % %         pcis(i,:) = stdp; % confidence intervals

        % simulate and plot results
        setup.plotEachSimCF = 1; %1;
        setup.simAllProfiles = 0; %1;
        setup.weightDataR = 1;
        setup.weightDataF = 1;
        setup.weightHaldane = 1;
        [error] = costfun_pH_FandR(xres,data,setup);
        setup.weightDataR = temp_wDR;
        setup.weightDataF = temp_wDF;
        setup.weightHaldane = temp_wH;
        setup.simAllProfiles = 0;
        setup.plotEachSimCF = 0;
        % calculating errors
        errorData(i) = sum(abs(error(1:end-7)));
        errorDataReverse(i) = sum(abs(error(1:26)));
        errorDataForward(i) = sum(abs(error(27:end-7)));
        errorHaldane(i) = sum(abs(error(end-6)));
        errorRegpars(i) = sum(abs(error(end-5,end)));
    end
    pvals_cell{o} = pvals;
    pcis_cell{o} = pcis;
    errorData_cell{o} = errorData;
    errorDataReverse_cell{o} = errorDataReverse;
    errorDataForward_cell{o} = errorDataForward;
    errorHaldane_cell{o} = errorHaldane;
    errorRegpars_cell{o} = errorRegpars;

    % %% parameter visualization + add confidence intervals
    sourceVm = setup.sourceVm;
    figure
    for i = 1:plength
        % plot parameter values
        subplot(4,3,i)
        plot(pHvals,pvals(:,i),'.-')
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
        titleName = setup.params{i};
        title(titleName);
        % plot errors
        if i == plength
            subplot(4,3,i+1)
            plot(pHvals,errorData,'.-')
            hold on
            plot(pHvals,errorHaldane,'.-')
            hold on
            plot(pHvals,errorRegpars,'.-')
            legend('error_{Data}','error_{Haldane}','error_{Regpars}','location','southoutside','orientation','horizontal')
            ylim([0 0.1])
        end
        % plot haldaner relationship
        if i == plength
            Keq_haldane = zeros(1,numpH);
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
                if setup.ode_pH == 'on'
                    H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
                else
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
            subplot(4,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated)
            hold on
            semilogy(pHvals,Keq_haldane_theory,'k+')
            legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
        end

    end
    suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')

end


% %% (2.2c) [Visualization all p-estimates] Pareto front
% parameter visualization + add confidence intervals
pvals = pvals_cell{10};

sourceVm = setup.sourceVm;
c = cool(10);
figure
for i = 1:plength
    
    % plot parameter values
    subplot(4,3,i)
    for o = 1:10
%     for o = 10
        pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    titleName = setup.params{i};
    title(titleName);
    
    % plot error Data
    if i == plength
        subplot(4,3,i+1)
        for o = 1:10
%         for o = 10
            errorData = errorData_cell{o};
            plot(pHvals,errorData,'.-','color',c(o,:))
            hold on
        end
        title('error_{Data}')
    end
    
    % plot haldane relationship
    if i == plength
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        Keq_haldane_estimated = zeros(1,numpH);
        for o = 1:10
%         for o = 10
            pvals = pvals_cell{o};
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
                if setup.ode_pH == 'on'
                    H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
                else
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            subplot(4,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated,'.-','color',c(o,:))
            hold on
            if o == 10
                semilogy(pHvals,Keq_haldane_theory,'k.','MarkerSize',10)
            end
        end
%         legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('Fitting the experimental data or the haldane relationship (k_{eq})')


% %% (2.2d) [visualization] Pareto front
sumErrorDataReverse = zeros(10,1);
sumErrorDataForward = zeros(10,1);
for o = 1:10
    sumErrorDataReverse(o) = sum(abs(errorDataReverse_cell{o}));
    sumErrorDataForward(o) = sum(abs(errorDataForward_cell{o}));
end

figure
plot(sumErrorDataForward, sumErrorDataReverse,'k.-')
xlabel('Error Data Forward (NADH)')
ylabel('Error Data Reverse (NADH)')
% adding labels (locations 3 to 6)
for i = 1:10 %3:6
    hold on
    txt = erase(sprintf('w_{DF}=%d, w_{DR}=%d',weightTest(i),1),".000000");
%     if i == 3, sumErrorHaldane(i) = sumErrorHaldane(i)-0.4; end
    text(sumErrorDataForward(i),sumErrorDataReverse(i),txt,'color','blue')
end


%% (2.3a) Regularization with selected setup
data.gapdhF = data_gapdhF;
data.gapdhR = data_gapdhR;
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode_pH = 'on';
setup.typeVm = 'specific';
typeVm = setup.typeVm;
switch typeVm
    case 'common'
    setup.params =  {'v_{maxFWD} [mM s^{-1}]';
                    'K_{gap} [mM]';
                    'K_{bpg} [mM]';
                    'K_{nad} [mM]';
                    'K_{nadh} [mM]';
                    'v_{maxREV} [mM s^{-1}]'};
    case 'specific'
    setup.params =  {'v_{maxFWD}^{REVdata} [mM s^{-1}]';
                    'K_{gap} [mM]';
                    'K_{bpg} [mM]';
                    'K_{nad} [mM]';
                    'K_{nadh} [mM]';
                    'v_{maxREV}^{REVdata} [mM s^{-1}]';
                    'v_{maxFWD}^{FWDdata} [mM s^{-1}]';
                    'v_{maxREV}^{FWDdata} [mM s^{-1}]'};
    otherwise
        disp('No parameters have been selected.');
end

setup.DFstudy = 4;
setup.costfun = 1;

pvals_cell = cell(10,1);
pcis_cell = cell(10,1);
errorData_cell = cell(10,1);
errorDataReverse_cell = cell(10,1);
errorDataForward_cell = cell(10,1);
errorHaldane_cell = cell(10,1);
errorRegpars_cell = cell(10,1);

% ntest = 10;
ntest = 20;
% weightTest = [0 1E-3 2E-3 5E-3 1E-2 2E-2 5E-2 1E-1 2E-1 1E0];
weightTest = [0 1E-5 2E-5 1E-4 2E-4 1E-3 2E-3 5E-3 1E-2 2E-2 5E-2 1E-1 2E-1 1E0 2E0 1E1 2E1 1E2 2E2 1E3];

for o = 1:ntest
    setup.plotResults = 0;
    setup.plotEachSim = 0;
    setup.plotEachSimCF = 0;
    setup.weightDataR = 1;
    setup.weightDataF = 1;
% % % %     setup.weightDataF = 0;
    setup.weightHaldane = 0;
    setup.selectedLambda = weightTest(o);
        temp_wDR = setup.weightDataR;
        temp_wDF = setup.weightDataF;
        temp_wH = setup.weightHaldane;
        temp_wL = setup.selectedLambda;

    ode_pH = setup.ode_pH;
    numpH = numpHtested;
    plength = length(setup.params);
    pvals = zeros(numpH,plength);
    pcis = zeros(numpH,plength);
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength); ub([1,6,7,8]) = 6;
    lb = -3*ones(1,plength); lb([1,6,7,8]) = -6;
    options = optimoptions('lsqnonlin','Display','iter');
    % options = optimoptions('lsqnonlin','Display','iter','FunctionTolerance',1e-10);

    errorData = zeros(numpH,1);
    errorDataReverse = zeros(numpH,1);
    errorDataForward = zeros(numpH,1);
    errorHaldane = zeros(numpH,1);
    errorRegpars = zeros(numpH,1);
    for i = 1:numpH 
    % for i = 1
        % setup here as well
        setup.plotResults = 0;
        setup.plotEachSim = 0;
        setup.plotEachSimCF = 0;
        setup.weightDataR = 1;
        setup.weightDataF = 0;
        setup.weightHaldane = 0;
        setup.selectedLambda = weightTest(o);
            temp_wDR = setup.weightDataR;
            temp_wDF = setup.weightDataF;
            temp_wH = setup.weightHaldane;
            temp_wL = setup.selectedLambda;    
        % inputs to be selected
        data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
        data.KeqPGK = setup.pH_Keq_pgk(i);
        setup.excessPGK = 1;
        data.gapdhR.NADH = data.gapdhR.conc_mean(i,:);
        data.gapdhR.Vprofs = data.gapdhR.RRs(i,:);
        data.gapdhR.tempTime = data.gapdhR.time(i,:);
        data.gapdhF.NADH = data.gapdhF.conc_mean(i,:);
        data.gapdhF.Vprofs = data.gapdhF.RRs(i,:);
        data.gapdhF.tempTime = data.gapdhF.time(i,:);
        data.i = i;
        % parameter estimation
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_FandR,x_temp,lb,ub,options,data,setup);
        t = toc;
        pvals(i,:) = xres;
        fprintf('Pest finished for pH #%d, time %d s\n',i,t);
        % confidence intervals estimated from covar./FIM. Only experimental
        % datapoins are considered for total N, and not regularization.
        lN = length(setup.DFstudy);
        switch lN
            case 1
                N = length(data.gapdhR.NADH{4});
            case 2
                N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3});
            case 4
                N = length(data.gapdhR.NADH{4}) + length(data.gapdhR.NADH{3}) + length(data.gapdhR.NADH{2}) + length(data.gapdhR.NADH{1});
            otherwise
                disp('No N has been selected');
        end
        Jacobian = full(Jacobian);  
        varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
        stdp = sqrt(diag(varp));
        pcis(i,:) = stdp; % confidence intervals

        % simulate and plot results
        setup.plotEachSimCF = 0; %1;
        setup.simAllProfiles = 0; %1;   
        setup.weightDataR = 1;
        setup.weightDataF = 1;
        setup.weightHaldane = 1;
        setup.selectedLambda = 1;
        [error] = costfun_pH_FandR(xres,data,setup);
        setup.weightDataR = temp_wDR;
        setup.weightDataF = temp_wDF;
        setup.weightHaldane = temp_wH;
        setup.selectedLambda = temp_wL;
        setup.simAllProfiles = 0;
        setup.plotEachSimCF = 0;
        % calculating errors
        switch typeVm
            case 'common'
                errorData(i) = sum(abs(error(1:end-7)));
                errorDataReverse(i) = sum(abs(error(1:26)));
                errorDataForward(i) = sum(abs(error(27:end-7)));
                errorHaldane(i) = sum(abs(error(end-6)));
                errorRegpars(i) = sum(abs(error(end-5,end)));
            case 'specific'
                errorData(i) = sum(abs(error(1:end-9)));
                errorDataReverse(i) = sum(abs(error(1:26)));
                errorDataForward(i) = sum(abs(error(27:end-9)));
                errorHaldane(i) = sum(abs(error(end-8)));
                errorRegpars(i) = sum(abs(error(end-7,end)));
            otherwise
                disp('No parameters have been selected.');
        end
    end
    pvals_cell{o} = pvals;
    pcis_cell{o} = pcis;
    errorData_cell{o} = errorData;
    errorDataReverse_cell{o} = errorDataReverse;
    errorDataForward_cell{o} = errorDataForward;
    errorHaldane_cell{o} = errorHaldane;
    errorRegpars_cell{o} = errorRegpars;
% end
end

%% (2.2c) [Visualization all p-estimates] Pareto front
% parameter visualization + add confidence intervals
% pvals = pvals_cell{10};

sourceVm = setup.sourceVm;
c = cool(ntest);
figure
for i = 1:plength
    
    % plot parameter values
    subplot(4,3,i)
    for o = 1:ntest
%     for o = 10
        pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    titleName = setup.params{i};
    title(titleName);
    
    % plot error Data
    if i == plength
        subplot(4,3,i+1)
        for o = 1:ntest
%         for o = 10
            errorData = errorData_cell{o};
            plot(pHvals,errorData,'.-','color',c(o,:))
            hold on
        end
        title('error_{Data}')
    end
    
    % plot haldane relationship
    if i == plength
        Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
        Keq_haldane_estimated = zeros(1,numpH);
        for o = 1:ntest
%         for o = 10
            pvals = pvals_cell{o};
            for j = 1:numpH
                data.i = j;
                switch sourceVm
                    case 'literature'
                        vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
                        vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
                    case 'experimentalSlopes'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
                    case 'experimentalSlopesFixed'
                        vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
                        vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
                    otherwise
                        disp('No source for vmax has been selected');
                end
                ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
                ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
                kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
                kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
                if setup.ode_pH == 'on'
                    H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
                else
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
                end
            end
            subplot(4,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated,'.-','color',c(o,:))
            hold on
            if o == ntest
                semilogy(pHvals,Keq_haldane_theory,'k.','MarkerSize',10)
            end
        end
%         legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('Fitting the experimental data or the haldane relationship (k_{eq})')

% %% plot error vs lambda values
sum_eDataF = zeros(ntest,1);
sum_eDataR = zeros(ntest,1);
sum_eParams = zeros(ntest,1);
for o = 1:ntest
% % % % for o = 10:ntest
    sum_eDataR(o) = sum(abs(errorDataReverse_cell{o}));
    sum_eDataF(o) = sum(abs(errorDataForward_cell{o}));
    sum_eParams(o) = sum(abs(errorRegpars_cell{o}));
end
sum_eData = sum_eDataR + sum_eDataF;


% figure
figure
subplot(2,2,1) % sum_eData
yyaxis left
% semilogx(weightTest,sum_eData,'.-')
loglog(weightTest,sum_eData,'.-')
ylabel('error_{Data}')
hold on
yyaxis right
% semilogx(weightTest,sum_eParams,'.-')
loglog(weightTest,sum_eParams,'.-')
ylabel('error_{Parameters}')
xlabel('regularization factor lambda')
title('')

subplot(2,2,3) % sum_eDataR
yyaxis left
% semilogx(weightTest,sum_eDataR,'.-')
loglog(weightTest,sum_eDataR,'.-')
ylabel('error_{Data}^{Reverse}')
hold on
yyaxis right
% semilogx(weightTest,sum_eParams,'.-')
loglog(weightTest,sum_eParams,'.-')
ylabel('error_{Parameters}')
xlabel('regularization factor lambda')
title('')

subplot(2,2,4) % sum_eDataF
yyaxis left
% semilogx(weightTest,sum_eDataF,'.-')
loglog(weightTest,sum_eDataF,'.-')
ylabel('error_{Data}^{Forward}')
hold on
yyaxis right
% semilogx(weightTest,sum_eParams,'.-')
loglog(weightTest,sum_eParams,'.-')
ylabel('error_{Parameters}')
xlabel('regularization factor lambda')
title('')

%% memoryDump
% 
% % parameter visualization + add confidence intervals
% for o = 1:ntest
%     pvals = pvals_cell{o};
%     pcis = pcis_cell{o};
%     errorData = errorData_cell{o};
%     errorDataReverse = errorDataReverse_cell{o};
%     errorDataForward = errorDataForward_cell{o};
%     errorHaldane = errorHaldane_cell{o};
%     errorRegpars = errorRegpars_cell{o};
%     
%     sourceVm = setup.sourceVm;
%     figure
%     for i = 1:plength
%         % plot parameter values
%         subplot(4,3,i)
%         plot(pHvals,pvals(:,i),'.-')
%     %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
%         titleName = setup.params{i};
%         title(titleName);
%         % plot errors
%         if i == plength
%             subplot(4,3,i+1)
%             plot(pHvals,errorData,'.-')
%             hold on
%             plot(pHvals,errorHaldane,'.-')
%             hold on
%             plot(pHvals,errorRegpars,'.-')
%             legend('error_{Data}','error_{Haldane}','error_{Regpars}','location','southoutside','orientation','horizontal')
%             ylim([0 0.1])
%         end
%         % plot haldaner relationship
%         if i == plength
%             Keq_haldane = zeros(1,numpH);
%             for j = 1:numpH
%                 data.i = j;
%                 switch sourceVm
%                     case 'literature'
%                         vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
%                         vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation          
%                     case 'experimentalSlopes'
%                         vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(data.i); % mM s^{-1}
%                         vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(data.i); % mM s^{-1}        
%                     case 'experimentalSlopesFixed'
%                         vmf = 10 .^ pvals(j,1) .* setup.exp_vmax_gapdhf(6); % mM s^{-1}
%                         vmr = 10 .^ pvals(j,6) .* setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                     otherwise
%                         disp('No source for vmax has been selected');
%                 end
%                 ks1 = 10 .^ pvals(j,2) .* 2.48; % mM
%                 ks2 = 10 .^ pvals(j,4) .* 2.92; %mM
%                 kp1 = 10 .^ pvals(j,3) .* 1.18; % mM
%                 kp2 = 10 .^ pvals(j,5) .* 0.022; % mM
%                 if setup.ode_pH == 'on'
%                     H_effect = 10^(setup.pH_vals(j) - setup.pH_vals(6));
%                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
%                 else
%                     Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
%                 end
%             end
%             Keq_haldane_theory = setup.pH_Keq_gapdh_eQ;
%             subplot(4,3,i+3)
%             semilogy(pHvals,Keq_haldane_estimated)
%             hold on
%             semilogy(pHvals,Keq_haldane_theory,'k+')
%             legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
%         end
% 
%     end
%     suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')
% end

