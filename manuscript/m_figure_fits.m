% % M_FIGURE_FITS.M
% tic
plotAllDFs = 1; % plots all dilution factors
% plotAllDFs = 0; % plots only DF1
% 
c_royalBlue = [65	105	225]/255; % royalblue
c_midnightblue = [25	25	112]/255; % midnightblue
c_CCCCCC = [204	204	204]/255; % #CCCCCC
c_E5E5E5 = [229 229 229]/255; % #E5E5E5
c_0f1076 = [15	16	118]/255; % #0f1076
c_chocolate = [210	105	30]/255; % (#e59400 temp orange)
% 
setup.c_royalBlue = c_royalBlue; %[65	105	225]/255; % royalblue
setup.c_midnightblue = c_midnightblue; %[25	25	112]/255; % midnightblue
setup.c_CCCCCC = c_CCCCCC; %[204	204	204]/255; % #CCCCCC
setup.c_E5E5E5 = c_E5E5E5; %[229 229 229]/255; % #E5E5E5
setup.c_0f1076 = c_0f1076; %[15	16	118]/255; % #0f1076
setup.c_chocolate = c_chocolate; %[210	105	30]/255; % (#e59400 temp orange)
% 
setup.only2DFs = 1;


%% PART 1. Vmax vs pH
enzymeList = {'hxk_parEst';... %1 hxk
    'pgi_parEst';... %2 pgi
    'pfk_parEst';... %3 pfk
    'ald_parEst';... %4 ald
    'tpi_parEst';... %5 tpi
    'gapdh_parEst';... %6 gapdh_fwd
    'gapdhr_parEst';... %7 gapdh_rev
    'pgm_parEst';... %8 pgm
    'eno_parEst';... %9 eno_kmfixed
    'pyk_parEst';... %10 pyk
    'pdc_parEst'}; %11 pdc
numEnz = length(enzymeList);
enzymeName = {'hxk';... %1 hxk
    'pgi';... %2 pgi
    'pfk';... %3 pfk
    'ald';... %4 ald
    'tpi';... %5 tpi
    'gapdh_{fwd}';... %6 gapdh_fwd
    'gapdh_{rev}';... %7 gapdh_rev
    'pgm';... %8 pgm
    'eno';... %9 eno_kmfixed
    'pyk';... %10 pyk
    'pdc'}; %11 pdc

colSteelBlue = [70/255 130/255 180/255]; % pH independet
colLightBlue = [173/255 216/255 203/255]; 
colLightSkyBlue = [135/255	206/255	250/255]; % pH dependent
col4 = [8,81,156]/255;
%%
if exist('mainFig','var')
    clf(1000)
end
%%
mainFig = figure(1000);
% for i =1:numEnz
for i = [7 8 9]
    % load and select the data
    load(enzymeList{i});
    if      i == 7, i2 = 4;
    elseif  i == 8, i2 = 7;
    elseif  i == 9, i2 = 1;
    end
    sp1 = subplot(3,3,i2); % select subplot
    tempName = ['output_',extractBefore(enzymeList{i},"_parEst")];
    eval(['dataset = ',tempName,';']); % select the data
    % 
    eb2 = plot(dataset.pHarray, dataset.vm_uChange,'o');%,...
    eb2.MarkerSize = 4.5;
    eb2.Color = [0 0 0];
    eb2.MarkerFaceColor = c_midnightblue;
    % 
    hold on
    plot(dataset.pHarray, dataset.vm_uChange, '-',...
        'LineWidth', 2, 'color', c_midnightblue);
    % 
    eb1 = plot(dataset.pHarray, dataset.Vmax_experimental,'o');%,...
    eb1.MarkerSize = 4.5;
    eb1.Color = [0 0 0];
    eb1.MarkerFaceColor = c_royalBlue;
    xticks([6 6.5 7 7.5 8])
    % 
    plot(dataset.pHarray, dataset.Vmax_experimental, '-',...
        'LineWidth', 2, 'color', c_royalBlue);
    hold on
    % readjust Y-axis
    ax = gca;
    ax.YLim(1) = 0;
    ax.YLim(2) = max([dataset.Vmax_experimental;dataset.vm_uChange])*1.1;
    % labels
    if i == 1
        ylabel_h = ylabel('Enzyme capacity (\mumol mg_{P}^{-1} min^{-1})','FontSize',14);
        ylabel_h.Position(2) = -1.175;
    end
    % xlabel
    if i == 11
        xlabel_h = xlabel('pH','FontSize',14);
        xlabel_h.Position(2) = -1;
    end
    xloc = (ax.XLim(2) - ax.XLim(1))*0.035 + ax.XLim(1);
    if((i == 6)||(i == 7))
        yloc = (ax.YLim(2) - ax.YLim(1))*0.85 + ax.YLim(1);
    else
        yloc = (ax.YLim(2) - ax.YLim(1))*0.87 + ax.YLim(1);
    end    
    [~,idx] = max(dataset.Vmax_experimental);
    pH_max = dataset.pHarray(idx);
    yloc2 = (ax.YLim(2) - ax.YLim(1))*0.05 + ax.YLim(1);
    % vertical line pointing at ideal pH
    line([pH_max pH_max], ax.YLim, 'Color', 'black', ...
        'LineStyle', '--', 'Linewidth', 2)
    plot(pH_max, dataset.Vmax_experimental(idx), ...
        'ko', 'MarkerSize', 12)
    plot(pH_max, dataset.vm_uChange(idx), ...
        'ko', 'MarkerSize', 12)

    % shading
    [l,p] = boundedline(dataset.pHarray, dataset.Vmax_experimental, dataset.stDev_experimental,...
        '-b*');
    outlinebounds(l,p);
    % 
    set(gca,'children',flipud(get(gca,'children')))
    
    % added labels in the end
    if i2 == 1
        hyl = ylabel('ENO','FontSize',12);
        text(7, 2 + (2 - 0) * 0.5, {'Enzyme capacity vs pH'; ...
            '(\mumol.min^{-1}.mg protein^{-1})'}, ...
            'FontSize', 12, ...
            'HorizontalAlignment', 'center')
    elseif i2 == 4
        hyl = ylabel('GAPDHR','FontSize',12);
    elseif i2 == 7
        hyl = ylabel('PGM','FontSize',12);
    end
    
    hyl.Position = hyl.Position + [-0.3 0 0];
    hyl.FontSize = 11;
    
    hold off
end
set(gcf,'color','w');


%% PART 2. data ENO
for recallENOstart = 1
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
            Vmax(i) = (conc_mean{i}(end) - conc_mean{i}(1)) ./ (time{i}(end) - time{i}(1)); 
        end
        % 
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
        % 
        pHvals = unique(import_eno.treatedData.pH_corrected);
    end
    % %% (0.1) Calculation of rates: moving window
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
    setup.plotOutput = 0;
    expRatesDetermination;

    % %% (1.1) Simple parameter fit. Parameter estimation
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
    optfun = @costfun_Kmfixed;
    plength = 14; % Kms (2) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
    options = optimset('Display','iter');

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
end
loadName = [setup.enzymeName, '_regularizationResults.mat'];
load(loadName);
selLambdaPos = 1; %12;

%% Simulation estimated parameters
setup.literatureValues = 0;
xres_selected = array_xres{selLambdaPos};
setup.colLightSkyBlue = colLightSkyBlue;
% 
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
if plotAllDFs == 1
    [error] = direct_vs_curvefitting_costfun_Kmfixed(xres_selected,data,setup);
elseif plotAllDFs == 0
    [error] = direct_vs_curvefitting_costfun_Kmfixed_oneDF(xres_selected,data,setup);
end
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
%% Simulation direct method vmax and literature km parameters
setup.literatureValues = 1;
setup.col4 = col4;
for i = 1:length(array_xres)
    array_xres{i} = zeros(size(array_xres{i}));
end
xres_selected = array_xres{selLambdaPos};

setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
if plotAllDFs == 1
%     [error] = direct_vs_curvefitting_costfun_Kmfixed(xres_selected,data,setup);
    [error] = direct_vs_curvefitting_costfun_Kmfixed_nolabel(xres_selected,data,setup);
elseif plotAllDFs == 0
    [error] = direct_vs_curvefitting_costfun_Kmfixed_oneDF(xres_selected,data,setup);
end
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;


%% PART 3. data GAPDHrev
for recallGAPDHrevstart = 1
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
%         % added
%         setup.saveOutput = 0;

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
            % Option 1. Vmax from the values obtained
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

    end
    % (0.1) Calculation of rates: moving window
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

    % % (1.1) Simple parameter fit. Parameter estimation
    % % % % setup.ode = 'vanHeerden2014';
    setup.ode = 'gapdh_rev_simplified';
    setup.sourceVm = 'experimentalSlopes';
    setup.ode_pH = 'on';
    setup.caseKm = 'pH_independent';
    % setup.caseKm = 'pH_dependent';

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
%     caseKm = setup.caseKm;
%     switch caseKm
%         case 'pH_independent'
%             % plength = 10; % Kms (4*0) + Vms (1) * numpH (10) % overly simplified 2020-09-30
            plength = 14; % Kms (4*1) + Vms (1) * numpH (10) % still keeping kms
%         case 'pH_dependent'
%             plength = 50; % ( Kms (4) + Vms (1) ) * numpH (10) % still keeping kms
%         otherwise
%             disp('Warning: no specification has been made on Km being pH dependent or independent');
%     end
    optfun = @costfun_Kmfixed;
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
    options = optimset('Display','iter');
    
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
end
loadName = [setup.enzymeName, '_regularizationResults.mat'];
load(loadName);
selLambdaPos = 1; %12;
%% Simulation estimated parameters
setup.literatureValues = 0;
xres_selected = array_xres{selLambdaPos};
setup.colLightSkyBlue = colLightSkyBlue;
% option that can be added for pfk
% xres_selected = array_xres{8}; %lambdalist based in 'ones', lam=0.1, loc=5.
%     xres_selected(10:12) = xres_selected(9);
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
if plotAllDFs == 1
    [error] = direct_vs_curvefitting_costfun_Kmfixed(xres_selected,data,setup);
elseif plotAllDFs == 0
    [error] = direct_vs_curvefitting_costfun_Kmfixed_oneDF(xres_selected,data,setup);
end
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
%% Simulation direct method vmax and literature km parameters
setup.literatureValues = 1;
setup.col4 = col4;
for i = 1:length(array_xres)
    array_xres{i} = zeros(size(array_xres{i}));
%     array_xres{i}(5:end) = zeros(size(array_xres{i}(5:end)));
end
xres_selected = array_xres{selLambdaPos};
%
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
if plotAllDFs == 1
%     [error] = direct_vs_curvefitting_costfun_Kmfixed(xres_selected,data,setup);
    [error] = direct_vs_curvefitting_costfun_Kmfixed_nolabel(xres_selected,data,setup);
elseif plotAllDFs == 0
    [error] = direct_vs_curvefitting_costfun_Kmfixed_oneDF(xres_selected,data,setup);
end
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;


%% PART 4. data PGM
safecopy_setup = setup;
clear setup
for recallPGMstart = 1
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
        setup.caseStudyPGM = 1;
        setup.caseStudyPYK = 0;
        setup.caseStudyTPI = 0;
        setup.caseStudyENO = 0;
        selectSetup_pH;

        load('expData.mat','expData');
        import_pgm = expData.pgm;

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

        pHarray = unique(import_pgm.treatedData.pH_corrected);
        for i = 1:numpHtested
            pHval = pHarray(i);
            tempID = find(import_pgm.treatedData.pH_corrected==pHval);
            pHTemp(:,i) = import_pgm.treatedData.pH_corrected(tempID);
            DFTemp(:,i) = import_pgm.treatedData.dilution_corrected(tempID);
            for j = 1:4
                abs_meanTemp{j,i} = import_pgm.treatedData.absorbance_mean{tempID(j)};
                abs_stdTemp{j,i} = import_pgm.treatedData.absorbance_std{tempID(j)};
                conc_meanTemp{j,i} = import_pgm.treatedData.concentration_mean{tempID(j)};
                conc_stdTemp{j,i} = import_pgm.treatedData.concentration_std{tempID(j)};
                timeTemp{j,i} = import_pgm.treatedData.time{tempID(j)};
                RRsTemp{j,i} = import_pgm.treatedData.reaction_rate{tempID(j)};
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
        temp1 = import_pgm.rawData.absorbance_corrected{4,4};
        temp2 = import_pgm.rawData.absorbance_corrected{5,4};
        temp3 = import_pgm.rawData.absorbance_corrected{6,4};
        data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
        data.raw.time = import_pgm.rawData.time{1};

        % 
        dps = length(NADH{1,1});
        endPoint = zeros(numpHtested,DFs);
        % locate the minimum
        for i = 1:DFs
            for j = 1:numpHtested
                endPoint(j,i) = min(NADH{j,i});
            end
        end
        % bringing the minimum to zero
        for i = 1:DFs
            for j = 1:numpHtested
                if(((i == 1)&&(j>=7))||((i == 2)&&(j>=8)&&(j<=9)))
                    for k = 1:dps
                        NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,4);
                    end
                else
                    for k = 1:dps
                        NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,i);
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
        % correcting the shifted dilution factors
        data.DF(1:2,1:4) = [16 8 4 2; 16 8 4 2];
        data.conc_mean = NADH;
        pHvals = unique(import_pgm.treatedData.pH_corrected);
    end

    % %% (0.1) Calculation of rates: moving window
    DF(1:2,1:4) = [16 8 4 2; 16 8 4 2];
    % intial things that could be in the setup
    minwindow = 7; % minimum size of the window
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

    % Experimental rates determination and plotting
    setup.plotOutput = 0;
    expRatesDetermination;

    % %% (1.1) Simple parameter fit. Parameter estimation
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
    optfun = @costfun_Kmfixed;
    plength = 14; % Kms (2) + Vms (1) * numpH (12)
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
    options = optimset('Display','iter');
    %
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
end
loadName = [setup.enzymeName, '_regularizationResults.mat'];
load(loadName);
selLambdaPos = 16; %1; %12;

% %% Simulation estimated parameters
% 
setup.c_royalBlue = c_royalBlue; %[65	105	225]/255; % royalblue
setup.c_midnightblue = c_midnightblue; %[25	25	112]/255; % midnightblue
setup.c_CCCCCC = c_CCCCCC; %[204	204	204]/255; % #CCCCCC
setup.c_E5E5E5 = c_E5E5E5; %[229 229 229]/255; % #E5E5E5
setup.c_0f1076 = c_0f1076; %[15	16	118]/255; % #0f1076
setup.c_chocolate = c_chocolate; %[210	105	30]/255; % (#e59400 temp orange)
% 
setup.only2DFs = 1;
%     
setup.literatureValues = 0;
xres_selected = array_xres{selLambdaPos};
setup.colLightSkyBlue = colLightSkyBlue;
% option that can be added for pfk
% xres_selected = array_xres{8}; %lambdalist based in 'ones', lam=0.1, loc=5.
%     xres_selected(10:12) = xres_selected(9);
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
if plotAllDFs == 1
    [error] = direct_vs_curvefitting_costfun_Kmfixed(xres_selected,data,setup);
elseif plotAllDFs == 0
    [error] = direct_vs_curvefitting_costfun_Kmfixed_oneDF(xres_selected,data,setup);
end
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
% %% Simulation direct method vmax and literature km parameters
setup.literatureValues = 1;
setup.col4 = col4;
for i = 1:length(array_xres)
    array_xres{i} = zeros(size(array_xres{i}));
end
xres_selected = array_xres{selLambdaPos};
%
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
if plotAllDFs == 1
%     [error] = direct_vs_curvefitting_costfun_Kmfixed(xres_selected,data,setup);
    [error] = direct_vs_curvefitting_costfun_Kmfixed_nolabel(xres_selected,data,setup);
elseif plotAllDFs == 0
    [error] = direct_vs_curvefitting_costfun_Kmfixed_oneDF(xres_selected,data,setup);
end
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;

% last edits
% set size
% set(1000, 'Position', [100 100 750 750])
set(1000, 'Position', [-1500 -150 750 750])

% size plots re organization
hfig = get(1000,'Children');
% sp4 -> hfig(9)
% sp5 -> hfig(4) 
% sp6 -> hfig(3) 
hfig(9).Position(4) = hfig(9).Position(3);
hfig(9).Position = hfig(9).Position + [0 -0.025 0 0];
hfig(4).Position(4) = hfig(4).Position(3);
hfig(4).Position = hfig(4).Position + [0 -0.025 0 0];
hfig(3).Position(4) = hfig(3).Position(3);
hfig(3).Position = hfig(3).Position + [0 -0.025 0 0];
% sp1 -> hfig(7) 
% sp2 -> hfig(6) 
% sp3 -> hfig(5) 
hfig(7).Position(4) = hfig(7).Position(3);
hfig(7).Position = hfig(7).Position + [0 -0.05 0 0];
hfig(6).Position(4) = hfig(6).Position(3);
hfig(6).Position = hfig(6).Position + [0 -0.05 0 0];
hfig(5).Position(4) = hfig(5).Position(3);
hfig(5).Position = hfig(5).Position + [0 -0.05 0 0];
% sp7 -> hfig(8) 
% sp8 -> hfig(2) 
% sp9 -> hfig(1) 
hfig(8).Position(4) = hfig(8).Position(3);
hfig(2).Position(4) = hfig(2).Position(3);
hfig(1).Position(4) = hfig(1).Position(3);

% adding legends
hL1 = legend(sp1.Children([4 9]), 'Direct method', 'Standard error');
hL1.Orientation = 'horizontal';
hL1.Box = 'off';
hL1.FontSize = 11;
hL1.Position = [0.03    0.04    0.3960    0.0340];
% 
hold on
% 
ah = gca;
ah2 = copyobj( ah, gcf);
hL2 = legend(ah2, sp1.Children([2 6]), 'Curve fitting', 'Maximum Vmax');
hL2.Orientation = 'horizontal';
hL2.Box = 'off';
hL2.FontSize = 11;
hL2.Position = [0.03    0.01    0.3960    0.0340];
% 
sp_PEP = mainFig.Children(3);
hL3 = legend(sp_PEP.Children([2 7 1]), 'Simulation direct method',...
                                        'Simulation curve fitting',...
                                        'Experimental data');
% 
hL3.Orientation = 'vertical';
hL3.Box = 'off';
hL3.FontSize = 11;
hL3.Position = [0.5    0.02    0.3960    0.0340];
% 

%% some points from the latest round of comments
% Changing a bit string names
% 
mainFig.Children(2).Children(5).String = '[NADH] . s^{-1}, pH 7.06';
mainFig.Children(4).Children(5).String = '[NADH] . s^{-1}, pH 7.06';
mainFig.Children(5).Children(5).String = '[NADH], pH 7.06';
% 
mainFig.Children(6).Children(9).String = '[NADH] . s^{-1}, pH 6.60';
mainFig.Children(7).Children(9).String = '[NADH], pH 6.60';
% 
mainFig.Children(8).Children(10).String = '[PEP] . s^{-1}, pH 7.90';
mainFig.Children(9).Children(10).String = '[PEP], pH 7.90';


%% temp annotation
dim = [.55 .085 .3 .0];
str = '...';
temp = annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'EdgeColor','none','FontSize',40);


%% needed stop not to get truncated... #matlabUselessSecrets
% save
setup = safecopy_setup ;
if setup.saveOutput == 1
    savefig(1000,'pHmanus_fits')
    % specs printing (method 3)
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters fits
end
% toc





