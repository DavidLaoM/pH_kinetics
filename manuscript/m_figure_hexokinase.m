% % M_FIGURE_HEXOKINASE.M
% 
% tic
if setup.fast_option == 0 % directly recall the workspace from here
    temp_setup = setup;
    pEst_hxk_keq_dependency
    setup = temp_setup;
elseif setup.fast_option == 1 % directly recall the workspace from here
    load('workspace_hexokinase.mat');
end


%%
blueTriplet = [0, 0.4470, 0.7410];
orangeTriplet = [0.8500, 0.3250, 0.0980];
greyTriplet = [0.3 0.3 0.3];
rotationValue = 40 * ones(1,18);
    rotationValue(1) = 15;
    rotationValue(2) = 20;
    rotationValue(3) = 25;
    
% % % % clf(1001), 
figure(1001)


%% SUBPLOT(1,3,1). PARAMETER ESTIMATES
    
% added for colors
c_royalBlue = [65	105	225]/255; % royalblue
c_midnightblue = [25	25	112]/255; % midnightblue
c_CCCCCC = [204	204	204]/255; % #CCCCCC
c_E5E5E5 = [229 229 229]/255; % #E5E5E5
c_0f1076 = [15	16	118]/255; % #0f1076
c_chocolate = [210	105	30]/255; % (#e59400 temp orange)

% 
spA = subplot(2,3,[1 2 4 5]);
plot(output_hxk_changing_pH.pHarray, output_hxk_changing_pH.vm_uChange, ...
    'o-', 'color', c_midnightblue, 'Linewidth', 1.2, 'MarkerSize', 5,...
    'MarkerFaceColor', c_midnightblue, 'LineWidth', 2)
hold on
plot(output_hxk_changing_pH.pHarray, output_hxk_constant_pH.vm_uChange, ...
    'o-', 'color', c_chocolate, 'Linewidth', 1.2, 'MarkerSize', 5,...
    'MarkerFaceColor', c_chocolate, 'LineWidth', 2)
% 
ax = gca;
ax.FontSize = 12; 
ylabel('Enzyme capacity (\mumol.min^{-1}.mg protein^{-1})','FontSize',12); 
xlabel('pH','FontSize',12);
% % % % title('Enzyme capacity vs pH range')


%% SUBPLOT(1,3,3). PROGRESSION CURVES
% clear
% added for colors
c_royalBlue = [65	105	225]/255; % royalblue
c_midnightblue = [25	25	112]/255; % midnightblue
c_CCCCCC = [204	204	204]/255; % #CCCCCC
c_E5E5E5 = [229 229 229]/255; % #E5E5E5
c_0f1076 = [15	16	118]/255; % #0f1076
c_chocolate = [210	105	30]/255; % (#e59400 temp orange)
%     
for recallHXKestimates = 1
    for step0 = 1
        % select specific case and recall data
        setup.caseStudyALD = 0;
        setup.caseStudyENO = 0;
        setup.caseStudyGAPDH = 0;
        setup.caseStudyGAPDHr = 0;
        setup.caseStudyHXK = 1;
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
        import_eno = expData.hxk;

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

        NADPH = blankCell;
        Vmax = blank;
        for i = 1:(DFs*numpHtested)
            tempDiff = conc_mean{i} - conc_mean{i}(1); % all stoichiometries are 1-to-1.
            NADPH{i} = conc_mean{i};
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
        data.chosenNADPHini = 0.05;
        temp1 = import_eno.rawData.absorbance_corrected{4,4};
        temp2 = import_eno.rawData.absorbance_corrected{5,4};
        temp3 = import_eno.rawData.absorbance_corrected{6,4};
        data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
        data.raw.time = import_eno.rawData.time{1};

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % Directly changing the concentration here, sicne the extinction
        % coefficient did not change.
        dps = length(NADPH{1,1});
        iniPointDF4 = zeros(numpHtested,DFs);
        for i = 1:DFs
            for j = 1:numpHtested
                iniPointDF4(j,i) = NADPH{j,1}(1);
            end
        end
        for i = 1:DFs
            for j = 1:numpHtested
                for k = 1:dps
                    NADPH{j,i}(k) = NADPH{j,i}(k) - iniPointDF4(j,1);
                end
            end
        end
        data.conc_mean = NADPH;
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pHvals = unique(import_eno.treatedData.pH_corrected);
    end
%     % %% (0.1) Experimental Vmax determination
%     tempSave = setup.plotOutput;
    setup.plotOutput = 0;
    % %% (0.1) Calculation of rates: moving window
    minwindow = 6; % minimum size of the window
    limRates = [0 5E-4]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.1];  %Ylims plot conc
    idxs_one_initial = ones(size(data.conc_mean));
    idxs_one_initial(1:2,:) = (20 - minwindow) * idxs_one_initial(1:2,:);
    idxs_one_initial(3:end,:) = (61 - minwindow) * idxs_one_initial(3:end,:);
    dp_start = ones(size(data.conc_mean));
    dp_start(:,1) = 6 * ones(size(dp_start(:,1)));
    % blank cell total length
    total_len = zeros(size(dp_start));
    % DFs considered
    DFarray = [1/8 1/4 1/2 1/1];
    % idxs2consider
    idxs2consider = ones(size(DF));
    % Experimental rates determination and plotting
    expRatesDetermination;
%     setup.plotOutput = tempSave;
    
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
    plength = 16; % Kms (4) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
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
saveName = [setup.enzymeName, '_regularizationResults.mat'];
load(saveName)
selLambdaPos = 1;%11;%1;%14;

% %% simulation changing keq with pH
xres_selected = array_xres{selLambdaPos};
setup.plotEachSimCF = 0; % 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 0; % 1;
% [error] = optfun(xres_selected,data,setup);
[~,simResult1,vObs1] = costfun_Kmfixed_simOutputNew(xres_selected,data,setup);
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;


%% simulation constant keq with pH
% 
nTests = length(setup.Keq_HXK);
lastSim_full = cell(1,nTests);
for o = 1:nTests
    % 
    temp_Keq_HXK = setup.Keq_HXK;
    temp_Keq_G6PDH = setup.Keq_G6PDH;
    for i = 1:12
        setup.Keq_HXK(i) = setup.Keq_HXK(o);
        setup.Keq_G6PDH(i) = setup.Keq_G6PDH(o);
    end
    % 
    xres_selected = array_xres{selLambdaPos};
    setup.plotEachSimCF = 0; % 1;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0; % 1;
    % [error] = optfun(xres_selected,data,setup);
    [~,simResult2,vObs2] = costfun_Kmfixed_simOutputNew(xres_selected,data,setup);
    setup.plotEachSimCF = 0;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0;
    % 
    setup.Keq_HXK = temp_Keq_HXK;
    setup.Keq_G6PDH = temp_Keq_G6PDH;
    % 
    lastSim_full{o}.simRes = simResult2;
    lastSim_full{o}.v = vObs2;
end


%%
figure(1001)
% subplot(1,3,3); % vmHXK
spB = subplot(2,3,3); % vmHXK
% 
x2 = [lastSim_full{o}.simRes.t(1:end)', fliplr(lastSim_full{o}.simRes.t(1:end)')];
inBetween = [lastSim_full{12}.simRes.y(1:end,7)', fliplr(lastSim_full{1}.simRes.y(1:end,7)')];
f1 = fill(x2, inBetween, 'g');
f1.FaceColor = [0.90 0.90 0.90];
f1.EdgeColor = 'none';
hold on
% 
plot(simResult1.t,simResult1.y(:,7),'-','LineWidth',2,'color', c_midnightblue)
plot(data.time{12,4}(linspace(1,61,16)), data.conc_mean{12,4}(linspace(1,61,16)), 'k.', 'MarkerSize', 10)
% % % % title('Progression curve (pH7.90) with different pH regimes')
ax = gca;
ax.FontSize = 12; 
ylabel('NAPDH (mM)','FontSize',12)
xlabel('Assay time (s)','FontSize',12)
% legend('Keq@7.90','kEQ@6.81','location','southeast')
% box on

%% SUBPLOT(1,3,2). VMAX VS KEQ
Keq_array_hxk = setup.Keq_HXK;
Keq_array_g6pdh = setup.Keq_G6PDH;
simResTotal = cell(1,length(Keq_array_hxk));
vObsTotal = cell(1,length(Keq_array_hxk));

for j = 1:length(Keq_array_hxk)
    %
    temp_Keq_HXK = setup.Keq_HXK;
    temp_Keq_G6PDH = setup.Keq_G6PDH;
    for i = 1:12
        setup.Keq_HXK(i) = setup.Keq_HXK(j);
        setup.Keq_G6PDH(i) = setup.Keq_G6PDH(j);
    end
    %
    xres_selected = array_xres{selLambdaPos};
    setup.plotEachSimCF = 0;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0;
    % [error] = optfun(xres_selected,data,setup);
    [~,simResultN,vObsN] = costfun_Kmfixed_simOutputNew(xres_selected,data,setup);
    setup.plotEachSimCF = 0;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0;
    %
    setup.Keq_HXK = temp_Keq_HXK;
    setup.Keq_G6PDH = temp_Keq_G6PDH;
    %
    simResTotal{j} = simResultN;
    vObsTotal{j} = vObsN;
end

%%
Vmax_array = zeros(size(Keq_array_hxk));
for j = 1:length(Keq_array_hxk)
    Vmax_array(j) = max(vObsTotal{j}) .* 60 .* 60 ./ setup.concProtein;
end

figure(1001), 
% subplot(1,3,2); % vmHXK
spC = subplot(2,3,6); % vmHXK
% yyaxis right
plot(Keq_array_hxk, Vmax_array, 'k-', 'LineWidth', 2)
ax = gca;
ax.FontSize = 12; 
xlabel('pH-independent K_{eq}','FontSize',12)
% ylabel('Enzyme capacity','FontSize',12)
ylabel({'Enzyme capacity';'(\mumol.min^{-1}.mg protein^{-1})'},'FontSize',12)
% 
set(1001,'color','w')


%%
spA = subplot(2,3,[1 2 4 5]);
    % axes and tick labels
    xlim([6 8])
    ylim([0 0.85])
    xticks([6 6.5 7 7.5 8])
    yticks([0 0.2 0.4 0.6 0.8])
spB = subplot(2,3,3);
    % axes and tick labels
%     set(gca, 'YAxisLocation', 'right')
    xlim([0 300])
    ylim([0.05 0.15])
    xticks([0 100 200 300])
    yticks([0.05 0.10 0.15])
spC = subplot(2,3,6);
    hold on
    % axes and tick labels
%     set(gca, 'YAxisLocation', 'right')
    ylim([0.42 0.52])
    xlim([0 8000])
    yticks([0.44 0.46 0.48 0.50])
    xticks([0 2000 4000 6000 8000])
    xticklabels({'0','2E3','4E3','6E3','8E3'})


%%
set(1001, 'Position', [100 100 800 500])
%%
spA = subplot(2,3,[1 2 4 5]);
    % text
    text(0.05 * [spA.XLim(2) - spA.XLim(1)] + spA.XLim(1), 0.9 * [spA.YLim(2) - spA.YLim(1)] + spA.YLim(1), ...
        'A', 'FontSize', 20)
spB = subplot(2,3,3);
    % text
    text(0.05 * [spB.XLim(2) - spB.XLim(1)] + spB.XLim(1), 0.9 * [spB.YLim(2) - spB.YLim(1)] + spB.YLim(1), ...
        'B', 'FontSize', 20)
spC = subplot(2,3,6);
    % text
    text(0.05 * [spC.XLim(2) - spC.XLim(1)] + spC.XLim(1), 0.9 * [spC.YLim(2) - spC.YLim(1)] + spC.YLim(1), ...
        'C', 'FontSize', 20)
    
    
%%
% 
% delete(hL1)
% 
hL1 = legend(spA.Children([2 3]), 'pH-dependent K_{eq}', 'pH-independent K_{eq}');
hL1.Orientation = 'vertical';
hL1.Box = 'off';
hL1.FontSize = 11;
hL1.Position = [0.25    0.20    0.3960    0.0340];


%% Latest edits (mail 2021 11 15)
spB.Children(2).MarkerSize = 15; % dots
c_classicBlue = [0 0.4470 0.7410];
spB.Children(3).Color = c_classicBlue; %line, previous [0.0980 0.0980 0.4392]
%
% delete(hL2)
% 
hL2 = legend(spB.Children([2 3 4]), 'Exp.data','pH-dep. K_{eq}', 'pH-indep. K_{eq}');
hL2.Orientation = 'vertical';
hL2.Box = 'off';
hL2.FontSize = 9;
% hL2.Position = [0.63    0.645    0.3960    0.0340];
hL2.Position = [0.685    0.645    0.3960    0.0340];
% 
spA.Position = [0.1300 0.1100 0.4942 0.8150];
spB.Position = [0.7500 0.5961 0.2100 0.3289];
spC.Position = [0.7500 0.1223 0.2100 0.3289];

%
%%
% save
if setup.saveOutput == 1
    savefig(1001,'pHmanus_hexokinase')
    % 
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters hexokinase
end
% toc


