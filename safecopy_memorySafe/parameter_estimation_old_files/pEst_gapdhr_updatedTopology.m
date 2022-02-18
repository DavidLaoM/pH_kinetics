% GAPDH Parameter estimation after the changes discussed with Bas an Laura:
% - Additional assays to test robutness of regularization results
% - Keq is no more a parameter to be estimated, but a pH-dependent constant
% - Haldane relationship used to calculate Vr and Vf
% - Kinetics from van Heerden 2014 are used.

% Files changed from previous sections
% 1 - 'odeGAPDHr.m' to 'odeGAPDH_vHeerden.m'
% 2 - (different ode to be called)

% Section order in this file
% (0) Setup and data load
% % % (1) System simulation and PSA


%% (0) Setup and data load
clear
set_paths_pHstudy;
dbstop if error

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
data.chosenNADini = 0.15;
temp1 = import_gapdhR.rawData.absorbance_corrected{4,4};
temp2 = import_gapdhR.rawData.absorbance_corrected{5,4};
temp3 = import_gapdhR.rawData.absorbance_corrected{6,4};
data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
data.raw.time = import_gapdhR.rawData.time{1};

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


%% (1) System simulation (@ph7) and PSA
setup.ode = 'vanHeerden2014';
% setup.ode = 'gapdhr_s_revMM';
setup.sourceVm = 'experimentalSlopes';
% setup.sourceVm = 'literature'; 
setup.nLinSpace = 21;
setup.PSArefval = 11;
setup.plotResults = 1;
setup.plotEachSim = 0;

%% % (1.1) Simple system simulation
% test run/structure with added vmf an vmr check for SS concentration
setup.ode = 'vanHeerden2014';
setup.ode_pH = 'off';

% inputs before simSys required before running
xtemp = zeros(6,1);
data.chosenDF = 1; % reference: DF = 1
data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(6); %setup.pH_Keq_gapdh(6); % reference from pH = 7.06
data.chosenKeqPGK = setup.pH_Keq_pgk(6); % reference from pH = 7.06
data.i = 6; % to select things inside the simSys file (datapsoint from the pH array).
setup.excessPGK = 1; % excess of PGK protein for faster linking reaction
setup.plotEachSim = 1;
% % vmf (x1) and vmr (x6) adjusted for equilibrium concentration:
% xtemp(1) = 1;
% xtemp(6) = -1;
% xtemp(1) = 3;
% xtemp(6) = 3;
% simulation
[testSim] = simSys(xtemp,data,setup);

%% (1.1b) Test with the effect of pH in the model
% inputs before simSys required before running
setup.ode = 'vanHeerden2014';
setup.ode_pH = 'off';
xtemp = zeros(6,1);
data.chosenDF = 1; % reference: DF = 1
data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(6); %setup.pH_Keq_gapdh(6); % reference from pH = 7.06
data.chosenKeqPGK = setup.pH_Keq_pgk(6); % reference from pH = 7.06
data.i = 6; % to select things inside the simSys file (datapsoint from the pH array).
% data.chosenpH = setup.pH_vals(6);
setup.excessPGK = 1; % excess of PGK protein for faster linking reaction
setup.plotEachSim = 1;
[testSim] = simSys(xtemp,data,setup);


% % % %%
% % % setup.plotEachSim = 0;
% % % setup.ode = 'gapdhr_s_revMM';
% % % xtemp(7) = 3;
% % % data.chosenLink = 1;
% % % % vals = [1E-4 3E-4 1E-3 3E-3 1E-2];
% % % % vals = [-1.5 -1 0 1 1.5]; %keq = 5.6E-3
% % % vals = linspace(-1.5, 1.5, 20);
% % % sims = cell(length(vals),1);
% % % % for i = 1
% % % for i = 1:length(vals)
% % % %     data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(i);
% % %     xtemp(1) = vals(i);
% % %     xtemp(6) = 2;
% % %     [testSim] = simSys(xtemp,data,setup);
% % %     sims{i} = testSim;
% % % end
% % % setup.ode = 'vanHeerden2014';
% % % setup.plotEachSim = 1;
% % % 
% % % figure
% % % for i = 1:length(vals)
% % %     plot(sims{i}.t,sims{i}.y(:,8),'color',[.5 .5 .5])
% % %     hold on
% % % end
% % % xlabel('time [s]')
% % % ylabel('NADH [mM]')
% % % %% using setup options with exact values
% % % setup.plotEachSim = 0;
% % % setup.ode = 'gapdhr_s_revMM';
% % % xtemp(7) = 3;
% % % data.chosenLink = 1;
% % % % vals = [1E-4 3E-4 1E-3 3E-3 1E-2];
% % % % vals = [-1.5 -1 0 1 1.5]; %keq = 5.6E-3
% % % % % % % vals = linspace(-1.5, 1.5, 20);
% % % sims = cell(length(setup.pH_Keq_gapdh_eQ),1);
% % % % for i = 1
% % % for i = 1:length(setup.pH_Keq_gapdh_eQ)
% % %     setup.chosenKeqgapdh = setup.pH_Keq_gapdh_eQ(i);
% % % %     data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(i);
% % %     xtemp(1) = vals(i);
% % %     xtemp(6) = 2;
% % %     [testSim] = simSys(xtemp,data,setup);
% % %     sims{i} = testSim;
% % % end
% % % setup.ode = 'vanHeerden2014';
% % % setup.plotEachSim = 1;
% % % 
% % % figure
% % % for i = 1:length(setup.pH_Keq_gapdh_eQ)
% % %     plot(sims{i}.t,sims{i}.y(:,8),'color',[.5 .5 .5])
% % %     hold on
% % % end
% % % xlabel('time [s]')
% % % ylabel('NADH [mM]')


%% % (1.2) PSA: influence of K in MA kinetics in the linking reaction
% test via changing 'setup.excessPGK'.
n = 21;
sval = logspace(-1,1,n);
c = cool(n);
xobs = 8;

allSims = cell(n,1);
figure(11)
for i = 1:n % values per parameter
    % inputs before simSys required before running
    xtemp = zeros(6,1);
    data.chosenDF = 1; % reference: DF = 1
    data.chosenKeqGAPDH = setup.pH_Keq_gapdh_eQ(6); %setup.pH_Keq_gapdh(6); % reference from pH = 7.06
    data.chosenKeqPGK = setup.pH_Keq_pgk(6); % reference from pH = 7.06
    data.i = 6;
    setup.excessPGK = sval(i); % excess of PGK protein for faster linking reaction
    setup.plotEachSim = 0;

    % simulation
    [tempSim] = simSys(xtemp,data,setup);
    allSims{i} = tempSim;
    
    % plotting results
    plot(tempSim.t,tempSim.y(:,xobs),'Color',c(i,:));
    hold on
    if i == n
        title('Change in the speed of the linking reaction')
        colormap(cool);
        ylabel('NADH [mM]')
        xlabel('time[s]')
        hold off
    end
end


%% % (1.3) PSA test: NADH, vappPGK and vappGAPDH plotted
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
        setup.excessPGK = 1; % excess of PGK protein for faster linking reaction
        setup.plotEachSim = 0;

        % simulation
        [tempSim] = simSys(xtemp,data,setup);
        
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
        P3G = tempSim.y(:,1);
        ATP = tempSim.y(:,2);
        BPG = tempSim.y(:,3);
        ADP = tempSim.y(:,4);
        NAD = tempSim.y(:,5);
        GAP = tempSim.y(:,6);
        PHOS = tempSim.y(:,7);
        NADH = tempSim.y(:,8);
        % calulate v (rateEquations)
        v_GAPDH = (-(p.TDH1_Vmr .* BPG .* NADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* NAD .* GAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./ p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
        v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq);
        tempSim.v = [v_GAPDH, v_PGK];

        % write down in output matrix
        allSims2{i,k} = tempSim;
    end
end

% Plotting results
figure(12)
for k = 1:np
    % plots NADH (3,6,[1:6])
    subplot(3,6,k)
    for i = 1:n
        % select specific data
        Tdata = allSims2{i,k}.t;
        Ydata = allSims2{i,k}.y;
        plot(Tdata,Ydata(:,8),'Color',[47/255 126/255 178/255])
        ylim([0 0.15])
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
        ylim([-0.08 0.02])
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
        ylim([-0.08 0.02])
        if k == 1
            ylabel('v_{PGK} [mM s^{-1}]','FontWeight','bold')
        else
            set(gca,'ytick',[])
            set(gca,'xtick',[])
        end
        hold on
    end
end
suptitle('PSA: Each parameter is changed 3 orders of magnitude up and down');


%% (1.4) Simulation of the effect of the dilution factor
xtemp = zeros(6,1);
data.chosenKeqGAPDH = setup.pH_Keq_gapdh(6); % reference from pH = 7.06
data.chosenKeqPGK = setup.pH_Keq_pgk(6); % reference from pH = 7.06
setup.excessPGK = 1; % excess of PGK protein for faster linking reaction
setup.plotEachSim = 0;

DF = [1 2 4 8];
DFsims = cell(length(DF),1);
figure(13)
for i = 1:4
    data.chosenDF = DF(i);
    [testSim] = simSys(xtemp,data,setup);
    DFsims = testSim;
    
    plot(testSim.t,testSim.y(:,xobs))
    hold on
end
legend('DF1','DF2','DF4','DF8')
xlabel('time [s]')
ylabel('NADH concentration [mM]')
title('Check that increasing dilution factor, reaction gets slowed down')

%% (2) Parameter estimation


%% (2.1) (pre-)Parameter estimation: simulation of all cases using the cost funtion
% In a way, this is a check that keq increases with time right now
setup.DFstudy = 4;
setup.costfun = 1;

setup.ode = 'vanHeerden2014';
setup.ode_pH = 'off';

setup.plotResults = 0;
setup.plotEachSim = 0;
setup.plotEachSimCF = 1;
setup.weightData = 0;
setup.weightHaldane = 0;
setup.selectedLambda = 0;

numpH = numpHtested;
for i = 1:numpH    
    % inputs to be selected
    x_temp = zeros(1,6);
    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(i);
    setup.excessPGK = 1;
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;
        
    [error] = costfun_pH_new(x_temp,data,setup);
end


%% (2.2) Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries
% setup.sourceVm = 'experimentalSlopes';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode = 'vanHeerden2014';
setup.ode_pH = 'on';

setup.DFstudy = 4;
setup.costfun = 1;

setup.plotResults = 0;
setup.plotEachSim = 0;
setup.plotEachSimCF = 0;
setup.weightData = 0; %0.1;
setup.weightHaldane = 1; %1;
setup.selectedLambda = 0;
    temp_wD = setup.weightData;
    temp_wH = setup.weightHaldane;

ode_pH = setup.ode_pH;
numpH = numpHtested;
plength = length(setup.params);
pvals = zeros(numpH,plength);
pcis = zeros(numpH,plength);
x_temp = zeros(1,plength);
ub = 3*ones(1,plength); ub([1,6]) = 6;
lb = -3*ones(1,plength); lb([1,6]) = -6;
% ub = 6*ones(1,plength);
% lb = -6*ones(1,plength);
options = optimset('Display','iter');

errorData = zeros(numpH,1);
errorHaldane = zeros(numpH,1);
errorRegpars = zeros(numpH,1);
for i = 1:numpH 
% for i = 1
    % inputs to be selected
    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(i);
    setup.excessPGK = 1;
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;
    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_new,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    % confidence intervals estimated from covar./FIM. Only experimental
    % datapoins are considered for total N, and not regularization.
    lN = length(setup.DFstudy);
    switch lN
        case 1
            N = length(data.NADH{4});
        case 2
            N = length(data.NADH{4}) + length(data.NADH{3});
        case 4
            N = length(data.NADH{4}) + length(data.NADH{3}) + length(data.NADH{2}) + length(data.NADH{1});
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
    setup.weightData = 1;
    setup.weightHaldane = 1;
    [error] = costfun_pH_new(xres,data,setup);
    setup.weightData = temp_wD;
    setup.weightHaldane = temp_wH;
    setup.simAllProfiles = 0;
    setup.plotEachSimCF = 0;
    % calculating errors
    errorData(i) = sum(abs(error(1:end-7)));
    errorHaldane(i) = sum(abs(error(end-6)));
    errorRegpars(i) = sum(abs(error(end-5,end)));
end

% %% parameter visualization + add confidence intervals
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


%% (2.3) Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries
% setup.sourceVm = 'experimentalSlopes';
setup.sourceVm = 'experimentalSlopesFixed';
setup.ode = 'vanHeerden2014';
setup.ode_pH = 'on';

setup.DFstudy = 4;
setup.costfun = 1;

pvals_cell = cell(10,1);
pcis_cell = cell(10,1);
errorData_cell = cell(10,1);
errorHaldane_cell = cell(10,1);
errorRegpars_cell = cell(10,1);
weightTest = [0 1E-2 2E-2 1E-1 2E-1 1E0 2E0 1E1 2E1 3E1 3E2];
for o = 1:10
    setup.plotResults = 0;
    setup.plotEachSim = 0;
    setup.plotEachSimCF = 0;
    % % % % setup.weightData = 0.1;%0; %1;
    % % % % setup.weightHaldane = 1; %0.1; %1;
    setup.weightData = 1;%0; %1;
    setup.weightHaldane = weightTest(o); %0.1; %1;
    setup.selectedLambda = 0;
        temp_wD = setup.weightData;
        temp_wH = setup.weightHaldane;

    numpH = numpHtested;
    plength = length(setup.params);
    pvals = zeros(numpH,plength);
    pcis = zeros(numpH,plength);
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength); ub([1,6]) = 6;
    lb = -3*ones(1,plength); lb([1,6]) = -6;
    % ub = 6*ones(1,plength);
    % lb = -6*ones(1,plength);
    options = optimset('Display','iter');

    errorData = zeros(numpH,1);
    errorHaldane = zeros(numpH,1);
    errorRegpars = zeros(numpH,1);
    for i = 1:numpH 
    % for i = 1
        % inputs to be selected
        data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
        data.KeqPGK = setup.pH_Keq_pgk(i);
        setup.excessPGK = 1;
        data.NADH = data.conc_mean(i,:);
        data.Vprofs = data.RRs(i,:);
        data.tempTime = data.time(i,:);
        data.i = i;
        % parameter estimation
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_new,x_temp,lb,ub,options,data,setup);
        t = toc;
        pvals(i,:) = xres;
        fprintf('Pest finished for pH #%d, time %d s\n',i,t);
        % confidence intervals estimated from covar./FIM. Only experimental
        % datapoins are considered for total N, and not regularization.
        lN = length(setup.DFstudy);
        switch lN
            case 1
                N = length(data.NADH{4});
            case 2
                N = length(data.NADH{4}) + length(data.NADH{3});
            case 4
                N = length(data.NADH{4}) + length(data.NADH{3}) + length(data.NADH{2}) + length(data.NADH{1});
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
        % %     
        setup.weightData = 1;
        setup.weightHaldane = 1;
        % % 
        [error] = costfun_pH_new(xres,data,setup);
        % %     
        setup.weightData = temp_wD;
        setup.weightHaldane = temp_wH;
        % %     
        setup.simAllProfiles = 0;
        setup.plotEachSimCF = 0;   

        % calculating errors
        errorData(i) = sum(abs(error(1:end-7)));
        errorHaldane(i) = sum(abs(error(end-6)));
        errorRegpars(i) = sum(abs(error(end-5,end)));
    end
    pvals_cell{o} = pvals;
    pcis_cell{o} = pcis;
    errorData_cell{o} = errorData;
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
                Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
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


% %% (2.3b) [Visualization all together] Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries

% parameter visualization + add confidence intervals
sourceVm = setup.sourceVm;
c = cool(10);
figure
for i = 1:plength
    
    % plot parameter values
    subplot(3,3,i)
    for o = 1:10
        pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
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
                Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
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


% %% (2.3c) [Pareto front visualization] Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries
sumErrorHaldane = zeros(10,1);
sumErrorData = zeros(10,1);
for o = 1:10
    sumErrorHaldane(o) = sum(abs(errorHaldane_cell{o}));
    sumErrorData(o) = sum(abs(errorData_cell{o}));
end

figure
plot(sumErrorHaldane,sumErrorData,'k.-')
xlabel('Error Haldane (k_{eq})')
ylabel('Error Data (NADH)')
% adding labels (locations 3 to 6)
for i = 3:6
    hold on
    txt = erase(sprintf('w_{H}=%d, w_{D}=%d',weightTest(i),1),".000000");
%     txt = sprintf('w_{H}=%d, w_{D}=%d',weightTest(i),1);
    if i == 3, sumErrorHaldane(i) = sumErrorHaldane(i)-0.4; end
    text(sumErrorHaldane(i),sumErrorData(i),txt,'color','blue')
end


%% (2.4) Understand if there is a theoretical way to fit Haldane and data at the same time
%


%% (2.5) Data regularization (lambda) without Haldane, if necessary
% setup.sourceVm = 'experimentalSlopes';
setup.sourceVm = 'experimentalSlopesFixed';

setup.DFstudy = 4;
setup.costfun = 1;

% the loop is run for the following
% ntest = 10;
ntest = 10 + 4;
% weightTest = [0                1E-2 2E-2      1E-1 2E-1 1E0 2E0 1E1 2E1 3E1];
weightTest = [0 1E-3 2E-3 5E-3 1E-2 2E-2 5E-2 1E-1 2E-1 1E0 2E0 1E1 2E1 3E1];

pvals_cell = cell(ntest,1);
pcis_cell = cell(ntest,1);
errorData_cell = cell(ntest,1);
errorHaldane_cell = cell(ntest,1);
errorRegpars_cell = cell(ntest,1);
for o = 1:ntest
    setup.plotResults = 0;
    setup.plotEachSim = 0;
    setup.plotEachSimCF = 0;
    setup.weightData = 1;%0; %1;
    setup.weightHaldane = 0; %0.1; %1;
    setup.selectedLambda = weightTest(o);
        temp_wD = setup.weightData;
        temp_wH = setup.weightHaldane;
        temp_wL = setup.selectedLambda;

    numpH = numpHtested;
    plength = length(setup.params);
    pvals = zeros(numpH,plength);
    pcis = zeros(numpH,plength);
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength); 
    ub([1,6]) = 6;
    lb = -3*ones(1,plength); 
    lb([1,6]) = -6;
    options = optimset('Display','iter');

    errorData = zeros(numpH,1);
    errorHaldane = zeros(numpH,1);
    errorRegpars = zeros(numpH,1);
    for i = 1:numpH 
    % for i = 1
        % inputs to be selected
        data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(i); %setup.pH_Keq_gapdh(i);
        data.KeqPGK = setup.pH_Keq_pgk(i);
        setup.excessPGK = 1;
        data.NADH = data.conc_mean(i,:);
        data.Vprofs = data.RRs(i,:);
        data.tempTime = data.time(i,:);
        data.i = i;
        % parameter estimation
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_new,x_temp,lb,ub,options,data,setup);
        t = toc;
        pvals(i,:) = xres;
        fprintf('Pest finished for pH #%d, time %d s\n',i,t);
        % confidence intervals estimated from covar./FIM. Only experimental
        % datapoins are considered for total N, and not regularization.
        lN = length(setup.DFstudy);
        switch lN
            case 1
                N = length(data.NADH{4});
            case 2
                N = length(data.NADH{4}) + length(data.NADH{3});
            case 4
                N = length(data.NADH{4}) + length(data.NADH{3}) + length(data.NADH{2}) + length(data.NADH{1});
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
        % %     
        setup.weightData = 1;
        setup.weightHaldane = 1;
        setup.selectedLambda = 1;
        % % 
        [error] = costfun_pH_new(xres,data,setup);
        % %     
        setup.weightData = temp_wD;
        setup.weightHaldane = temp_wH;
        setup.selectedLambda = temp_wL;
        % %     
        setup.simAllProfiles = 0;
        setup.plotEachSimCF = 0;   

        % calculating errors
        errorData(i) = sum(abs(error(1:end-7)));
        errorHaldane(i) = sum(abs(error(end-6)));
        errorRegpars(i) = sum(abs(error(end-5,end)));
    end
    pvals_cell{o} = pvals;
    pcis_cell{o} = pcis;
    errorData_cell{o} = errorData;
    errorHaldane_cell{o} = errorHaldane;
    errorRegpars_cell{o} = errorRegpars;

    % %% parameter visualization + add confidence intervals
%     for o = 1:ntest
%     pvals = pvals_cell{o};
%     pcis = pcis_cell{o};
%     for i = 1:numpH
        
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
                Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
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

%% temp plot with CIs
for o = 1:ntest
    pvals = pvals_cell{o};
    pcis = pcis_cell{o};
        
    sourceVm = setup.sourceVm;
    figure
        
        for i = 1:plength
            
            % plot parameter values
            subplot(3,3,i)
%             plot(pHvals,pvals(:,i),'.-')
            errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
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
                    Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
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


%%

% %% (2.5b) [Visualization all together] Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries

% parameter visualization + add confidence intervals
sourceVm = setup.sourceVm;
% c = cool(ntest);
c = cool(ntest-9);
c = [zeros(9,3);c];
figure
for i = 1:plength
    
    % plot parameter values
    subplot(3,3,i)
    for o = 1:ntest
% % % %     for o = 10:ntest
        pvals = pvals_cell{o};
        plot(pHvals,pvals(:,i),'.-','color',c(o,:))
        hold on
    %     errorbar(pHvals,pvals(:,i),pcis(:,i),'.-')
    end
    titleName = setup.params{i};
    title(titleName);
    
    % plot error Data
    if i == plength
        subplot(3,3,i+1)
        for o = 1:ntest
% % % %         for o = 10:ntest
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
% % % %         for o = 10:ntest
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
                Keq_haldane_estimated(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
            end
            subplot(3,3,i+3)
            semilogy(pHvals,Keq_haldane_estimated,'.-','color',c(o,:))
            hold on
            if o == ntest
                semilogy(pHvals,Keq_haldane_theory,'k.','MarkerSize',10)
            end
        end
%         legend('K_{eq,estimated}','K_{eq,haldane}','location','southoutside','orientation','horizontal')
    end
    
end
suptitle('Fitting the experimental data or regularizing parameters')

% %% plot error vs lambda values
sum_eData = zeros(ntest,1);
sum_eParams = zeros(ntest,1);
for o = 1:ntest
% % % % for o = 10:ntest
    sum_eData(o) = sum(abs(errorData_cell{o}));
    sum_eParams(o) = sum(abs(errorRegpars_cell{o}));
end

figure
yyaxis left
% plot(weightTest,sum_eData)
semilogx(weightTest,sum_eData,'.-')
ylabel('error_{Data}')
hold on
yyaxis right
% plot(weightTest,sum_eParams)
semilogx(weightTest,sum_eParams,'.-')
ylabel('error_{Parameters}')
xlabel('regularization factor lambda')



% % % % %% (2.5c) [Pareto front visualization] Test loop: Parameter estimation: unconstrained paremeter estimation with -3 +3 parameter boundaries
% % % % sumErrorHaldane = zeros(ntest,1);
% % % % sumErrorData = zeros(ntest,1);
% % % % for o = 1:ntest
% % % %     sumErrorHaldane(o) = sum(abs(errorHaldane_cell{o}));
% % % %     sumErrorData(o) = sum(abs(errorData_cell{o}));
% % % % end
% % % % 
% % % % figure
% % % % plot(sumErrorHaldane,sumErrorData,'k.-')
% % % % xlabel('Error Haldane (k_{eq})')
% % % % ylabel('Error Data (NADH)')
% % % % % adding labels (locations 3 to 6)
% % % % for i = 3:6
% % % %     hold on
% % % %     txt = erase(sprintf('w_{H}=%d, w_{D}=%d',weightTest(i),1),".000000");
% % % % %     txt = sprintf('w_{H}=%d, w_{D}=%d',weightTest(i),1);
% % % %     if i == 3, sumErrorHaldane(i) = sumErrorHaldane(i)-0.4; end
% % % %     text(sumErrorHaldane(i),sumErrorData(i),txt,'color','blue')
% % % % end






%% (2.3) Parameter Estimation: MultiStart test the first pH value
% Because for the previous case just one (boundaries 6, pH 7) hit the spot.
setup.DFstudy = 4;
setup.costfun = 1;

setup.plotResults = 0;
setup.plotEachSim = 0;
setup.plotEachSimCF = 0;
setup.weightData = 1;
setup.weightHaldane = 0;
setup.selectedLambda = 0;

numpH = numpHtested;
plength = length(setup.params);
pvals = zeros(numpH,plength);
x_temp = zeros(1,plength);
% ub = 3*ones(1,plength);
% lb = -3*ones(1,plength);
ub = 6*ones(1,plength);
lb = -6*ones(1,plength);
% options = optimset('Display','iter');
nMS = 5;


% for i = 1:numpH 
for i = 7
    % inputs to be selected
    data.KeqGAPDH = setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(i);
    setup.excessPGK = 1;
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;

    % parameter estimation (using MultiStart)
    model = @(x_temp)costfun_pH_new(x_temp,data,setup);
    problem = createOptimProblem('lsqnonlin', 'objective', model, ...
        'xdata', data, 'ydata', data, 'x0', x_temp, 'lb', lb, ...
        'ub', ub);
    ms = MultiStart('Display','iter');
    tic
    [b,fval,exitflag,output,solutions] = run(ms, problem, nMS);
    t = toc;
    fprintf('MS Pest finished for pH #%d, time %d s\n',i,t);
    pvals(i,:) = b;
%     res.b = b;
%     res.fval = fval;
%     res.exitflag = exitflag;
%     res.output = output;
%     res.solutions = solutions;
%     res.t = t;
%     MSresult{i} = res;

    % simulate and plot results
    setup.plotEachSimCF = 1;
    [error] = costfun_pH_new(b,data,setup);
    setup.plotEachSimCF = 0;
end








%     data.Vmaxs = data.Vmax(i,:);
%     data.NADH = data.conc_mean(i,:);
%     data.Vprofs = data.RRs(i,:);
%     data.tempTime = data.time(i,:);

    model = @(x_temp)costfun_pH(x_temp,data,setup);
    problem = createOptimProblem('lsqnonlin', 'objective', model, ...
        'xdata', data, 'ydata', data, 'x0', x_temp, 'lb', lb, ...
        'ub', ub);
    ms = MultiStart('Display','iter');
    tic
    [b,fval,exitflag,output,solutions] = run(ms, problem, nMS);
    t = toc;
    fprintf('MS Pest finished for pH #%d, time %d s\n',i,t);




%% parameter estimation. No extra constraint
data.pH_Keq_gapdh = setup.pH_Keq_gapdh;
data.pH_Keq_pgk = setup.pH_Keq_pgk;
data.pH_vmr_gapdh = [0.0017    0.0017    0.0019    0.0021    0.0016    0.0014    0.0011    0.0010    0.0010    0.0009];
data.pH_vmf_gapdh = [1.8780e-04 8.6996e-05 1.2221e-04 1.4983e-04 2.3061e-04 4.2531e-04 5.8826e-04 8.1749e-04 9.3141e-04 0.0011];

numpH = numpHtested;
setup.DFstudy = 4;
setup.costfun = 1;
setup.weightHaldane = 0;

setup.ode = 'vanHeerden2014';
setup.plotResults = 0;
setup.plotEachSim = 0;

plength = length(setup.params)-1;
pvals = zeros(numpH,plength);

x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
lb = -1*ones(1,plength);
setup.selectedLambda = 0;
options = optimset('Display','iter');

for i = 1:numpH
    x_temp = zeros(1,6);
    % select required data
    data.Vmr = data.pH_vmr_gapdh(i);
    data.Vmf = data.pH_vmf_gapdh(i);
    data.KeqGAPDH = data.pH_Keq_gapdh(i);
    data.KeqPGK = data.pH_Keq_pgk(i);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;
%     [error] = costfun_pH_new(x_temp,data,setup);
    
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_new,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
end

% %%
DFstudy = setup.DFstudy;
setup.plotEachSim = 0;
figure
for i = 1:numpH
    x_temp = zeros(1,6);
    % select required data
    data.Vmr = data.pH_vmr_gapdh(i);
    data.Vmf = data.pH_vmf_gapdh(i);
    data.KeqGAPDH = data.pH_Keq_gapdh(i);
    data.KeqPGK = data.pH_Keq_pgk(i);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;
    x_temp1 = pvals(i,:);
    simResultAll = cell(1:4);
% % % %     for j = DFstudy
    for j = 1:4
            % recall vmax for the specific value and simulate
            data.chosenVmf = data.Vmf/data.DF(1,j);
            data.chosenVmr = data.Vmr/data.DF(1,j);
            data.chosenKeqGAPDH = data.KeqGAPDH;
            data.chosenKeqPGK = data.KeqPGK;
            data.chosenNADini = data.NADH{j}(1);
            data.chosenDF = data.DF(1,j);
            % simulations ODEs + calculation fluxes
            [simResult] = simSys(x_temp1,data,setup);
            simResultAll{j} = simResult;
    end
    
    subplot(3,4,i)
    for j = 1:4
% % % %         simRes = simResult;
        simRes = simResultAll{j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{data.i,j},data.conc_mean{data.i,j},'k+')
        hold on
    end
    titleName = erase(sprintf('pH = %d', pHvals(i)),"0000e+00");
    title(titleName);
end
suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')

% %%
figure
for i = 1:plength
    subplot(3,3,i)
    plot(pHvals,pvals(:,i),'.-')
    titleName = setup.params{i};
    title(titleName);
    
    if i == plength
        Keq_haldane = zeros(1,numpH);
        for j = 1:numpH
            vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
            vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation
            kp1 = 10 .^ pvals(j,2) .* 2.48; % mM
            kp2 = 10 .^ pvals(j,4) .* 2.92; %mM
            ks1 = 10 .^ pvals(j,3) .* 1.18; % mM
            ks2 = 10 .^ pvals(j,5) .* 0.022; % mM
            Keq_haldane(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
        end
        Keq_theory = data.pH_Keq_gapdh;
        subplot(3,3,i+1)
        semilogy(pHvals,Keq_theory)
        hold on
        semilogy(pHvals,Keq_haldane)
        legend('K_{eq, Theory}','K_{eq, Haldane}')
    end
    
end
suptitle('vanHeerden 2014 kinetics. NADH fit. No Haldane Constraint')


%% parameter estimation. Added haldane extra constraint
data.pH_Keq_gapdh = setup.pH_Keq_gapdh;
data.pH_Keq_pgk = setup.pH_Keq_pgk;
data.pH_vmr_gapdh = [0.0017    0.0017    0.0019    0.0021    0.0016    0.0014    0.0011    0.0010    0.0010    0.0009];
data.pH_vmf_gapdh = [1.8780e-04 8.6996e-05 1.2221e-04 1.4983e-04 2.3061e-04 4.2531e-04 5.8826e-04 8.1749e-04 9.3141e-04 0.0011];

numpH = numpHtested;
setup.DFstudy = 4;
setup.costfun = 1;
setup.weightHaldane = 0.1;

setup.ode = 'vanHeerden2014';
setup.plotResults = 0;
setup.plotEachSim = 0;

plength = length(setup.params)-1;
pvals = zeros(numpH,plength);

x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
lb = -1*ones(1,plength);
setup.selectedLambda = 0;
options = optimset('Display','iter');

for i = 1:numpH
    x_temp = zeros(1,6);
    % select required data
    data.Vmr = data.pH_vmr_gapdh(i);
    data.Vmf = data.pH_vmf_gapdh(i);
    data.KeqGAPDH = data.pH_Keq_gapdh(i);
    data.KeqPGK = data.pH_Keq_pgk(i);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;
%     [error] = costfun_pH_new(x_temp,data,setup);
    
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_new,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
end

% %%
DFstudy = setup.DFstudy;
setup.plotEachSim = 0;
figure
for i = 1:numpH
    x_temp = zeros(1,6);
    % select required data
    data.Vmr = data.pH_vmr_gapdh(i);
    data.Vmf = data.pH_vmf_gapdh(i);
    data.KeqGAPDH = data.pH_Keq_gapdh(i);
    data.KeqPGK = data.pH_Keq_pgk(i);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    data.i = i;
    x_temp1 = pvals(i,:);
    simResultAll = cell(1:4);
% % % %     for j = DFstudy
    for j = 1:4
            % recall vmax for the specific value and simulate
            data.chosenVmf = data.Vmf/data.DF(1,j);
            data.chosenVmr = data.Vmr/data.DF(1,j);
            data.chosenKeqGAPDH = data.KeqGAPDH;
            data.chosenKeqPGK = data.KeqPGK;
            data.chosenNADini = data.NADH{j}(1);
            data.chosenDF = data.DF(1,j);
            % simulations ODEs + calculation fluxes
            [simResult] = simSys(x_temp1,data,setup);
            simResultAll{j} = simResult;
    end
    
    subplot(3,4,i)
    for j = 1:4
% % % %         simRes = simResult;
        simRes = simResultAll{j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{data.i,j},data.conc_mean{data.i,j},'k+')
        hold on
    end
    titleName = erase(sprintf('pH = %d', pHvals(i)),"0000e+00");
    title(titleName);
end
suptitle('vanHeerden 2014 kinetics. NADH fit. Added Haldane Constraint')

% %%
figure
for i = 1:plength
    subplot(3,3,i)
    plot(pHvals,pvals(:,i),'.-')
    titleName = setup.params{i};
    title(titleName);
    
    if i == plength
        Keq_haldane = zeros(1,numpH);
        for j = 1:numpH
            vmf = 10 .^ pvals(j,1) .* 1184.52/60; % mM s^{-1} % old implementation
            vmr = 10 .^ pvals(j,6) .* 6549.8/60; % mM s^{-1} % old implementation
            kp1 = 10 .^ pvals(j,2) .* 2.48; % mM
            kp2 = 10 .^ pvals(j,4) .* 2.92; %mM
            ks1 = 10 .^ pvals(j,3) .* 1.18; % mM
            ks2 = 10 .^ pvals(j,5) .* 0.022; % mM
            Keq_haldane(j) =  (vmf * kp1 * kp2) / (vmr * ks1 * ks2);
        end
        Keq_theory = data.pH_Keq_gapdh;
        subplot(3,3,i+1)
        semilogy(pHvals,Keq_theory)
        hold on
        semilogy(pHvals,Keq_haldane)
        legend('K_{eq, Theory}','K_{eq, Haldane}')
    end
end
suptitle('vanHeerden 2014 kinetics. NADH fit. Added Haldane Constraint')


%% memoryDump

% % One case works if
% 
% ub = 6*ones(1,plength);
% lb = -6*ones(1,plength);
% 
% pHvals(7)
% ans =
%     7.2900
% pvals(7,:)
%   Columns 1 through 5
%     0.2542    0.6243    1.8228    1.3603    0.4749
%   Column 6
%    -1.5722
