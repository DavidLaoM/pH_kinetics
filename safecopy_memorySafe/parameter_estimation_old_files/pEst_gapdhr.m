% pEst_gapdhR.m
% Estimate parameters for the gapdh_reverse reaction.

% 0 - Initial setup and data import DONE
% 1 - PSA DONE
% 2 - Parameter estimation: non-regularized -
% 2a - Simple case. Non regularized. DF1. CI and covariance
% 2b - Simple case. Non regularized. DF12. CI and covariance
% 2c - Simple case. Non regularized. DF1234. CI and covariance
% 2d - Study the cost function. Regularization. DF1234

% 2c - Global optimization. MultiStart
% 2e - Selected case. Uncertainty analysis. PLA
% 2f - Selected case. Uncertainty analysis. Anything else?

% 3 - Parameter estimation: regularized
% 3b - Uncertainty analysis: regularized
% 4 - Final visualization


%% (0) Setup
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
end
suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
suptitle(suptitleName);


%% (1) PSA
% Variables that have little influence in the output can already be seen.
setup.nLinSpace = 21;
setup.PSArefval = 11;
setup.plotResults = 1;
[PSAresult] = PSA(setup,data);

if setup.saveOutput == 1
    saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_PSAsims.mat'];
    saveFigName1 = ['results\', setup.enzymeName, '\', setup.enzymeName, '_PSA.fig'];
    %saving
    save(saveName,'PSAresult');
    savefig(1,saveFigName1);
end


%% (2) Parameter estimation. Non regularized
setup.ode = 'gapdhr_s_revMM';
setup.plotEachSim = 0;

% %% (2a) Simple case. Non regularized. DF1.
% parameters estimated changed from the power idea to ease the calculation
% of the confidence intervals. Seemed it was harder to understand.
setup.DFstudy = 4;
setup.costfun = 1;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1 = zeros(numpH,plength);
pcis_DF1 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
% % % %     [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pHcheck,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF1(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF1(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF1(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end

% % % % % % % % 
% % % % figure
% % % % for i = 1:plength
% % % %     subplot(3,3,i)
% % % %     plot(pHvals,pvals_DF1(:,i))
% % % %     hold on
% % % % %     plot(pHvals,pvals_DF12(:,i))
% % % % %     hold on
% % % % %     plot(pHvals,pvals_DF1234(:,i))
% % % % %     hold on
% % % %     ylim([-1 1])
% % % %     title(setup.params{i})
% % % %     if i == plength
% % % % %         hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
% % % %         hL = legend('pvalsDF1');
% % % %         newPosition = [0.7 0.1 0.2 0.2];
% % % %         newUnits = 'normalized';
% % % %         set(hL,'Position', newPosition,'Units', newUnits);
% % % %         hold off
% % % %     end
% % % % end
% % % % suptitle('Plot parameter estimates vs pH. Non-regularized')
% % % % 
% % % % figure
% % % % for i = 1:plength
% % % %     subplot(3,3,i)
% % % %     errorbar(pHvals,pvals_DF1(:,i),pcis_DF1(:,i))
% % % %     hold on
% % % % %     errorbar(pHvals,pvals_DF12(:,i),pcis_DF12(:,i))
% % % % %     hold on
% % % % %     errorbar(pHvals,pvals_DF1234(:,i),pcis_DF1234(:,i))
% % % % %     hold on
% % % %     x = [6 8];
% % % %     y = [0 0];
% % % %     line(x,y,'Color','black','LineStyle','--')
% % % %     hold on
% % % % %     ylim([-1 1])
% % % %     title(setup.params{i})
% % % %     if i == plength
% % % % %         hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
% % % %         hL = legend('pvalsDF1');
% % % %         newPosition = [0.7 0.1 0.2 0.2];
% % % %         newUnits = 'normalized';
% % % %         set(hL,'Position', newPosition,'Units', newUnits);
% % % %         hold off
% % % %     end
% % % % end
% % % % suptitle('Errorbar parameter estimates vs pH. Non-regularized')
% % % % % % % % 


% % % % 
% %% (2b) Study the cost function. Test use of DF1 and DF2
setup.DFstudy = [3 4];
setup.costfun = 2;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF12 = zeros(numpH,plength);
pcis_DF12 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF12(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF12(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF12(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
    end
    
%     figure
%     subplot(1,2,1)
%     for j = 1:DFs
%         simRes = simulations{i,j};
%         plot(simRes.t,simRes.y(:,8),'-')
%         hold on
%         plot(data.time{i,j},data.conc_mean{i,j},'k+')
%     end
%     title('NADH')
%     subplot(1,2,2)
%     for j = 1:DFs
%         simRes = simulations{i,j};
%         simRRs = simRates{i,j};
%         plot(simRes.t,simRRs,'-')
%         hold on
%         plot(data.time{i,j},-data.RRs{i,j},'k+')
%     end
%     title('v_{GAPDHr}')   
%     suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (2c) Study the cost function. Test use of DF1, DF2, DF3 and DF4
setup.DFstudy = [1 2 3 4];
setup.costfun = 3;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1234 = zeros(numpH,plength);
pcis_DF1234 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF1234(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF1234(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF12(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
    end
    
%     figure
%     subplot(1,2,1)
%     for j = 1:DFs
%         simRes = simulations{i,j};
%         plot(simRes.t,simRes.y(:,8),'-')
%         hold on
%         plot(data.time{i,j},data.conc_mean{i,j},'k+')
%     end
%     title('NADH')
%     subplot(1,2,2)
%     for j = 1:DFs
%         simRes = simulations{i,j};
%         simRRs = simRates{i,j};
%         plot(simRes.t,simRRs,'-')
%         hold on
%         plot(data.time{i,j},-data.RRs{i,j},'k+')
%     end
%     title('v_{GAPDHr}')   
%     suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (2a,2b,2c) Visualization
figure
for i = 1:plength
    subplot(3,3,i)
    plot(pHvals,pvals_DF1(:,i))
    hold on
    plot(pHvals,pvals_DF12(:,i))
    hold on
    plot(pHvals,pvals_DF1234(:,i))
    hold on
    ylim([-1 1])
    title(setup.params{i})
    if i == plength
        hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Plot parameter estimates vs pH. Non-regularized')

figure
for i = 1:plength
    subplot(3,3,i)
    errorbar(pHvals,pvals_DF1(:,i),pcis_DF1(:,i))
    hold on
    errorbar(pHvals,pvals_DF12(:,i),pcis_DF12(:,i))
    hold on
    errorbar(pHvals,pvals_DF1234(:,i),pcis_DF1234(:,i))
    hold on
    x = [6 8];
    y = [0 0];
    line(x,y,'Color','black','LineStyle','--')
    hold on
%     ylim([-1 1])
    title(setup.params{i})
    if i == plength
        hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Errorbar parameter estimates vs pH. Non-regularized')


% %%
% % % enzName = setup.enzymeName;
% % % saveFigName31b = ['results\', enzName,'\', enzName, '_pEst_lsq_6p_plot.fig'];
% % % saveFigName32b = ['results\', enzName,'\', enzName, '_pEst_lsq_6p_errorbar.fig'];
% % % savefig(31,saveFigName31b)
% % % savefig(32,saveFigName32b)
% % saveFigName31a = ['results\', enzName,'\', enzName, '_pEst_lsq_7p_plot.fig'];
% % saveFigName32a = ['results\', enzName,'\', enzName, '_pEst_lsq_7p_errorbar.fig'];
% % savefig(31,saveFigName31a)
% % savefig(32,saveFigName32a)
% saveFigName32c = ['results\', enzName,'\', enzName, '_pEst_lsq_7p_non_regularized_plot.fig'];
% saveFigName33c = ['results\', enzName,'\', enzName, '_pEst_lsq_7p_non_regularized_errorbar.fig'];
% savefig(32,saveFigName32c)
% savefig(33,saveFigName33c)


%% (2d) Study the cost function. Test use of DF1 and DF2. Added regularization
setup.DFstudy = [1 2 3 4];
setup.costfun = 3;
setup.lambdaArray = [1E-3 1E-2 1E-1 1E0 1E1 1E2 1E3];

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1234 = zeros(numpH,plength);
pcis_DF1234 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.lambdaArray = [1E-4 2E-4 1E-3 2E-3 1E-2 2E-2 1E-1 2E-1 1E0 2E0 1E1 2E1 1E2 2E2];

llam = length(setup.lambdaArray);
errorData = zeros(size(setup.lambdaArray));
errorParameters = zeros(size(setup.lambdaArray));
errorData_cell = cell(1,numpH);
errorParameters_cell = cell(1,numpH);
% parameter estiamtion for each pH value
% for i = 1
for i = 1:numpH
    errorData_cell{i} = zeros(size(errorData));
    errorParameters_cell{i} = zeros(size(errorParameters));
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
    
    % parameter estimation with selected lambda value
%     tic
    for j = 1:llam
        setup.selectedLambda = setup.lambdaArray(j);
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
        pvals_DF1234(i,:) = xres;
        errorData(j) = sum(residual(1:end-7));
        errorParameters(j) = abs(sum(residual(end-6:end))/setup.selectedLambda);
        fprintf('Iteration lambda#%d done\n',j);
    end
    errorData_cell{i} = errorData;
    errorParameters_cell{i} = errorParameters;
    fprintf('Iteration pH#%d done\n',i);
%     t = toc;
%     fprintf('Pest done for pH #%d all lambda values, time %d s\n',i,t);
end

% add line @0.1
figure
for i = 1:numpH
    subplot(3,4,i)
    temp_errorData = errorData_cell{i};
    temp_errorParameters = errorParameters_cell{i};
    yyaxis left
    % plot(setup.lambdaArray,errorData)
    semilogx(setup.lambdaArray,temp_errorData)
    ax1 = gca;
    % semilogx(setup.lambdaArray,errorData,'-o')
    hold on
    yyaxis right
    % plot(setup.lambdaArray,errorParameters)
    semilogx(setup.lambdaArray,temp_errorParameters)
    ax2 = gca;
    % semilogx(setup.lambdaArray,errorParameters,'-o')
    hold on
    mins = [ax1.YLim(1), ax2.YLim(1)];
    maxs = [ax1.YLim(2), ax2.YLim(2)];
    line([0.01 0.01],[min(mins) max(maxs)],'Color','black','LineStyle','--')
% % % %     title(sprintf('regulatization for pH#%d',setup.params{i})); % pH vals
    if i == plength
        hL = legend('errorData','errorParameters','lambda = 0.01');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end

% %% saving
% % enzName = setup.enzymeName;
% % saveResname = ['results\', enzName,'\', enzName, '_pEst_regularization_pHall.mat'];
% % save(saveResname,'errorParameters_cell','errorData_cell');
% % saveFigName = ['results\', enzName,'\', enzName, '_pEst_regularization_pHall.fig'];
% % savefig(1,saveFigName)


%% (3) Parameter estimation. Regularized

% %% (3a) Simple case. Regularized. DF1.
% parameters estimated changed from the power idea to ease the calculation
% of the confidence intervals. Seemed it was harder to understand.
setup.DFstudy = 4;
setup.costfun = 1;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1 = zeros(numpH,plength);
pcis_DF1 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0.01;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF1(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF1(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF1(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (3b) Study the cost function. Test use of DF1 and DF2. Regularized
setup.DFstudy = [3 4];
setup.costfun = 2;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF12 = zeros(numpH,plength);
pcis_DF12 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0.01;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF12(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF12(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF12(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (3c) Study the cost function. Test use of DF1, DF2, DF3 and DF4. Regularized
setup.DFstudy = [1 2 3 4];
setup.costfun = 3;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1234 = zeros(numpH,plength);
pcis_DF1234 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0.01;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF1234(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF1234(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF12(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (3a,3b,3c) Visualization
figure
for i = 1:plength
    subplot(3,3,i)
    plot(pHvals,pvals_DF1(:,i))
    hold on
    plot(pHvals,pvals_DF12(:,i))
    hold on
    plot(pHvals,pvals_DF1234(:,i))
    hold on
    ylim([-1 1])
    title(setup.params{i})
    if i == plength
        hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
% % % %         hL = legend('pvalsDF1');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Plot parameter estimates vs pH. Regularized')

figure
for i = 1:plength
    subplot(3,3,i)
    errorbar(pHvals,pvals_DF1(:,i),pcis_DF1(:,i))
    hold on
    errorbar(pHvals,pvals_DF12(:,i),pcis_DF12(:,i))
    hold on
    errorbar(pHvals,pvals_DF1234(:,i),pcis_DF1234(:,i))
    hold on
    x = [6 8];
    y = [0 0];
    line(x,y,'Color','black','LineStyle','--')
    hold on
    ylim([-1 1])
%     ylim([-1 1])
    title(setup.params{i})
    if i == plength
        hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
% % % %         hL = legend('pvalsDF1');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Errorbar parameter estimates vs pH. Regularized')

% % saving
% saveFigName31c = ['results\', enzName,'\', enzName, '_pEst_lsq_7p_regularized_plot.fig'];
% saveFigName32c = ['results\', enzName,'\', enzName, '_pEst_lsq_7p_regularized_errorbar.fig'];
% savefig(31,saveFigName31c)
% savefig(32,saveFigName32c)


%% (4) Global optimization. MultiStart
% setup.DFstudy = [3 4];
% setup.costfun = 2;
setup.DFstudy = 4;
setup.costfun = 1;

setup.MSruns = 500;
setup.figWindow = 20;
setup.plotResults = 1;
setup.limsMS = 1;
setup.selectedLambda = 0.01;

[MSresult_DF1_500] = MS_pH(setup,data);
% %%
% save('MSresult_DF1_500','MSresult_DF1_500')
%% firs plots. deppendencies

MStempRes = MSresult_DF1_500{1}.solutions;
lM = length(MStempRes);
lP = length(setup.paramRefVals);
msParray = zeros(lP,lM);
for i = 1:lP
    for j = 1:lM
        msParray(i,j) = MStempRes(j).X(i);
    end
end
% %%
figure
for i = 1:(lP*lP)
    i1 = mod(i,lP);
    i2 = fix(i/lP)+1;
    if i1 == 0
        i1 = lP;
        i2 = fix(i/lP);
    end

    disp(i1); disp(i2); disp(i1+(i2-1)*lP); 
    subplot(lP,lP,i1+(i2-1)*lP)
    
    for j = 1:lM
        plot(msParray(i1,j),msParray(i2,j),'.','MarkerSize',3,'color','blue')
        hold on
    end
    
end

figure('name','Correlation estimated parameters','units','normalized','outerposition',[0 0 1 1])
cell77={{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']};{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']};{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']};{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']};{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']};{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']};{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']},{['-g']}};
C = {cell77};
[h,labelfontsize] = subplotplus(C);

% % %% to cell subidx numbers and location
% for i = 1:length(h)
%     subidx = i;
%     set(gcf,'CurrentAxes',h(subidx));
%     text(0.5,0.5,num2str(i))
% end

% %% content. Parameter values
for i = 1:(lP*lP)

    i1 = mod(i,lP);
    i2 = fix(i/lP)+1;
    if i1 == 0
        i1 = lP;
        i2 = fix(i/lP);
    end
    
%     subplot(lP,lP,i1+(i2-1)*lP)
    subidx = i1+(i2-1)*lP;
    set(gcf,'CurrentAxes',h(subidx));
    for j = 1:lM
        plot(msParray(i1,j),msParray(i2,j),'.','MarkerSize',3,'color','blue')
        hold on
    end
    % 
    set(gca,'xaxisLocation','top')
% %     set(gca,'xtick',[])
% %     set(gca,'ytick',[])
    
    if i == 1
        ylabel(setup.params{i}(1:7))
        xlabel(setup.params{i}(1:7))
    elseif i == 2
        xlabel(setup.params{i}(1:7))
    elseif i == 3
        xlabel(setup.params{i}(1:7))
    elseif i == 4
        xlabel(setup.params{i}(1:7))
    elseif i == 5
        xlabel(setup.params{i}(1:8))
    elseif i == 6
        xlabel(setup.params{i}(1:7))
    elseif i == 7
        xlabel('p_{Link}')
    elseif i == 8
        ylabel(setup.params{2}(1:7))
    elseif i == 15
        ylabel(setup.params{3}(1:7))
    elseif i == 22
        ylabel(setup.params{4}(1:7))
    elseif i == 29
        ylabel(setup.params{5}(1:8))
    elseif i == 36
        ylabel(setup.params{6}(1:7))
    elseif i == 43
        ylabel('p_{Link}')
    end
%     if i == (26 || 27 || 32 || 34 || 39 || 40)
    if ((i == 26) || (i == 27) || (i == 32) || (i == 34) || (i == 39) || (i == 40))
        set(gca,'Color',[0.9 0.9 0.9])
    end
end

% %% saving
enzName = setup.enzymeName;
saveResname = ['results\', enzName,'\', enzName, '_pEst_MSresult_DF1_500_regularized.mat'];
save(saveResname,'MSresult_DF1_500');
% enzName = setup.enzymeName;
% % saveResname = ['results\', enzName,'\', enzName, '_pEst_MSresult_DF12_regularized.mat'];
% % save(saveResname,'MSresult_DF12');
% saveFigName = ['results\', enzName,'\', enzName, '_pEst_lsq_7p_MS_pDeppendencies.fig'];
% savefig(1,saveFigName)


%% coefficient intervals vs 
% 
[l1,l2] = size(MSresult_DF1_500);
msParray = cell(1,l1);
for k = 1:l1
    MStempRes = MSresult_DF1_500{k}.solutions;
    lM = length(MStempRes);
    lP = length(setup.paramRefVals);
    msParray{k} = zeros(lP,lM);
    for i = 1:lP
        for j = 1:lM
            msParray{k}(i,j) = MStempRes(j).X(i);
        end
    end
end
% %%
% save('gapdhr_pEst_MSresult_DF1_500_regularized_reorg.mat','msParray');

%% example analysed for pH#1 (pH=6.19)
% ntest = 10;
% limval = [2 5 10 20 40 80 160 320 440 500];
ntest = 500;
limval = 1:ntest;
mstemp = msParray{1};
pmean = zeros(1,ntest);
pstd = zeros(1,ntest);
for j = 1:7
    for i = 1:ntest
        pmean(i) = mean(mstemp(j,1:limval(i)));
        pstd(i) = std(mstemp(j,1:limval(i)));
    end
    figure

    subplot(121)
    plot(limval,pmean)
    title({'Evolution of mean for parameter 1 vs';
    'number of samples taken into account';
    'Estimation for pH = 6.19'});

    subplot(122)
    plot(limval,pstd)
    title({'Evolution of std for parameter 1 vs';
    'number of samples taken into account';
    'Estimation for pH = 6.19'});

end

%% (5) Selected case. Uncertainty analysis. PLA
dbstop if error
setup.limsPLA = 3;
setup.step = 0.05; %0.025; % 0.01;
setup.caseStudy.parameters = 1:7;
setup.selectedLambda = 0.01;
setup.plotResults = 0;

setup.DFstudy = [3 4];
setup.costfun = 2;

runPLA_pH;

% % saving results
% saveResname = ['results\', enzName,'\', enzName, '_pEst_PLA_DF12.mat'];
% save(saveResname,'plPar_cell','plRes_cell');

%%

% figure
% for j = 1:length(setup.params)
%     subplot(3,3,j)
%     for i = 1:numpHtested
%         plPar = plPar_cell{j};
%         plRes = plRes_cell{j};
%         plot(plPar, plRes)
% %         plot(pProfile.pPar.names{j},pProfile.pRes.names{j})
%         hold on
%         if i == numpHtested
%             title(setup.params{j})
%         end
%     end
%     xlim([-3 3])
%     if j == length(setup.params)
%         hL = legend('pH1','pH2','pH3','pH4','pH5','pH6','pH7','pH8','pH9','pH10');
%         newPosition = [0.7 0.1 0.2 0.2];
%         newUnits = 'normalized';
%         set(hL,'Position', newPosition,'Units', newUnits);
%     end
%     
% end
% suptitle('PLA GAPDH_{reverse}')
% % savefig


figure
for j = 1:length(setup.params)
    subplot(3,3,j)
    for i = 1:numpHtested
        importname = sprintf('GAPDH_{reverse}_PLA_pH%d.mat',i);
        load(importname)
        plPar = pProfile.pPar.names{j};
        plRes = pProfile.pRes.names{j};
        
        [val, idx] = max(plRes(1:(end-20)));
%         disp(val); disp(idx);
        if (idx ~= 1)&&((j == 2)||(j == 3)||(j == 4)||(j == 5)||(j == 7))
%             if (idx ~= 124)
                plRes(idx-2) = plRes(idx-1);
                plRes(idx-1) = plRes(idx-1);
                plRes(idx) = plRes(idx-1);
                plRes(idx+1) = plRes(idx-1);
                plRes(idx+2) = plRes(idx-1);
%             end
        end
%         [val, idx] = max(plRes(40:75));  
%         if((j == 2)||(j == 3)||(j == 4)||(j == 5)||(j == 7))
%             plRes(idx-2) = plRes(idx-1);
%             plRes(idx-1) = plRes(idx-1);
%             plRes(idx) = plRes(idx-1);
%             plRes(idx+1) = plRes(idx-1);
%             plRes(idx+2) = plRes(idx-1);   
%         end
% %         [val, idx] = min(plRes);
% %         plRes(idx-1) = plRes(idx);
% %         plRes(idx) = plRes(idx);
% %         plRes(idx+1) = plRes(idx);

        if j == 2
            [val, idx] = max(plRes(50:70));
%             plRes(50+idx-2) = plRes(50+idx-1);
%             plRes(50+idx-1) = plRes(50+idx-1);
%             plRes(50+idx) = plRes(50+idx-1);
%             plRes(50+idx+1) = plRes(50+idx-1);
%             plRes(50+idx+2) = plRes(50+idx-1); 
            plRes(50+idx-2) = plRes(50+idx);
            plRes(50+idx-1) = plRes(50+idx);
            plRes(50+idx) = plRes(50+idx);
            plRes(50+idx+1) = plRes(50+idx);
            plRes(50+idx+2) = plRes(50+idx);             
        elseif j == 3
            [val, idx] = max(plRes(50:70));
            plRes(50+idx-2) = plRes(50+idx);
            plRes(50+idx-1) = plRes(50+idx);
            plRes(50+idx) = plRes(50+idx);
            plRes(50+idx+1) = plRes(50+idx);
            plRes(50+idx+2) = plRes(50+idx);
        elseif j == 4    
            [val, idx] = max(plRes(50:70)); 
            plRes(50+idx-2) = plRes(50+idx);
            plRes(50+idx-1) = plRes(50+idx);
            plRes(50+idx) = plRes(50+idx);
            plRes(50+idx+1) = plRes(50+idx);
            plRes(50+idx+2) = plRes(50+idx);       
        elseif j == 5    
            [val, idx] = max(plRes(50:90)); 
            plRes(50+idx-2) = plRes(50+idx);
            plRes(50+idx-1) = plRes(50+idx);
            plRes(50+idx) = plRes(50+idx);
            plRes(50+idx+1) = plRes(50+idx);
            plRes(50+idx+2) = plRes(50+idx);        
        elseif j == 7
            [val, idx] = max(plRes(50:70));
            plRes(50+idx-2) = plRes(50+idx);
            plRes(50+idx-1) = plRes(50+idx);
            plRes(50+idx) = plRes(50+idx);
            plRes(50+idx+1) = plRes(50+idx);
            plRes(50+idx+2) = plRes(50+idx);
        end
        
        plot(plPar, plRes)

%         plot(pProfile.pPar.names{j},pProfile.pRes.names{j})
        hold on
        if i == numpHtested
            title(setup.params{j})
        end
    end
    xlim([-3 3])
    if j == length(setup.params)
        hL = legend('pH1','pH2','pH3','pH4','pH5','pH6','pH7','pH8','pH9','pH10');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
    end
    
end
suptitle('PLA GAPDH_{reverse}')
% savefig


%% (100) After meeting 2020-04-20: Evaluate the thing of the bias
% We'll evaluate if the results given after the regularization bias are
% robust. For this:
% (1) Fix vm and estimate parameters. See if results can be reproduced.
%       - compare to data fit
%       - see also that Km valeus do not go too crazy
% (2) Regularize NADH to a different value


%% (101) Fix vm and re estimate parameters
% Mostly copied from section (2). Adjustments in the parameter value below
% in this code:
%   -> 'simSys.m' and in 'calcRates.m': call simple vm
%   -> 'pEst_gapdhr.m' and 'costfun_pH.m' fix vm and divide by DF.
setup.constantVm = 1;
% setup.constantVm = 0;
setup.costfun2 = 0;

% %% (2a) Simple case. Non regularized. DF1.
% parameters estimated changed from the power idea to ease the calculation
% of the confidence intervals. Seemed it was harder to understand.
setup.DFstudy = 4;
setup.costfun = 1;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1 = zeros(numpH,plength);
pcis_DF1 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0;
% setup.selectedLambda = 0.01;
% setup.selectedLambda = 0.0000001;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF1(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF1(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
errorArrayDF1 = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        % check only one vmax
        if setup.constantVm == 1
            data.chosenVmax = 0.0021/DF(i,j); % fixing at the maximum
%             data.chosenVmax = 0.0009/DF(i,j); % fixing at the minimum
%             data.chosenVmax = 0.0015/DF(i,j); % middle value
        else
        end
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF1(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
        x_temp2 = xres; errorArrayDF1{i,j} = costfun_pH(x_temp2,data,setup);
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (2b) Study the cost function. Test use of DF1 and DF2
setup.DFstudy = [3 4];
setup.costfun = 2;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF12 = zeros(numpH,plength);
pcis_DF12 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0;
% setup.selectedLambda = 0.01;
% setup.selectedLambda = 0.0000001;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF12(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF12(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
errorArrayDF12 = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        % check only one vmax
        if setup.constantVm == 1
            data.chosenVmax = 0.0021/DF(i,j); % fixing at the maximum
%             data.chosenVmax = 0.0009/DF(i,j); % fixing at the minimum
%             data.chosenVmax = 0.0015/DF(i,j); % middle value
        else
        end
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF12(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
        x_temp2 = xres; errorArrayDF12{i,j} = costfun_pH(x_temp2,data,setup);
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (2c) Study the cost function. Test use of DF1, DF2, DF3 and DF4
setup.DFstudy = [1 2 3 4];
setup.costfun = 3;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1234 = zeros(numpH,plength);
pcis_DF1234 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
setup.selectedLambda = 0;
% setup.selectedLambda = 0.01;
% setup.selectedLambda = 0.0000001;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF1234(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF1234(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
errorArrayDF1248 = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        % check only one vmax
        if setup.constantVm == 1
            data.chosenVmax = 0.0021/DF(i,j); % fixing at the maximum
%             data.chosenVmax = 0.0009/DF(i,j); % fixing at the minimum
%             data.chosenVmax = 0.0015/DF(i,j); % middle value
        else
        end
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF12(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
        x_temp2 = xres; errorArrayDF1248{i,j} = costfun_pH(x_temp2,data,setup);
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (2a,2b,2c) Visualization
figure
for i = 1:plength
    subplot(3,3,i)
    plot(pHvals,pvals_DF1(:,i))
    hold on
    plot(pHvals,pvals_DF12(:,i))
    hold on
    plot(pHvals,pvals_DF1234(:,i))
    hold on
    ylim([-1 1])
    title(setup.params{i})
    if i == plength
        hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Plot parameter estimates vs pH. Non-regularized')

figure
for i = 1:plength
    subplot(3,3,i)
    errorbar(pHvals,pvals_DF1(:,i),pcis_DF1(:,i))
    hold on
    errorbar(pHvals,pvals_DF12(:,i),pcis_DF12(:,i))
    hold on
    errorbar(pHvals,pvals_DF1234(:,i),pcis_DF1234(:,i))
    hold on
    x = [6 8];
    y = [0 0];
    line(x,y,'Color','black','LineStyle','--')
    hold on
%     ylim([-1 1])
    title(setup.params{i})
    if i == plength
        hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Errorbar parameter estimates vs pH. Non-regularized')


% %%
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorNonReg.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorReg.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorLimitUp.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorLimitUp_check.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorLimitIntermediate.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorLimitIntermediate_check.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorLimitDown.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorLimitDown_check.mat'];
% % save(saveName,'errorArrayDF1','errorArrayDF12','errorArrayDF1248');

%% Display all results

% recall data
load('gapdhr_errorNonReg.mat');
errNonReg.DF1 = errorArrayDF1;
errNonReg.DF12 = errorArrayDF12;
errNonReg.DF1248 = errorArrayDF1248;

load('gapdhr_errorReg.mat');
errReg.DF1 = errorArrayDF1;
errReg.DF12 = errorArrayDF12;
errReg.DF1248 = errorArrayDF1248;

load('gapdhr_errorLimitUp.mat');
% load('gapdhr_errorLimitUp_check.mat');
errLimUp.DF1 = errorArrayDF1;
errLimUp.DF12 = errorArrayDF12;
errLimUp.DF1248 = errorArrayDF1248;

load('gapdhr_errorLimitIntermediate.mat');
% load('gapdhr_errorLimitIntermediate_check.mat');
errLimInt.DF1 = errorArrayDF1;
errLimInt.DF12 = errorArrayDF12;
errLimInt.DF1248 = errorArrayDF1248;

load('gapdhr_errorLimitDown.mat');
% load('gapdhr_errorLimitDown_check.mat');
errLimDown.DF1 = errorArrayDF1;
errLimDown.DF12 = errorArrayDF12;
errLimDown.DF1248 = errorArrayDF1248;

% re-structure data
errDF1_nonReg = zeros(10,1);
errDF1_reg = zeros(10,1);
errDF1_up = zeros(10,1);
errDF1_int = zeros(10,1);
errDF1_low = zeros(10,1);

errDF12_nonReg = zeros(10,1);
errDF12_reg = zeros(10,1);
errDF12_up = zeros(10,1);
errDF12_int = zeros(10,1);
errDF12_low = zeros(10,1);

errDF1248_nonReg = zeros(10,1);
errDF1248_reg = zeros(10,1);
errDF1248_up = zeros(10,1);
errDF1248_int = zeros(10,1);
errDF1248_low = zeros(10,1);

for i = 1:10
%     errDF1_nonReg(i) = abs(sum(errNonReg.DF1{i,4}));
%     errDF1_reg(i) = abs(sum(errReg.DF1{i,4}));
%     errDF1_up(i) = abs(sum(errLimUp.DF1{i,4}));
%     errDF1_int(i) = abs(sum(errLimInt.DF1{i,4}));
%     errDF1_low(i) = abs(sum(errLimDown.DF1{i,4}));
%     
%     errDF12_nonReg(i) = abs(sum(errNonReg.DF12{i,4}));
%     errDF12_reg(i) = abs(sum(errReg.DF12{i,4}));
%     errDF12_up(i) = abs(sum(errLimUp.DF12{i,4}));
%     errDF12_int(i) = abs(sum(errLimInt.DF12{i,4}));
%     errDF12_low(i) = abs(sum(errLimDown.DF12{i,4}));
%     
%     errDF1248_nonReg(i) = abs(sum(errNonReg.DF1248{i,4}));
%     errDF1248_reg(i) = abs(sum(errReg.DF1248{i,4}));
%     errDF1248_up(i) = abs(sum(errLimUp.DF1248{i,4}));
%     errDF1248_int(i) = abs(sum(errLimInt.DF1248{i,4}));
%     errDF1248_low(i) = abs(sum(errLimDown.DF1248{i,4}));
    
    errDF1_nonReg(i) = sum(abs(errNonReg.DF1{i,4}(1:26)));
    errDF1_reg(i) = sum(abs(errReg.DF1{i,4}(1:26)));
    errDF1_up(i) = sum(abs(errLimUp.DF1{i,4}(1:26)));
    errDF1_int(i) = sum(abs(errLimInt.DF1{i,4}(1:26)));
    errDF1_low(i) = sum(abs(errLimDown.DF1{i,4}(1:26)));
    
    errDF12_nonReg(i) = sum(abs(errNonReg.DF12{i,4}(1:52)));
    errDF12_reg(i) = sum(abs(errReg.DF12{i,4}(1:52)));
    errDF12_up(i) = sum(abs(errLimUp.DF12{i,4}(1:52)));
    errDF12_int(i) = sum(abs(errLimInt.DF12{i,4}(1:52)));
    errDF12_low(i) = sum(abs(errLimDown.DF12{i,4}(1:52)));
    
    errDF1248_nonReg(i) = sum(abs(errNonReg.DF1248{i,4}(1:104)));
    errDF1248_reg(i) = sum(abs(errReg.DF1248{i,4}(1:104)));
    errDF1248_up(i) = sum(abs(errLimUp.DF1248{i,4}(1:104)));
    errDF1248_int(i) = sum(abs(errLimInt.DF1248{i,4}(1:104)));
    errDF1248_low(i) = sum(abs(errLimDown.DF1248{i,4}(1:104)));
end

% %%
% plotting
figure

% error DF1
subplot(2,2,1)
plot(pHvals, errDF1_nonReg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_reg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_up,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_int,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_low,'.-','Linewidth',1.5,'MarkerSize',10), hold on
% legend('errDF1_{nonReg}','errDF1_{reg}','errDF1_{up}','errDF1_{int}','errDF1_{low}','Location','SouthOutside','Orientation','Horizontal')
legend('errDF1_{nonReg}','errDF1_{reg}','errDF1_{up}','errDF1_{int}','errDF1_{low}','Location','EastOutside')
title('Error vs pH: DF1')

% error DF12
subplot(2,2,2)
plot(pHvals, errDF12_nonReg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_reg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_up,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_int,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_low,'.-','Linewidth',1.5,'MarkerSize',10), hold on
% legend('errDF12_{nonReg}','errDF12_{reg}','errDF12_{up}','errDF12_{int}','errDF12_{low}','Location','SouthOutside','Orientation','Horizontal')
legend('errDF12_{nonReg}','errDF12_{reg}','errDF12_{up}','errDF12_{int}','errDF12_{low}','Location','EastOutside')
title('Error vs pH: DF12')

% error DF1248
subplot(2,2,3)
plot(pHvals, errDF1248_nonReg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_reg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_up,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_int,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_low,'.-','Linewidth',1.5,'MarkerSize',10), hold on
% legend('errDF1248_{nonReg}','errDF1248_{reg}','errDF1248_{up}','errDF1248_{int}','errDF1248_{low}','Location','SouthOutside','Orientation','Horizontal')
legend('errDF1248_{nonReg}','errDF1248_{reg}','errDF1248_{up}','errDF1248_{int}','errDF1248_{low}','Location','EastOutside')
title('Error vs pH: DF1248')

% Total error
data2plot = [sum(errDF1_nonReg) sum(errDF1_reg) sum(errDF1_up) sum(errDF1_int) sum(errDF1_low);...
            sum(errDF12_nonReg) sum(errDF12_reg) sum(errDF12_up) sum(errDF12_int) sum(errDF12_low);...
            sum(errDF1248_nonReg) sum(errDF1248_reg)  sum(errDF1248_up) sum(errDF1248_int) sum(errDF1248_low)];
clabels = categorical({'DF1','DF12','DF1248'});
subplot(2,2,4)
bar(clabels,data2plot)
% bar(data2plot)
% legend('err_{nonReg}','err_{reg}','err_{up}','err_{int}','err_{low}','Location','SouthOutside','Orientation','Horizontal')
legend('err_{nonReg}','err_{reg}','err_{up}','err_{int}','err_{low}','Location','EastOutside')
title('overall error different routines')

% hold on
% suptitle('Errors obtained with the different methods')

% %%
% saveFigName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_overallError.fig'];
% savefig(1,saveFigName);


%% (201) Regularize NADH to a different value
% Mostly copied from section (101). Knadh is regularized to values -2 -1 +1
% +2. To begin with, lambda is also regularized to 0.01. Same error
% comparison/measurement as previous section (101).
% - Where to add the change? 
setup.constantVm = 0;
setup.costfun2 = 1;

% %% (2a) Simple case. Non regularized. DF1.
% parameters estimated changed from the power idea to ease the calculation
% of the confidence intervals. Seemed it was harder to understand.
setup.DFstudy = 4;
setup.costfun = 1;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1 = zeros(numpH,plength);
pcis_DF1 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
% setup.selectedLambda = 0;
setup.selectedLambda = 0.01;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF1(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF1(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
errorArrayDF1 = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        % check only one vmax
        if setup.constantVm == 1
            data.chosenVmax = 0.0021/DF(i,j); % fixing at the maximum
            data.chosenVmax = 0.0009/DF(i,j); % fixing at the minimum
            data.chosenVmax = 0.0015/DF(i,j); % middle value
        else
        end
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF1(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
        x_temp2 = xres; errorArrayDF1{i,j} = costfun_pH(x_temp2,data,setup);
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (2b) Study the cost function. Test use of DF1 and DF2
setup.DFstudy = [3 4];
setup.costfun = 2;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF12 = zeros(numpH,plength);
pcis_DF12 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
% setup.selectedLambda = 0;
setup.selectedLambda = 0.01;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF12(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF12(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
errorArrayDF12 = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        % check only one vmax
        if setup.constantVm == 1
            data.chosenVmax = 0.0021/DF(i,j); % fixing at the maximum
            data.chosenVmax = 0.0009/DF(i,j); % fixing at the minimum
            data.chosenVmax = 0.0015/DF(i,j); % middle value
        else
        end
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF12(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
        x_temp2 = xres; errorArrayDF12{i,j} = costfun_pH(x_temp2,data,setup);
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (2c) Study the cost function. Test use of DF1, DF2, DF3 and DF4
setup.DFstudy = [1 2 3 4];
setup.costfun = 3;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 10).
plength = length(setup.params);
% % % % plength = length(setup.params)-1;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
ub(end) = 10;
lb = -1*ones(1,plength);
lb(end) = -10;
numpH = numpHtested;
pvals_DF1234 = zeros(numpH,plength);
pcis_DF1234 = zeros(numpH,plength);
DFstudy = setup.DFstudy;
% setup.selectedLambda = 0;
setup.selectedLambda = 0.01;

% parameter estiamtion for each pH value
for i = 1:numpH
    % select required data
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    options = optimset('Display','off');
%     options = optimset('Display','iter');
%     options = optimset('Display','iter','PlotFcns',{@optimplotx,@optimplotfval});

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH,x_temp,lb,ub,options,data,setup);
    t = toc;
    pvals_DF1234(i,:) = xres;
    fprintf('Pest finished for pH #%d, time %d s\n',i,t);
    
    %Estimated parameter variance, covariance matrix and confidence
    %intervals
    lN = length(DFstudy);
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
%     stdp = 100*stdp'./p; %[%]
    pcis_DF1234(i,:) = stdp; % confidence intervals
%     disp([' Keq:    ', num2str(xres(1)), ' +/- ', num2str(stdp(1))]);
%     disp([' Kgap:   ', num2str(xres(2)), ' +/- ', num2str(stdp(2))]);
%     disp([' Kbpg:   ', num2str(xres(3)), ' +/- ', num2str(stdp(3))]);
%     disp([' Knad:   ', num2str(xres(4)), ' +/- ', num2str(stdp(4))]);
%     disp([' Knadh:  ', num2str(xres(5)), ' +/- ', num2str(stdp(5))]);
%     disp([' Vm:     ', num2str(xres(6)), ' +/- ', num2str(stdp(6))]);
%     disp([' linkR:  ', num2str(xres(7)), ' +/- ', num2str(stdp(7))]);
%     % Something like this to get the confidence intervals % p95 =[pmin-psigma*tcr; pmin+psigma*tcr]; %+-95%    
%     figure
%     heatmap(varp);
   
%     % intermediate display of results (if needed)
%     data.chosenNADini = data.NADH{4}(1);
%     data.chosenLink = data.DF(1,4);
%     data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,4);
%     [simResult] = simSys(xres,data,setup);
%     x_temp = xres;
%     [vObs,~] = calcRates(x_temp,simResult,data,setup);    
% 
%     figure,
% 
%     subplot(1,2,1)
%     plot(simResult.t,simResult.y(:,8),'-')
%     hold on
%     plot(data.time{i,4},data.conc_mean{i,4},'r+')
%     title('NADH')
% 
%     subplot(1,2,2)
%     plot(simResult.t,vObs,'-')
%     hold on
%     plot(data.time{i,4},-data.RRs{i,4},'r+')
%     title('v_{GAPDHr}')   
%     
%     drawnow()
    
end

% full system simulation
simulations = cell(numpHtested,DFs);
simRates = cell(numpHtested,DFs);
errorArrayDF1248 = cell(numpHtested,DFs);
for i = 1:numpHtested
    for j = 1:DFs       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
        % check only one vmax
        if setup.constantVm == 1
            data.chosenVmax = 0.0021/DF(i,j); % fixing at the maximum
            data.chosenVmax = 0.0009/DF(i,j); % fixing at the minimum
            data.chosenVmax = 0.0015/DF(i,j); % middle value
        else
        end
        
%         xres = MSresult_temp{i}.b;
        xres = pvals_DF12(i,:);
        
        [simResult] = simSys(xres,data,setup);
        simulations{i,j} = simResult;
        [vObs,~] = calcRates(xres,simResult,data,setup);   
        simRates{i,j} = vObs;
        x_temp2 = xres; errorArrayDF1248{i,j} = costfun_pH(x_temp2,data,setup);
    end
    
    figure
    subplot(1,2,1)
    for j = 1:DFs
        simRes = simulations{i,j};
        plot(simRes.t,simRes.y(:,8),'-')
        hold on
        plot(data.time{i,j},data.conc_mean{i,j},'k+')
    end
    title('NADH')
    subplot(1,2,2)
    for j = 1:DFs
        simRes = simulations{i,j};
        simRRs = simRates{i,j};
        plot(simRes.t,simRRs,'-')
        hold on
        plot(data.time{i,j},-data.RRs{i,j},'k+')
    end
    title('v_{GAPDHr}')   
    suptitle(erase(sprintf('pH = %d',pHarray(i)),"0000e+00"));

end


% %% (2a,2b,2c) Visualization
figure
for i = 1:plength
    subplot(3,3,i)
    plot(pHvals,pvals_DF1(:,i))
    hold on
    plot(pHvals,pvals_DF12(:,i))
    hold on
    plot(pHvals,pvals_DF1234(:,i))
    hold on
    ylim([-1 1])
    title(setup.params{i})
    if i == plength
        hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Plot parameter estimates vs pH. Non-regularized')

figure
for i = 1:plength
    subplot(3,3,i)
    errorbar(pHvals,pvals_DF1(:,i),pcis_DF1(:,i))
    hold on
    errorbar(pHvals,pvals_DF12(:,i),pcis_DF12(:,i))
    hold on
    errorbar(pHvals,pvals_DF1234(:,i),pcis_DF1234(:,i))
    hold on
    x = [6 8];
    y = [0 0];
    line(x,y,'Color','black','LineStyle','--')
    hold on
%     ylim([-1 1])
    title(setup.params{i})
    if i == plength
        hL = legend('pvalsDF1','pvalsDF12','pvalsDF1234');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Errorbar parameter estimates vs pH. Non-regularized')


% %%
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorReg_Knadh_0.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorReg_Knadh_Reg_n1.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorReg_Knadh_Reg_n.5.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorReg_Knadh_Reg_p.5.mat'];
% % saveName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_errorReg_Knadh_Reg_p1.mat'];
% % save(saveName,'errorArrayDF1','errorArrayDF12','errorArrayDF1248');
% % close(1:30)

%% Display all results

% recall data
load('gapdhr_errorReg.mat');
errReg.DF1 = errorArrayDF1;
errReg.DF12 = errorArrayDF12;
errReg.DF1248 = errorArrayDF1248;

load('gapdhr_errorReg_Knadh_0.mat');
errReg_0.DF1 = errorArrayDF1;
errReg_0.DF12 = errorArrayDF12;
errReg_0.DF1248 = errorArrayDF1248;

load('gapdhr_errorReg_Knadh_Reg_n1.mat');
errReg_n1.DF1 = errorArrayDF1;
errReg_n1.DF12 = errorArrayDF12;
errReg_n1.DF1248 = errorArrayDF1248;

load('gapdhr_errorReg_Knadh_Reg_n.5.mat');
errReg_n_5.DF1 = errorArrayDF1;
errReg_n_5.DF12 = errorArrayDF12;
errReg_n_5.DF1248 = errorArrayDF1248;

load('gapdhr_errorReg_Knadh_Reg_p.5.mat');
errReg_p_5.DF1 = errorArrayDF1;
errReg_p_5.DF12 = errorArrayDF12;
errReg_p_5.DF1248 = errorArrayDF1248;

load('gapdhr_errorReg_Knadh_Reg_p1.mat');
errReg_p1.DF1 = errorArrayDF1;
errReg_p1.DF12 = errorArrayDF12;
errReg_p1.DF1248 = errorArrayDF1248;

% re-structure data
errDF1_Reg = zeros(10,1);
errDF1_0 = zeros(10,1);
errDF1_n1 = zeros(10,1);
errDF1_n_5 = zeros(10,1);
errDF1_p_5 = zeros(10,1);
errDF1_p1 = zeros(10,1);

errDF12_Reg = zeros(10,1);
errDF12_0 = zeros(10,1);
errDF12_n1 = zeros(10,1);
errDF12_n_5 = zeros(10,1);
errDF12_p_5 = zeros(10,1);
errDF12_p1 = zeros(10,1);

errDF1248_Reg = zeros(10,1);
errDF1248_0 = zeros(10,1);
errDF1248_n1 = zeros(10,1);
errDF1248_n_5 = zeros(10,1);
errDF1248_p_5 = zeros(10,1);
errDF1248_p1 = zeros(10,1);

for i = 1:10
%     errDF1_Reg(i) = abs(sum(errReg.DF1{i,4}));
%     errDF1_0(i) = abs(sum(errReg_0.DF1{i,4}));
%     errDF1_n1(i) = abs(sum(errReg_n1.DF1{i,4}));
%     errDF1_n_5(i) = abs(sum(errReg_n_5.DF1{i,4}));
%     errDF1_p_5(i) = abs(sum(errReg_p_5.DF1{i,4}));
%     errDF1_p1(i) = abs(sum(errReg_p1.DF1{i,4}));
%     
%     errDF12_Reg(i) = abs(sum(errReg.DF12{i,4}));
%     errDF12_0(i) = abs(sum(errReg_0.DF12{i,4}));
%     errDF12_n1(i) = abs(sum(errReg_n1.DF12{i,4}));
%     errDF12_n_5(i) = abs(sum(errReg_n_5.DF12{i,4}));
%     errDF12_p_5(i) = abs(sum(errReg_p_5.DF12{i,4}));
%     errDF12_p1(i) = abs(sum(errReg_p1.DF12{i,4})); 
%     
%     errDF1248_Reg(i) = abs(sum(errReg.DF1248{i,4}));
%     errDF1248_0(i) = abs(sum(errReg_0.DF1248{i,4}));
%     errDF1248_n1(i) = abs(sum(errReg_n1.DF1248{i,4}));
%     errDF1248_n_5(i) = abs(sum(errReg_n_5.DF1248{i,4}));
%     errDF1248_p_5(i) = abs(sum(errReg_p_5.DF1248{i,4}));
%     errDF1248_p1(i) = abs(sum(errReg_p1.DF1248{i,4}));
    errDF1_Reg(i) = sum(abs(errReg.DF1{i,4}(1:26)));
    errDF1_0(i) = sum(abs(errReg_0.DF1{i,4}(1:26)));
    errDF1_n1(i) = sum(abs(errReg_n1.DF1{i,4}(1:26)));
    errDF1_n_5(i) = sum(abs(errReg_n_5.DF1{i,4}(1:26)));
    errDF1_p_5(i) = sum(abs(errReg_p_5.DF1{i,4}(1:26)));
    errDF1_p1(i) = sum(abs(errReg_p1.DF1{i,4}(1:26)));
    
    errDF12_Reg(i) = sum(abs(errReg.DF12{i,4}(1:52)));
    errDF12_0(i) = sum(abs(errReg_0.DF12{i,4}(1:52)));
    errDF12_n1(i) = sum(abs(errReg_n1.DF12{i,4}(1:52)));
    errDF12_n_5(i) = sum(abs(errReg_n_5.DF12{i,4}(1:52)));
    errDF12_p_5(i) = sum(abs(errReg_p_5.DF12{i,4}(1:52)));
    errDF12_p1(i) = sum(abs(errReg_p1.DF12{i,4}(1:52))); 
    
    errDF1248_Reg(i) = sum(abs(errReg.DF1248{i,4}(1:104)));
    errDF1248_0(i) = sum(abs(errReg_0.DF1248{i,4}(1:104)));
    errDF1248_n1(i) = sum(abs(errReg_n1.DF1248{i,4}(1:104)));
    errDF1248_n_5(i) = sum(abs(errReg_n_5.DF1248{i,4}(1:104)));
    errDF1248_p_5(i) = sum(abs(errReg_p_5.DF1248{i,4}(1:104)));
    errDF1248_p1(i) = sum(abs(errReg_p1.DF1248{i,4}(1:104)));
end


% %%
% plotting
figure

% error DF1
subplot(2,2,1)
plot(pHvals, errDF1_Reg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_0,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_n1,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_n_5,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_p_5,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1_p1,'.-','Linewidth',1.5,'MarkerSize',10), hold on
legend('errDF1_{Reg}','errDF1_{0}','errDF1_{n1}','errDF1_{n.5}','errDF1_{p.5}','errDF1_{p1}','Location','EastOutside')
title('Error vs pH: DF1')

% error DF12
subplot(2,2,2)
plot(pHvals, errDF12_Reg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_0,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_n1,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_n_5,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_p_5,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF12_p1,'.-','Linewidth',1.5,'MarkerSize',10), hold on
legend('errDF12_{Reg}','errDF12_{0}','errDF12_{n1}','errDF12_{n.5}','errDF12_{p.5}','errDF12_{p1}','Location','EastOutside')
title('Error vs pH: DF12')

% error DF1248
subplot(2,2,3)
plot(pHvals, errDF1248_Reg,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_0,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_n1,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_n_5,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_p_5,'.-','Linewidth',1.5,'MarkerSize',10), hold on
plot(pHvals, errDF1248_p1,'.-','Linewidth',1.5,'MarkerSize',10), hold on
legend('errDF1248_{Reg}','errDF1248_{0}','errDF1248_{n1}','errDF1248_{n.5}','errDF1248_{p.5}','errDF1248_{p1}','Location','EastOutside')
title('Error vs pH: DF1248')

% Total error
data2plot = [sum(errDF1_Reg) sum(errDF1_0) sum(errDF1_n1) sum(errDF1_n_5) sum(errDF1_p_5) sum(errDF1_p1);...
            sum(errDF12_Reg) sum(errDF12_0) sum(errDF12_n1) sum(errDF12_n_5) sum(errDF12_p_5) sum(errDF12_p1);...
            sum(errDF1248_Reg) sum(errDF1248_0) sum(errDF1248_n1) sum(errDF1248_n_5) sum(errDF1248_p_5) sum(errDF1248_p1)];
clabels = categorical({'DF1','DF12','DF1248'});
subplot(2,2,4)
bar(clabels,data2plot)
% bar(data2plot)
% legend('err_{nonReg}','err_{reg}','err_{up}','err_{int}','err_{low}','Location','SouthOutside','Orientation','Horizontal')
legend('errDF1248_{Reg}','errDF1248_{0}','errDF1248_{n1}','errDF1248_{n.5}','errDF1248_{p.5}','errDF1248_{p1}','Location','EastOutside')
title('overall error different routines')

% hold on
% suptitle('Errors obtained with the different methods')


%%
saveFigName = ['results\', setup.enzymeName, '\', setup.enzymeName, '_overallError_KnadhStudy.fig'];
savefig(33,saveFigName);


%% (301) Big scale optimization. Unique Km estimated for the entire profile. Keq, vm and linkReaction change with pH
% The parameter array will be increased so that an independent parameter is
% present for each pH calue for keq, vm and linkR.
% first we start with no regularization.

% %% (2a) Simple case. Non regularized. DF1.
% parameters estimated changed from the power idea to ease the calculation
% of the confidence intervals. Seemed it was harder to understand.
setup.DFstudy = 4;
setup.costfun = 1;
setup.selectedLambda = 0;
setup.numpHtested = numpHtested;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 34).
% The array of parameters is increased to a length of 100. keq, vm and linkR
% are picked specifically for each case
plength = 34;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
lb = -1*ones(1,plength);
ub([7,10,13,16,19,22,25,28,31,34]) = 10;
lb([7,10,13,16,19,22,25,28,31,34]) = -10;
options = optimset('Display','iter');

tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_Kmconstant,x_temp,lb,ub,options,data,setup);
t = toc;
%Estimated parameter variance, covariance matrix and confidence
%intervals
%     N = ...;
%     Jacobian = full(Jacobian);  
%     varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
%     stdp = sqrt(diag(varp));
%     pcis_DF1(i,:) = stdp; % confidence intervals

% simulation to get the overall error
error_sum = resnorm;
error_pH1 = residual(1:26);
error_pH2 = residual(27:52);
error_pH3 = residual(53:78);
error_pH4 = residual(79:104);
error_pH5 = residual(105:130);
error_pH6 = residual(131:156);
error_pH7 = residual(157:182);
error_pH8 = residual(183:208);
error_pH9 = residual(209:234);
error_pH10 = residual(235:260);

error_Kmconstant = cell(10,1);
error_Kmconstant{1} = error_pH1;
error_Kmconstant{2} = error_pH2;
error_Kmconstant{3} = error_pH3;
error_Kmconstant{4} = error_pH4;
error_Kmconstant{5} = error_pH5;
error_Kmconstant{6} = error_pH6;
error_Kmconstant{7} = error_pH7;
error_Kmconstant{8} = error_pH8;
error_Kmconstant{9} = error_pH9;
error_Kmconstant{10} = error_pH10;

save('gapdhr_error_Kmconstant.mat','error_sum','error_Kmconstant');


%% visualization
% reorganise (not regularized)
pReg = [0.678	0.0081	-0.0452	0.0085	-0.2049	0.1972	0.0095;
        0.7315	0.0082	-0.0426	0.0079	-0.1976	0.1887	-0.0068;
        0.7435	0.0063	-0.0422	0.0062	-0.1605	0.1947	-0.0097;
        0.7244	0.003	-0.0362	0.0036	-0.1046	0.1823	0.0001;
        0.8479	-0.0084	-0.0192	-0.0149	0.2264	0.157	0.0002;
        0.8873	-0.0158	-0.0255	-0.0388	0.4656	0.1992	0.0056;
        0.8476	-0.0188	-0.0201	-0.0507	0.5604	0.2203	0.0091;
        0.8292	-0.0205	-0.0207	-0.0554	0.584	0.2182	-0.0085;
        0.7975	-0.0237	-0.0186	-0.0644	0.6265	0.2097	-0.0094;
        0.7132	-0.0172	-0.0146	-0.0423	0.5111	0.189	-0.01];

load('E:\tempWork\pHdata\results\gapdhr\gapdhr_xres_full.mat');
pKmfixed = [...
        xres(5) xres(1) xres(2) xres(3) xres(4) xres(6) xres(7);
        xres(8) xres(1) xres(2) xres(3) xres(4) xres(9) xres(10);
        xres(11) xres(1) xres(2) xres(3) xres(4) xres(12) xres(13);
        xres(14) xres(1) xres(2) xres(3) xres(4) xres(15) xres(16);
        xres(17) xres(1) xres(2) xres(3) xres(4) xres(18) xres(19);
        xres(20) xres(1) xres(2) xres(3) xres(4) xres(21) xres(22);
        xres(23) xres(1) xres(2) xres(3) xres(4) xres(24) xres(25);
        xres(26) xres(1) xres(2) xres(3) xres(4) xres(27) xres(28);
        xres(29) xres(1) xres(2) xres(3) xres(4) xres(30) xres(31);
        xres(32) xres(1) xres(2) xres(3) xres(4) xres(33) xres(34)];

load('E:\tempWork\pHdata\results\gapdhr\gapdhr_errorReg.mat');
errorReg = cell(10,1);
errorReg{1} = errorArrayDF1{1,4}(1:26);
errorReg{2} = errorArrayDF1{2,4}(1:26);
errorReg{3} = errorArrayDF1{3,4}(1:26);
errorReg{4} = errorArrayDF1{4,4}(1:26);
errorReg{5} = errorArrayDF1{5,4}(1:26);
errorReg{6} = errorArrayDF1{6,4}(1:26);
errorReg{7} = errorArrayDF1{7,4}(1:26);
errorReg{8} = errorArrayDF1{8,4}(1:26);
errorReg{9} = errorArrayDF1{9,4}(1:26);
errorReg{10} = errorArrayDF1{10,4}(1:26);


% Option A. Call from the saving (something went unexpected. Did not fin
% reason in the costfun_pH_Kmconstant yet)
load('E:\tempWork\pHdata\results\gapdhr\gapdhr_error_Kmconstant.mat');
errorKmfixed = error_Kmconstant;
% errorcalc =
%     0.0297
%     0.0328
%     0.0191
%     0.0106
%     0.0087
%     0.0281
%     0.0375
%     0.0439
%     0.0518
%     0.0442
% % % Option B. Use costfun_pH
% % setup.DFstudy = 4;
% % setup.costfun = 1;
% % errorKmfixed = cell(10,1);
% % for i = 1:10
% %     data.Vmaxs = data.Vmax(i,:);
% %     data.NADH = data.conc_mean(i,:);
% %     data.Vprofs = data.RRs(i,:);
% %     data.tempTime = data.time(i,:);
% %     x_temp3 = pKmfixed(i,:);
% %     temp = costfun_pH(x_temp3,data,setup);
% %     errorKmfixed{i} = temp(1:26);
% % end

% compare resulting parameters
figure
for i = 1:7
    subplot(3,3,i)
    plot(pHvals,pReg(:,i))
    hold on
    plot(pHvals,pKmfixed(:,i)) 
    hold on
    ylim([-1 1])
    title(setup.params{i})
    if i == 7
        hL = legend('pReg','pKmfixed');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Plot parameter estimates vs pH. Km is fixed')

% compare error
errorRegSummed = zeros(10,1);
for i = 1:10
    errorRegSummed(i) = sum(abs(errorReg{i}));
end
errorKmfixedSummed = zeros(10,1);
for i = 1:10
    errorKmfixedSummed(i) = sum(abs(errorKmfixed{i}));
end

figure
plot(pHvals,errorRegSummed,'.-')
hold on
plot(pHvals,errorKmfixedSummed,'.-')
legend('errorReg','errorKmfixed')


%% (302) Big scale optimization. 
% Same as (301) but regularization factor tested 0.01

% The parameter array will be increased so that an independent parameter is
% present for each pH calue for keq, vm and linkR.
% first we start with no regularization.

% %% (2a) Simple case. Non regularized. DF1.
% parameters estimated changed from the power idea to ease the calculation
% of the confidence intervals. Seemed it was harder to understand.
setup.DFstudy = 4;
setup.costfun = 1;
% setup.selectedLambda = 0;
setup.selectedLambda = 0.01;
setup.numpHtested = numpHtested;

% parameter estimation setup. Boundaries are set to 1 order of magnitude
% (positive or negative) except for the linkin reaction (to 34).
% The array of parameters is increased to a length of 100. keq, vm and linkR
% are picked specifically for each case
plength = 34;
x_temp = zeros(1,plength);
ub = 1*ones(1,plength);
lb = -1*ones(1,plength);
ub([7,10,13,16,19,22,25,28,31,34]) = 10;
lb([7,10,13,16,19,22,25,28,31,34]) = -10;
options = optimset('Display','iter');

tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(@costfun_pH_Kmconstant,x_temp,lb,ub,options,data,setup);
t = toc;
%Estimated parameter variance, covariance matrix and confidence
%intervals
%     N = ...;
%     Jacobian = full(Jacobian);  
%     varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
%     stdp = sqrt(diag(varp));
%     pcis_DF1(i,:) = stdp; % confidence intervals

% simulation to get the overall error
error_sum = resnorm;
error_pH1 = residual(1:26);
error_pH2 = residual(27:52);
error_pH3 = residual(53:78);
error_pH4 = residual(79:104);
error_pH5 = residual(105:130);
error_pH6 = residual(131:156);
error_pH7 = residual(157:182);
error_pH8 = residual(183:208);
error_pH9 = residual(209:234);
error_pH10 = residual(235:260);

error_Kmconstant_reg = cell(10,1);
error_Kmconstant_reg{1} = error_pH1;
error_Kmconstant_reg{2} = error_pH2;
error_Kmconstant_reg{3} = error_pH3;
error_Kmconstant_reg{4} = error_pH4;
error_Kmconstant_reg{5} = error_pH5;
error_Kmconstant_reg{6} = error_pH6;
error_Kmconstant_reg{7} = error_pH7;
error_Kmconstant_reg{8} = error_pH8;
error_Kmconstant_reg{9} = error_pH9;
error_Kmconstant_reg{10} = error_pH10;

% %%
% save('gapdhr_xres_full_reg.mat','xres');
% save('gapdhr_error_Kmconstant_reg.mat','error_sum','error_Kmconstant');


%% visualization
% reorganise (regularized)
pReg = [0.678	0.0081	-0.0452	0.0085	-0.2049	0.1972	0.0095;
        0.7315	0.0082	-0.0426	0.0079	-0.1976	0.1887	-0.0068;
        0.7435	0.0063	-0.0422	0.0062	-0.1605	0.1947	-0.0097;
        0.7244	0.003	-0.0362	0.0036	-0.1046	0.1823	0.0001;
        0.8479	-0.0084	-0.0192	-0.0149	0.2264	0.157	0.0002;
        0.8873	-0.0158	-0.0255	-0.0388	0.4656	0.1992	0.0056;
        0.8476	-0.0188	-0.0201	-0.0507	0.5604	0.2203	0.0091;
        0.8292	-0.0205	-0.0207	-0.0554	0.584	0.2182	-0.0085;
        0.7975	-0.0237	-0.0186	-0.0644	0.6265	0.2097	-0.0094;
        0.7132	-0.0172	-0.0146	-0.0423	0.5111	0.189	-0.01];

load('E:\tempWork\pHdata\results\gapdhr\gapdhr_xres_full_reg.mat');
pKmfixed = [...
        xres(5) xres(1) xres(2) xres(3) xres(4) xres(6) xres(7);
        xres(8) xres(1) xres(2) xres(3) xres(4) xres(9) xres(10);
        xres(11) xres(1) xres(2) xres(3) xres(4) xres(12) xres(13);
        xres(14) xres(1) xres(2) xres(3) xres(4) xres(15) xres(16);
        xres(17) xres(1) xres(2) xres(3) xres(4) xres(18) xres(19);
        xres(20) xres(1) xres(2) xres(3) xres(4) xres(21) xres(22);
        xres(23) xres(1) xres(2) xres(3) xres(4) xres(24) xres(25);
        xres(26) xres(1) xres(2) xres(3) xres(4) xres(27) xres(28);
        xres(29) xres(1) xres(2) xres(3) xres(4) xres(30) xres(31);
        xres(32) xres(1) xres(2) xres(3) xres(4) xres(33) xres(34)];

load('E:\tempWork\pHdata\results\gapdhr\gapdhr_errorReg.mat');
errorReg = cell(10,1);
errorReg{1} = errorArrayDF1{1,4}(1:26);
errorReg{2} = errorArrayDF1{2,4}(1:26);
errorReg{3} = errorArrayDF1{3,4}(1:26);
errorReg{4} = errorArrayDF1{4,4}(1:26);
errorReg{5} = errorArrayDF1{5,4}(1:26);
errorReg{6} = errorArrayDF1{6,4}(1:26);
errorReg{7} = errorArrayDF1{7,4}(1:26);
errorReg{8} = errorArrayDF1{8,4}(1:26);
errorReg{9} = errorArrayDF1{9,4}(1:26);
errorReg{10} = errorArrayDF1{10,4}(1:26);


% Option A. Call from the saving (something went unexpected. Did not fin
% reason in the costfun_pH_Kmconstant yet)
load('E:\tempWork\pHdata\results\gapdhr\gapdhr_error_Kmconstant_reg.mat');
errorKmfixed = error_Kmconstant;
% errorcalc =
%     0.0297
%     0.0328
%     0.0191
%     0.0106
%     0.0087
%     0.0281
%     0.0375
%     0.0439
%     0.0518
%     0.0442
% % Option B. Use costfun_pH
% setup.DFstudy = 4;
% setup.costfun = 1;
% errorKmfixed = cell(10,1);
% for i = 1:10
%     data.Vmaxs = data.Vmax(i,:);
%     data.NADH = data.conc_mean(i,:);
%     data.Vprofs = data.RRs(i,:);
%     data.tempTime = data.time(i,:);
%     x_temp3 = pKmfixed(i,:);
%     temp = costfun_pH(x_temp3,data,setup);
%     errorKmfixed{i} = temp(1:26);
% end



% compare resulting parameters
figure
for i = 1:7
    subplot(3,3,i)
    plot(pHvals,pReg(:,i))
    hold on
    plot(pHvals,pKmfixed(:,i)) 
    hold on
    ylim([-1 1])
    title(setup.params{i})
    if i == 7
        hL = legend('pReg','pKmfixed');
        newPosition = [0.7 0.1 0.2 0.2];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);
        hold off
    end
end
suptitle('Plot parameter estimates vs pH. Km is fixed')

% compare error
errorRegSummed = zeros(10,1);
for i = 1:10
    errorRegSummed(i) = sum(abs(errorReg{i}));
end
errorKmfixedSummed = zeros(10,1);
for i = 1:10
    errorKmfixedSummed(i) = sum(abs(errorKmfixed{i}));
end

figure
plot(pHvals,errorRegSummed,'.-')
hold on
plot(pHvals,errorKmfixedSummed,'.-')
legend('errorReg','errorKmfixed')


%% (303) recalculate the error directly, using the same costfunction
% (1) costfun_pH get parameters reg
% (2) costfun_pH parameters from the km constant, nonreg
% (3) costfun_pH parameters from the km constant, reg

% (1) reorganise (regularized)
pReg = [0.678	0.0081	-0.0452	0.0085	-0.2049	0.1972	0.0095;
        0.7315	0.0082	-0.0426	0.0079	-0.1976	0.1887	-0.0068;
        0.7435	0.0063	-0.0422	0.0062	-0.1605	0.1947	-0.0097;
        0.7244	0.003	-0.0362	0.0036	-0.1046	0.1823	0.0001;
        0.8479	-0.0084	-0.0192	-0.0149	0.2264	0.157	0.0002;
        0.8873	-0.0158	-0.0255	-0.0388	0.4656	0.1992	0.0056;
        0.8476	-0.0188	-0.0201	-0.0507	0.5604	0.2203	0.0091;
        0.8292	-0.0205	-0.0207	-0.0554	0.584	0.2182	-0.0085;
        0.7975	-0.0237	-0.0186	-0.0644	0.6265	0.2097	-0.0094;
        0.7132	-0.0172	-0.0146	-0.0423	0.5111	0.189	-0.01];

% (2) pKmfixed_nonreg
load('E:\tempWork\pHdata\results\gapdhr\gapdhr_xres_full.mat');
pKmfixed_nonreg = [...
        xres(5) xres(1) xres(2) xres(3) xres(4) xres(6) xres(7);
        xres(8) xres(1) xres(2) xres(3) xres(4) xres(9) xres(10);
        xres(11) xres(1) xres(2) xres(3) xres(4) xres(12) xres(13);
        xres(14) xres(1) xres(2) xres(3) xres(4) xres(15) xres(16);
        xres(17) xres(1) xres(2) xres(3) xres(4) xres(18) xres(19);
        xres(20) xres(1) xres(2) xres(3) xres(4) xres(21) xres(22);
        xres(23) xres(1) xres(2) xres(3) xres(4) xres(24) xres(25);
        xres(26) xres(1) xres(2) xres(3) xres(4) xres(27) xres(28);
        xres(29) xres(1) xres(2) xres(3) xres(4) xres(30) xres(31);
        xres(32) xres(1) xres(2) xres(3) xres(4) xres(33) xres(34)];

% (3) pKmfixed_reg
load('E:\tempWork\pHdata\results\gapdhr\gapdhr_xres_full_reg.mat');
pKmfixed_reg = [...
        xres(5) xres(1) xres(2) xres(3) xres(4) xres(6) xres(7);
        xres(8) xres(1) xres(2) xres(3) xres(4) xres(9) xres(10);
        xres(11) xres(1) xres(2) xres(3) xres(4) xres(12) xres(13);
        xres(14) xres(1) xres(2) xres(3) xres(4) xres(15) xres(16);
        xres(17) xres(1) xres(2) xres(3) xres(4) xres(18) xres(19);
        xres(20) xres(1) xres(2) xres(3) xres(4) xres(21) xres(22);
        xres(23) xres(1) xres(2) xres(3) xres(4) xres(24) xres(25);
        xres(26) xres(1) xres(2) xres(3) xres(4) xres(27) xres(28);
        xres(29) xres(1) xres(2) xres(3) xres(4) xres(30) xres(31);
        xres(32) xres(1) xres(2) xres(3) xres(4) xres(33) xres(34)];

% simulations
setup.DFstudy = 4;
setup.costfun = 1;
error_reg = cell(10,1);
errorKmfixed_nonReg = cell(10,1);
errorKmfixed_Reg = cell(10,1);
for i = 1:10
    data.Vmaxs = data.Vmax(i,:);
    data.NADH = data.conc_mean(i,:);
    data.Vprofs = data.RRs(i,:);
    data.tempTime = data.time(i,:);
    
    x_temp1 = pReg(i,:);
    x_temp2 = pKmfixed_nonreg(i,:);
    x_temp3 = pKmfixed_reg(i,:);
    
    temp1 = costfun_pH(x_temp1,data,setup);
    temp2 = costfun_pH(x_temp2,data,setup);
    temp3 = costfun_pH(x_temp3,data,setup);
    
    error_reg{i} = temp1(1:26);
    errorKmfixed_nonReg{i} = temp2(1:26);
    errorKmfixed_Reg{i} = temp3(1:26);    
end


% compare error
error_Reg_Summed = zeros(10,1);
errorKmfixed_nonReg_Summed = zeros(10,1);
errorKmfixed_Reg_Summed = zeros(10,1);
for i = 1:10
    error_Reg_Summed(i) = sum(abs(error_reg{i}));
    errorKmfixed_nonReg_Summed(i) = sum(abs(errorKmfixed_nonReg{i}));
    errorKmfixed_Reg_Summed(i) = sum(abs(errorKmfixed_Reg{i}));
end

figure
plot(pHvals,error_Reg_Summed,'.-')
hold on
plot(pHvals,errorKmfixed_nonReg_Summed,'.-')
hold on
plot(pHvals,errorKmfixed_nonReg_Summed,'.-')
legend('errorReg','errorKmfixed_{nonReg}','errorKmfixed_{Reg}','location','southeast')


%%
setup.DFstudy = 4;
setup.costfun = 1;

error_reg = cell(10,1);
errorKmfixed_nonReg = cell(10,1);
errorKmfixed_Reg = cell(10,1);
for i = 1:10
    for j = 4       
        data.chosenVmax = data.Vmax(i,4)/DF(i,j);
        data.chosenNADini = data.conc_mean{i,j}(1);
        data.chosenLink = data.DF(i,j);
%         data.tempTime = data.time(i,:);
        
        x_temp1 = pReg(i,:);
        x_temp2 = pKmfixed_nonreg(i,:);
        x_temp3 = pKmfixed_reg(i,:);

        temp1 = costfun_pH(x_temp1,data,setup);
        temp2 = costfun_pH(x_temp2,data,setup);
        temp3 = costfun_pH(x_temp3,data,setup);

        error_reg{i} = temp1(1:26);
        errorKmfixed_nonReg{i} = temp2(1:26);
        errorKmfixed_Reg{i} = temp3(1:26);
    end
end

% compare error
error_Reg_Summed = zeros(10,1);
errorKmfixed_nonReg_Summed = zeros(10,1);
errorKmfixed_Reg_Summed = zeros(10,1);
for i = 1:10
    error_Reg_Summed(i) = sum(abs(error_reg{i}));
    errorKmfixed_nonReg_Summed(i) = sum(abs(errorKmfixed_nonReg{i}));
    errorKmfixed_Reg_Summed(i) = sum(abs(errorKmfixed_Reg{i}));
end

figure
plot(pHvals,error_Reg_Summed,'.-')
hold on
plot(pHvals,errorKmfixed_nonReg_Summed,'.-')
hold on
plot(pHvals,errorKmfixed_nonReg_Summed,'.-')
legend('errorReg','errorKmfixed_{nonReg}','errorKmfixed_{Reg}','location','northeast')







% 