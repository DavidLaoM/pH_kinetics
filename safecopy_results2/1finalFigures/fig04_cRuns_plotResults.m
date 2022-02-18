%1 Get data
%2 Make array together
%3 Plot
%4 Save array
%5 Save plot

%% 0
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
pHvals = unique(import_gapdhR.treatedData.pH_corrected);


% % % % %% 1 get data
% % % % % vf
% % % % load('output_vf_run1.mat'); output_vf_run1 = output; clear output
% % % % load('output_vf_run2.mat'); output_vf_run2 = output; clear output
% % % % load('output_vf_run3.mat'); output_vf_run3 = output; clear output
% % % % load('output_vf_run4.mat'); output_vf_run4 = output; clear output
% % % % load('output_vf_run5.mat'); output_vf_run5 = output; clear output
% % % % load('output_vf_run6.mat'); output_vf_run6 = output; clear output
% % % % 
% % % % %kf
% % % % load('output_kf_run1.mat'); output_kf_run1 = output; clear output
% % % % load('output_kf_run2.mat'); output_kf_run2 = output; clear output
% % % % load('output_kf_run3.mat'); output_kf_run3 = output; clear output
% % % % load('output_kf_run4.mat'); output_kf_run4 = output; clear output
% % % % load('output_kf_run5.mat'); output_kf_run5 = output; clear output
% % % % load('output_kf_run6.mat'); output_kf_run6 = output; clear output
% % % % 
% % % % 
% % % % %% 2 make array together
% % % % % vf
% % % % output_vf.xres = [output_vf_run1.xres; output_vf_run2.xres; output_vf_run3.xres; output_vf_run4.xres; output_vf_run5.xres; output_vf_run6.xres];
% % % % output_vf.weightTest = [output_vf_run1.weightTest, output_vf_run2.weightTest, output_vf_run3.weightTest, output_vf_run4.weightTest, output_vf_run5.weightTest, output_vf_run6.weightTest]';
% % % % output_vf.errorData = [output_vf_run1.errorData; output_vf_run2.errorData; output_vf_run3.errorData; output_vf_run4.errorData; output_vf_run5.errorData; output_vf_run6.errorData];
% % % % output_vf.errorHaldane = [output_vf_run1.errorHaldane; output_vf_run2.errorHaldane; output_vf_run3.errorHaldane; output_vf_run4.errorHaldane; output_vf_run5.errorHaldane; output_vf_run6.errorHaldane];
% % % % output_vf.errorRegpars = [output_vf_run1.errorRegpars; output_vf_run2.errorRegpars; output_vf_run3.errorRegpars; output_vf_run4.errorRegpars; output_vf_run5.errorRegpars; output_vf_run6.errorRegpars];
% % % % output_vf.eData = [output_vf_run1.eData; output_vf_run2.eData; output_vf_run3.eData; output_vf_run4.eData; output_vf_run5.eData; output_vf_run6.eData];
% % % % output_vf.eHaldane = [output_vf_run1.eHaldane; output_vf_run2.eHaldane; output_vf_run3.eHaldane; output_vf_run4.eHaldane; output_vf_run5.eHaldane; output_vf_run6.eHaldane];
% % % % clear output_vf_run1 output_vf_run2 output_vf_run3 output_vf_run4 output_vf_run5 output_vf_run6
% % % % 
% % % % %kf
% % % % output_kf.xres = [output_kf_run1.xres; output_kf_run2.xres; output_kf_run3.xres; output_kf_run4.xres; output_kf_run5.xres; output_kf_run6.xres];
% % % % output_kf.weightTest = [output_kf_run1.weightTest, output_kf_run2.weightTest, output_kf_run3.weightTest, output_kf_run4.weightTest, output_kf_run5.weightTest, output_kf_run6.weightTest]';
% % % % output_kf.errorData = [output_kf_run1.errorData; output_kf_run2.errorData; output_kf_run3.errorData; output_kf_run4.errorData; output_kf_run5.errorData; output_kf_run6.errorData];
% % % % output_kf.errorHaldane = [output_kf_run1.errorHaldane; output_kf_run2.errorHaldane; output_kf_run3.errorHaldane; output_kf_run4.errorHaldane; output_kf_run5.errorHaldane; output_kf_run6.errorHaldane];
% % % % output_kf.errorRegpars = [output_kf_run1.errorRegpars; output_kf_run2.errorRegpars; output_kf_run3.errorRegpars; output_kf_run4.errorRegpars; output_kf_run5.errorRegpars; output_kf_run6.errorRegpars];
% % % % output_kf.eData = [output_kf_run1.eData; output_kf_run2.eData; output_kf_run3.eData; output_kf_run4.eData; output_kf_run5.eData; output_kf_run6.eData];
% % % % output_kf.eHaldane = [output_kf_run1.eHaldane; output_kf_run2.eHaldane; output_kf_run3.eHaldane; output_kf_run4.eHaldane; output_kf_run5.eHaldane; output_kf_run6.eHaldane];
% % % % clear output_kf_run1 output_kf_run2 output_kf_run3 output_kf_run4 output_kf_run5 output_kf_run6
% % % % 
% % % % 
% % % % %% 1 vf
% % % % setup.sourceVm = 'experimentalSlopesFixed';
% % % % setup.ode = 'vanHeerden2014';
% % % % setup.ode_pH = 'on';
% % % % ode_pH = setup.ode_pH;
% % % % setup.numpHtested = numpHtested;
% % % % 
% % % % % setup.DFstudy = 4;
% % % % setup.DFstudy = [1 2 3 4]; 
% % % % setup.costfun = 1; 
% % % % 
% % % % setup.plotResults = 0;
% % % % setup.plotEachSim = 0;
% % % % setup.plotEachSimCF = 0;
% % % % 
% % % % % Km variable, Vm fixed
% % % % optfun = @costfun_Vmfixed;
% % % % % xres = output_vf.xres{16};
% % % % xres = output_vf.xres{1};
% % % % 
% % % % % simulation
% % % % setup.weightData = 1;
% % % % setup.weightDataEsp = ones(1,setup.numpHtested);
% % % % setup.weightHaldane = 1;
% % % % setup.selectedLambda = 1;
% % % % 
% % % % setup.plotEachSimCF = 1;
% % % % setup.simAllProfiles = 1;%0;
% % % % 
% % % % [error] = optfun(xres,data,setup);
% % % % 
% % % % 
% % % % %% 1 kf
% % % % setup.sourceVm = 'experimentalSlopesFixed';
% % % % setup.ode = 'vanHeerden2014';
% % % % setup.ode_pH = 'on';
% % % % ode_pH = setup.ode_pH;
% % % % setup.numpHtested = numpHtested;
% % % %  
% % % % setup.DFstudy = 4;
% % % % setup.costfun = 1;
% % % % 
% % % % setup.plotResults = 0;
% % % % setup.plotEachSim = 0;
% % % % setup.plotEachSimCF = 0;
% % % % 
% % % % % Vm variable, Km fixed
% % % % optfun = @costfun_Kmfixed;
% % % % % xres = output_kf.xres{12};
% % % % 
% % % % % simulation
% % % % setup.weightData = 1;
% % % % setup.weightDataEsp = ones(1,setup.numpHtested);
% % % % setup.weightHaldane = 1;
% % % % setup.selectedLambda = 0;
% % % % 
% % % % setup.plotEachSimCF = 1;
% % % % setup.simAllProfiles = 0;
% % % % 
% % % % [error] = optfun(xres,data,setup);


%% 3 plot
% only for vf
load('fullSimRes_1.mat'); fullSimRes_1 = fullSimRes;
load('fullSimRes_16.mat'); fullSimRes_16 = fullSimRes;

figure(41)
%1
for i = 1:4
    plot(fullSimRes_1{i}.t, fullSimRes_1{i}.y(:,8),'-','color','blue')
    hold on
    plot(fullSimRes_16{i}.t, fullSimRes_16{i}.y(:,8),'-','color','black')
    hold on
    plot(data.time{8,i}, data.conc_mean{8,i},'k+')
end
%2
p = get(gca, 'Position');

legend('wH = 0','wH = 1E3','experimental data','location','southwest')
xlabel('time [s]')
ylabel('NADH concentration [mM]')

hold on
rectangle('Position',[250 0 50 0.01],'LineStyle','--')

h = axes('Parent', gcf, 'Position', [p(1)+.49 p(2)+.59 p(3)-.5 p(4)-.6]);
%3
for i = 1:4
    plot(fullSimRes_1{i}.t, fullSimRes_1{i}.y(:,8),'-','color','blue')
    hold on
    plot(fullSimRes_16{i}.t, fullSimRes_16{i}.y(:,8),'-','color','black')
    hold on
    plot(data.time{8,i}, data.conc_mean{8,i},'k+')
end
%4
set(h, 'XTick', [], 'YTick', []);
set(h, 'Xlim', [250 300], 'Ylim', [0 0.005]);
hold off

%% 4 Save array
% 


%% 5 Save plot
savefig(41,'fig04_dataFitLost_vf.fig');


%% memoryDump
% eArray = zeros(36,1);
% for i = 1:36
%     xres = output_vf.xres{i};
%     [error] = optfun(xres,data,setup);
%     eArray(i) = sum(abs(error(1:260)));
% end
% 
% output_vf.eData - eArray
% % %%
% xres = output_vf.xres{1};
% [error1] = optfun(xres,data,setup);
% xres = output_vf.xres{16};
% [error16] = optfun(xres,data,setup);
%     
% e1 = [sum(abs(error1(1:26))),...
%     sum(abs(error1(27:52))),...
%     sum(abs(error1(53:78))),...
%     sum(abs(error1(79:104))),...
%     sum(abs(error1(105:130))),...
%     sum(abs(error1(131:156))),...
%     sum(abs(error1(157:182))),...
%     sum(abs(error1(183:208))),...
%     sum(abs(error1(209:234))),...
%     sum(abs(error1(235:260)))];
%     
% e16 = [sum(abs(error16(1:26))),...
%     sum(abs(error16(27:52))),...
%     sum(abs(error16(53:78))),...
%     sum(abs(error16(79:104))),...
%     sum(abs(error16(105:130))),...
%     sum(abs(error16(131:156))),...
%     sum(abs(error16(157:182))),...
%     sum(abs(error16(183:208))),...
%     sum(abs(error16(209:234))),...
%     sum(abs(error16(235:260)))];
% 
% figure
% plot(pHarray,e1,'.-')
% hold on
% plot(pHarray,e16,'.-')
% legend('e1','e16')
% xlabel('pH')
% ylabel('error')
