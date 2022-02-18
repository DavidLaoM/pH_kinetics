% Kms are assumed constant
% Keq from the eQuilibrator

% % Contents
% AP 1  Estimation at pH #6
% AP 2  Estimation for all the other pH values
% AP 3  Regularization for the selected pH value

%% (0) Setup and data load
clear
set_paths_pHstudy;
dbstop if error
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
    % added
    setup.saveOutput = 0;

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


%% (1.1) Simple parameter fit. Parameter estimation. pH#6
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';

setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
setup.costfun = 1;

setup.weightData = 1;
setup.weightDataEsp = ones(1,12);
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% % Km fixed
% optfun = @costfun_Kmfixed;
% plength = 14; % Kms (2) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
% x_temp = zeros(1,plength);
% ub = 3*ones(1,plength);
% lb = -3*ones(1,plength);
% options = optimset('Display','iter');

% Vm fixed
optfun = @costfun_Vmfixed;
setup.startVm = 12;
plength = 25; % Kms (2) * numpH (12) + Vms (1) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');
% %%
% % testing the costfunction
% [error] = optfun(x_temp,data,setup);



%% (2.1) Simple parameter fit. Parameter estimation. All pH#
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
xres_cell = cell(numpHtested,1);
% loop to estimate each pH value
for i = 1:numpHtested
    fprintf('testing pH#%f\n',i);
    % setup
    setup.startVm = i;
    % parameter estimation
    tic
    [xres_temp,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;
    xres_cell{i} = xres_temp;
end
%% (2.2) plotting 1 by 1 + overall error
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
error_array = zeros(numpHtested,1);
for i = 1:numpHtested
    xres_temp = xres_cell{i};
    [error] = optfun(xres_temp,data,setup);
    error_array(i) = sum(abs(error));
end
%%
% figure(1001)
% semilogy(setup.fullpHarray, error_array, '.-')
% hold on
% line([6 8],[1E-1 1E-1],'color','blue','LineStyle','--')
% xlabel('pH value')
% ylabel('overall error simulation vs experimental data')
% title({'Overall error obtained at each pH value';'\color{blue} error with km_{fixed} and vm_{variable} below 0.1.'})
% ylim([1E-2 1E4])
% %%
% set(1001,'color','w');
% savefig(1001,'ENO_vmfixed_checkAllpHs.fig')
% save('ENO_vmfixed_testAllpH.mat','xres_cell','error_array')


%% (2) check with ENO km fixed

%% (2.1) Simple parameter fit. Parameter estimation
setup.ode = 'vanHeerden2014';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';

setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
setup.costfun = 1;

setup.weightData = 1;
setup.weightDataEsp = ones(1,12);
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% Km fixed
optfun = @costfun_Kmfixed;
plength = 14; % Kms (2) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');


%% (2.2) parameter estimation
tic
[xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t = toc;

% %% (1.2) Simple parameter fit. Results Visualization
setup.plotEachSimCF = 1;
setup.plotEachSim = 0;
setup.simAllProfiles = 1;
[error_kmfixed] = optfun(xres,data,setup);
sum_error_kmfixed = sum(abs(error_kmfixed));
% sum(abs(error))

%% (2.3) plotting
figure(1001)
semilogy(setup.fullpHarray, error_array, '.-')
hold on
line([6 8],[sum_error_kmfixed sum_error_kmfixed],'color','blue','LineStyle','--')
xlabel('pH value')
ylabel('overall error simulation vs experimental data')
title({'Overall error obtained at each pH value';'\color{blue} error with km_{fixed} and vm_{variable} below 0.1.'})
ylim([1E0 1E4])
%%
set(1001,'color','w');
savefig(1001,'ENO_vmfixed_checkAllpHs.fig')
save('ENO_vmfixed_testAllpH.mat','xres_cell','error_array','xres','error_kmfixed','sum_error_kmfixed')

