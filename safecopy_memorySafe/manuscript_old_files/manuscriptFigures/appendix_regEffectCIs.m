% % APPENDIX_REGEFFECTCIS
% Appendix to show the regularization effect on the confidence intervals,
% both of the parameter estimates and the simulations


%% SECTION 1
% Reproduce the regularization plot for aldolase and write it down in
% subplot(3,3,2).
enzymeNames = {'ald'};
selLamPosArray = 17;
figure(202)
for j = 1:length(enzymeNames)
    tempName = [enzymeNames{j},'_regularizationResults.mat']; % get the name
    load(tempName);% load the regularization results
    p1 = subplot(3,3,2);% subplot
    selLambdaPos = selLamPosArray(j);% selLamPos
    regularizationSimple;% run regularization simple
    %p1.YLim = [0 p1.YLim(2)];
    tempText = [enzymeNames{j}, ', lam=', sprintf('%d',lambdalist(selLambdaPos))];
    tempText2 = erase(tempText,".000000");
%     tempText = [enzymeNames{j}, ', lam=', erase(sprintf('%d',lambdalist(selLambdaPos),".000000"))];
    text(1E0,p1.YLim(2)*1.1,tempText2);
end
% suptitle('Regularization plots for all the enzymes')
set(202,'color','white')

% From the plot, we'll select lambda1 = 5E-5 (low), lambda2 = 1E-1 (mid), 
% lambda3 = 5E-1 (high).
lambda1 = 5E-5;
lambda2 = 1E-1;
lambda3 = 5E-1;


%% SECTION 2
% Estimating parameters for the 3 lambda values and saving the result.

% recalling sections from 'standard' code
for case2 = 1
    %clear, close all
    set_paths_pHstudy;
    dbstop if error
    for step0 = 1
        % select specific case and recall data
        setup.caseStudyALD = 1;
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
        setup.caseStudyENO = 0;
        selectSetup_pH;
        % added
        setup.saveOutput = 0;

        load('expData.mat','expData');
        import_ald = expData.ald;

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

        pHarray = unique(import_ald.treatedData.pH_corrected);
        for i = 1:numpHtested
            pHval = pHarray(i);
            tempID = find(import_ald.treatedData.pH_corrected==pHval);
            pHTemp(:,i) = import_ald.treatedData.pH_corrected(tempID);
            DFTemp(:,i) = import_ald.treatedData.dilution_corrected(tempID);
            for j = 1:4
                abs_meanTemp{j,i} = import_ald.treatedData.absorbance_mean{tempID(j)};
                abs_stdTemp{j,i} = import_ald.treatedData.absorbance_std{tempID(j)};
                conc_meanTemp{j,i} = import_ald.treatedData.concentration_mean{tempID(j)};
                conc_stdTemp{j,i} = import_ald.treatedData.concentration_std{tempID(j)};
                timeTemp{j,i} = import_ald.treatedData.time{tempID(j)};
                RRsTemp{j,i} = import_ald.treatedData.reaction_rate{tempID(j)};
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
        temp1 = import_ald.rawData.absorbance_corrected{4,4};
        temp2 = import_ald.rawData.absorbance_corrected{5,4};
        temp3 = import_ald.rawData.absorbance_corrected{6,4};
        data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
        data.raw.time = import_ald.rawData.time{1};

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % Directly changing the concentration here, sicne the extinction
        % coefficient did not change.
        dps = length(NADH{1,1});
        endPoint = zeros(numpHtested,DFs);
        for i = 1:DFs
            for j = 1:numpHtested
    %             endPoint(j,i) = min(NADH{j,i});
                endPoint(j,i) = min([NADH{j,1}' NADH{j,2}' NADH{j,3}' NADH{j,4}']);
            end
        end
        for i = 1:DFs
            for j = 1:numpHtested
                for k = 1:dps
                    NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,i);
                end
            end
        end
        data.conc_mean = NADH;
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        % % % % %     Addition ALDolase to delete the first part of the profile
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % NADH concentration
        NADH2 = cell(size(NADH));
        for i = 1:DFs
            for j = 1:numpHtested
                if j == 1
                    NADH2{j,i} = NADH{j,i}(25:end);
                else
                    NADH2{j,i} = NADH{j,i}(8:end);
                end
            end
        end
        % bring time to zero start
        for i = 1:DFs
            for j = 1:numpHtested
                for k = 1:dps
                    if j == 1
                        time{j,i}(k) = time{j,i}(k) - 120;
                    else
                        time{j,i}(k) = time{j,i}(k) - 35;
                    end
                end
            end
        end
        % GPD reaction rate
        RRs2 = cell(size(RRs));
        for i = 1:DFs
            for j = 1:numpHtested
                if j == 1
                    RRs2{j,i} = RRs{j,i}(25:end);
                else
                    RRs2{j,i} = RRs{j,i}(8:end);
                end
            end
        end
        data.RRs = RRs2;
        % time
        time2 = cell(size(time));
        for i = 1:DFs
            for j = 1:numpHtested
                if j == 1
                    time2{j,i} = time{j,i}(25:end);
                else
                    time2{j,i} = time{j,i}(8:end);
                end
            end
        end
        clear NADH, NADH = NADH2; clear NADH2
        data.conc_mean = NADH;    
        clear time, time = time2; clear time2
        data.time = time;
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        pHvals = unique(import_ald.treatedData.pH_corrected);
        % visualize: check calculations made
        hs = figure('units','normalized','outerposition',[0 0 1 1]);
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
        end
        suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
        suptitle(suptitleName);
    %     if setup.plotOutput == 1
    %         saveFigName = ['results\', setup.enzymeName,'\', setup.enzymeName, '_concentrations_basezero.fig'];
    %         savefig(hs,saveFigName);
    %     end

    %     figure
    %     plot(pHvals, Vmax(:,4),'.-')
    %     title('Starting estimate: Vmax [mM s-1] vs pH')    
    end
    
    
    
    
end
for case2 = 2
    % %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 6; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
    % select start point (this needs manual selection deppending on previous plots)
    dp_start = 2 * ones(size(data.conc_mean));
    dp_start(:,1) = 10 * ones(size(dp_start(:,1)));
    dp_start(:,2) = 5 * ones(size(dp_start(:,2)));
    dp_start(:,3) = 3 * ones(size(dp_start(:,3)));
    dp_start(2,4) = 3 * ones(size(dp_start(2,4)));
    % blank cell total length
    total_len = zeros(size(dp_start));
    % DFs considered
    DFarray = [1/8 1/4 1/2 1/1];
    % idxs2consider
    idxs2consider = [1 1 1 0;
                    1 1 0 0;
                    1 1 0 0;
                    1 1 1 0;
                    1 1 1 0;
                    1 1 1 1;
                    1 1 1 1;
                    1 1 1 1];

    % Experimental rates determination and plotting
    expRatesDetermination;
    close all
end
for case2 = 3
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
    plength = 13; % Kms (3) + Vms (1) * numpH (8) + vm linking (2) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 3*ones(1,plength);
    lb = -3*ones(1,plength);
    options = optimset('Display','iter');
end


%% Second try
middleLambdas = [5E-3,1E-2,2E-2,5E-2,1E-1,2E-1,3E-1,4E-1];


% simulation display off
setup.plotEachSimCF = 0;
setup.plotEachSim = 0;
setup.simAllProfiles = 0;
% no lambda selected
setup.selectedLambda = 0; % by now just testing
tic
[xres_nonreg,resnorm_nonreg,residual_nonreg,~,~,~,Jacobian_nonreg] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t_nonreg = toc;
    [error_nonreg] = optfun(xres_nonreg,data,setup);
    N_nonreg = length(error_nonreg);
    Jacobian_nonreg = full(Jacobian_nonreg);  
    varp_nonreg = resnorm_nonreg*inv(Jacobian_nonreg'*Jacobian_nonreg)/N_nonreg; % covariance matrix
    stdp_nonreg = sqrt(diag(varp_nonreg));

% % selected lambda value
% setup.selectedLambda = 1E-1; % by now just testing
% tic
% [xres_middlereg,resnorm_middlereg,residual_middlereg,~,~,~,Jacobian_middlereg] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
% t_middlereg = toc;
%     [error_middlereg] = optfun(xres_middlereg,data,setup);
%     N_middlereg = length(error_middlereg);
%     Jacobian_middlereg = full(Jacobian_middlereg);  
%     varp_middlereg = resnorm_middlereg*inv(Jacobian_middlereg'*Jacobian_middlereg)/N_middlereg; % covariance matrix
%     stdp_middlereg = sqrt(diag(varp_middlereg));

for casetemp = 1
    % 1
    % selected lambda value
    setup.selectedLambda = middleLambdas(1); % by now just testing
    tic
    [xres_middlereg1,resnorm_middlereg1,residual_middlereg1,~,~,~,Jacobian_middlereg1] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t_middlereg1 = toc;
        [error_middlereg1] = optfun(xres_middlereg1,data,setup);
        N_middlereg1 = length(error_middlereg1);
        Jacobian_middlereg1 = full(Jacobian_middlereg1);  
        varp_middlereg1 = resnorm_middlereg1*inv(Jacobian_middlereg1'*Jacobian_middlereg1)/N_middlereg1; % covariance matrix
        stdp_middlereg1 = sqrt(diag(varp_middlereg1));
    % 1
    % selected lambda value
    setup.selectedLambda = middleLambdas(2); % by now just testing
    tic
    [xres_middlereg2,resnorm_middlereg2,residual_middlereg2,~,~,~,Jacobian_middlereg2] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t_middlereg2 = toc;
        [error_middlereg2] = optfun(xres_middlereg2,data,setup);
        N_middlereg2 = length(error_middlereg2);
        Jacobian_middlereg2 = full(Jacobian_middlereg2);  
        varp_middlereg2 = resnorm_middlereg2*inv(Jacobian_middlereg2'*Jacobian_middlereg2)/N_middlereg2; % covariance matrix
        stdp_middlereg2 = sqrt(diag(varp_middlereg2));
    % 1
    % selected lambda value
    setup.selectedLambda = middleLambdas(3); % by now just testing
    tic
    [xres_middlereg3,resnorm_middlereg3,residual_middlereg3,~,~,~,Jacobian_middlereg3] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t_middlereg3 = toc;
        [error_middlereg3] = optfun(xres_middlereg3,data,setup);
        N_middlereg3 = length(error_middlereg3);
        Jacobian_middlereg3 = full(Jacobian_middlereg3);  
        varp_middlereg3 = resnorm_middlereg3*inv(Jacobian_middlereg3'*Jacobian_middlereg3)/N_middlereg3; % covariance matrix
        stdp_middlereg3 = sqrt(diag(varp_middlereg3));
    % 1
    % selected lambda value
    setup.selectedLambda = middleLambdas(4); % by now just testing
    tic
    [xres_middlereg4,resnorm_middlereg4,residual_middlereg4,~,~,~,Jacobian_middlereg4] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t_middlereg4 = toc;
        [error_middlereg4] = optfun(xres_middlereg4,data,setup);
        N_middlereg4 = length(error_middlereg4);
        Jacobian_middlereg4 = full(Jacobian_middlereg4);  
        varp_middlereg4 = resnorm_middlereg4*inv(Jacobian_middlereg4'*Jacobian_middlereg4)/N_middlereg4; % covariance matrix
        stdp_middlereg4 = sqrt(diag(varp_middlereg4));
    % 1
    % selected lambda value
    setup.selectedLambda = middleLambdas(5); % by now just testing
    tic
    [xres_middlereg5,resnorm_middlereg5,residual_middlereg5,~,~,~,Jacobian_middlereg5] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t_middlereg5 = toc;
        [error_middlereg5] = optfun(xres_middlereg5,data,setup);
        N_middlereg5 = length(error_middlereg5);
        Jacobian_middlereg5 = full(Jacobian_middlereg5);  
        varp_middlereg5 = resnorm_middlereg5*inv(Jacobian_middlereg5'*Jacobian_middlereg5)/N_middlereg5; % covariance matrix
        stdp_middlereg5 = sqrt(diag(varp_middlereg5));
    % 1
    % selected lambda value
    setup.selectedLambda = middleLambdas(6); % by now just testing
    tic
    [xres_middlereg6,resnorm_middlereg6,residual_middlereg6,~,~,~,Jacobian_middlereg6] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t_middlereg6 = toc;
        [error_middlereg6] = optfun(xres_middlereg6,data,setup);
        N_middlereg6 = length(error_middlereg6);
        Jacobian_middlereg6 = full(Jacobian_middlereg6);  
        varp_middlereg6 = resnorm_middlereg6*inv(Jacobian_middlereg6'*Jacobian_middlereg6)/N_middlereg6; % covariance matrix
        stdp_middlereg6 = sqrt(diag(varp_middlereg6));
    % 1
    % selected lambda value
    setup.selectedLambda = middleLambdas(7); % by now just testing
    tic
    [xres_middlereg7,resnorm_middlereg7,residual_middlereg7,~,~,~,Jacobian_middlereg7] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t_middlereg7 = toc;
        [error_middlereg7] = optfun(xres_middlereg7,data,setup);
        N_middlereg7 = length(error_middlereg7);
        Jacobian_middlereg7 = full(Jacobian_middlereg7);  
        varp_middlereg7 = resnorm_middlereg7*inv(Jacobian_middlereg7'*Jacobian_middlereg7)/N_middlereg7; % covariance matrix
        stdp_middlereg7 = sqrt(diag(varp_middlereg7));
    % 1
    % selected lambda value
    setup.selectedLambda = middleLambdas(8); % by now just testing
    tic
    [xres_middlereg8,resnorm_middlereg8,residual_middlereg8,~,~,~,Jacobian_middlereg8] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t_middlereg8 = toc;
        [error_middlereg8] = optfun(xres_middlereg8,data,setup);
        N_middlereg8 = length(error_middlereg8);
        Jacobian_middlereg8 = full(Jacobian_middlereg8);  
        varp_middlereg8 = resnorm_middlereg8*inv(Jacobian_middlereg8'*Jacobian_middlereg8)/N_middlereg8; % covariance matrix
        stdp_middlereg8 = sqrt(diag(varp_middlereg8));
        
end
% selected lambda value
setup.selectedLambda = 5E-1; % by now just testing
tic
[xres_reg,resnorm_reg,residual_reg,~,~,~,Jacobian_reg] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
t_reg = toc;
    [error_reg] = optfun(xres_reg,data,setup);
    N_reg = length(error_reg);
    Jacobian_reg = full(Jacobian_reg);  
    varp_reg = resnorm_reg*inv(Jacobian_reg'*Jacobian_reg)/N_reg; % covariance matrix
    stdp_reg = sqrt(diag(varp_reg));
    
% % saving
% saveName2 = ['results/',setup.enzymeName,'/',setup.enzymeName, '_regularizationONOFF.mat'];
% save(saveName2,'xres_nonreg','resnorm_nonreg','residual_nonreg','Jacobian_nonreg','t_nonreg',...
%     'xres_reg','resnorm_reg','residual_reg','Jacobian_reg','t_reg',...
%     'error_nonreg','N_nonreg','Jacobian_nonreg','varp_nonreg','stdp_nonreg',...
%     'error_reg','N_reg','Jacobian_reg','varp_reg','stdp_reg');

pNames = {'Kf16bp';...
    'Kglyceral3p';...
    'Kdhap';...
    'Vm_pH6_32';...
    'Vm_pH6_81';...
    'Vm_pH7_06';...
    'Vm_pH7_29';...
    'Vm_pH7_51';...
    'Vm_pH7_68';...
    'Vm_pH7_81';...
    'Vm_pH7_90';...
    'blank';...
    'blank'};
% T = table(pNames,stdp_nonreg,stdp_middlereg,stdp_reg)
T = table(pNames,stdp_nonreg,stdp_middlereg1,stdp_middlereg2,stdp_middlereg3,stdp_middlereg4,stdp_middlereg5,stdp_middlereg6,stdp_middlereg7,stdp_middlereg8,stdp_reg)


%%
% figure,
% bar([stdp_nonreg(4:11), stdp_middlereg1(4:11), stdp_middlereg2(4:11), stdp_middlereg3(4:11), stdp_middlereg4(4:11), stdp_middlereg5(4:11), stdp_middlereg6(4:11), stdp_middlereg7(4:11), stdp_middlereg8(4:11), stdp_reg(4:11)])
% % legend(pNames(4:11))
% legend('stdp_nonreg','stdp_middlereg1','stdp_middlereg2','stdp_middlereg3','stdp_middlereg4','stdp_middlereg5','stdp_middlereg6','stdp_middlereg7','stdp_middlereg8','stdp_reg')

% Select cases 'stdp_nonreg', 'stdp_middlereg1' and 'stdp_reg'.


%% create noised dataset
stdp_nonreg = ones(size(stdp_nonreg)); % artificial change since it was 'Inf'
n = 500;
% x_temp
x_tempOri_nonreg = zeros(n,length(xres_nonreg));
x_tempOri_middlereg1 = zeros(n,length(xres_middlereg1));
x_tempOri_reg = zeros(n,length(xres_reg));
% std
stdpOri_nonreg = zeros(n,length(stdp_nonreg));
stdpOri_middlereg1 = zeros(n,length(stdp_middlereg1));
stdpOri_reg = zeros(n,length(stdp_reg));
for i = 1:n
    % x_temp
    x_tempOri_nonreg(i,:) = xres_nonreg;
    x_tempOri_middlereg1(i,:) = xres_middlereg1;
    x_tempOri_reg(i,:) = xres_reg;
    % std
    stdpOri_nonreg(i,:) = stdp_nonreg;
    stdpOri_middlereg1(i,:) = stdp_middlereg1;
    stdpOri_reg(i,:) = stdp_reg;
end
% create randomness
rng(0)
% r = randn(n,length(stdpOri_nonreg));
[~,tempLen] = size(stdpOri_nonreg);
r = randn(n,length(tempLen));
% add randomness
x_tempOri_nonreg = x_tempOri_nonreg + stdpOri_nonreg .* r;
x_tempOri_middlereg1 = x_tempOri_middlereg1 + stdpOri_middlereg1 .* r;
x_tempOri_reg = x_tempOri_reg + stdpOri_reg .* r;

% add on the latest position, the correct parameter
x_tempOri_nonreg = [x_tempOri_nonreg; xres_nonreg];
x_tempOri_middlereg1 = [x_tempOri_middlereg1; xres_middlereg1];
x_tempOri_reg = [x_tempOri_reg; xres_reg];


% %% simulate
%
[~,FullSim_nonreg] = simSysALD_multipleParameterSets(x_tempOri_nonreg,data,setup);
disp('non_reg done');
%
[~,FullSim_middlereg1] = simSysALD_multipleParameterSets(x_tempOri_middlereg1,data,setup);
disp('middlereg1 done');
%
[~,FullSim_reg] = simSysALD_multipleParameterSets(x_tempOri_reg,data,setup);
disp('reg done');

% %%
% %% plot (simple + heatmap)
enzymeNames = {'ald'};
selLamPosArray = 17;
figure(202)
subplot(4,3,2)
for j = 1:length(enzymeNames)
    tempName = [enzymeNames{j},'_regularizationResults.mat']; % get the name
    load(tempName);% load the regularization results
    p1 = subplot(3,3,2);% subplot
    selLambdaPos = selLamPosArray(j);% selLamPos
    regularizationSimple;% run regularization simple
end
hold on

%
% recall data
load('ald_parEst.mat');
pHvals = data.pH(:,1);
% nonreg
tempVals = data.Vmax(:,4)' .* 10 .^ xres_nonreg(4:11);
nonreg_pVals = ( tempVals') * (60 * 60 / setup.concProtein);
tempVals = data.Vmax(:,4) .* (10 .^ stdp_nonreg(4:11) - 1);
nonreg_stdev = ( tempVals ) * (60 * 60 / setup.concProtein);
% middlereg1
tempVals = data.Vmax(:,4)' .* 10 .^ xres_middlereg1(4:11);
middlereg1_pVals = ( tempVals') * (60 * 60 / setup.concProtein);
tempVals = data.Vmax(:,4) .* (10 .^ stdp_middlereg1(4:11) - 1);
middlereg1_stdev = ( tempVals ) * (60 * 60 / setup.concProtein);
% reg
tempVals = data.Vmax(:,4)' .* 10 .^ xres_reg(4:11);
reg_pVals = ( tempVals') * (60 * 60 / setup.concProtein);
tempVals = data.Vmax(:,4) .* ( 10 .^ stdp_reg(4:11) - 1);
reg_stdev = ( tempVals ) * (60 * 60 / setup.concProtein);

%
% nonreg pEstimates
subplot(4,3,4)
errorbar(pHvals, nonreg_pVals, nonreg_stdev,'k.-','Linewidth',1.2,'MarkerSize',10)
ylabel({'Vmax estimated';'[umol mgP^{-1} min^{-1}]'});
xlabel('pH');
ylim([0 3])

% middlereg1 pEstimates
subplot(4,3,5)
errorbar(pHvals, middlereg1_pVals, middlereg1_stdev,'k.-','Linewidth',1.2,'MarkerSize',10)
% ylabel('Vmax estimated [umol mgP^{-1} min^{-1}]');
xlabel('pH');
ylim([0 3])

% reg pEstimates
subplot(4,3,6)
errorbar(pHvals, reg_pVals, reg_stdev,'k.-','Linewidth',1.2,'MarkerSize',10)
% ylabel('Vmax estimated [umol mgP^{-1} min^{-1}]');
xlabel('pH');
ylim([0 3])

% nonreg line plot
% plotting the case of pH7.51 (#5)
subplot(4,3,7)
for i = 1:4
    % experimental data
    plot(data.time{5,i}, data.conc_mean{5,i},'k.','MarkerSize',4)
    hold on
    % noised simulations
    len = length(FullSim_nonreg) - 1;
    for j = 1:len
        plot(FullSim_nonreg{j}{i,5}.t, FullSim_nonreg{j}{i,5}.y(:,5),'-','LineWidth',1,'color',[0.85 0.85 1])
        hold on
    end
    % real estimated parameters
%     e
    plot(FullSim_nonreg{end}{i,5}.t, FullSim_nonreg{end}{i,5}.y(:,5),'-','LineWidth',2,'color',[0.5 0.5 1])
    hold on
end
ylabel({'NADH concentration';'[mM]'}) 
xlabel('assay time [s]')
xlim([0 400])

% middlereg1 line plot
% plotting the case of pH7.51 (#5)
subplot(4,3,8)
for i = 1:4
    % experimental data
    plot(data.time{5,i}, data.conc_mean{5,i},'k.','MarkerSize',4)
    hold on
    % noised simulations
    len = length(FullSim_middlereg1) - 1;
    for j = 1:len
        plot(FullSim_middlereg1{j}{i,5}.t, FullSim_middlereg1{j}{i,5}.y(:,5),'-','LineWidth',1,'color',[0.85 0.85 1])
        hold on
    end
    % real estimated parameters
%     e
    plot(FullSim_middlereg1{end}{i,5}.t, FullSim_middlereg1{end}{i,5}.y(:,5),'-','LineWidth',2,'color',[0.5 0.5 1])
    hold on
end
% ylabel('NADH concentration [mM]') 
xlabel('assay time [s]')
xlim([0 400])

% reg line plot
% plotting the case of pH7.51 (#5)
subplot(4,3,9)
for i = 1:4
    % experimental data
    plot(data.time{5,i}, data.conc_mean{5,i},'k.','MarkerSize',4)
    hold on
    % noised simulations
    len = length(FullSim_reg) - 1;
    for j = 1:len
        plot(FullSim_reg{j}{i,5}.t, FullSim_reg{j}{i,5}.y(:,5),'-','LineWidth',1,'color',[0.85 0.85 1])
        hold on
    end
    % real estimated parameters
%     e
    plot(FullSim_reg{end}{i,5}.t, FullSim_reg{end}{i,5}.y(:,5),'-','LineWidth',2,'color',[0.5 0.5 1])
    hold on
end
% ylabel('NADH concentration [mM]') 
xlabel('assay time [s]')

%white color
set(202,'color','white')
xlim([0 400])


% %% plot (heatmap)
sp10 = subplot(4,3,10);
sp11 = subplot(4,3,11);
sp12 = subplot(4,3,12);

% nonreg
Xarray = [];
Yarray = [];
for o = 1:n
    Xarray = [Xarray; FullSim_nonreg{o}{4,5}.t];
    Yarray = [Yarray; FullSim_nonreg{o}{4,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_nonreg{o}{3,5}.t];
    Yarray = [Yarray; FullSim_nonreg{o}{3,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_nonreg{o}{2,5}.t];
    Yarray = [Yarray; FullSim_nonreg{o}{2,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_nonreg{o}{1,5}.t];
    Yarray = [Yarray; FullSim_nonreg{o}{1,5}.y(:,5)];
end
[outfile,hplot10] = heatscatter_modified(Xarray, Yarray, [], 'test.png',[],[],'.',0,0,[],[],[]);
copyobj(hplot10,sp10)
close(203)

% middlereg1
Xarray = [];
Yarray = [];
for o = 1:n
    Xarray = [Xarray; FullSim_middlereg1{o}{4,5}.t];
    Yarray = [Yarray; FullSim_middlereg1{o}{4,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_middlereg1{o}{3,5}.t];
    Yarray = [Yarray; FullSim_middlereg1{o}{3,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_middlereg1{o}{2,5}.t];
    Yarray = [Yarray; FullSim_middlereg1{o}{2,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_middlereg1{o}{1,5}.t];
    Yarray = [Yarray; FullSim_middlereg1{o}{1,5}.y(:,5)];
end
[outfile,hplot11] = heatscatter_modified(Xarray, Yarray, [], 'test.png',[],[],'.',0,0,[],[],[]);
copyobj(hplot11,sp11)
close(203)
% reg
Xarray = [];
Yarray = [];
for o = 1:n
    Xarray = [Xarray; FullSim_reg{o}{4,5}.t];
    Yarray = [Yarray; FullSim_reg{o}{4,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_reg{o}{3,5}.t];
    Yarray = [Yarray; FullSim_reg{o}{3,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_reg{o}{2,5}.t];
    Yarray = [Yarray; FullSim_reg{o}{2,5}.y(:,5)];
end
for o = 1:n
    Xarray = [Xarray; FullSim_reg{o}{1,5}.t];
    Yarray = [Yarray; FullSim_reg{o}{1,5}.y(:,5)];
end
[outfile,hplot12] = heatscatter_modified(Xarray, Yarray, [], 'test.png',[],[],'.',0,0,[],[],[]);
copyobj(hplot12,sp12)
close(203)

 
% %%
% 
%  figure;
%  x = linspace(-2*pi,2*pi,500);
%  y = sin(x);
%  hCurve = plot(x,y,'r');
%  
%  figure;
%  hSub1 = subplot(2,1,1);
%  hSub2 = subplot(2,1,2);
%  
%  copyobj(hCurve,hSub1);
 

%% memoryDump
% % % % savefig(201,'results/manuscriptFigures/appendixes_allRegPlots.fig');

% %%
% lambdalist = [5E-5, 1E-1, 5E-1];
% array_xres = cell(1,length(lambdalist));
% array_eData = cell(1,length(lambdalist));
% array_eParams = cell(1,length(lambdalist));
% array_fullError = cell(1,length(lambdalist));
% cell_resnorm = cell(1,length(lambdalist));
% cell_residual = cell(1,length(lambdalist));
% cell_Jacobian = cell(1,length(lambdalist));
% for i = 1:length(lambdalist)
%     fprintf('pEst for lambda=%d\n',lambdalist(i));
%     setup.selectedLambda = lambdalist(i);
%     tic
%     [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
%     t = toc;
%     setup.selectedLambda = 1;
%     [error] = optfun(xres,data,setup);
%     % seting in output arrays
%     array_xres{i} = xres;
%     array_eData{i} = error(1:end-16);
%     array_eParams{i} = error(end-15:end);
%     array_fullError{i} = error;
%     cell_resnorm{i} = resnorm;
%     cell_residual{i} = residual;
%     cell_Jacobian{i} = Jacobian;
% end
% 
% 
% %%
% N_cell = cell(1,length(lambdalist));
% Jacobian_cell = cell(1,length(lambdalist));
% varp_cell = cell(1,length(lambdalist));
% stdp_cell = cell(1,length(lambdalist));
% for i = 1:length(lambdalist)
%     % recall
%     error = array_fullError{i};
%     Jacobian = cell_Jacobian{i};
%     resnorm = cell_resnorm{i};
%     
%     % confidence intervals
%     N = length(error);
%     Jacobian = full(Jacobian);  
%     varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
%     stdp = sqrt(diag(varp));
%     
%     % reassign
%     N_cell{i} = N;
%     Jacobian_cell{i} = Jacobian;
%     varp_cell{i} = varp;
%     stdp_cell{i} = stdp;
% end
% 
% %% 
% cell_vm = cell(1,length(lambdalist));
% cell_vm_up = cell(1,length(lambdalist));
% cell_vm_down = cell(1,length(lambdalist));
% for j = 1:length(lambdalist)
%     % recall
%     xres_selected = array_xres{j};
%     vm = zeros(numpHtested,1);
%     vm_up = zeros(numpHtested,1);
%     vm_down = zeros(numpHtested,1);
%     % calculate
%     for i = 1:numpHtested
%         vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+3);
%         vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+3)+stdp(i+3)); %up
%         vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+3)-stdp(i+3)); %down
%     end
%     % reassign
%     cell_vm{j} = vm;
%     cell_vm_up{j} = vm_up;
%     cell_vm_down{j} = vm_down;
% end

