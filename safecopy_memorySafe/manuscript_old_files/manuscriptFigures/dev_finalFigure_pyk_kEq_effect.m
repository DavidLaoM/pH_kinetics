% % finalFigure_pyk_kEq_effect.m
% 0. Setup
% 1. Estimation normal keq
% 2. Estimation constant Keq
% 3. plot comparison parameter values


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
    setup.caseStudyPYK = 1;
    setup.caseStudyTPI = 0;
    setup.caseStudyENO = 0;
    selectSetup_pH;
    % added
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_pyk = expData.pyk;
    
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
    
    pHarray = unique(import_pyk.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_pyk.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_pyk.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_pyk.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_pyk.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_pyk.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_pyk.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_pyk.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_pyk.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_pyk.treatedData.reaction_rate{tempID(j)};
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
    temp1 = import_pyk.rawData.absorbance_corrected{4,4};
    temp2 = import_pyk.rawData.absorbance_corrected{5,4};
    temp3 = import_pyk.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_pyk.rawData.time{1};
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Directly changing the concentration here, since the extinction
    % coefficient did not change.
    endPoint = zeros(numpHtested,DFs);
    for i = 1:DFs
        for j = 1:numpHtested
%             endPoint(j,i) = min(NADH{j,i});
            endPoint(j,i) = min([NADH{j,1}' NADH{j,2}' NADH{j,3}' NADH{j,4}']);
        end
    end
    for i = 1:DFs
        for j = 1:numpHtested
            dps = length(NADH{j,i});
            for k = 1:dps
                NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,i);
            end
        end
    end
    % making the last values increasing, zero
    for i = 1:DFs
        for j = 1:numpHtested
            [~,I] = min(NADH{j,i});
            % get the idx of the minimum, and tehen from the next idx to
            % end, make all zero
            for k = I:length(NADH{j,i})
                NADH{j,i}(k) = NADH{j,i}(I);
            end
        end
    end
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % % % % %     Addition PGI to delete the first part of the profile
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % NADH concentration 
        % pHs 1-3 | time 250 | dps 51
        % pHs 4 | time 200 | dps 41
        % pHs 5 | time 175 | dps 36
        % pHs 6 | time 125 | dps 26
        % pHs 7 | time 100 | dps 21
        % pHs 8-12 | time 75 | dps 16
    NADH2 = cell(size(NADH));
    for i = 1:DFs
        for j = 1:numpHtested
            NADH2{j,i} = NADH{j,i}(2:end);
        end
    end
    % bring time to zero start
    for i = 1:DFs
        for j = 1:numpHtested
            dps = length(time{j,i});
            for k = 1:dps
                time{j,i}(k) = time{j,i}(k) - 5;
            end
        end
    end
    % GPD reaction rate
    RRs2 = cell(size(RRs));
    for i = 1:DFs
        for j = 1:numpHtested
            RRs2{j,i} = RRs{j,i}(2:end);
        end
    end
    data.RRs = RRs2;
    % time
    time2 = cell(size(time));
    for i = 1:DFs
        for j = 1:numpHtested
            time2{j,i} = time{j,i}(2:end);
        end
    end
    clear NADH, NADH = NADH2; clear NADH2
    data.conc_mean = NADH;    
    clear time, time = time2; clear time2
    data.time = time;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    pHvals = unique(import_pyk.treatedData.pH_corrected);
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
    end
    suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
    suptitle(suptitleName);
    
    figure
    plot(pHvals, Vmax(:,4),'.-')
    title('Starting estimate: Vmax [mM s-1] vs pH')    
end


%% (1) Simple parameter fit. Parameter estimation
for caseChangeKeq = 1
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
    setup.weightDataEsp = ones(1,numpHtested);
    setup.weightHaldane = 0; % off for this case (Keq fixed)
    setup.selectedLambda = 1E-5; % by now just testing

    % Km fixed
    optfun = @costfun_Kmfixed;
    plength = 18; % others (6) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 1*ones(1,plength);
    lb = -1*ones(1,plength);
    options = optimset('Display','iter');

    % %%
    % % testing the costfunction
    % setup.plotEachSimCF = 1;
    % [error] = optfun(x_temp,data,setup);
    % setup.plotEachSimCF = 0;
    % %%

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;

    % calc error
    [error] = optfun(xres,data,setup);
    for case2 = 1 % confidence intervals
        N = length(error);
        Jacobian = full(Jacobian);  
        varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
        stdp = sqrt(diag(varp));
    end

    % calc and plot parameter values
    xres_selected = xres; %lambdalist based in 'ones', lam=0.1, loc=5.
    % xres_selected = array_xres{14}; %lambdalist based in 'ones', lam=0.1, loc=5.

    kadp = ones(numpHtested,1);
    katp = ones(numpHtested,1);
    kfbp = ones(numpHtested,1);
    kpep = ones(numpHtested,1);
    L = ones(numpHtested,1);
    nHill = ones(numpHtested,1);
    vm = zeros(numpHtested,1);
    for i = 1:numpHtested
        kadp(i) = 0.243 * 10.^xres_selected(1); %mM
        katp(i) = 9.3 * 10.^xres_selected(2); %mM
        kfbp(i) = 0.2 * 10.^xres_selected(3); %mM
        kpep(i) = 0.281 * 10.^xres_selected(4); %mM
        L(i) = 60000 * 10.^xres_selected(5);
        nHill(i) = 4 * 10.^xres_selected(6);
        vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+6);
    end
    keq_pyk = setup.Keq_PYK;
    keq_ldh = setup.Keq_LDH;

    % limits
    kadp_up = ones(numpHtested,1);
    katp_up = ones(numpHtested,1);
    kfbp_up = ones(numpHtested,1);
    kpep_up = ones(numpHtested,1);
    L_up = ones(numpHtested,1);
    nHill_up = ones(numpHtested,1);
    vm_up = zeros(numpHtested,1);

    kadp_down = ones(numpHtested,1);
    katp_down = ones(numpHtested,1);
    kfbp_down = ones(numpHtested,1);
    kpep_down = ones(numpHtested,1);
    L_down = ones(numpHtested,1);
    nHill_down = ones(numpHtested,1);
    vm_down = zeros(numpHtested,1);
    for i = 1:numpHtested
        % up
        kadp_up(i) = 0.243 * 10.^(xres_selected(1)+stdp(1)); %mM
        katp_up(i) = 9.3 * 10.^(xres_selected(2)+stdp(2)); %mM
        kfbp_up(i) = 0.2 * 10.^(xres_selected(3)+stdp(3)); %mM
        kpep_up(i) = 0.281 * 10.^(xres_selected(4)+stdp(4)); %mM
        L_up(i) = 60000 * 10.^(xres_selected(5)+stdp(5)); %mM
        nHill_up(i) = 4 * 10.^(xres_selected(6)+stdp(6)); %mM
        vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+6)+stdp(i+6)); %mM
        % down
        kadp_down(i) = 0.243 * 10.^(xres_selected(1)-stdp(1)); %mM
        katp_down(i) = 9.3 * 10.^(xres_selected(2)-stdp(2)); %mM
        kfbp_down(i) = 0.2 * 10.^(xres_selected(3)-stdp(3)); %mM
        kpep_down(i) = 0.281 * 10.^(xres_selected(4)-stdp(4)); %mM
        L_down(i) = 60000 * 10.^(xres_selected(5)-stdp(5)); %mM
        nHill_down(i) = 4 * 10.^(xres_selected(6)-stdp(6)); %mM
        vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+6)-stdp(i+6)); %mM
    end


    % %% (3.3) Study on paarameter values: plotting
    vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
    vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
    vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;

    figure(105)
    % figure(106)

    subplot(4,3,1) % vm
    %     plot(pHarray,vm_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    %     plot(pHarray,vm_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    %     plot(pHarray,vm,'.-','color','black')
    %     title('v_{m} [mM]')
        plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_uChange,'.-','color','black')
        title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})
    subplot(4,3,2) % L
        plot(pHarray,L_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,L_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,L,'.-','color','black')
        title('k_{L} [mM]')
    subplot(4,3,3) % nHill
        plot(pHarray,nHill_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,nHill_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,nHill,'.-','color','black')
        title('k_{nHill} [mM]')

    subplot(4,3,5) % kadp
        plot(pHarray,kadp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kadp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kadp,'.-','color','black')
        title('k_{adp} [mM]')
    subplot(4,3,6) % katp
        plot(pHarray,katp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,katp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,katp,'.-','color','black')
        title('k_{atp} [mM]')
    subplot(4,3,7) % kfbp
        plot(pHarray,kfbp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kfbp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kfbp,'.-','color','black')
        title('k_{fbp} [mM]')
    subplot(4,3,8) % kpep
        plot(pHarray,kpep_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kpep_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kpep,'.-','color','black')
        title('k_{pep} [mM]')

    subplot(4,3,10) % keq_pyk
        plot(pHarray,keq_pyk,'.-','color','black')
        title('k_{eq.PYK} [mM]')
    subplot(4,3,11) % keq_ldh
        plot(pHarray,keq_ldh,'.-','color','black')
        title('k_{eq.LDH} [mM]')


    subplot(4,3,4)
    plot(pHarray,xres_selected(7:18)','.-','color','black')
    title('vm reference')

    suptitle('parameter estimates vs pH')
    
    
        output_pyk_changeKeq.xres_selected = xres_selected;

        output_pyk_changeKeq.pHarray = pHarray;

        output_pyk_changeKeq.kadp = kadp;% = ones(numpHtested,1);
        output_pyk_changeKeq.katp = katp;% = ones(numpHtested,1);
        output_pyk_changeKeq.kfbp = kfbp;% = ones(numpHtested,1);
        output_pyk_changeKeq.kpep = kpep;% = ones(numpHtested,1);
        output_pyk_changeKeq.L = L;% = ones(numpHtested,1);
        output_pyk_changeKeq.nHill = nHill;% = ones(numpHtested,1);
        output_pyk_changeKeq.vm = vm;% = zeros(numpHtested,1);
        output_pyk_changeKeq.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

        output_pyk_changeKeq.kadp_up = kadp_up;% = ones(numpHtested,1);
        output_pyk_changeKeq.katp_up = katp_up;% = ones(numpHtested,1);
        output_pyk_changeKeq.kfbp_up = kfbp_up;% = ones(numpHtested,1);
        output_pyk_changeKeq.kpep_up = kpep_up;% = ones(numpHtested,1);
        output_pyk_changeKeq.L_up = L_up;% = ones(numpHtested,1);
        output_pyk_changeKeq.nHill_up = nHill_up;% = ones(numpHtested,1);
        output_pyk_changeKeq.vm_up = vm_up;% = zeros(numpHtested,1);
        output_pyk_changeKeq.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

        output_pyk_changeKeq.kadp_down = kadp_down;% = ones(numpHtested,1);
        output_pyk_changeKeq.katp_down = katp_down;% = ones(numpHtested,1);
        output_pyk_changeKeq.kfbp_down = kfbp_down;% = ones(numpHtested,1);
        output_pyk_changeKeq.kpep_down = kpep_down;% = ones(numpHtested,1);
        output_pyk_changeKeq.L_down = L_down;% = ones(numpHtested,1);
        output_pyk_changeKeq.nHill_down = nHill_down;% = ones(numpHtested,1);
        output_pyk_changeKeq.vm_down = vm_down;% = zeros(numpHtested,1);
        output_pyk_changeKeq.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

        output_pyk_changeKeq.keq_pyk = keq_pyk;% = setup.Keq_PYK;
        output_pyk_changeKeq.keq_ldh = keq_ldh;% = setup.Keq_LDH;
    
    
end



%% (1) Simple parameter fit. Parameter estimation
temp_Keq_PYK = setup.Keq_PYK;
temp_Keq_LDH = setup.Keq_LDH;
for i = 1:12
    setup.Keq_PYK(i) = setup.Keq_PYK(6);
    setup.Keq_LDH(i) = setup.Keq_PYK(6);
end

for caseChangeKeq = 1
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
    setup.weightDataEsp = ones(1,numpHtested);
    setup.weightHaldane = 0; % off for this case (Keq fixed)
    setup.selectedLambda = 1E-5; % by now just testing

    % Km fixed
    optfun = @costfun_Kmfixed;
    plength = 18; % others (6) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
    x_temp = zeros(1,plength);
    ub = 1*ones(1,plength);
    lb = -1*ones(1,plength);
    options = optimset('Display','iter');

    % %%
    % % testing the costfunction
    % setup.plotEachSimCF = 1;
    % [error] = optfun(x_temp,data,setup);
    % setup.plotEachSimCF = 0;
    % %%

    % parameter estimation
    tic
    [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
    t = toc;

    % calc error
    [error] = optfun(xres,data,setup);
    for case2 = 1 % confidence intervals
        N = length(error);
        Jacobian = full(Jacobian);  
        varp = resnorm*inv(Jacobian'*Jacobian)/N; % covariance matrix
        stdp = sqrt(diag(varp));
    end

    % calc and plot parameter values
    xres_selected = xres; %lambdalist based in 'ones', lam=0.1, loc=5.
    % xres_selected = array_xres{14}; %lambdalist based in 'ones', lam=0.1, loc=5.

    kadp = ones(numpHtested,1);
    katp = ones(numpHtested,1);
    kfbp = ones(numpHtested,1);
    kpep = ones(numpHtested,1);
    L = ones(numpHtested,1);
    nHill = ones(numpHtested,1);
    vm = zeros(numpHtested,1);
    for i = 1:numpHtested
        kadp(i) = 0.243 * 10.^xres_selected(1); %mM
        katp(i) = 9.3 * 10.^xres_selected(2); %mM
        kfbp(i) = 0.2 * 10.^xres_selected(3); %mM
        kpep(i) = 0.281 * 10.^xres_selected(4); %mM
        L(i) = 60000 * 10.^xres_selected(5);
        nHill(i) = 4 * 10.^xres_selected(6);
        vm(i) = data.Vmax(i,4) * 10.^xres_selected(i+6);
    end
    keq_pyk = setup.Keq_PYK;
    keq_ldh = setup.Keq_LDH;

    % limits
    kadp_up = ones(numpHtested,1);
    katp_up = ones(numpHtested,1);
    kfbp_up = ones(numpHtested,1);
    kpep_up = ones(numpHtested,1);
    L_up = ones(numpHtested,1);
    nHill_up = ones(numpHtested,1);
    vm_up = zeros(numpHtested,1);

    kadp_down = ones(numpHtested,1);
    katp_down = ones(numpHtested,1);
    kfbp_down = ones(numpHtested,1);
    kpep_down = ones(numpHtested,1);
    L_down = ones(numpHtested,1);
    nHill_down = ones(numpHtested,1);
    vm_down = zeros(numpHtested,1);
    for i = 1:numpHtested
        % up
        kadp_up(i) = 0.243 * 10.^(xres_selected(1)+stdp(1)); %mM
        katp_up(i) = 9.3 * 10.^(xres_selected(2)+stdp(2)); %mM
        kfbp_up(i) = 0.2 * 10.^(xres_selected(3)+stdp(3)); %mM
        kpep_up(i) = 0.281 * 10.^(xres_selected(4)+stdp(4)); %mM
        L_up(i) = 60000 * 10.^(xres_selected(5)+stdp(5)); %mM
        nHill_up(i) = 4 * 10.^(xres_selected(6)+stdp(6)); %mM
        vm_up(i) = data.Vmax(i,4) * 10.^(xres_selected(i+6)+stdp(i+6)); %mM
        % down
        kadp_down(i) = 0.243 * 10.^(xres_selected(1)-stdp(1)); %mM
        katp_down(i) = 9.3 * 10.^(xres_selected(2)-stdp(2)); %mM
        kfbp_down(i) = 0.2 * 10.^(xres_selected(3)-stdp(3)); %mM
        kpep_down(i) = 0.281 * 10.^(xres_selected(4)-stdp(4)); %mM
        L_down(i) = 60000 * 10.^(xres_selected(5)-stdp(5)); %mM
        nHill_down(i) = 4 * 10.^(xres_selected(6)-stdp(6)); %mM
        vm_down(i) = data.Vmax(i,4) * 10.^(xres_selected(i+6)-stdp(i+6)); %mM
    end


    % %% (3.3) Study on paarameter values: plotting
    vm_up_uChange = vm_up .* 60 .* 60 ./ setup.concProtein;
    vm_down_uChange = vm_down .* 60 .* 60 ./ setup.concProtein;
    vm_uChange = vm .* 60 .* 60 ./ setup.concProtein;

    figure(105)
    % figure(106)

    subplot(4,3,1) % vm
    %     plot(pHarray,vm_up,'.-','color',[0.5 0.5 0.5]), hold on, 
    %     plot(pHarray,vm_down,'.-','color',[0.5 0.5 0.5]), hold on, 
    %     plot(pHarray,vm,'.-','color','black')
    %     title('v_{m} [mM]')
        plot(pHarray,vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,vm_uChange,'.-','color','black')
        title({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})
    subplot(4,3,2) % L
        plot(pHarray,L_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,L_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,L,'.-','color','black')
        title('k_{L} [mM]')
    subplot(4,3,3) % nHill
        plot(pHarray,nHill_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,nHill_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,nHill,'.-','color','black')
        title('k_{nHill} [mM]')

    subplot(4,3,5) % kadp
        plot(pHarray,kadp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kadp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kadp,'.-','color','black')
        title('k_{adp} [mM]')
    subplot(4,3,6) % katp
        plot(pHarray,katp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,katp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,katp,'.-','color','black')
        title('k_{atp} [mM]')
    subplot(4,3,7) % kfbp
        plot(pHarray,kfbp_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kfbp_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kfbp,'.-','color','black')
        title('k_{fbp} [mM]')
    subplot(4,3,8) % kpep
        plot(pHarray,kpep_up,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kpep_down,'.-','color',[0.5 0.5 0.5]), hold on, 
        plot(pHarray,kpep,'.-','color','black')
        title('k_{pep} [mM]')

    subplot(4,3,10) % keq_pyk
        plot(pHarray,keq_pyk,'.-','color','black')
        title('k_{eq.PYK} [mM]')
    subplot(4,3,11) % keq_ldh
        plot(pHarray,keq_ldh,'.-','color','black')
        title('k_{eq.LDH} [mM]')


    subplot(4,3,4)
    plot(pHarray,xres_selected(7:18)','.-','color','black')
    title('vm reference')

    suptitle('parameter estimates vs pH')
    
    
        output_pyk_constantKeq.xres_selected = xres_selected;

        output_pyk_constantKeq.pHarray = pHarray;

        output_pyk_constantKeq.kadp = kadp;% = ones(numpHtested,1);
        output_pyk_constantKeq.katp = katp;% = ones(numpHtested,1);
        output_pyk_constantKeq.kfbp = kfbp;% = ones(numpHtested,1);
        output_pyk_constantKeq.kpep = kpep;% = ones(numpHtested,1);
        output_pyk_constantKeq.L = L;% = ones(numpHtested,1);
        output_pyk_constantKeq.nHill = nHill;% = ones(numpHtested,1);
        output_pyk_constantKeq.vm = vm;% = zeros(numpHtested,1);
        output_pyk_constantKeq.vm_uChange = vm_uChange;% = zeros(numpHtested,1);

        output_pyk_constantKeq.kadp_up = kadp_up;% = ones(numpHtested,1);
        output_pyk_constantKeq.katp_up = katp_up;% = ones(numpHtested,1);
        output_pyk_constantKeq.kfbp_up = kfbp_up;% = ones(numpHtested,1);
        output_pyk_constantKeq.kpep_up = kpep_up;% = ones(numpHtested,1);
        output_pyk_constantKeq.L_up = L_up;% = ones(numpHtested,1);
        output_pyk_constantKeq.nHill_up = nHill_up;% = ones(numpHtested,1);
        output_pyk_constantKeq.vm_up = vm_up;% = zeros(numpHtested,1);
        output_pyk_constantKeq.vm_up_uChange = vm_up_uChange;% = zeros(numpHtested,1);

        output_pyk_constantKeq.kadp_down = kadp_down;% = ones(numpHtested,1);
        output_pyk_constantKeq.katp_down = katp_down;% = ones(numpHtested,1);
        output_pyk_constantKeq.kfbp_down = kfbp_down;% = ones(numpHtested,1);
        output_pyk_constantKeq.kpep_down = kpep_down;% = ones(numpHtested,1);
        output_pyk_constantKeq.L_down = L_down;% = ones(numpHtested,1);
        output_pyk_constantKeq.nHill_down = nHill_down;% = ones(numpHtested,1);
        output_pyk_constantKeq.vm_down = vm_down;% = zeros(numpHtested,1);
        output_pyk_constantKeq.vm_down_uChange = vm_down_uChange;% = zeros(numpHtested,1);

        output_pyk_constantKeq.keq_pyk = keq_pyk;% = setup.Keq_PYK;
        output_pyk_constantKeq.keq_ldh = keq_ldh;% = setup.Keq_LDH;
    
    
end

setup.Keq_PYK = temp_Keq_PYK;
setup.Keq_LDH = temp_Keq_LDH;

% % % % %% plotting both
% % % % % changes are actually really small.
% % % % figure, 
% % % % plot(output_pyk_changeKeq.pHarray, output_pyk_changeKeq.vm_uChange)
% % % % hold on
% % % % plot(output_pyk_constantKeq.pHarray, output_pyk_constantKeq.vm_uChange)
% % % % legend('Keq_{changing}','Keq_{constant}','location','southeast')
% % % % xlabel('pH')
% % % % ylabel({'v_{m} [umol_{NADH} mg_{P}^{-1} min^{-1}]';'not normalized'})

%% section for part 2
% PSA to see the effect of the Keq.

keqvalsTested = [1E-5 1E-4 1E-3 1E-2 1E-1 1E0 1E1 1E2 1E3 1E4 1E5 1E6 1E7 1E8 1E9];
tempResult = cell(15,1);

for k = 1:length(keqvalsTested)
    temp_Keq_PYK = setup.Keq_PYK;
    temp_Keq_LDH = setup.Keq_LDH;
    for q = 1:12
        setup.Keq_PYK(q) = keqvalsTested(k);
        setup.Keq_LDH(q) = keqvalsTested(k);
    end
    setup.plotEachSimCF = 0;
    setup.plotEachSim = 0;
    setup.simAllProfiles = 0;
    
    for tempFigure = 1
        %UNTITLED4 Summary of this function goes here
        %   Detailed explanation goes here
        %     x_temp(1) = Kgap
        %     x_temp(2) = Kbpg
        %     x_temp(3) = Knad
        %     x_temp(4) = Knadh
        %     x_temp([5:6]) = {Vmf, Vmr}, pH#1
        %     x_temp([7:8]) = {Vmf, Vmr}, pH#2
        %     x_temp([9:10]) = {Vmf, Vmr}, pH#3
        %     x_temp([11:12]) = {Vmf, Vmr}, pH#4
        %     x_temp([13:14]) = {Vmf, Vmr}, pH#5
        %     x_temp([15:16]) = {Vmf, Vmr}, pH#6
        %     x_temp([17:18]) = {Vmf, Vmr}, pH#7
        %     x_temp([19:20]) = {Vmf, Vmr}, pH#8
        %     x_temp([21:22]) = {Vmf, Vmr}, pH#9
        %     x_temp([23:24]) = {Vmf, Vmr}, pH#10
        enzyme = setup.enzymeName;
        DFs = setup.DFactorsTotal;
        DFstudy = setup.DFstudy; % default is the lowest dilution factor (usually DF1, location 4)
        obsMet = setup.observableMetabolite;
        costfun = setup.costfun; % default value is 1. No regularization
        lambda = setup.selectedLambda;
        selectedLambda = setup.selectedLambda;
        numpH = setup.numpHtested;
        sourceVm = setup.sourceVm;
        ode_pH = setup.ode_pH;
        wD = setup.weightData;
        wDesp = setup.weightDataEsp;
        wH = setup.weightHaldane;
        wL = setup.selectedLambda;
        plotEachSimCF = setup.plotEachSimCF;
        simAllProfiles = setup.simAllProfiles;


        switch enzyme

            case 'gapdhr'
                simNADH = cell(DFs,numpH);
                expNADH = cell(DFs,numpH);
                simGAPDHr = cell(DFs,numpH);
                expGAPDHr = cell(DFs,numpH);

                % simulations loop for each pH value
                for j = 1:numpH
                    % select required data
                    data.Vmaxs = data.Vmax(j,:);
                    data.NADH = data.conc_mean(j,:);
                    data.Vprofs = data.RRs(j,:);
                    data.tempTime = data.time(j,:);
                    % inputs to be selected
                    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ; %setup.pH_Keq_gapdh(i);
                    data.KeqPGK = setup.pH_Keq_pgk;
                    data.chosenKeqGAPDH = data.KeqGAPDH(j);
                    data.chosenKeqPGK = data.KeqPGK(j);

                    % selecting the right parameters
                    for temp11 = 1
                    switch j
                        case 1
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(5);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(6);
                        case 2
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(7);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(8);
                        case 3
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(9);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(10);
                        case 4
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(11);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(12);
                        case 5
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(13);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(14);
                        case 6
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(15);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(16);
                        case 7
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(17);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(18);
                        case 8
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(19);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(20);
                        case 9
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(21);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(22);
                        case 10
                            xassay = zeros(1,6);
                            xassay(1) = x_temp(23);
                            xassay(2) = x_temp(1);
                            xassay(3) = x_temp(2);
                            xassay(4) = x_temp(3);
                            xassay(5) = x_temp(4);
                            xassay(6) = x_temp(24);
                        otherwise
                            disp('Something went wrong is selecting the pH value');
                    end
                    end

                    % simulations
                    for i = DFstudy
                        % recall vmax for the specific value and simulate
                        data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,i); % vmax from the highest DF is taken and then divided
                        data.chosenLink = data.DF(1,i);
                        data.chosenNADini = data.NADH{i}(1);
                        data.chosenDF = data.DF(j,i);
                        setup.excessPGK = 1;

                        data.NADH = data.conc_mean(j,:);
                        data.Vprofs = data.RRs(j,:);
                        data.tempTime = data.time(j,:);                
                        data.i = j;

                        % simulate metabolites
                        [simResult] = simSys(xassay,data,setup);
                        % calculation of reaction rates
                        [vObs,~] = calcRates(xassay,simResult,data,setup);   
                        % cost function (NADHexp + vGAPDHr)
                        simTime = simResult.t;
                        simMet = simResult.y(:,obsMet);
                        simRate = vObs;

                        simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                        simGAPDHr{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                        expNADH{i,j} = data.NADH{i};
                        expGAPDHr{i,j} = -data.Vprofs{i};
                    end
                    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                    if plotEachSimCF == 1
                        if simAllProfiles == 0
                            figure
                            subplot(1,2,1)
                            for i = DFstudy
                                simRes = simResult;
                                plot(simRes.t,simRes.y(:,8),'-')
                                hold on
                                plot(data.time{j,i},data.conc_mean{j,i},'k+')
                            end
                            title('NADH')
                            subplot(1,2,2)
                            for i = DFstudy
                                simRRs = vObs;
                                plot(simRes.t,simRRs,'-')
                                hold on
                                plot(data.time{j,i},-data.RRs{j,i},'k+')
                            end
                            title('v_{GAPDHr}')       
                            suptitle(erase(sprintf('pH = %d',setup.pH_vals(j)),"0000e+00"));
                        elseif simAllProfiles == 1
                            figure
                            subplot(1,2,1)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADH{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.y(:,8),'-')
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k+')
                                hold on
                            end
                            title('NADH')
                            subplot(1,2,2)
                            for i = DFstudy
                                plot(data.tempTime{i}, simGAPDHr{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.v, '-')
                                hold on
                                plot(data.tempTime{i}, expGAPDHr{i,j},'k+')
                                hold on
                            end
                            title('v_{apparent.GAPDHr}')
                            suptitle(erase(sprintf('pH = %d',setup.pH_vals(j)),"0000e+00"));
        %                     suptitle(erase(sprintf('pH = %d',pHarray(data.i)),"0000e+00"));

        %                     figure
        %                     subplot(1,2,1)
        %                     for i = DFstudy
        %             %             plot(data.tempTime{j}, simNADH{j},'-')
        %                         plot(simRes{i}.t, simRes{i}.y(:,8),'-')
        %                         hold on
        %                         plot(data.tempTime{i}, expNADH{i},'k+')
        %                         hold on
        %                     end
        %                     title('NADH')
        %                     subplot(1,2,2)
        %                     for i = DFstudy
        %             %             plot(data.tempTime{j}, simGAPDHr{j},'-')
        %                         plot(simRes{i}.t, simRes{i}.v, '-')
        %                         hold on
        %                         plot(data.tempTime{i}, expGAPDHr{i},'k+')
        %                         hold on
        %                     end
        %                     title('v_{apparent.GAPDHr}')
        %                     suptitle(erase(sprintf('pH = %d',pHarray(data.i)),"0000e+00"));
                        end
                    end
                    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                end

                % calculation cost function
                switch costfun
                    case 1 % DF1
                            errorNADH1 = simNADH{4,1} - expNADH{4,1};
                            errorNADH2 = simNADH{4,2} - expNADH{4,2};
                            errorNADH3 = simNADH{4,3} - expNADH{4,3};
                            errorNADH4 = simNADH{4,4} - expNADH{4,4};
                            errorNADH5 = simNADH{4,5} - expNADH{4,5};
                            errorNADH6 = simNADH{4,6} - expNADH{4,6};
                            errorNADH7 = simNADH{4,7} - expNADH{4,7};
                            errorNADH8 = simNADH{4,8} - expNADH{4,8};
                            errorNADH9 = simNADH{4,9} - expNADH{4,9};
                            errorNADH10 = simNADH{4,10} - expNADH{4,10};
                            errorNADH = [...
                                wDesp(1)*errorNADH1;
                                wDesp(2)*errorNADH2;
                                wDesp(3)*errorNADH3;
                                wDesp(4)*errorNADH4;
                                wDesp(5)*errorNADH5;
                                wDesp(6)*errorNADH6;
                                wDesp(7)*errorNADH7;
                                wDesp(8)*errorNADH8;
                                wDesp(9)*errorNADH9;
                                wDesp(10)*errorNADH10];
                            for temp1 = 1                    
                                Keq = setup.pH_Keq_gapdh_eQ; %[]
        %                         Keq = data.chosenKeqGAPDH; %[]
                                switch sourceVm
                                    case 'experimentalSlopesFixed'
                                        vmf1 = 10.^x_temp(5).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr1 = 10.^x_temp(6).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf2 = 10.^x_temp(7).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr2 = 10.^x_temp(8).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf3 = 10.^x_temp(9).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr3 = 10.^x_temp(10).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf4 = 10.^x_temp(11).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr4 = 10.^x_temp(12).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf5 = 10.^x_temp(13).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr5 = 10.^x_temp(14).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf6 = 10.^x_temp(15).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr6 = 10.^x_temp(16).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf7 = 10.^x_temp(17).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr7 = 10.^x_temp(18).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf8 = 10.^x_temp(19).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr8 = 10.^x_temp(20).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf9 = 10.^x_temp(21).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr9 = 10.^x_temp(22).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                        vmf10 = 10.^x_temp(23).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                        vmr10 = 10.^x_temp(24).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                                    otherwise
                                        disp('No source for vmax has been selected');
                                end
                                ks1 = 10 .^ x_temp(1) .* 2.48; % mM %k_gap
                                ks2 = 10 .^ x_temp(3) .* 2.92; %mM %k_nad
                                kp1 = 10 .^ x_temp(2) .* 1.18; % mM %k_bpg
                                kp2 = 10 .^ x_temp(4) .* 0.022; % mM %k_nadh
                                switch ode_pH
                                    case 'on'
                                        H_effect = zeros(1:10);
                                        for j = 1:numpH
                                            H_effect(j) = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                                        end
                                        errorHaldane1 = (Keq(1) - (vmf1 * kp1 * kp2 * H_effect(1)) / (vmr1 * ks1 * ks2) );
                                        errorHaldane2 = (Keq(2) - (vmf2 * kp1 * kp2 * H_effect(2)) / (vmr2 * ks1 * ks2) );
                                        errorHaldane3 = (Keq(3) - (vmf3 * kp1 * kp2 * H_effect(3)) / (vmr3 * ks1 * ks2) );
                                        errorHaldane4 = (Keq(4) - (vmf4 * kp1 * kp2 * H_effect(4)) / (vmr4 * ks1 * ks2) );
                                        errorHaldane5 = (Keq(5) - (vmf5 * kp1 * kp2 * H_effect(5)) / (vmr5 * ks1 * ks2) );
                                        errorHaldane6 = (Keq(6) - (vmf6 * kp1 * kp2 * H_effect(6)) / (vmr6 * ks1 * ks2) );
                                        errorHaldane7 = (Keq(7) - (vmf7 * kp1 * kp2 * H_effect(7)) / (vmr7 * ks1 * ks2) );
                                        errorHaldane8 = (Keq(8) - (vmf8 * kp1 * kp2 * H_effect(8)) / (vmr8 * ks1 * ks2) );
                                        errorHaldane9 = (Keq(9) - (vmf9 * kp1 * kp2 * H_effect(9)) / (vmr9 * ks1 * ks2) );
                                        errorHaldane10 =(Keq(10) - (vmf10 * kp1 * kp2 * H_effect(10)) / (vmr10 * ks1 * ks2) );
                                        errorHaldane = [...
                                            errorHaldane1;
                                            errorHaldane2;
                                            errorHaldane3;
                                            errorHaldane4;
                                            errorHaldane5;
                                            errorHaldane6;
                                            errorHaldane7;
                                            errorHaldane8;
                                            errorHaldane9;
                                            errorHaldane10];
                                %         errorHaldane = sum(abs(errorHaldane));
                                    otherwise
                                % % % %         errorHaldane = wH * (Keq - (vmf * kp1 * kp2) / (vmr * ks1 * ks2) );
                                        disp('No source for vmax has been selected');
                                end
                            end
                            errorReg = lambda * x_temp';

                    otherwise
                        disp('No specific cost function has been appointed');
                end
                error = [
                    wD * errorNADH;
                    wH * errorHaldane;
                    wL * errorReg;
                    ];        
        %         disp(lambda);

            case 'eno'
                % simulations loop for each pH value
                for j = 1:numpH
        %         for j = 11:12
                    % select required data
                    data.Vmaxs = data.Vmax(j,:);
                    data.PEP = data.conc_mean(j,:);
                    data.Vprofs = data.RRs(j,:);
                    data.tempTime = data.time(j,:);
                    % inputs to be selected
                    data.chosenKeq = setup.keq(j);              
                    data.i = j;
                    % selecting the right parameters
                    for temp11 = 1
                    switch j
                        case 1
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                        case 2
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(4);
                        case 3
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(5);
                        case 4
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(6);
                        case 5
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(7);
                        case 6
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(8);
                        case 7
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(9);
                        case 8
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(10);
                        case 9
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(11);
                        case 10
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(12);
                        case 11
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(13);
                        case 12
                            xassay = zeros(1,3);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(14);
                        otherwise
                            disp('Something went wrong is selecting the pH value');
                    end
                    end

                    % simulations
        %             simRes = cell(numpH,DFstudy);
                    for i = DFstudy
                        % recall vmax for the specific value and simulate
                        data.chosenDF = data.DF(j,i);
                        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
        %                 data.chosenLink = data.DF(1,i);
                        data.chosenPEPini = data.PEP{i}(1);
        %                 setup.excessPGK = 1;

        %                 data.PEP = data.conc_mean(j,:);
        %                 data.Vprofs = data.RRs(j,:);
        %                 data.tempTime = data.time(j,:);                
        %                 data.i = j;

                        % simulate metabolites
                        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                        % calculation of reaction rates
                        [vObs,~] = calcRates(xassay,simResult,data,setup);   
                        % cost function (NADHexp + vGAPDHr)
                        simTime = simResult.t;
                        simMet = simResult.y(:,obsMet);
                        simRate = vObs;

                        simPEP{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                        simENO{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                        expPEP{i,j} = data.PEP{i};
                        expENO{i,j} = data.Vprofs{i};
                    end
                end

                for j = 1:numpH
                    if plotEachSimCF == 1
                        if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                            if j == 1
                                h101 = figure(101);
                            end
        %                     set(101);
                            set(0,'CurrentFigure', h101)
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
                                ylim([0 1.2])
                            end
                            if j == numpH
                                suptitle('PEP concentration [mM] vs assay time [s]')
                            end
        %                 end
        %             end
        %         end
                            if j == 1
                                h102 = figure(102);
                            end
        %                     set(102)
                            set(0,'CurrentFigure', h102)
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
                                ylim([0 1.5E-3])
                            end
                            if j == numpH
                                suptitle('ENO reaction rate [mM s^{-1}] vs assay time [s]')
                            end                    
        %                 end
        %             end
        %         end

        %                     figure
        %                     subplot(1,2,1)
        %                     for i = DFstudy
        %                         simRes = simResult;
        %                         plot(simRes.t,simRes.y(:,obsMet),'-')
        %                         hold on
        %                         plot(data.time{j,i},data.conc_mean{j,i},'k+')
        %                     end
        %                     title('PEP')
        %                     subplot(1,2,2)
        %                     for i = DFstudy
        %                         simRRs = vObs;
        %                         plot(simRes.t,simRRs,'-')
        %                         hold on
        %                         plot(data.time{j,i},-data.RRs{j,i},'k+')
        %                     end
        %                     title('v_{ENO}')       
        %                     suptitle(erase(sprintf('pH = %d',setup.pH_vals(j)),"0000e+00"));
                        elseif((simAllProfiles == 0)&&(setup.plotEachSim == 1))
                            figure
                            subplot(1,2,1)
                            for i = DFstudy
                                plot(data.tempTime{i}, simPEP{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.y(:,8),'-')
                                hold on
                                plot(data.tempTime{i}, expPEP{i,j},'k+')
                                hold on
                            end
                            title('PEP')
                            subplot(1,2,2)
                            for i = DFstudy
                                plot(data.tempTime{i}, simENO{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.v, '-')
                                hold on
                                plot(data.tempTime{i}, expENO{i,j},'k+')
                                hold on
                            end
                            title('v_{apparent.ENO}')
                            % % % % ONLY FOR ENO
                            suptitle(erase(sprintf('pH = %d',setup.fullpHarray(j)),"0000e+00"));
                            % % % % ONLY FOR ENO
                        end
                    end
                end

                % calculation cost function
                switch costfun
                    case 1 % DF1
                            errorHaldane = 0;
                            errorPEP1 = simPEP{4,1} - expPEP{4,1};
                            errorPEP2 = simPEP{4,2} - expPEP{4,2};
                            errorPEP3 = simPEP{4,3} - expPEP{4,3};
                            errorPEP4 = simPEP{4,4} - expPEP{4,4};
                            errorPEP5 = simPEP{4,5} - expPEP{4,5};
                            errorPEP6 = simPEP{4,6} - expPEP{4,6};
                            errorPEP7 = simPEP{4,7} - expPEP{4,7};
                            errorPEP8 = simPEP{4,8} - expPEP{4,8};
                            errorPEP9 = simPEP{4,9} - expPEP{4,9};
                            errorPEP10 = simPEP{4,10} - expPEP{4,10};
                            errorPEP11 = simPEP{4,11} - expPEP{4,11};
                            errorPEP12 = simPEP{4,12} - expPEP{4,12};
                                errorPEP1_2 = simPEP{3,1} - expPEP{3,1};
                                errorPEP2_2 = simPEP{3,2} - expPEP{3,2};
                                errorPEP3_2 = simPEP{3,3} - expPEP{3,3};
                                errorPEP4_2 = simPEP{3,4} - expPEP{3,4};
                                errorPEP5_2 = simPEP{3,5} - expPEP{3,5};
                                errorPEP6_2 = simPEP{3,6} - expPEP{3,6};
                                errorPEP7_2 = simPEP{3,7} - expPEP{3,7};
                                errorPEP8_2 = simPEP{3,8} - expPEP{3,8};
                                errorPEP9_2 = simPEP{3,9} - expPEP{3,9};
                                errorPEP10_2 = simPEP{3,10} - expPEP{3,10};
                                errorPEP11_2 = simPEP{3,11} - expPEP{3,11};
                                errorPEP12_2 = simPEP{3,12} - expPEP{3,12};
                            errorPEP = [...
                                wDesp(1)*errorPEP1; wDesp(1)*errorPEP1_2;
                                wDesp(2)*errorPEP2; wDesp(1)*errorPEP2_2;
                                wDesp(3)*errorPEP3; wDesp(1)*errorPEP3_2;
                                wDesp(4)*errorPEP4; wDesp(1)*errorPEP4_2;
                                wDesp(5)*errorPEP5; wDesp(1)*errorPEP5_2;
                                wDesp(6)*errorPEP6; wDesp(1)*errorPEP6_2;
                                wDesp(7)*errorPEP7; wDesp(1)*errorPEP7_2;
                                wDesp(8)*errorPEP8; wDesp(1)*errorPEP8_2;
                                wDesp(9)*errorPEP9; wDesp(1)*errorPEP9_2;
                                wDesp(12)*errorPEP12; wDesp(12)*errorPEP12_2;
                                wDesp(11)*errorPEP11; wDesp(10)*errorPEP11_2;
                                wDesp(10)*errorPEP10; wDesp(11)*errorPEP10_2];
                            errorReg = lambda * x_temp';

                    otherwise
                        disp('No specific cost function has been appointed');
                end
                error = [
                    wD * errorPEP;
                    wH * errorHaldane;
                    wL * errorReg;
                    ];        
        %         disp(lambda); 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            case 'hxk'
                % simulations loop for each pH value
                for j = 1:numpH
        %         for j = 11:12
                    % select required data
                    data.Vmaxs = data.Vmax(j,:);
                    data.NADPH = data.conc_mean(j,:);
                    data.Vprofs = data.RRs(j,:);
                    data.tempTime = data.time(j,:);
                    % inputs to be selected
        %             data.chosenKeq = setup.keq(j);   
                    data.chosenKeq_HXK = setup.Keq_HXK(j);
                    data.chosenKeq_G6PDH = setup.Keq_G6PDH(j);
                    data.i = j;
                    % selecting the right parameters
                    for temp11 = 1
                    switch j
                        case 1
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                        case 2
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(6);
                        case 3
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(7);
                        case 4
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(8);
                        case 5
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(9);
                        case 6
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(10);
                        case 7
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(11);
                        case 8
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(12);
                        case 9
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(13);
                        case 10
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(14);
                        case 11
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(15);
                        case 12
                            xassay = zeros(1,5);
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(16);
                        otherwise
                            disp('Something went wrong is selecting the pH value');
                    end
                    end

                    % simulations
        %             simRes = cell(numpH,DFstudy);
                    for i = DFstudy
                        % recall vmax for the specific value and simulate
                        data.chosenDF = data.DF(j,i);
                        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
        %                 data.chosenLink = data.DF(1,i);
                        data.chosenNADPHini = data.NADPH{i}(1);
        %                 setup.excessPGK = 1;
        %                 data.PEP = data.conc_mean(j,:);
        %                 data.Vprofs = data.RRs(j,:);
        %                 data.tempTime = data.time(j,:);                
        %                 data.i = j;

                        % simulate metabolites
                        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                        % calculation of reaction rates
                        [vObs,~] = calcRates(xassay,simResult,data,setup);   
                        % cost function (NADHexp + vGAPDHr)
                        simTime = simResult.t;
                        simMet = simResult.y(:,obsMet);
                        simRate = vObs;

                        simNADPH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                        simG6PDH{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                        expNADPH{i,j} = data.NADPH{i};
                        expG6PDH{i,j} = data.Vprofs{i};
                    end
                end

                for j = 1:numpH
                    if plotEachSimCF == 1
                        if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                            if j == 1
                                h101 = figure(101);
                            end
        %                     set(101);
                            set(0,'CurrentFigure', h101);
        %                     set(gcf,'color','w');
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADPH{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, expNADPH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                            end
                            if j == numpH
                                suptitle('NADPH concentration [mM] vs assay time [s]')
                            end
        %                 end
        %             end
        %         end
                            if j == 1
                                h102 = figure(102);
                            end
        %                     set(102)
                            set(0,'CurrentFigure', h102);
        %                     set(gcf,'color','w');
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simG6PDH{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, expG6PDH{i,j},'k.','MarkerSize',4)
                                ylim([0 5E-4])
                            end
                            if j == numpH
                                suptitle('G6PDH reaction rate [mM s^{-1}] vs assay time [s]')
                            end                    
        %                 end
        %             end
        %         end

        %                     figure
        %                     subplot(1,2,1)
        %                     for i = DFstudy
        %                         simRes = simResult;
        %                         plot(simRes.t,simRes.y(:,obsMet),'-')
        %                         hold on
        %                         plot(data.time{j,i},data.conc_mean{j,i},'k+')
        %                     end
        %                     title('PEP')
        %                     subplot(1,2,2)
        %                     for i = DFstudy
        %                         simRRs = vObs;
        %                         plot(simRes.t,simRRs,'-')
        %                         hold on
        %                         plot(data.time{j,i},-data.RRs{j,i},'k+')
        %                     end
        %                     title('v_{ENO}')       
        %                     suptitle(erase(sprintf('pH = %d',setup.pH_vals(j)),"0000e+00"));
                        elseif((simAllProfiles == 0)&&(setup.plotEachSim == 1))
                            figure
                            subplot(1,2,1)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADPH{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.y(:,8),'-')
                                hold on
                                plot(data.tempTime{i}, expNADPH{i,j},'k+')
                                hold on
                            end
                            title('PEP')
                            subplot(1,2,2)
                            for i = DFstudy
                                plot(data.tempTime{i}, simG6PDH{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.v, '-')
                                hold on
                                plot(data.tempTime{i}, expG6PDH{i,j},'k+')
                                hold on
                            end
                            title('v_{apparent.G6PDH}')
                            % % % % ONLY FOR ENO
                            suptitle(erase(sprintf('pH = %d',setup.fullpHarray(j)),"0000e+00"));
                            % % % % ONLY FOR ENO
                        end
                    end
                end

                % calculation cost function
                switch costfun
                    case 1 % DF1
                            errorHaldane = 0;
                            errorNADPH1 = simNADPH{4,1} - expNADPH{4,1};
                            errorNADPH2 = simNADPH{4,2} - expNADPH{4,2};
                            errorNADPH3 = simNADPH{4,3} - expNADPH{4,3};
                            errorNADPH4 = simNADPH{4,4} - expNADPH{4,4};
                            errorNADPH5 = simNADPH{4,5} - expNADPH{4,5};
                            errorNADPH6 = simNADPH{4,6} - expNADPH{4,6};
                            errorNADPH7 = simNADPH{4,7} - expNADPH{4,7};
                            errorNADPH8 = simNADPH{4,8} - expNADPH{4,8};
                            errorNADPH9 = simNADPH{4,9} - expNADPH{4,9};
                            errorNADPH10 = simNADPH{4,10} - expNADPH{4,10};
                            errorNADPH11 = simNADPH{4,11} - expNADPH{4,11};
                            errorNADPH12 = simNADPH{4,12} - expNADPH{4,12};
                                errorNADPH1_2 = simNADPH{3,1} - expNADPH{3,1};
                                errorNADPH2_2 = simNADPH{3,2} - expNADPH{3,2};
                                errorNADPH3_2 = simNADPH{3,3} - expNADPH{3,3};
                                errorNADPH4_2 = simNADPH{3,4} - expNADPH{3,4};
                                errorNADPH5_2 = simNADPH{3,5} - expNADPH{3,5};
                                errorNADPH6_2 = simNADPH{3,6} - expNADPH{3,6};
                                errorNADPH7_2 = simNADPH{3,7} - expNADPH{3,7};
                                errorNADPH8_2 = simNADPH{3,8} - expNADPH{3,8};
                                errorNADPH9_2 = simNADPH{3,9} - expNADPH{3,9};
                                errorNADPH10_2 = simNADPH{3,10} - expNADPH{3,10};
                                errorNADPH11_2 = simNADPH{3,11} - expNADPH{3,11};
                                errorNADPH12_2 = simNADPH{3,12} - expNADPH{3,12};
                            errorNADPH = [...
                                wDesp(1)*errorNADPH1; wDesp(1)*errorNADPH1_2;
                                wDesp(2)*errorNADPH2; wDesp(2)*errorNADPH2_2;
                                wDesp(3)*errorNADPH3; wDesp(3)*errorNADPH3_2;
                                wDesp(4)*errorNADPH4; wDesp(4)*errorNADPH4_2;
                                wDesp(5)*errorNADPH5; wDesp(5)*errorNADPH5_2;
                                wDesp(6)*errorNADPH6; wDesp(6)*errorNADPH6_2;
                                wDesp(7)*errorNADPH7; wDesp(7)*errorNADPH7_2;
                                wDesp(8)*errorNADPH8; wDesp(8)*errorNADPH8_2;
                                wDesp(9)*errorNADPH9; wDesp(9)*errorNADPH9_2;
                                wDesp(12)*errorNADPH12; wDesp(12)*errorNADPH12_2;
                                wDesp(11)*errorNADPH11; wDesp(11)*errorNADPH11_2;
                                wDesp(10)*errorNADPH10; wDesp(10)*errorNADPH10_2];
                            errorReg = lambda * x_temp';

                    otherwise
                        disp('No specific cost function has been appointed');
                end
                error = [
                    wD * errorNADPH;
                    wH * errorHaldane;
                    wL * errorReg;
                    ];        
        %         disp(lambda); 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            case 'ald'
                % simulations loop for each pH value
                for j = 1:numpH
        %         for j = 11:12
                    % select required data
                    data.Vmaxs = data.Vmax(j,:);
                    data.NADH = data.conc_mean(j,:);
                    data.Vprofs = data.RRs(j,:);
                    data.tempTime = data.time(j,:);
                    % inputs to be selected
        %             data.chosenKeq = setup.keq(j);   
                    data.chosenKeq_FBA = setup.Keq_FBA(j);% = [1.0E-3 7.9E-4 7.3E-4 7E-4 6.8E-4 6.7E-4 6.6E-4 6.5E-4];  %dir+
                    data.chosenKeq_TPI = setup.Keq_TPI(j);% = [1/(8.31) 1/(8.97) 1/(9.16) 1/(9.26) 1/(9.33) 1/(9.36) 1/(9.38) 1/(9.39)];  %dir-
                    data.chosenKeq_GPD = setup.Keq_GPD(j);% = [1/(4.2E-6) 1/(1.5E-5) 1/(2.7E-5) 1/(4.7E-5) 1/(7.9E-5) 1/(1.2E-4) 1/(1.6E-4) 1/(2.0E-4) ]; %dir-
                    data.i = j;
                    % selecting the right parameters
                    for temp11 = 1
                        xassay = zeros(1,4);
                    switch j
                        case 1
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(12);
                            xassay(6) = x_temp(13);
                        case 2
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(5);
                            xassay(5) = x_temp(12);
                            xassay(6) = x_temp(13);
                        case 3
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(6);
                            xassay(5) = x_temp(12);
                            xassay(6) = x_temp(13);
                        case 4
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(7);
                            xassay(5) = x_temp(12);
                            xassay(6) = x_temp(13);
                        case 5
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(8);
                            xassay(5) = x_temp(12);
                            xassay(6) = x_temp(13);
                        case 6
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(9);
                            xassay(5) = x_temp(12);
                            xassay(6) = x_temp(13);
                        case 7
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(10);
                            xassay(5) = x_temp(12);
                            xassay(6) = x_temp(13);
                        case 8
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(11);
                            xassay(5) = x_temp(12);
                            xassay(6) = x_temp(13);
                        otherwise
                            disp('Something went wrong is selecting the pH value');
                    end
                    end

                    % simulations
        %             simRes = cell(numpH,DFstudy);
                    for i = DFstudy
                        % recall vmax for the specific value and simulate
                        data.chosenDF = data.DF(j,i);
                        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
        %                 data.chosenLink = data.DF(1,i);
                        data.chosenNADHini = data.NADH{i}(1);
        %                 setup.excessPGK = 1;
        %                 data.PEP = data.conc_mean(j,:);
        %                 data.Vprofs = data.RRs(j,:);
        %                 data.tempTime = data.time(j,:);                
        %                 data.i = j;

                        % simulate metabolites
                        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                        % calculation of reaction rates
                        [vObs,~] = calcRates(xassay,simResult,data,setup);   
                        % cost function (NADHexp + vGAPDHr)
                        simTime = simResult.t;
                        simMet = simResult.y(:,obsMet);
                        simRate = vObs(:,2);

                        simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                        simGPD{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                        expNADH{i,j} = data.NADH{i};
                        expGPD{i,j} = data.Vprofs{i};
        %                 % diplay
        %                 for temp50 = 1
        %                     figure
        %                     for o = 1:6
        %                         subplot(3,3,o)
        %                         plot(simResult.t, simResult.y(:,o))
        %                         title(setup.PSAmets{o})
        %                     end
        %                     subplot(3,3,7)
        %                     plot(simResult.t, vObs(:,1))
        %                     title('v_{ALD}')
        %                     subplot(3,3,8)
        %                     plot(simResult.t, vObs(:,2))
        %                     title('v_{GPD}')
        %                     subplot(3,3,9)
        %                     plot(simResult.t, vObs(:,3))
        %                     title('v_{TPI}')
        %                 end
                    end
                end

                for j = 1:numpH
                    data.tempTime = data.time(j,:);
                    if plotEachSimCF == 1
                        if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                            if j == 1
                                h101 = figure(101);
                            end
        %                     set(101);
                            set(0,'CurrentFigure', h101);
        %                     set(gcf,'color','w');
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                            end
                            if j == numpH
                                suptitle('NADH concentration [mM] vs assay time [s]')
                            end
        %                 end
        %             end
        %         end
                            if j == 1
                                h102 = figure(102);
                            end
        %                     set(102)
                            set(0,'CurrentFigure', h102);
        %                     set(gcf,'color','w');
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
        %                         ylim([0 5E-4])
                            end
                            if j == numpH
                                suptitle('GPD reaction rate [mM s^{-1}] vs assay time [s]')
                            end                    
        %                 end
        %             end
        %         end

        %                     figure
        %                     subplot(1,2,1)
        %                     for i = DFstudy
        %                         simRes = simResult;
        %                         plot(simRes.t,simRes.y(:,obsMet),'-')
        %                         hold on
        %                         plot(data.time{j,i},data.conc_mean{j,i},'k+')
        %                     end
        %                     title('PEP')
        %                     subplot(1,2,2)
        %                     for i = DFstudy
        %                         simRRs = vObs;
        %                         plot(simRes.t,simRRs,'-')
        %                         hold on
        %                         plot(data.time{j,i},-data.RRs{j,i},'k+')
        %                     end
        %                     title('v_{ENO}')       
        %                     suptitle(erase(sprintf('pH = %d',setup.pH_vals(j)),"0000e+00"));
                        elseif((simAllProfiles == 0)&&(setup.plotEachSim == 1))
                            figure
                            subplot(1,2,1)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADH{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.y(:,8),'-')
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k+')
                                hold on
                            end
                            title('NADH')
                            subplot(1,2,2)
                            for i = DFstudy
                                plot(data.tempTime{i}, simGPD{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.v, '-')
                                hold on
                                plot(data.tempTime{i}, expGPD{i,j},'k+')
                                hold on
                            end
                            title('v_{apparent.GPD}')
                            % % % % ONLY FOR ENO
                            suptitle(erase(sprintf('pH = %d',setup.fullpHarray(j)),"0000e+00"));
                            % % % % ONLY FOR ENO
                        end
                    end
                end

                % calculation cost function
                switch costfun
                    case 1 % DF1
                            errorHaldane = 0;
                            errorNADH1 = simNADH{4,1} - expNADH{4,1};
                            errorNADH2 = simNADH{4,2} - expNADH{4,2};
                            errorNADH3 = simNADH{4,3} - expNADH{4,3};
                            errorNADH4 = simNADH{4,4} - expNADH{4,4};
                            errorNADH5 = simNADH{4,5} - expNADH{4,5};
                            errorNADH6 = simNADH{4,6} - expNADH{4,6};
                            errorNADH7 = simNADH{4,7} - expNADH{4,7};
                            errorNADH8 = simNADH{4,8} - expNADH{4,8};
                                errorNADH1_2 = simNADH{3,1} - expNADH{3,1};
                                errorNADH2_2 = simNADH{3,2} - expNADH{3,2};
                                errorNADH3_2 = simNADH{3,3} - expNADH{3,3};
                                errorNADH4_2 = simNADH{3,4} - expNADH{3,4};
                                errorNADH5_2 = simNADH{3,5} - expNADH{3,5};
                                errorNADH6_2 = simNADH{3,6} - expNADH{3,6};
                                errorNADH7_2 = simNADH{3,7} - expNADH{3,7};
                                errorNADH8_2 = simNADH{3,8} - expNADH{3,8};
                            errorNADH = [...
                                wDesp(1)*errorNADH1; wDesp(1)*errorNADH1_2; 
                                wDesp(2)*errorNADH2; wDesp(2)*errorNADH2_2;
                                wDesp(3)*errorNADH3; wDesp(3)*errorNADH3_2;
                                wDesp(4)*errorNADH4; wDesp(4)*errorNADH4_2;
                                wDesp(5)*errorNADH5; wDesp(5)*errorNADH5_2;
                                wDesp(6)*errorNADH6; wDesp(6)*errorNADH6_2;
                                wDesp(7)*errorNADH7; wDesp(7)*errorNADH7_2;
                                wDesp(8)*errorNADH8; wDesp(8)*errorNADH8_2];
                            errorReg = lambda * x_temp';

                    otherwise
                        disp('No specific cost function has been appointed');
                end
                error = [
                    wD * errorNADH;
                    wH * errorHaldane;
                    wL * errorReg;
                    ];        
        %         disp(lambda); 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            case 'pyk'
                % simulations loop for each pH value
                for j = 1:numpH
        %         for j = 11:12
                    % select required data
                    data.Vmaxs = data.Vmax(j,:);
        %             data.Vmaxs = data.Vmax(8,:);
                    data.NADH = data.conc_mean(j,:);
                    data.Vprofs = data.RRs(j,:);
                    data.tempTime = data.time(j,:);
                    % inputs to be selected
        %             data.chosenKeq = setup.keq(j); 
                    data.chosenKeq_PYK = setup.Keq_PYK(j);% = [1/(3.1E-6) 1/(3.5E-6) 1/(3.9E-6) 1/(4.5E-6) 1/(6.5E-6) 1/(9.7E-6) 1/(1.6E-5) 1/(2.5E-5) 1/(4.1E-5) 1/(5.9E-5) 1/(7.9E-5) 1/(9.6E-5)];  %dir-
                    data.chosenKeq_LDH = setup.Keq_LDH(j);%
                    data.i = j;
                    % selecting the right parameters
                    for temp11 = 1
                        xassay = zeros(1,7);
                    switch j
                        case 1
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                        case 2
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(8);
                        case 3
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(9);
                        case 4
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(10);
                        case 5
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(11);
                        case 6
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(12);
                        case 7
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(13);
                        case 8
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(14);
                        case 9
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(15);
                        case 10
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(16);
                        case 11
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(17);
                        case 12
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(18);
                        otherwise
                            disp('Something went wrong is selecting the pH value');
                    end
                    end

                    % simulations
        %             simRes = cell(numpH,DFstudy);
        %             j, data.Vmaxs,
                    for i = DFstudy
                        % recall vmax for the specific value and simulate
                        data.chosenDF = data.DF(j,i);
                        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
                        data.chosenNADHini = data.NADH{i}(1);                
                        % simulate metabolites
        % % % %                 disp(j);
                        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                        tempResult{k} = simResult;
                        % calculation of reaction rates
                        [vObs,~] = calcRates(xassay,simResult,data,setup);   
                        % cost function (NADHexp + vGAPDHr)
                        simTime = simResult.t;
                        simMet = simResult.y(:,obsMet);
                        simRate = vObs(:,2);

                        simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                        simLDH{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                        expNADH{i,j} = data.NADH{i};
                        expLDH{i,j} = data.Vprofs{i};
        %                 % diplay
        %                 for temp50 = 1
        %                     figure
        %                     for o = 1:6
        %                         subplot(3,3,o)
        %                         plot(simResult.t, simResult.y(:,o))
        %                         title(setup.PSAmets{o})
        %                     end
        %                     subplot(3,3,7)
        %                     plot(simResult.t, vObs(:,1))
        %                     title('v_{ALD}')
        %                     subplot(3,3,8)
        %                     plot(simResult.t, vObs(:,2))
        %                     title('v_{GPD}')
        %                     subplot(3,3,9)
        %                     plot(simResult.t, vObs(:,3))
        %                     title('v_{TPI}')
        %                 end
                    end
                end

                for j = 1:numpH
                    data.tempTime = data.time(j,:);
                    if plotEachSimCF == 1
                        if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                            if j == 1
                                h101 = figure(101);
                            end
        %                     set(101);
                            set(0,'CurrentFigure', h101);
        %                     set(gcf,'color','w');
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                            end
                            if j == numpH
                                suptitle('NADH concentration [mM] vs assay time [s]')
                            end
        %                 end
        %             end
        %         end
                            if j == 1
                                h102 = figure(102);
                            end
        %                     set(102)
                            set(0,'CurrentFigure', h102);
        %                     set(gcf,'color','w');
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simLDH{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, -expLDH{i,j},'k.','MarkerSize',4)
        %                         ylim([0 5E-4])
                            end
                            if j == numpH
                                suptitle('LDH reaction rate [mM s^{-1}] vs assay time [s]')
                            end                    
        %                 end
        %             end
        %         end

        %                     figure
        %                     subplot(1,2,1)
        %                     for i = DFstudy
        %                         simRes = simResult;
        %                         plot(simRes.t,simRes.y(:,obsMet),'-')
        %                         hold on
        %                         plot(data.time{j,i},data.conc_mean{j,i},'k+')
        %                     end
        %                     title('PEP')
        %                     subplot(1,2,2)
        %                     for i = DFstudy
        %                         simRRs = vObs;
        %                         plot(simRes.t,simRRs,'-')
        %                         hold on
        %                         plot(data.time{j,i},-data.RRs{j,i},'k+')
        %                     end
        %                     title('v_{ENO}')       
        %                     suptitle(erase(sprintf('pH = %d',setup.pH_vals(j)),"0000e+00"));
                        elseif((simAllProfiles == 0)&&(setup.plotEachSim == 1))
                            figure
                            subplot(1,2,1)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADH{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.y(:,8),'-')
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k+')
                                hold on
                            end
                            title('NADH')
                            subplot(1,2,2)
                            for i = DFstudy
                                plot(data.tempTime{i}, simLDH{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.v, '-')
                                hold on
                                plot(data.tempTime{i}, -expLDH{i,j},'k+')
                                hold on
                            end
                            title('v_{apparent.LDH}')
                            % % % % ONLY FOR ENO
                            suptitle(erase(sprintf('pH = %d',setup.fullpHarray(j)),"0000e+00"));
                            % % % % ONLY FOR ENO
                        end
                    end
                end

                % calculation cost function
                switch costfun
                    case 1 % DF1
                            errorHaldane = 0;
                            errorNADH1 = simNADH{4,1} - expNADH{4,1};
                            errorNADH2 = simNADH{4,2} - expNADH{4,2};
                            errorNADH3 = simNADH{4,3} - expNADH{4,3};
                            errorNADH4 = simNADH{4,4} - expNADH{4,4};
                            errorNADH5 = simNADH{4,5} - expNADH{4,5};
                            errorNADH6 = simNADH{4,6} - expNADH{4,6};
                            errorNADH7 = simNADH{4,7} - expNADH{4,7};
                            errorNADH8 = simNADH{4,8} - expNADH{4,8};
                            errorNADH9 = simNADH{4,9} - expNADH{4,9};
                            errorNADH10 = simNADH{4,10} - expNADH{4,10};
                            errorNADH11 = simNADH{4,11} - expNADH{4,11};
                            errorNADH12 = simNADH{4,12} - expNADH{4,12};
                                errorNADH1_2 = simNADH{3,1} - expNADH{3,1};
                                errorNADH2_2 = simNADH{3,2} - expNADH{3,2};
                                errorNADH3_2 = simNADH{3,3} - expNADH{3,3};
                                errorNADH4_2 = simNADH{3,4} - expNADH{3,4};
                                errorNADH5_2 = simNADH{3,5} - expNADH{3,5};
                                errorNADH6_2 = simNADH{3,6} - expNADH{3,6};
                                errorNADH7_2 = simNADH{3,7} - expNADH{3,7};
                                errorNADH8_2 = simNADH{3,8} - expNADH{3,8};
                                errorNADH9_2 = simNADH{3,9} - expNADH{3,9};
                                errorNADH10_2 = simNADH{3,10} - expNADH{3,10};
                                errorNADH11_2 = simNADH{3,11} - expNADH{3,11};
                                errorNADH12_2 = simNADH{3,12} - expNADH{3,12};
                            errorNADH = [...
                                wDesp(1)*errorNADH1; wDesp(1)*errorNADH1_2; 
                                wDesp(2)*errorNADH2; wDesp(2)*errorNADH2_2;
                                wDesp(3)*errorNADH3; wDesp(3)*errorNADH3_2;
                                wDesp(4)*errorNADH4; wDesp(4)*errorNADH4_2;
                                wDesp(5)*errorNADH5; wDesp(5)*errorNADH5_2;
                                wDesp(6)*errorNADH6; wDesp(6)*errorNADH6_2;
                                wDesp(7)*errorNADH7; wDesp(7)*errorNADH7_2;
                                wDesp(8)*errorNADH8; wDesp(8)*errorNADH8_2;
                                wDesp(9)*errorNADH9; wDesp(9)*errorNADH9_2;
                                wDesp(10)*errorNADH10; wDesp(10)*errorNADH10_2;
                                wDesp(11)*errorNADH11; wDesp(11)*errorNADH11_2;
                                wDesp(12)*errorNADH12; wDesp(12)*errorNADH12_2];
                            errorReg = lambda * x_temp';
                    otherwise
                        disp('No specific cost function has been appointed');
                end
                error = [
                    wD * errorNADH;
                    wH * errorHaldane;
                    wL * errorReg;
                    ];        
        %         disp(lambda); 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            case 'pgi'
                % simulations loop for each pH value
                for j = 1:numpH
        %         for j = 11:12
                    % select required data
                    data.Vmaxs = data.Vmax(j,:);
                    data.NADH = data.conc_mean(j,:);
                    data.Vprofs = data.RRs(j,:);
                    data.tempTime = data.time(j,:);
                    % inputs to be selected
        %             data.chosenKeq = setup.keq(j); 
                    data.chosenKeq_PGI = setup.Keq_PGI(j);%
                    data.chosenKeq_PFK = setup.Keq_PFK(j);%
                    data.chosenKeq_FBA = setup.Keq_FBA(j);%
                    data.chosenKeq_TPI = setup.Keq_TPI(j);%
                    data.chosenKeq_GPD = setup.Keq_GPD(j);%
                    data.i = j;
                    % selecting the right parameters
                    for temp11 = 1
                        xassay = zeros(1,3);
        % % % %                 xassay = zeros(1,4);
                    switch j
                        case 1
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            % % % % xassay(4) = x_temp(15);
                        case 2
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(4);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(16);
                        case 3
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(5);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(17);
                        case 4
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(6);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(18);
                        case 5
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(7);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(19);
                        case 6
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(8);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(20);
                        case 7
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(9);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(21);
                        case 8
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(10);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(22);
                        case 9
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(11);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(23);
                        case 10
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(12);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(24);
                        case 11
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(13);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(25);
                        case 12
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(14);
        % % % %                     % % % % xassay(4) = x_temp(15);
                            % % % % xassay(4) = x_temp(26);
                        otherwise
                            disp('Something went wrong is selecting the pH value');
                    end
                    end

                    % simulations
        %             simRes = cell(numpH,DFstudy);
                    for i = DFstudy
                        % recall vmax for the specific value and simulate
                        data.chosenDF = data.DF(j,i);
                        data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
                        data.chosenNADHini = data.NADH{i}(1);                
                        % simulate metabolites
                        [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                        % calculation of reaction rates
                        [vObs,~] = calcRates(xassay,simResult,data,setup);   
                        % cost function (NADHexp + vGAPDHr)
                        simTime = simResult.t;
                        simMet = simResult.y(:,obsMet);
                        simRate = vObs(:,5);

                        simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                        simGPD{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                        expNADH{i,j} = data.NADH{i};
                        expGPD{i,j} = data.Vprofs{i};
        %                 % diplay
        %                 for temp50 = 1
        %                     figure
        %                     for o = 1:6
        %                         subplot(3,3,o)
        %                         plot(simResult.t, simResult.y(:,o))
        %                         title(setup.PSAmets{o})
        %                     end
        %                     subplot(3,3,7)
        %                     plot(simResult.t, vObs(:,1))
        %                     title('v_{ALD}')
        %                     subplot(3,3,8)
        %                     plot(simResult.t, vObs(:,2))
        %                     title('v_{GPD}')
        %                     subplot(3,3,9)
        %                     plot(simResult.t, vObs(:,3))
        %                     title('v_{TPI}')
        %                 end
                    end
                end

                for j = 1:numpH
                    data.tempTime = data.time(j,:);
                    if plotEachSimCF == 1
                        if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                            if j == 1
                                h101 = figure(101);
                            end
        %                     set(101);
                            set(0,'CurrentFigure', h101);
        %                     set(gcf,'color','w');
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                                xlim([0 1000])
                            end
                            if j == numpH
                                suptitle('NADH concentration [mM] vs assay time [s]')
                            end
        %                 end
        %             end
        %         end
                            if j == 1
                                h102 = figure(102);
                            end
        %                     set(102)
                            set(0,'CurrentFigure', h102);
        %                     set(gcf,'color','w');
                            subplot(3,4,j)
                            for i = DFstudy
                                plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2)
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
        %                         ylim([0 5E-4])
                            end
                            xlim([0 1000])
                            if j == numpH
                                suptitle('GPD reaction rate [mM s^{-1}] vs assay time [s]')
                            end                    
        %                 end
        %             end
        %         end

        %                     figure
        %                     subplot(1,2,1)
        %                     for i = DFstudy
        %                         simRes = simResult;
        %                         plot(simRes.t,simRes.y(:,obsMet),'-')
        %                         hold on
        %                         plot(data.time{j,i},data.conc_mean{j,i},'k+')
        %                     end
        %                     title('PEP')
        %                     subplot(1,2,2)
        %                     for i = DFstudy
        %                         simRRs = vObs;
        %                         plot(simRes.t,simRRs,'-')
        %                         hold on
        %                         plot(data.time{j,i},-data.RRs{j,i},'k+')
        %                     end
        %                     title('v_{ENO}')       
        %                     suptitle(erase(sprintf('pH = %d',setup.pH_vals(j)),"0000e+00"));
                        elseif((simAllProfiles == 0)&&(setup.plotEachSim == 1))
                            figure
                            subplot(1,2,1)
                            for i = DFstudy
                                plot(data.tempTime{i}, simNADH{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.y(:,8),'-')
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k+')
                                hold on
                            end
                            title('NADH')
                            xlim([0 1000])
                            ylim([0 0.15])
                            subplot(1,2,2)
                            for i = DFstudy
                                plot(data.tempTime{i}, simGPD{i,j},'-')
        %                         plot(simRes{j}.t, simRes{j}.v, '-')
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k+')
                                hold on
                            end
                            title('v_{apparent.GPD}')
                            xlim([0 1000])
                            % % % % ONLY FOR ENO
                            suptitle(erase(sprintf('pH = %d',setup.fullpHarray(j)),"0000e+00"));
                            % % % % ONLY FOR ENO
                        end
                    end
                end

                % calculation cost function
                switch costfun
                    case 1 % DF1
                            errorHaldane = 0;
                            errorNADH1 = simNADH{4,1} - expNADH{4,1};
                            errorNADH2 = simNADH{4,2} - expNADH{4,2};
                            errorNADH3 = simNADH{4,3} - expNADH{4,3};
                            errorNADH4 = simNADH{4,4} - expNADH{4,4};
                            errorNADH5 = simNADH{4,5} - expNADH{4,5};
                            errorNADH6 = simNADH{4,6} - expNADH{4,6};
                            errorNADH7 = simNADH{4,7} - expNADH{4,7};
                            errorNADH8 = simNADH{4,8} - expNADH{4,8};
                            errorNADH9 = simNADH{4,9} - expNADH{4,9};
                            errorNADH10 = simNADH{4,10} - expNADH{4,10};
                            errorNADH11 = simNADH{4,11} - expNADH{4,11};
                            errorNADH12 = simNADH{4,12} - expNADH{4,12};
        %                         errorNADH1_2 = simNADH{3,1} - expNADH{3,1};
        %                         errorNADH2_2 = simNADH{3,2} - expNADH{3,2};
        %                         errorNADH3_2 = simNADH{3,3} - expNADH{3,3};
        %                         errorNADH4_2 = simNADH{3,4} - expNADH{3,4};
        %                         errorNADH5_2 = simNADH{3,5} - expNADH{3,5};
        %                         errorNADH6_2 = simNADH{3,6} - expNADH{3,6};
        %                         errorNADH7_2 = simNADH{3,7} - expNADH{3,7};
        %                         errorNADH8_2 = simNADH{3,8} - expNADH{3,8};
        %                         errorNADH9_2 = simNADH{3,9} - expNADH{3,9};
        %                         errorNADH10_2 = simNADH{3,10} - expNADH{3,10};
        %                         errorNADH11_2 = simNADH{3,11} - expNADH{3,11};
        %                         errorNADH12_2 = simNADH{3,12} - expNADH{3,12};
                            errorNADH = [...
                                wDesp(1)*errorNADH1; %wDesp(1)*errorNADH1_2; 
                                wDesp(2)*errorNADH2; %wDesp(2)*errorNADH2_2;
                                wDesp(3)*errorNADH3; %wDesp(3)*errorNADH3_2;
                                wDesp(4)*errorNADH4; %wDesp(4)*errorNADH4_2;
                                wDesp(5)*errorNADH5; %wDesp(5)*errorNADH5_2;
                                wDesp(6)*errorNADH6; %wDesp(6)*errorNADH6_2;
                                wDesp(7)*errorNADH7; %wDesp(7)*errorNADH7_2;
                                wDesp(8)*errorNADH8; %wDesp(8)*errorNADH8_2;
                                wDesp(9)*errorNADH9; %wDesp(9)*errorNADH9_2;
                                wDesp(10)*errorNADH10; %wDesp(10)*errorNADH10_2;
                                wDesp(11)*errorNADH11; %wDesp(11)*errorNADH11_2;
                                wDesp(12)*errorNADH12]; %wDesp(12)*errorNADH12_2];
                            errorReg = lambda * x_temp';
                    otherwise
                        disp('No specific cost function has been appointed');
                end
                error = [
                    wD * errorNADH;
                    wH * errorHaldane;
                    wL * errorReg;
                    ];        
        %         disp(lambda); 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            otherwise
                disp('No enzyme has been selected in the cost function file');

        end

    end
  
    setup.Keq_PYK = temp_Keq_PYK;
    setup.Keq_LDH = temp_Keq_LDH;
end

save('new_tempRes_figure_Keq_pyk.mat');

% % % % %%
% % % % load('tempRes_figure_Keq_pyk.mat');
% % % % legNames = cell(length(tempResult),1);
% % % % c = cool(length(tempResult));
% % % % c(8,:) = [0 0 0];
% % % % figure(104)
% % % % for i = 1:length(tempResult)
% % % %     legNames{i} = mat2str(keqvalsTested(i));
% % % %     plot(tempResult{i}.t,tempResult{i}.y(:,2),'.-','color',c(i,:))
% % % %     hold on
% % % % end
% % % % % colormap('cool')
% % % % lgd = legend(legNames,'location','EastOutside','Orientation','vertical');
% % % % title({'The effect of equilibrium onstant just becomes apparent';'when values reach 1E2 or below'})
% % % % ylabel('NADH [mM] @ pH7.90')
% % % % xlabel('time [s]')
% % % % title(lgd,'K_{Eq} values tested')