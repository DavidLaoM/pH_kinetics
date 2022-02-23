% % MS_FIGURE_KEQ_PH_DEPENDENCY.M
% 
temp_setup = setup;
if setup.fast_option == 0 
    % %% (1/4) recall ALD
    for recall_ALD = 1
%         clear, close all
%         set_paths_pHstudy;
%         dbstop if error
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


        % %% (0.1) Experimental Vmax determination
        setup.plotOutput = 0;
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
        % % % % plength = 13; % Kms (3) + Vms (1) * numpH (8) + vm linking (2) + (Keq is fixed to experimental data)
        plength = 11; % Kms (3) + Vms (1) * numpH (8) + (Keq is fixed to experimental data)
        x_temp = zeros(1,plength);
        ub = 3*ones(1,plength);
        lb = -3*ones(1,plength);
        options = optimset('Display','iter');

        % %% (2.1) Parameter estimation with regularization
        % parameter estimation
        lambdalist = 1;
        parameterEstimation_lambdalist;
        selLambdaPos = 1;
        % 
        xres_Keq_pH_dependent = array_xres{selLambdaPos};

        % %% (2) Second run: keq pH independent (constant
        for i = 1:length(setup.Keq_FBA)
            setup.Keq_FBA(i) = setup.Keq_FBA(6);
            setup.Keq_TPI(i) = setup.Keq_TPI(6);
            setup.Keq_GPD(i) = setup.Keq_GPD(6);
        end
    %     %%
        parameterEstimation_lambdalist;
        % 
        xres_Keq_pH_independent = array_xres{selLambdaPos};

        % %% changing units
        vm_keq_pH_dependent = zeros(numpHtested,1);
        vm_keq_pH_independent = zeros(numpHtested,1);
        for i = 1:numpHtested
            vm_keq_pH_dependent(i) = data.Vmax(i,4) * 10.^xres_Keq_pH_dependent(i+3);
            vm_keq_pH_independent(i) = data.Vmax(i,4) * 10.^xres_Keq_pH_independent(i+3);
        end
        vm_keq_pH_dependent_uChange = vm_keq_pH_dependent .* 60 .* 60 ./ setup.concProtein;
        vm_keq_pH_independent_uChange = vm_keq_pH_independent .* 60 .* 60 ./ setup.concProtein;

        % % %%
        % figure
        % plot(pHarray, vm_keq_pH_dependent_uChange)
        % hold on
        % plot(pHarray, vm_keq_pH_independent_uChange)
        % legend('pH-dependent','pH-independent')
        % to save
        % % % % save('updated_keqconst_ald2.mat', 'updated_keqconst_ald')
    end
    updated_keqconst_ald.vm_keq_pH_dependent_uChange = vm_keq_pH_dependent_uChange;
    updated_keqconst_ald.vm_keq_pH_independent_uChange = vm_keq_pH_independent_uChange;
    updated_keqconst_ald.pH = pH;

    % %% (2/4) recall PDC
    for tempRecall = 1
%         clear, close all
%         set_paths_pHstudy;
%         dbstop if error
        for step0 = 1
            % select specific case and recall data
            setup.caseStudyALD = 0;
            setup.caseStudyENO = 0;
            setup.caseStudyGAPDH = 0;
            setup.caseStudyGAPDHr = 0;
            setup.caseStudyHXK = 0;
            setup.caseStudyPDC = 1;
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
            import_pdc = expData.pdc;

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

            pHarray = unique(import_pdc.treatedData.pH_corrected);
            for i = 1:numpHtested
                pHval = pHarray(i);
                tempID = find(import_pdc.treatedData.pH_corrected==pHval);
                pHTemp(:,i) = import_pdc.treatedData.pH_corrected(tempID);
                DFTemp(:,i) = import_pdc.treatedData.dilution_corrected(tempID);
                for j = 1:4
                    abs_meanTemp{j,i} = import_pdc.treatedData.absorbance_mean{tempID(j)};
                    abs_stdTemp{j,i} = import_pdc.treatedData.absorbance_std{tempID(j)};
                    conc_meanTemp{j,i} = import_pdc.treatedData.concentration_mean{tempID(j)};
                    conc_stdTemp{j,i} = import_pdc.treatedData.concentration_std{tempID(j)};
                    timeTemp{j,i} = import_pdc.treatedData.time{tempID(j)};
                    RRsTemp{j,i} = import_pdc.treatedData.reaction_rate{tempID(j)};
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
            temp1 = import_pdc.rawData.absorbance_corrected{4,4};
            temp2 = import_pdc.rawData.absorbance_corrected{5,4};
            temp3 = import_pdc.rawData.absorbance_corrected{6,4};
            data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
            data.raw.time = import_pdc.rawData.time{1};

                % (1) Correct for minimum value
                % (2) Bring the minimum to zero (apply to all)
                % (3) In principle, use the 3 first dilution rates
                % (4) Watch out with the dilution factors (first 2 cases are
                % reversed)

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            % Adjusted to PDC

            % relocating well DF2,#7,8
            dps = length(NADH{1,1});
            for i = 3
                for j = 7:8
                        for k = 1:dps
                            NADH{j,i}(k) = NADH{j,i}(k) - (NADH{j,3}(end)-NADH{j,4}(end));
                        end
                end
            end

            % Directly changing the concentration here, since the extinction
            % coefficient did not change.
            endPoint = zeros(numpHtested,DFs);
            % locate the minimum
            for i = 1:DFs
                for j = 1:numpHtested
                    endPoint(j,i) = min(NADH{j,i});
        %             endPoint(j,i) = min([NADH{j,1}' NADH{j,2}' NADH{j,3}' NADH{j,4}']);
                end
            end
            % bringing the minimum to zero
            for i = 1:DFs
                for j = 1:numpHtested
                    if((i==4)||(i==3)||((i==2)&&((j>=2)&&(j<=7))))
                        for k = 1:dps
                            NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,i);
                        end
                    else
                        for k = 1:dps
                            NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,4);
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

            data.conc_mean = NADH;
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

            pHvals = unique(import_pdc.treatedData.pH_corrected);
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

        % % % %     figure
        % % % %     plot(pHvals, Vmax(:,4),'.-')
        % % % %     title('Starting estimate: Vmax [mM s-1] vs pH')    
        end


        % %% (0.1) Experimental Vmax determination
        setup.plotOutput = 0;
        % correcting DF, specific for pdc case
        data.DF(7:8,3) = [4;4]; 
        DF(7:8,3) = [4;4];
        % %% (0.1) Calculation of rates: moving window
            % intial things that could be in the setup
            minwindow = 6; % minimum size of the window
            limRates = [0 2E-3]; %Ylims plot vmaxs
            limR2 = [0 1]; %Ylims plot R2
            limcConc = [0 0.15];  %Ylims plot conc
        % select start point (this needs manual selection deppending on previous plots)
        dp_start = ones(size(data.conc_mean));
        dp_start(:,1) = 10 * ones(size(dp_start(:,1)));
        dp_start(:,2) = 10 * ones(size(dp_start(:,2)));
        dp_start(:,3) = 5 * ones(size(dp_start(:,3)));
        dp_start(:,4) = 2 * ones(size(dp_start(:,4)));
        % blank cell total length
        total_len = zeros(size(dp_start));
        % DFs considered
        DFarray = [1/8 1/4 1/2 1/1];
        % idxs2consider
        idxs2consider = [0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 0 1 1;
                        0 0 1 1;
                        0 0 1 1;
                        0 0 1 1;
                        0 0 1 1;
                        0 0 1 1];

        % Experimental rates determination and plotting
        expRatesDetermination;

        % %% (1.1) Simple parameter fit. Parameter estimation
        setup.ode = 'vanHeerden2014';
        setup.sourceVm = 'experimentalSlopes';
        setup.ode_pH = 'on';

        setup.plotResults = 0;
        setup.plotEachSimCF = 0;
        setup.simAllProfiles = 0;
        setup.plotEachSim = 1;

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
        % % % % ub = 3*ones(1,plength);
        % % % % lb = -3*ones(1,plength);
        ub = 1*ones(1,plength);
        lb = -1*ones(1,plength);
        options = optimset('Display','iter');
        % %%

        % %% (2.1) Parameter estimation with regularization
        % parameter estimation
        lambdalist = 1;
        parameterEstimation_lambdalist;

        % %% (2.2) Regularization. Results Visualization
        selLambdaPos = 1;%15;%1;%17;%14;
        regularizationVisualization;

        %
        xres_Keq_pH_dependent = array_xres{selLambdaPos};

        %
        for i = 1:length(setup.Keq_ADH)
            setup.Keq_ADH(i) = setup.Keq_ADH(6);
        end

        % %% parameter estimation pH independent
        tic
        [xres,resnorm,residual,~,~,~,Jacobian] = lsqnonlin(optfun,x_temp,lb,ub,options,data,setup);
        t = toc;
        % 
        xres_Keq_pH_independent = xres;

        % parameter values
        vm_keq_pH_dependent = zeros(numpHtested,1);
        vm_keq_pH_independent = zeros(numpHtested,1);
        for i = 1:numpHtested
            vm_keq_pH_dependent(i) = data.Vmax(i,4) * 10.^xres_Keq_pH_dependent(i+2);
            vm_keq_pH_independent(i) = data.Vmax(i,4) * 10.^xres_Keq_pH_independent(i+2);
        end
        vm_keq_pH_dependent_uChange = vm_keq_pH_dependent .* 60 .* 60 ./ setup.concProtein;
        vm_keq_pH_independent_uChange = vm_keq_pH_independent .* 60 .* 60 ./ setup.concProtein;

        % %%
        % figure
        % plot(pHarray, vm_keq_pH_dependent_uChange)
        % hold on
        % plot(pHarray, vm_keq_pH_independent_uChange)
        % legend('pH-dependent','pH-independent')
        % to save
        % % % % save('updated_keqconst_pdc2.mat', 'updated_keqconst_pdc')
    end
    updated_keqconst_pdc.vm_keq_pH_dependent_uChange = vm_keq_pH_dependent_uChange;
    updated_keqconst_pdc.vm_keq_pH_independent_uChange = vm_keq_pH_independent_uChange;
    updated_keqconst_pdc.pH = pH;   

    % %% (3/4) recall PFK
    for tempRecall = 1
%     clear, close all
%     set_paths_pHstudy;
%     dbstop if error
    for step0 = 1
            % select specific case and recall data
            setup.caseStudyALD = 0;
            setup.caseStudyENO = 0;
            setup.caseStudyGAPDH = 0;
            setup.caseStudyGAPDHr = 0;
            setup.caseStudyHXK = 0;
            setup.caseStudyPDC = 0;
            setup.caseStudyPFK = 1;
            setup.caseStudyPGI = 0;
            setup.caseStudyPGM = 0;
            setup.caseStudyPYK = 0;
            setup.caseStudyTPI = 0;
            setup.caseStudyENO = 0;
            selectSetup_pH;
            % added
            setup.saveOutput = 0;

            load('expData.mat','expData');
            import_pfk = expData.pfk;

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

            pHarray = unique(import_pfk.treatedData.pH_corrected);
            for i = 1:numpHtested
                pHval = pHarray(i);
                tempID = find(import_pfk.treatedData.pH_corrected==pHval);
                pHTemp(:,i) = import_pfk.treatedData.pH_corrected(tempID);
                DFTemp(:,i) = import_pfk.treatedData.dilution_corrected(tempID);
                for j = 1:4
                    abs_meanTemp{j,i} = import_pfk.treatedData.absorbance_mean{tempID(j)};
                    abs_stdTemp{j,i} = import_pfk.treatedData.absorbance_std{tempID(j)};
                    conc_meanTemp{j,i} = import_pfk.treatedData.concentration_mean{tempID(j)};
                    conc_stdTemp{j,i} = import_pfk.treatedData.concentration_std{tempID(j)};
                    timeTemp{j,i} = import_pfk.treatedData.time{tempID(j)};
                    RRsTemp{j,i} = import_pfk.treatedData.reaction_rate{tempID(j)};
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
            temp1 = import_pfk.rawData.absorbance_corrected{4,4};
            temp2 = import_pfk.rawData.absorbance_corrected{5,4};
            temp3 = import_pfk.rawData.absorbance_corrected{6,4};
            data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
            data.raw.time = import_pfk.rawData.time{1};

                % (1) Correct for minimum value
                % (2) Bring the minimum to zero (apply to all)
                % (3) In principle, use the 3 first dilution rates
                % (4) Watch out with the dilution factors (first 2 cases are
                % reversed)

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            % Adjusted to PGM
            % Directly changing the concentration here, since the extinction
            % coefficient did not change.
            dpsarray = zeros(numpHtested,DFs);
            for i = 1:DFs
                for j = 1:numpHtested
                    dpsarray(j,i) = length(NADH{j,i});
                end
            end
            endPoint = zeros(numpHtested,DFs);
            endLocation = zeros(numpHtested,DFs);
            % locate the minimum
            for i = 1:DFs
                for j = 1:numpHtested
                    [endPoint(j,i),endLocation(j,i)] = min(NADH{j,i});
                end
            end
            % bringing the minimum to zero
            for i = 1:DFs
                for j = 1:numpHtested
                    dps = dpsarray(j,i);
                    if ((i==3)&&(j>=3)&&(j<=9))
                        for k = 1:dps
                            NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,3);
                        end
                    else
                        for k = 1:dps
                            NADH{j,i}(k) = NADH{j,i}(k) - endPoint(j,4);
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
                    dps = dpsarray(j,i);
                    for k = (tempidx+1):dps
                        NADH{j,i}(k) = NADH{j,i}(tempidx);
                    end
                end
            end

            data.conc_mean = NADH;
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            takenTime = [700, 700, 550, 500, 450, 450, 350, 350, 350, 350, 350, 350];
            cuttingPoints = [141, 141, 111, 101, 91, 91, 71, 71, 71, 71, 71, 71];
            for i = 1:DFs
                for j = 1:numpHtested
                    startVal = cuttingPoints(j);
                    data.abs_mean{j,i} = data.abs_mean{j,i}(startVal:end);
                    data.abs_std{j,i} = data.abs_std{j,i}(startVal:end);
                    data.conc_mean{j,i} = data.conc_mean{j,i}(startVal:end);
                    data.conc_std{j,i} = data.conc_std{j,i}(startVal:end);
                    data.time{j,i} = data.time{j,i}(startVal:end) - takenTime(j);
                    data.RRs{j,i} = data.RRs{j,i}(startVal:end);
                    time{j,i} = time{j,i}(startVal:end) - takenTime(j);
                    NADH{j,i} = NADH{j,i}(startVal:end);
                end
            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

            pHvals = unique(import_pfk.treatedData.pH_corrected);
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
                if setup.caseStudyPFK == 1
                    ylim([0 0.12])
                    xlim([0 2000])
                end
            end
            suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
            suptitle(suptitleName);

        %     figure
        %     plot(pHvals, Vmax(:,4),'.-')
        %     title('Starting estimate: Vmax [mM s-1] vs pH')    
    end


    %     %% (0.1) Experimental Vmax determination
        setup.plotOutput = 0;
        % %% (0.1) Calculation of rates: moving window
            % intial things that could be in the setup
            minwindow = 70; % minimum size of the window
            limRates = [0 1.5E-4]; %Ylims plot vmaxs
            limR2 = [0 1]; %Ylims plot R2
            limcConc = [0 0.15];  %Ylims plot conc
        % select start point (this needs manual selection deppending on previous plots)
        dp_start = 6 * ones(size(data.conc_mean));
        % blank cell total length
        total_len = zeros(size(dp_start));
        % DFs considered
        DFarray_pre = 1./data.DF;
        DFarray = DFarray_pre(1,:);
        % idxs2consider
        idxs2consider = [0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1;
                        0 1 1 1];

        % Experimental rates determination and plotting
        expRatesDetermination;

%         % %% (0.3) saving
%         if setup.plotOutput == 1
%             save(['results/',setup.enzymeName,'/',setup.enzymeName, '_initial_variables.mat'],'Vmax_mw_opt_corr','idxs2consider','DF','pH');
%             set(1,'color','white'), savefig(1,['results/',setup.enzymeName,'/',setup.enzymeName, '_concentrations_basezero.fig']);
%             set(2,'color','white'), savefig(2,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_df.fig']);
%             set(3,'color','white'), savefig(3,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_movingWindow.fig']);
%             set(4,'color','white'), savefig(4,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_R2_vs_movingWindow.fig']);
%             set(5,'color','white'), savefig(5,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_iniPoint_vs_movingWindow.fig']);
%             set(6,'color','white'), savefig(6,['results/',setup.enzymeName,'/',setup.enzymeName, '_experimental_vmax_vs_pH.fig']);
%         end

    %     %% (1.1) Simple parameter fit. Parameter estimation
        setup.ode = 'vanHeerden2014';
        setup.sourceVm = 'experimentalSlopes';
        setup.ode_pH = 'on';

        setup.plotResults = 0;
        setup.plotEachSimCF = 0;
        setup.simAllProfiles = 0;
        setup.plotEachSim = 1;

        setup.numpHtested = numpHtested;
        setup.DFstudy = 1:4;
        setup.costfun = 3;

        setup.weightData = 1;
        setup.weightDataEsp = idxs2consider;
        % % % % idxs2consider2 = [0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1;
        % % % %                 0 0 0 1];
        % % % % setup.weightDataEsp = idxs2consider2;
        setup.weightHaldane = 0; % off for this case (Keq fixed)
        setup.selectedLambda = 0; % by now just testing

        % three stypes of topologies and study
        % % % % setup.problemStudy = 'onlyVmax';
        setup.problemStudy = 'fullKinetics_paramsPartFixed';
        % setup.problemStudy = 'fullKinetics_paramsAllFlexible';
        % prepare setup based of 'problemStudy'
        problemStudy = setup.problemStudy;
        optfun = @costfun_Kmfixed;
        switch problemStudy
            case 'onlyVmax'
                plength = 12; % only the vmax change with pH. Other params OFF
            case 'fullKinetics_paramsPartFixed'
                plength = 25; % only the vmax change with pH. Other params ON
            case 'fullKinetics_paramsAllFlexible'
                plength = 168; % All parameters change with pH. Other params ON
            otherwise
                disp('Warning: problem study (reaction kinetics) have not been selected');
        end
        x_temp = zeros(1,plength);
        ub = 3*ones(1,plength);
        lb = -3*ones(1,plength);
        %     ub(1:13) = zeros;
        %     lb(1:13) = zeros; 
        % options = optimset('Display','iter');
        options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4);
        % % plength = 25; % vms (12) + other params (12)
        % % x_temp = zeros(1,plength);
        % plength = 12; % vms (12) + other params off
        % % x_temp = -0.45*ones(1,plength);

        % 
        lambdalist = 0;
        parameterEstimation_lambdalist;

        % %% (2.2) Regularization. Results Visualization
        selLambdaPos = 1;%1;%25; %8;
    %     regularizationVisualization;

        % %% (3.1) Study on paarameter values: estimation with the lambda value
        % parameter values
        % xres_selected = xres; %lambdalist based in 'ones', lam=0.1, loc=5.
        % xres_selected = array_xres{15}; %lambdalist based in 'ones', lam=0.1, loc=5.
        xres_selected_pHdependent = array_xres{selLambdaPos};

        % %% (2.1) Parameter estimation pH independent
        for i = 1:length(setup.Keq_FBA)
            setup.Keq_FBA(i) = setup.Keq_FBA(6);
            setup.Keq_TPI(i) = setup.Keq_TPI(6);
            setup.Keq_GPD(i) = setup.Keq_GPD(6);
            setup.Keq_PFK(i) = setup.Keq_PFK(6);
        end
        parameterEstimation_lambdalist;
        xres_selected_pHindependent = array_xres{selLambdaPos};

        % 
        vm_pHdependent = zeros(numpHtested,1);
        vm_pHindependent = zeros(numpHtested,1);
        for i = 1:numpHtested
            vm_pHdependent(i) = data.Vmax(i,4) * 10.^xres_selected_pHdependent(i+13);
            vm_pHindependent(i) = data.Vmax(i,4) * 10.^xres_selected_pHindependent(i+13);
        end
        keq_fba = setup.Keq_FBA;
        keq_gpd = setup.Keq_GPD;
        keq_tpi = setup.Keq_TPI;
        keq_pfk = setup.Keq_PFK;
        vm_pHdependent_uChange = vm_pHdependent .* 60 .* 60 ./ setup.concProtein;
        vm_pHindependent_uChange = vm_pHindependent .* 60 .* 60 ./ setup.concProtein;

        % %%
        % figure
        % plot(pHarray, vm_pHdependent_uChange)
        % hold on
        % plot(pHarray, vm_pHindependent_uChange)
        % legend('pH-dependent','pH-independent')
        % 
        % %% to save
        % % % % save('updated_keqconst_pfk2.mat', 'updated_keqconst_pfk')

    end
    updated_keqconst_pfk.vm_keq_pH_dependent_uChange = vm_pHdependent_uChange;
    updated_keqconst_pfk.vm_keq_pH_independent_uChange = vm_pHindependent_uChange;
    updated_keqconst_pfk.pH = pH;

    % %% (4/4) recall PYK
    for tempRecall = 1
%         clear, close all
%         set_paths_pHstudy;
%         dbstop if error
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

        % % % %     figure
        % % % %     plot(pHvals, Vmax(:,4),'.-')
        % % % %     title('Starting estimate: Vmax [mM s-1] vs pH')    
        end

    %     %% (0.1) Experimental Vmax determination
        setup.plotOutput = 0;
        % %% (0.1) Calculation of rates: moving window
            % intial things that could be in the setup
            minwindow = 7; % minimum size of the window
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

    %     %% (1.1) Simple parameter fit. Parameter estimation
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
        plength = 18; % others (6) + Vms (1) * numpH (12) + (Keq is fixed to experimental data)
        x_temp = zeros(1,plength);
        ub = 1*ones(1,plength);
        lb = -1*ones(1,plength);
        options = optimset('Display','iter');
        % %%

        % %% (2.1) Parameter estimation with regularization
        lambdalist = 0;
        parameterEstimation_lambdalist;
        % %% (2.2) Regularization. Results Visualization
        selLambdaPos = 1;%13;%17;%16;
        % regularizationVisualization;
        % %% (3.2) Study on paarameter values: recalculation
        xres_selected_pHdependent = array_xres{selLambdaPos};

        % %% pH-independent
        for i = 1:length(setup.Keq_PYK)
            setup.Keq_PYK(i) = setup.Keq_PYK(6);
            setup.Keq_LDH(i) = setup.Keq_LDH(6);
        end
        % 
        parameterEstimation_lambdalist;
        xres_selected_pHindependent = array_xres{selLambdaPos};

        % %% units
        vm_pHdependent = zeros(numpHtested,1);
        vm_pHindependent = zeros(numpHtested,1);
        for i = 1:numpHtested
            vm_pHdependent(i) = data.Vmax(i,4) * 10.^xres_selected_pHdependent(i+6);
            vm_pHindependent(i) = data.Vmax(i,4) * 10.^xres_selected_pHindependent(i+6);
        end

        % %% (3.3) Study on paarameter values: plotting
        vm_pHdependent_uChange = vm_pHdependent .* 60 .* 60 ./ setup.concProtein;
        vm_pHindependent_uChange = vm_pHindependent .* 60 .* 60 ./ setup.concProtein;

        % %%
        % figure
        % plot(pHarray, vm_pHdependent_uChange)
        % hold on
        % plot(pHarray, vm_pHindependent_uChange)
        % legend('pH-dependent','pH-independent')
        % 
        % %% to save
        % % % % save('updated_keqconst_pyk2.mat', 'updated_keqconst_pyk')

    end
    updated_keqconst_pyk.vm_keq_pH_dependent_uChange = vm_pHdependent_uChange;
    updated_keqconst_pyk.vm_keq_pH_independent_uChange = vm_pHindependent_uChange;
    updated_keqconst_pyk.pH = pH;

elseif setup.fast_option == 1 % directly recall the workspace from here
    load('workspace_keqconst_ald.mat')
    load('workspace_keqconst_pdc.mat')
    load('workspace_keqconst_pfk.mat')
    load('workspace_keqconst_pyk.mat')
end

%% plots
% c_royalBlue = [65	105	225]/255; % royalblue
c_midnightblue = [25	25	112]/255; % midnightblue
% c_CCCCCC = [204	204	204]/255; % #CCCCCC
% c_E5E5E5 = [229 229 229]/255; % #E5E5E5
% c_0f1076 = [15	16	118]/255; % #0f1076
c_chocolate = [210	105	30]/255; % (#e59400 temp orange)

% fill with right colors
% check which ones different values
% check why? possibly older version

fig_h101 = figure(101);
% ALD
sp1 = subplot(2,2,1);
plot(updated_keqconst_ald.pH, updated_keqconst_ald.vm_keq_pH_dependent_uChange,...% '.-',...
    'o-', 'color', c_midnightblue, 'Linewidth', 1.2, 'MarkerSize', 5,...
    'MarkerFaceColor', c_midnightblue, 'LineWidth', 2)
%     'MarkerSize',10,'LineWidth',2,'color',c_midnightblue)
hold on
plot(updated_keqconst_ald.pH, updated_keqconst_ald.vm_keq_pH_independent_uChange,...% '.-',...
    'o-', 'color', c_chocolate, 'Linewidth', 1.2, 'MarkerSize', 3,...
    'MarkerFaceColor', c_chocolate, 'LineWidth', 2)
%     'MarkerSize',10,'LineWidth',2,'color',c_chocolate)
% legend('Keq pH-dependent', 'Keq pH-independent')
title('ALD')
ylim([0 3]), yticks([0 1 2])
xlim([6 8]), xticks([6.0 6.5 7.0 7.5 8.0])
xlabel('pH')
ylab = ylabel('Enzyme capacity (\mumol min^{-1} mg protein^{-1})');
hold off
% PDC
sp2 = subplot(2,2,2);
plot(updated_keqconst_pdc.pH, updated_keqconst_pdc.vm_keq_pH_dependent_uChange,...% '.-',...
    'o-', 'color', c_midnightblue, 'Linewidth', 1.2, 'MarkerSize', 5,...
    'MarkerFaceColor', c_midnightblue, 'LineWidth', 2)
%     'MarkerSize',10,'LineWidth',2,'color',c_midnightblue)
hold on
plot(updated_keqconst_pdc.pH, updated_keqconst_pdc.vm_keq_pH_independent_uChange,...% '.-',...
    'o-', 'color', c_chocolate, 'Linewidth', 1.2, 'MarkerSize', 3,...
    'MarkerFaceColor', c_chocolate, 'LineWidth', 2)
%     'MarkerSize',10,'LineWidth',2,'color',c_chocolate)
% legend('Keq pH-dependent', 'Keq pH-independent')
title('PDC')
ylim([0 3]), yticks([0 1 2])
xlim([6 8]), xticks([6.0 6.5 7.0 7.5 8.0])
xlabel('pH')
hold off
% PFK
sp3 = subplot(2,2,3);
plot(updated_keqconst_pfk.pH, updated_keqconst_pfk.vm_keq_pH_dependent_uChange,...% '.-',...
    'o-', 'color', c_midnightblue, 'Linewidth', 1.2, 'MarkerSize', 5,...
    'MarkerFaceColor', c_midnightblue, 'LineWidth', 2)
%     'MarkerSize',10,'LineWidth',2,'color',c_midnightblue)
hold on
plot(updated_keqconst_pfk.pH, updated_keqconst_pfk.vm_keq_pH_independent_uChange,...% '.-',...
    'o-', 'color', c_chocolate, 'Linewidth', 1.2, 'MarkerSize', 3,...
    'MarkerFaceColor', c_chocolate, 'LineWidth', 2)
%     'MarkerSize',10,'LineWidth',2,'color',c_chocolate)
% legend('Keq pH-dependent', 'Keq pH-independent')
title('PFK')
ylim([0 0.6]), yticks([0 0.2 0.4])
xlim([6 8]), xticks([6.0 6.5 7.0 7.5 8.0])
xlabel('pH')
hold off
% PYK
sp4 = subplot(2,2,4);
% plot(updated_keqconst_pyk.pH, updated_keqconst_pyk.vm_keq_pH_dependent_uChange, '.-',...
%     'MarkerSize',10,'LineWidth',2,'color',c_midnightblue)
% hold on
% plot(updated_keqconst_pyk.pH, updated_keqconst_pyk.vm_keq_pH_independent_uChange, '.-',...
%     'MarkerSize',10,'LineWidth',2,'color',c_chocolate)
plot(updated_keqconst_pyk.pH, updated_keqconst_pyk.vm_keq_pH_dependent_uChange,...
    'o-', 'color', c_midnightblue, 'Linewidth', 1.2, 'MarkerSize', 5,...
    'MarkerFaceColor', c_midnightblue, 'LineWidth', 2)
%     '.-',...
%     'MarkerSize',10,'LineWidth',2,'color',c_midnightblue)
hold on
plot(updated_keqconst_pyk.pH, updated_keqconst_pyk.vm_keq_pH_independent_uChange,...
    'o-', 'color', c_chocolate, 'Linewidth', 1.2, 'MarkerSize', 3,...
    'MarkerFaceColor', c_chocolate, 'LineWidth', 2)

% legend('Keq pH-dependent', 'Keq pH-independent')
title('PYK')
ylim([0 8]), yticks([0 2 4 6 8])
xlim([6 8]), xticks([6.0 6.5 7.0 7.5 8.0])
xlabel('pH')
hold off

% sp size
sp1.FontSize = 12;
sp2.FontSize = 12;
sp3.FontSize = 12;
sp4.FontSize = 12;
%
% delete(hL1)
hL1 = legend(sp1.Children([1 5]),'pH-independent K_{eq}','pH-dependent K_{eq}');
% hL1 = legend(sp1.Children([1 2 3 4 5 6 7 8]), 'pH-dependent Keq', 'pH-independent Keq');
hL1.Orientation = 'horizontal';
hL1.Box = 'off';
hL1.FontSize = 11;
hL1.Position = [0.3    0.03    0.3960    0.0340];
% ylab.Position = ylab.Position + [5.8491 1.5000 -1.0000];
ylab.Position = [5.6 -0.75 -1.0000];
% set position
sp1.Position = [0.1300    0.5838+0.025    0.3347    0.3412-0.025];
sp2.Position = [0.5703    0.5838+0.025    0.3347    0.3412-0.025];
sp3.Position = [0.1300    0.1100+0.05    0.3347    0.3412-0.025];
sp4.Position = [0.5703    0.1100+0.05    0.3347    0.3412-0.025];  
% 
set(101, 'Position', [100 100 750 750],'color','w')


%% needed stop not to get truncated... #matlabUselessSecrets
% save
setup = temp_setup;
if setup.saveOutput == 1
    savefig(101, '1appendix_keq_otherEnzymes')
    % specs printing (method 3)
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters 1appendix_keq_otherEnzymes
    print -dpng -painters 1appendix_keq_otherEnzymes
end



% %%
% plot(output_hxk_changing_pH.pHarray, output_hxk_changing_pH.vm_uChange, ...
%     'o-', 'color', c_midnightblue, 'Linewidth', 1.2, 'MarkerSize', 5,...
%     'MarkerFaceColor', c_midnightblue, 'LineWidth', 2)
% hold on
% plot(output_hxk_changing_pH.pHarray, output_hxk_constant_pH.vm_uChange, ...
%     'o-', 'color', c_chocolate, 'Linewidth', 1.2, 'MarkerSize', 5,...
%     'MarkerFaceColor', c_chocolate, 'LineWidth', 2)






