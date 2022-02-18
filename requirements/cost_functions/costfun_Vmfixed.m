function [error] = costfun_Vmfixed(x_temp,data,setup)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%     x_temp(1) = Vmf
%     x_temp(2) = Vmr
%     x_temp([3:6]) = {Kgap, Kbpg, Knad, Knadh}, pH#1
%     x_temp([7:10]) = {Kgap, Kbpg, Knad, Knadh}, pH#2
%     x_temp([11:14]) = {Kgap, Kbpg, Knad, Knadh}, pH#3
%     x_temp([15:18]) = {Kgap, Kbpg, Knad, Knadh}, pH#4
%     x_temp([19:22]) = {Kgap, Kbpg, Knad, Knadh}, pH#5
%     x_temp([23:26]) = {Kgap, Kbpg, Knad, Knadh}, pH#6
%     x_temp([27:30]) = {Kgap, Kbpg, Knad, Knadh}, pH#7
%     x_temp([31:34]) = {Kgap, Kbpg, Knad, Knadh}, pH#8
%     x_temp([35:38]) = {Kgap, Kbpg, Knad, Knadh}, pH#9
%     x_temp([39:42]) = {Kgap, Kbpg, Knad, Knadh}, pH#10

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
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(3);
                    xassay(3) = x_temp(4);
                    xassay(4) = x_temp(5);
                    xassay(5) = x_temp(6);
                    xassay(6) = x_temp(2);
                case 2
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(7);
                    xassay(3) = x_temp(8);
                    xassay(4) = x_temp(9);
                    xassay(5) = x_temp(10);
                    xassay(6) = x_temp(2);
                case 3
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(11);
                    xassay(3) = x_temp(12);
                    xassay(4) = x_temp(13);
                    xassay(5) = x_temp(14);
                    xassay(6) = x_temp(2);
                case 4
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(15);
                    xassay(3) = x_temp(16);
                    xassay(4) = x_temp(17);
                    xassay(5) = x_temp(18);
                    xassay(6) = x_temp(2);
                case 5
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(19);
                    xassay(3) = x_temp(20);
                    xassay(4) = x_temp(21);
                    xassay(5) = x_temp(22);
                    xassay(6) = x_temp(2);
                case 6
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(23);
                    xassay(3) = x_temp(24);
                    xassay(4) = x_temp(25);
                    xassay(5) = x_temp(26);
                    xassay(6) = x_temp(2);
                case 7
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(27);
                    xassay(3) = x_temp(28);
                    xassay(4) = x_temp(29);
                    xassay(5) = x_temp(30);
                    xassay(6) = x_temp(2);
                case 8
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(31);
                    xassay(3) = x_temp(32);
                    xassay(4) = x_temp(33);
                    xassay(5) = x_temp(34);
                    xassay(6) = x_temp(2);
                case 9
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(35);
                    xassay(3) = x_temp(36);
                    xassay(4) = x_temp(37);
                    xassay(5) = x_temp(38);
                    xassay(6) = x_temp(2);
                case 10
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(39);
                    xassay(3) = x_temp(40);
                    xassay(4) = x_temp(41);
                    xassay(5) = x_temp(42);
                    xassay(6) = x_temp(2);
                otherwise
                    disp('Something went wrong is selecting the pH value');
            end
            end
            
            fullSimRes = cell(4,1);
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
                
                simResult.v = vObs;
                fullSimRes{i} = simResult;
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
                        simRes = fullSimRes{i};
                        plot(simRes.t, simRes.y(:,8),'-')
                        hold on
                        plot(data.time{j,i}, data.conc_mean{j,i},'k+')
                        hold on
                    end
                    title('NADH')
                    subplot(1,2,2)
                    for i = DFstudy
                        simRes = fullSimRes{i};
                        plot(simRes.t, simRes.v, '-')
                        hold on
                        plot(data.tempTime{i}, expGAPDHr{i},'k+')
                        hold on
                    end
                    title('v_{apparent.GAPDHr}')
                    suptitle(erase(sprintf('pH = %d',setup.pH_vals(j)),"0000e+00"));
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
                                % x_temp([3:6]) = {Kgap, Kbpg, Knad, Knadh}, pH#1
                                ks1_1 = 10 .^ x_temp(3) .* 2.48; % mM %k_gap
                                ks2_1 = 10 .^ x_temp(5) .* 2.92; %mM %k_nad
                                kp1_1 = 10 .^ x_temp(4) .* 1.18; % mM %k_bpg
                                kp2_1 = 10 .^ x_temp(6) .* 0.022; % mM %k_nadh
                                % x_temp([7:10]) = {Kgap, Kbpg, Knad, Knadh}, pH#2                                
                                ks1_2 = 10 .^ x_temp(7) .* 2.48; % mM %k_gap
                                ks2_2 = 10 .^ x_temp(9) .* 2.92; %mM %k_nad
                                kp1_2 = 10 .^ x_temp(8) .* 1.18; % mM %k_bpg
                                kp2_2 = 10 .^ x_temp(10) .* 0.022; % mM %k_nadh
                                % x_temp([11:14]) = {Kgap, Kbpg, Knad, Knadh}, pH#3
                                ks1_3 = 10 .^ x_temp(11) .* 2.48; % mM %k_gap
                                ks2_3 = 10 .^ x_temp(13) .* 2.92; %mM %k_nad
                                kp1_3 = 10 .^ x_temp(12) .* 1.18; % mM %k_bpg
                                kp2_3 = 10 .^ x_temp(14) .* 0.022; % mM %k_nadh
                                % x_temp([15:18]) = {Kgap, Kbpg, Knad, Knadh}, pH#4                                
                                ks1_4 = 10 .^ x_temp(15) .* 2.48; % mM %k_gap
                                ks2_4 = 10 .^ x_temp(17) .* 2.92; %mM %k_nad
                                kp1_4 = 10 .^ x_temp(16) .* 1.18; % mM %k_bpg
                                kp2_4 = 10 .^ x_temp(18) .* 0.022; % mM %k_nadh
                                % x_temp([19:22]) = {Kgap, Kbpg, Knad, Knadh}, pH#5                                
                                ks1_5 = 10 .^ x_temp(19) .* 2.48; % mM %k_gap
                                ks2_5 = 10 .^ x_temp(21) .* 2.92; %mM %k_nad
                                kp1_5 = 10 .^ x_temp(20) .* 1.18; % mM %k_bpg
                                kp2_5 = 10 .^ x_temp(22) .* 0.022; % mM %k_nadh
                                % x_temp([23:26]) = {Kgap, Kbpg, Knad, Knadh}, pH#6                                
                                ks1_6 = 10 .^ x_temp(23) .* 2.48; % mM %k_gap
                                ks2_6 = 10 .^ x_temp(25) .* 2.92; %mM %k_nad
                                kp1_6 = 10 .^ x_temp(24) .* 1.18; % mM %k_bpg
                                kp2_6 = 10 .^ x_temp(26) .* 0.022; % mM %k_nadh
                                % x_temp([27:30]) = {Kgap, Kbpg, Knad, Knadh}, pH#7                                
                                ks1_7 = 10 .^ x_temp(27) .* 2.48; % mM %k_gap
                                ks2_7 = 10 .^ x_temp(29) .* 2.92; %mM %k_nad
                                kp1_7 = 10 .^ x_temp(28) .* 1.18; % mM %k_bpg
                                kp2_7 = 10 .^ x_temp(30) .* 0.022; % mM %k_nadh
                                % x_temp([31:34]) = {Kgap, Kbpg, Knad, Knadh}, pH#8                                
                                ks1_8 = 10 .^ x_temp(31) .* 2.48; % mM %k_gap
                                ks2_8 = 10 .^ x_temp(33) .* 2.92; %mM %k_nad
                                kp1_8 = 10 .^ x_temp(32) .* 1.18; % mM %k_bpg
                                kp2_8 = 10 .^ x_temp(34) .* 0.022; % mM %k_nadh
                                % x_temp([35:38]) = {Kgap, Kbpg, Knad, Knadh}, pH#9                                
                                ks1_9 = 10 .^ x_temp(35) .* 2.48; % mM %k_gap
                                ks2_9 = 10 .^ x_temp(37) .* 2.92; %mM %k_nad
                                kp1_9 = 10 .^ x_temp(36) .* 1.18; % mM %k_bpg
                                kp2_9 = 10 .^ x_temp(38) .* 0.022; % mM %k_nadh
                                % x_temp([39:42]) = {Kgap, Kbpg, Knad, Knadh}, pH#10                                
                                ks1_10 = 10 .^ x_temp(39) .* 2.48; % mM %k_gap
                                ks2_10 = 10 .^ x_temp(41) .* 2.92; %mM %k_nad
                                kp1_10 = 10 .^ x_temp(40) .* 1.18; % mM %k_bpg
                                kp2_10 = 10 .^ x_temp(42) .* 0.022; % mM %k_nadh
                            otherwise
                                disp('No source for kms has been selected');
                        end
                        vmf = 10.^x_temp(1).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                        vmr = 10.^x_temp(2).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
                        switch ode_pH
                            case 'on'
                                H_effect = zeros(1:10);
                                for j = 1:numpH
                                    H_effect(j) = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                                end
                                errorHaldane1 = (Keq(1) - (vmf * kp1_1 * kp2_1 * H_effect(1)) / (vmr * ks1_1 * ks2_1) );
                                errorHaldane2 = (Keq(2) - (vmf * kp1_2 * kp2_2 * H_effect(2)) / (vmr * ks1_2 * ks2_2) );
                                errorHaldane3 = (Keq(3) - (vmf * kp1_3 * kp2_3 * H_effect(3)) / (vmr * ks1_3 * ks2_3) );
                                errorHaldane4 = (Keq(4) - (vmf * kp1_4 * kp2_4 * H_effect(4)) / (vmr * ks1_4 * ks2_4) );
                                errorHaldane5 = (Keq(5) - (vmf * kp1_5 * kp2_5 * H_effect(5)) / (vmr * ks1_5 * ks2_5) );
                                errorHaldane6 = (Keq(6) - (vmf * kp1_6 * kp2_6 * H_effect(6)) / (vmr * ks1_6 * ks2_6) );
                                errorHaldane7 = (Keq(7) - (vmf * kp1_7 * kp2_7 * H_effect(7)) / (vmr * ks1_7 * ks2_7) );
                                errorHaldane8 = (Keq(8) - (vmf * kp1_8 * kp2_8 * H_effect(8)) / (vmr * ks1_8 * ks2_8) );
                                errorHaldane9 = (Keq(9) - (vmf * kp1_9 * kp2_9 * H_effect(9)) / (vmr * ks1_9 * ks2_9) );
                                errorHaldane10 = (Keq(10) - (vmf * kp1_10 * kp2_10 * H_effect(10)) / (vmr * ks1_10 * ks2_10) );
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

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    case 'eno'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = setup.fixedValues;
% % % %             data.Vmaxs = data.Vmax(j,:);
% % % % % % % %             data.Vmaxs = data.Vmax(setup.startVm,:); % In case of Vmfixed, we have only one value for Vm here (took the 6th data point)
% % % % % % % %             data.Vmaxs = data.Vmax(setup.startVm,:).*setup.factMultiply; % In case of Vmfixed, we have only one value for Vm here (took the 6th data point)
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
                    xassay(1) = x_temp(2);
                    xassay(2) = x_temp(14);
                    xassay(3) = x_temp(1);
                case 2
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(3);
                    xassay(2) = x_temp(15);
                    xassay(3) = x_temp(1);
                case 3
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(4);
                    xassay(2) = x_temp(16);
                    xassay(3) = x_temp(1);
                case 4
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(5);
                    xassay(2) = x_temp(17);
                    xassay(3) = x_temp(1);
                case 5
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(6);
                    xassay(2) = x_temp(18);
                    xassay(3) = x_temp(1);
                case 6
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(7);
                    xassay(2) = x_temp(19);
                    xassay(3) = x_temp(1);
                case 7
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(8);
                    xassay(2) = x_temp(20);
                    xassay(3) = x_temp(1);
                case 8
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(9);
                    xassay(2) = x_temp(21);
                    xassay(3) = x_temp(1);
                case 9
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(10);
                    xassay(2) = x_temp(22);
                    xassay(3) = x_temp(1);
                case 10
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(11);
                    xassay(2) = x_temp(23);
                    xassay(3) = x_temp(1);
                case 11
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(12);
                    xassay(2) = x_temp(24);
                    xassay(3) = x_temp(1);
                case 12
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(13);
                    xassay(2) = x_temp(25);
                    xassay(3) = x_temp(1);
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
%                 % 2020-08-10: intercalcate for 'prettier' plot
%                 selectedVal = 1:4:121;                
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
%                         h101 = figure(101);
                        num1 = setup.startVm + 100;
                        h101 = figure(num1);
                    end
%                     set(101);
                    set(0,'CurrentFigure', h101)
                    subplot(3,4,j)
                    for i = DFstudy
                        plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2)
                        hold on
                        plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
%                         plot(data.tempTime{i}(selectedVal), expPEP{i,j}(selectedVal),'k.','MarkerSize',4)
                        ylim([0 1.2])
                    end
                    if j == numpH
                        suptitle('PEP concentration [mM] vs assay time [s]')
                    end
%                 end
%             end
%         end
                    if j == 1
%                         h102 = figure(102);
%                         h101 = figure(101);
                        num2 = setup.startVm + 200;
                        h102 = figure(num2);
                    end
%                     set(102)
                    set(0,'CurrentFigure', h102)
                    subplot(3,4,j)
                    for i = DFstudy
                        plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2)
                        hold on
                        plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
%                         plot(data.tempTime{i}(selectedVal), expENO{i,j}(selectedVal),'k.','MarkerSize',4)
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
                    
            case 3 % selected DFs
                    errorHaldane = 0;
                    % DF1
                    errorPEP1_df1 = simPEP{4,1} - expPEP{4,1};
                    errorPEP2_df1 = simPEP{4,2} - expPEP{4,2};
                    errorPEP3_df1 = simPEP{4,3} - expPEP{4,3};
                    errorPEP4_df1 = simPEP{4,4} - expPEP{4,4};
                    errorPEP5_df1 = simPEP{4,5} - expPEP{4,5};
                    errorPEP6_df1 = simPEP{4,6} - expPEP{4,6};
                    errorPEP7_df1 = simPEP{4,7} - expPEP{4,7};
                    errorPEP8_df1 = simPEP{4,8} - expPEP{4,8};
                    errorPEP9_df1 = simPEP{4,9} - expPEP{4,9};
                    errorPEP10_df1 = simPEP{4,10} - expPEP{4,10};
                    errorPEP11_df1 = simPEP{4,11} - expPEP{4,11};
                    errorPEP12_df1 = simPEP{4,12} - expPEP{4,12};
                        % DF2
                        errorPEP1_df2 = simPEP{3,1} - expPEP{3,1};
                        errorPEP2_df2 = simPEP{3,2} - expPEP{3,2};
                        errorPEP3_df2 = simPEP{3,3} - expPEP{3,3};
                        errorPEP4_df2 = simPEP{3,4} - expPEP{3,4};
                        errorPEP5_df2 = simPEP{3,5} - expPEP{3,5};
                        errorPEP6_df2 = simPEP{3,6} - expPEP{3,6};
                        errorPEP7_df2 = simPEP{3,7} - expPEP{3,7};
                        errorPEP8_df2 = simPEP{3,8} - expPEP{3,8};
                        errorPEP9_df2 = simPEP{3,9} - expPEP{3,9};
                        errorPEP10_df2 = simPEP{3,10} - expPEP{3,10};
                        errorPEP11_df2 = simPEP{3,11} - expPEP{3,11};
                        errorPEP12_df2 = simPEP{3,12} - expPEP{3,12};
                            % DF4
                            errorPEP1_df4 = simPEP{2,1} - expPEP{2,1};
                            errorPEP2_df4 = simPEP{2,2} - expPEP{2,2};
                            errorPEP3_df4 = simPEP{2,3} - expPEP{2,3};
                            errorPEP4_df4 = simPEP{2,4} - expPEP{2,4};
                            errorPEP5_df4 = simPEP{2,5} - expPEP{2,5};
                            errorPEP6_df4 = simPEP{2,6} - expPEP{2,6};
                            errorPEP7_df4 = simPEP{2,7} - expPEP{2,7};
                            errorPEP8_df4 = simPEP{2,8} - expPEP{2,8};
                            errorPEP9_df4 = simPEP{2,9} - expPEP{2,9};
                            errorPEP10_df4 = simPEP{2,10} - expPEP{2,10};
                            errorPEP11_df4 = simPEP{2,11} - expPEP{2,11};
                            errorPEP12_df4 = simPEP{2,12} - expPEP{2,12};
                                % DF8
                                errorPEP1_df8 = simPEP{1,1} - expPEP{1,1};
                                errorPEP2_df8 = simPEP{1,2} - expPEP{1,2};
                                errorPEP3_df8 = simPEP{1,3} - expPEP{1,3};
                                errorPEP4_df8 = simPEP{1,4} - expPEP{1,4};
                                errorPEP5_df8 = simPEP{1,5} - expPEP{1,5};
                                errorPEP6_df8 = simPEP{1,6} - expPEP{1,6};
                                errorPEP7_df8 = simPEP{1,7} - expPEP{1,7};
                                errorPEP8_df8 = simPEP{1,8} - expPEP{1,8};
                                errorPEP9_df8 = simPEP{1,9} - expPEP{1,9};
                                errorPEP10_df8 = simPEP{1,10} - expPEP{1,10};
                                errorPEP11_df8 = simPEP{1,11} - expPEP{1,11};
                                errorPEP12_df8 = simPEP{1,12} - expPEP{1,12};
                    wDesp_t = wDesp';
                    errorPEP = [...
                        wDesp_t(1,1)*errorPEP1_df8; wDesp_t(2,1)*errorPEP1_df4; wDesp_t(3,1)*errorPEP1_df2; wDesp_t(4,1)*errorPEP1_df1; 
                        wDesp_t(1,2)*errorPEP2_df8; wDesp_t(2,2)*errorPEP2_df4; wDesp_t(3,2)*errorPEP2_df2; wDesp_t(4,2)*errorPEP2_df1;
                        wDesp_t(1,3)*errorPEP3_df8; wDesp_t(2,3)*errorPEP3_df4; wDesp_t(3,3)*errorPEP3_df2; wDesp_t(4,3)*errorPEP3_df1;
                        wDesp_t(1,4)*errorPEP4_df8; wDesp_t(2,4)*errorPEP4_df4; wDesp_t(3,4)*errorPEP4_df2; wDesp_t(4,4)*errorPEP4_df1;
                        wDesp_t(1,5)*errorPEP5_df8; wDesp_t(2,5)*errorPEP5_df4; wDesp_t(3,5)*errorPEP5_df2; wDesp_t(4,5)*errorPEP5_df1;
                        wDesp_t(1,6)*errorPEP6_df8; wDesp_t(2,6)*errorPEP6_df4; wDesp_t(3,6)*errorPEP6_df2; wDesp_t(4,6)*errorPEP6_df1;
                        wDesp_t(1,7)*errorPEP7_df8; wDesp_t(2,7)*errorPEP7_df4; wDesp_t(3,7)*errorPEP7_df2; wDesp_t(4,7)*errorPEP7_df1;
                        wDesp_t(1,8)*errorPEP8_df8; wDesp_t(2,8)*errorPEP8_df4; wDesp_t(3,8)*errorPEP8_df2; wDesp_t(4,8)*errorPEP8_df1;
                        wDesp_t(1,9)*errorPEP9_df8; wDesp_t(2,9)*errorPEP9_df4; wDesp_t(3,9)*errorPEP9_df2; wDesp_t(4,9)*errorPEP9_df1;
                        wDesp_t(1,10)*errorPEP10_df8; wDesp_t(2,10)*errorPEP10_df4; wDesp_t(3,10)*errorPEP10_df2; wDesp_t(4,10)*errorPEP10_df1;
                        wDesp_t(1,11)*errorPEP11_df8; wDesp_t(2,11)*errorPEP11_df4; wDesp_t(3,11)*errorPEP11_df2; wDesp_t(4,11)*errorPEP11_df1;
                        wDesp_t(1,12)*errorPEP12_df8; wDesp_t(2,12)*errorPEP12_df4; wDesp_t(3,12)*errorPEP12_df2; wDesp_t(4,12)*errorPEP12_df1];
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';  
                     
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
    case 'pgi'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = setup.fixedValues;
% % % %             data.Vmaxs = data.Vmax(j,:);
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
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(2);
                    xassay(2) = x_temp(14);
                    xassay(3) = x_temp(1);
                case 2
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(3);
                    xassay(2) = x_temp(15);
                    xassay(3) = x_temp(1);
                case 3
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(4);
                    xassay(2) = x_temp(16);
                    xassay(3) = x_temp(1);
                case 4
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(5);
                    xassay(2) = x_temp(17);
                    xassay(3) = x_temp(1);
                case 5
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(6);
                    xassay(2) = x_temp(18);
                    xassay(3) = x_temp(1);
                case 6
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(7);
                    xassay(2) = x_temp(19);
                    xassay(3) = x_temp(1);
                case 7
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(8);
                    xassay(2) = x_temp(20);
                    xassay(3) = x_temp(1);
                case 8
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(9);
                    xassay(2) = x_temp(21);
                    xassay(3) = x_temp(1);
                case 9
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(10);
                    xassay(2) = x_temp(22);
                    xassay(3) = x_temp(1);
                case 10
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(11);
                    xassay(2) = x_temp(23);
                    xassay(3) = x_temp(1);
                case 11
                    %xassay = zeros(1,3);
                    xassay(1) = x_temp(12);
                    xassay(2) = x_temp(24);
                    xassay(3) = x_temp(1);
                case 12
                    %xassay = zeros(1,3); 
                    xassay(1) = x_temp(13);
                    xassay(2) = x_temp(25);
                    xassay(3) = x_temp(1);
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
%                 disp(data.Vmaxs(1,4));
                data.chosenNADHini = data.NADH{i}(1);                
                % simulate metabolites
                [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                % calculation of reaction rates
                [vObs,~] = calcRates(xassay,simResult,data,setup);   
                % cost function (NADHexp + vGAPDHr)
                simTime = simResult.t;
                simMet = simResult.y(:,obsMet);
                simMet2 = simResult.y(:,1);
                simRate = vObs(:,5);

                simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                simGPD{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                expNADH{i,j} = data.NADH{i};
                expGPD{i,j} = data.Vprofs{i};
                simG6P{i,j} = interp1(simTime,simMet2,data.tempTime{i},'pchip');
                
                
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
                    if j == 1
                        h103 = figure(103);
                    end
%                     set(101);
                    set(0,'CurrentFigure', h103);
%                     set(gcf,'color','w');
                    subplot(3,4,j)
                    for i = DFstudy
                        plot(data.tempTime{i}, simG6P{i,j},'-','LineWidth',2)
                        ylim([4 5])
                        xlim([0 1000])
                    end
                    if j == numpH
                        suptitle('G6P concentration [mM] vs assay time [s]')
                    end
                    
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
    case 'hxk'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = setup.fixedValues;
% % % %             data.Vmaxs = data.Vmax(j,:);
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
                xassay = zeros(1,5);
                xassay(1) = x_temp(1+j);
                xassay(2) = x_temp(13+j);
                xassay(3) = x_temp(25+j);
                xassay(4) = x_temp(37+j);
                xassay(5) = x_temp(1);
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
            data.Vmaxs = setup.fixedValues;
%             data.Vmaxs = data.Vmax(j,:);
%                 data.Vmaxs(4) = mean([data.Vmaxs(1)*8, data.Vmaxs(2)*4]); % testing df4_8 case as starting vm
%                 mean([data.Vmax(:,1)*8 data.Vmax(:,2)*4],2)
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
                    xassay(1) = x_temp(1+j);
                    xassay(2) = x_temp(13+j);
                    xassay(3) = x_temp(26+j);
                    xassay(4) = x_temp(1);
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
                        plot(data.tempTime{i}, -expGPD{i,j},'k+')
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
            case 1 % DF1_2
                    errorHaldane = 0;
                    % DF1
                    errorNADH1 = simNADH{4,1} - expNADH{4,1};
                    errorNADH2 = simNADH{4,2} - expNADH{4,2};
                    errorNADH3 = simNADH{4,3} - expNADH{4,3};
                    errorNADH4 = simNADH{4,4} - expNADH{4,4};
                    errorNADH5 = simNADH{4,5} - expNADH{4,5};
                    errorNADH6 = simNADH{4,6} - expNADH{4,6};
                    errorNADH7 = simNADH{4,7} - expNADH{4,7};
                    errorNADH8 = simNADH{4,8} - expNADH{4,8};
                        % DF2
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
                     
            case 2 % DF4_8
                    errorHaldane = 0;
                    % DF4
                    errorNADH1 = simNADH{2,1} - expNADH{2,1};
                    errorNADH2 = simNADH{2,2} - expNADH{2,2};
                    errorNADH3 = simNADH{2,3} - expNADH{2,3};
                    errorNADH4 = simNADH{2,4} - expNADH{2,4};
                    errorNADH5 = simNADH{2,5} - expNADH{2,5};
                    errorNADH6 = simNADH{2,6} - expNADH{2,6};
                    errorNADH7 = simNADH{2,7} - expNADH{2,7};
                    errorNADH8 = simNADH{2,8} - expNADH{2,8};
                        % DF8
                        errorNADH1_2 = simNADH{1,1} - expNADH{1,1};
                        errorNADH2_2 = simNADH{1,2} - expNADH{1,2};
                        errorNADH3_2 = simNADH{1,3} - expNADH{1,3};
                        errorNADH4_2 = simNADH{1,4} - expNADH{1,4};
                        errorNADH5_2 = simNADH{1,5} - expNADH{1,5};
                        errorNADH6_2 = simNADH{1,6} - expNADH{1,6};
                        errorNADH7_2 = simNADH{1,7} - expNADH{1,7};
                        errorNADH8_2 = simNADH{1,8} - expNADH{1,8};
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
    otherwise
        disp('No enzyme has been selected in the cost function file');
        
end

end

% %% memoryDump
% eD1 = sum(abs(errorNADH1));
% eD2 = sum(abs(errorNADH2));
% eD3 = sum(abs(errorNADH3));
% eD4 = sum(abs(errorNADH4));
% eD5 = sum(abs(errorNADH5));
% eD6 = sum(abs(errorNADH6));
% eD7 = sum(abs(errorNADH7));
% eD8 = sum(abs(errorNADH8));
% eD9 = sum(abs(errorNADH9));
% eD10 = sum(abs(errorNADH10));
% eDarray = [eD1, eD2, eD3, eD4, eD5, eD6, eD7, eD8, eD9, eD10];
% 
% figure
% plot(setup.pH_vals, eDarray,'-o')
% xlabel('pH value')
% ylabel('errorData')

