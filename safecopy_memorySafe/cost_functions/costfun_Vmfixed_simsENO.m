function [error,simPEP,simENO,expPEP,expENO] = costfun_Vmfixed_simsENO(x_temp,data,setup)
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

