function [error] = costfun_allVars(x_temp,data,setup)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% 
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
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                case 2
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(7);
                    xassay(2) = x_temp(8);
                    xassay(3) = x_temp(9);
                    xassay(4) = x_temp(10);
                    xassay(5) = x_temp(11);
                    xassay(6) = x_temp(12);
                case 3
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(13);
                    xassay(2) = x_temp(14);
                    xassay(3) = x_temp(15);
                    xassay(4) = x_temp(16);
                    xassay(5) = x_temp(17);
                    xassay(6) = x_temp(18);
                case 4
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(19);
                    xassay(2) = x_temp(20);
                    xassay(3) = x_temp(21);
                    xassay(4) = x_temp(22);
                    xassay(5) = x_temp(23);
                    xassay(6) = x_temp(24);
                case 5
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(25);
                    xassay(2) = x_temp(26);
                    xassay(3) = x_temp(27);
                    xassay(4) = x_temp(28);
                    xassay(5) = x_temp(29);
                    xassay(6) = x_temp(30);
                case 6
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(31);
                    xassay(2) = x_temp(32);
                    xassay(3) = x_temp(33);
                    xassay(4) = x_temp(34);
                    xassay(5) = x_temp(35);
                    xassay(6) = x_temp(36);
                case 7
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(37);
                    xassay(2) = x_temp(38);
                    xassay(3) = x_temp(39);
                    xassay(4) = x_temp(40);
                    xassay(5) = x_temp(41);
                    xassay(6) = x_temp(42);
                case 8
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(43);
                    xassay(2) = x_temp(44);
                    xassay(3) = x_temp(45);
                    xassay(4) = x_temp(46);
                    xassay(5) = x_temp(47);
                    xassay(6) = x_temp(48);
                case 9
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(49);
                    xassay(2) = x_temp(50);
                    xassay(3) = x_temp(51);
                    xassay(4) = x_temp(52);
                    xassay(5) = x_temp(53);
                    xassay(6) = x_temp(54);
                case 10
                    xassay = zeros(1,6);
                    xassay(1) = x_temp(55);
                    xassay(2) = x_temp(56);
                    xassay(3) = x_temp(57);
                    xassay(4) = x_temp(58);
                    xassay(5) = x_temp(59);
                    xassay(6) = x_temp(60);
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
                                H_effect = 10^(setup.pH_vals(data.i) - setup.pH_vals(6));
                                errorHaldane1 = (Keq(1) - (vmf1 * kp1 * kp2 * H_effect) / (vmr1 * ks1 * ks2) );
                                errorHaldane2 = (Keq(2) - (vmf2 * kp1 * kp2 * H_effect) / (vmr2 * ks1 * ks2) );
                                errorHaldane3 = (Keq(3) - (vmf3 * kp1 * kp2 * H_effect) / (vmr3 * ks1 * ks2) );
                                errorHaldane4 = (Keq(4) - (vmf4 * kp1 * kp2 * H_effect) / (vmr4 * ks1 * ks2) );
                                errorHaldane5 = (Keq(5) - (vmf5 * kp1 * kp2 * H_effect) / (vmr5 * ks1 * ks2) );
                                errorHaldane6 = (Keq(6) - (vmf6 * kp1 * kp2 * H_effect) / (vmr6 * ks1 * ks2) );
                                errorHaldane7 = (Keq(7) - (vmf7 * kp1 * kp2 * H_effect) / (vmr7 * ks1 * ks2) );
                                errorHaldane8 = (Keq(8) - (vmf8 * kp1 * kp2 * H_effect) / (vmr8 * ks1 * ks2) );
                                errorHaldane9 = (Keq(9) - (vmf9 * kp1 * kp2 * H_effect) / (vmr9 * ks1 * ks2) );
                                errorHaldane10 =(Keq(10) - (vmf10 * kp1 * kp2 * H_effect) / (vmr10 * ks1 * ks2) );
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

