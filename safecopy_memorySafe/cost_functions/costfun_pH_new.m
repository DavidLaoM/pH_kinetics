function [error] = costfun_pH_new(x_temp,data,setup)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
enzyme = setup.enzymeName;
DFs = setup.DFactorsTotal;
DFstudy = setup.DFstudy; % default is the lowest dilution factor (usually DF1, location 4)
obsMet = setup.observableMetabolite;
costfun = setup.costfun; % default value is 1. No regularization
lambda = setup.selectedLambda;
plotEachSimCF = setup.plotEachSimCF;
pHarray = data.pH(:,1);
sourceVm = setup.sourceVm;
wD = setup.weightData;
wH = setup.weightHaldane;
ode_pH = setup.ode_pH;
% optional simulation and plotting of all the datasets
if((setup.plotEachSimCF == 1)&&(setup.simAllProfiles == 1))
    DFstudy = [1 2 3 4];
end
% x_temp(2:5) = zeros;
% x_temp([1,6]) = zeros;

switch enzyme
    
    case 'gapdhr'
        
        % simulations
        simNADH = cell(DFs,1);
        expNADH = cell(DFs,1);
        simGAPDHr = cell(DFs,1);
        expGAPDHr = cell(DFs,1);
        simRes = cell(1,length(DFstudy));
        for i = DFstudy
            % recall vmax for the specific value and simulate
% % % %             data.chosenVmf = data.Vmf/data.DF(1,i);
% % % %             data.chosenVmr = data.Vmr/data.DF(1,i);
            data.chosenDF = data.DF(1,i);
            data.chosenKeqGAPDH = data.KeqGAPDH;
            data.chosenKeqPGK = data.KeqPGK;
            data.chosenNADini = data.NADH{i}(1);
%             disp(data.chosenNADini);
            
            % simulations ODEs + calculation fluxes
            [simResult] = simSys(x_temp,data,setup);
            [vObs,~] = calcRates(x_temp,simResult,data,setup);   
            
            % cost function (NADHexp + vGAPDHr)
            simTime = simResult.t;
            simMet = simResult.y(:,obsMet);
            simRate = vObs;
            simRes{i} = simResult;
            simRes{i}.v = simRate;
            simNADH{i} = interp1(simTime,simMet,data.tempTime{i},'pchip');
            simGAPDHr{i} = interp1(simTime,simRate,data.tempTime{i},'pchip');
            expNADH{i} = data.NADH{i};
            expGAPDHr{i} = -data.Vprofs{i};
        end
        
        % create the cost function
        switch costfun
            case 1 % DF1
                    errorNADH = wD * (simNADH{4} - expNADH{4});
                        Keq = data.chosenKeqGAPDH; %[]
                        switch sourceVm
                            case 'literature'
                                vmf = 10 .^ x_temp(1) .* 1184.52/60; % mM s^{-1}
                                vmr = 10 .^ x_temp(6) .* 6549.8/60; % mM s^{-1}            
                            case 'experimentalSlopes'
                                vmf = 10.^x_temp(1).*setup.exp_vmax_gapdhf(data.i);% mM s^{-1}
                                vmr = 10.^x_temp(6).*setup.exp_vmax_gapdhr(data.i); % mM s^{-1} %.*data.chosenVmax 
                            case 'experimentalSlopesFixed'
                                vmf = 10.^x_temp(1).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
                                vmr = 10.^x_temp(6).*setup.exp_vmax_gapdhr(6); % mM s^{-1} %.*data.chosenVmax
                            otherwise
                                disp('No source for vmax has been selected');
                        end
                        ks1 = 10 .^ x_temp(2) .* 2.48; % mM %k_gap
                        ks2 = 10 .^ x_temp(4) .* 2.92; %mM %k_nad
                        kp1 = 10 .^ x_temp(3) .* 1.18; % mM %k_bpg
                        kp2 = 10 .^ x_temp(5) .* 0.022; % mM %k_nadh
                        switch ode_pH
                            case 'on'
                                H_effect = 10^(setup.pH_vals(data.i) - setup.pH_vals(6));
%                                 Keq_haldane_estimated(j) =  (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2);
                                errorHaldane = wH * (Keq - (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2) );
                            otherwise
                                errorHaldane = wH * (Keq - (vmf * kp1 * kp2) / (vmr * ks1 * ks2) );
                        end
                    errorReg = lambda * x_temp';
            case 2 % DF1,2
%                     errorNADH1 = simNADH{4} - expNADH{4};
%                     errorNADH2 = simNADH{3} - expNADH{3};
%                     errorNADH = [errorNADH1; errorNADH2];
%                     errorGAPDHr = 0;
%                     errorKeq = 0;
%                     errorVmax = 0;
%                     errorLinking = 0;
%                     errorKm = 0; 
%                     errorReg = lambda * x_temp';
%                     errorHaldane = [1 1]; % dummy to cause error (vertcat) and then adjust with final structure
            case 3 % DF1,2,4,8 (all)
%                     errorNADH1 = simNADH{4} - expNADH{4};
%                     errorNADH2 = simNADH{3} - expNADH{3};
%                     errorNADH3 = simNADH{2} - expNADH{2};
%                     errorNADH4 = simNADH{1} - expNADH{1};
%                     errorNADH = [errorNADH1; errorNADH2; errorNADH3; errorNADH4];
%                     errorGAPDHr = 0;
%                     errorKeq = 0;
%                     errorVmax = 0;
%                     errorLinking = 0;
%                     errorKm = 0;  
%                     errorReg = lambda * x_temp';
%                     errorHaldane = [1 1]; % dummy to cause error (vertcat) and then adjust with final structure
            otherwise
                disp('No specific cost function has been appointed');
        end
        error = [
            errorNADH;
            errorHaldane;
            errorReg;
            ];        
%         disp(lambda);
    otherwise
        disp('No enzyme has been selected in the cost function file');        
end

if plotEachSimCF == 1
    if setup.simAllProfiles == 0
        figure
        subplot(1,2,1)
        for j = DFstudy
            simRes = simResult;
            plot(simRes.t,simRes.y(:,8),'-')
            hold on
            plot(data.time{data.i,j},data.conc_mean{data.i,j},'k+')
        end
        title('NADH')
        subplot(1,2,2)
        for j = DFstudy
            simRRs = vObs;
            plot(simRes.t,simRRs,'-')
            hold on
            plot(data.time{i,j},-data.RRs{i,j},'k+')
        end
        title('v_{GAPDHr}')       
        suptitle(erase(sprintf('pH = %d',pHarray(data.i)),"0000e+00"));
    elseif setup.simAllProfiles == 1
        figure
        subplot(1,2,1)
        for j = DFstudy
%             plot(data.tempTime{j}, simNADH{j},'-')
            plot(simRes{j}.t, simRes{j}.y(:,8),'-')
            hold on
            plot(data.tempTime{j}, expNADH{j},'k+')
            hold on
        end
        title('NADH')
        subplot(1,2,2)
        for j = DFstudy
%             plot(data.tempTime{j}, simGAPDHr{j},'-')
            plot(simRes{j}.t, simRes{j}.v, '-')
            hold on
            plot(data.tempTime{j}, expGAPDHr{j},'k+')
            hold on
        end
        title('v_{apparent.GAPDHr}')
        suptitle(erase(sprintf('pH = %d',pHarray(data.i)),"0000e+00"));
    end

end

end

% %%
% 
% errorNADH = simNADH{4} - expNADH{4};
% errorReg = lambda * x_temp';
% errorHaldane = wH * 1/10 * (Keq - (vmf * kp1 * kp2) / (vmr * ks1 * ks2) );
                    
                    
                    