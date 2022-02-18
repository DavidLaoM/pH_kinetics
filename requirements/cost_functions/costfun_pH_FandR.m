function [error] = costfun_pH_FandR(x_temp,data,setup)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
enzyme = setup.enzymeName;
DFs = setup.DFactorsTotal;
DFstudy = setup.DFstudy; % default is the lowest dilution factor (usually DF1, location 4)
obsMet = setup.observableMetabolite;
costfun = setup.costfun; % default value is 1. No regularization
lambda = setup.selectedLambda;
plotEachSimCF = setup.plotEachSimCF;
pHarray = data.gapdhR.pH(:,1);
sourceVm = setup.sourceVm;
wDR = setup.weightDataR;
wDF = setup.weightDataF;
wH = setup.weightHaldane;
ode_pH = setup.ode_pH;
typeVm = setup.typeVm;

% % % % wH = 0;
% % % % disp(lambda);
% % % % disp(wDR);
% % % % disp(wDF);
% % % % disp(wH);

% optional simulation and plotting of all the datasets
if((setup.plotEachSimCF == 1)&&(setup.simAllProfiles == 1))
    DFstudy = [1 2 3 4];
end
% x_temp(2:5) = zeros;
% x_temp([1,6]) = zeros;
% x_temp([1,6,7,8]) = zeros;
% x_temp([1,6,7,8]) = zeros; x_temp(1) = 3; x_temp(6) = 0.7;
% x_temp(1) = 1.8620; x_temp(6) = 1.9253; % case fixed vms
% % % % x_temp(1) = 3; x_temp(6) = 0.7; % case fixed vms
% % % % x_temp(1) = 0; x_temp(6) = 0;

% % Fixing results from the previous round:
% % x_temp(1) = 1.8620;
% x_temp(2) = -0.9592;
% x_temp(3) = -0.4083;
% x_temp(4) = -0.9628;
% x_temp(5) = 0.9151;
% % x_temp(6) = 1.9253;
% x_temp(7) = -0.6412;
% x_temp(8) = -0.0939;

% % using both gapdhF and gapdhR data simultaneously
% x_temp(1) = 1.6934;
% % x_temp(2) = -0.9992;
% % x_temp(3) = -0.1930;
% % x_temp(4) = -0.9222;
% % x_temp(5) = 0.9769;
% x_temp(6) = 2.1853;
% x_temp(7) = -0.5059;
% x_temp(8) = 1.7793;

% % % % using full kinetics with Keq
% % % x_temp(1) = 3.4869;
% x_temp(1) = 2;
% % x_temp(2) = -0.5192;
% % x_temp(3) = 0.1867;
% % x_temp(4) = -0.5119;
% % x_temp(5) = -0.1666;
% % % % x_temp(6) = 3.4291;
% x_temp(7) = 0;
% x_temp(8) = 0;

% % % ode_pH = 'on'
% x_temp(1) = 1.9694;
% % x_temp(2) = -0.9122;
% % x_temp(3) = 0.0510;
% % x_temp(4) = -0.8527;
% % x_temp(5) = 0.9665;
% x_temp(6) = 2.3869;
% % % x_temp(7) = -1.1284;
% % % x_temp(8) = -1.1742;

% % % ode_pH = 'on'
% % x_temp(1) = 1.9684;
% x_temp(2) = -0.9993;
% x_temp(3) = 0.3829;
% x_temp(4) = -0.9773;
% x_temp(5) = 0.9914;
% % x_temp(6) = 2.7449;
% x_temp(7) = -1.1372;
% x_temp(8) = -2.7234;

% % % % % % % DF12, FWD+REV, one_pH = 'on'
% % % % x_temp(1) = 2.7481;   
% % % % % x_temp(2) = -0.4878;   
% % % % % x_temp(3) = -0.2746;  
% % % % % x_temp(4) = -0.5066;    
% % % % % x_temp(5) = 0.6720;  
% % % % x_temp(6) = 1.9146;   
% % % % x_temp(7) = -0.2367;    
% % % % x_temp(8) = 2.8005;

switch enzyme
    
    case 'gapdhr'
        
        % simulations
        % reverse
        rsimNADH = cell(DFs,1);
        rexpNADH = cell(DFs,1);
        rsimGAPDHr = cell(DFs,1);
        rexpGAPDHr = cell(DFs,1);
        % forward
        fsimNADH = cell(DFs,1);
        fexpNADH = cell(DFs,1);
        fsimGAPDHr = cell(DFs,1);
        fexpGAPDHr = cell(DFs,1);

        simRes = cell(1,length(DFstudy));
        for i = DFstudy
            % recall vmax for the specific value and simulate
% % % %             data.chosenVmf = data.Vmf/data.DF(1,i);
% % % %             data.chosenVmr = data.Vmr/data.DF(1,i);
            data.chosenDF = data.gapdhR.DF(1,i);
            data.chosenKeqGAPDH = data.KeqGAPDH;
            data.chosenKeqPGK = data.KeqPGK;
            data.gapdhR.chosenNADini = data.gapdhR.NADH{i}(1);
            data.gapdhF.chosenNADini = data.gapdhF.NADH{i}(1);

            % simulations ODEs + calculation fluxes
            [simResult] = simSys_FandR(x_temp,data,setup);
            [vObs,~] = calcRates_FandR(x_temp,simResult,data,setup);   
            
            % cost function (NADHexp + vGAPDHr)
            simTime = simResult.t;
            simMetR = simResult.y(:,obsMet);
            simMetF = simResult.y(:,2*obsMet);
            simRateR = vObs(:,1);
            simRateF = vObs(:,2);
            simRes{i} = simResult;
            simRes{i}.v = vObs;
            
            % reverse
            rsimNADH{i} = interp1(simTime,simMetR,data.gapdhR.tempTime{i},'pchip');
            rexpNADH{i} = data.gapdhR.NADH{i};
            rsimGAPDHr{i} = interp1(simTime,simRateR,data.gapdhR.tempTime{i},'pchip');
            rexpGAPDHr{i} = data.gapdhR.Vprofs{i};
            % forward
            fsimNADH{i} = interp1(simTime,simMetF,data.gapdhF.tempTime{i},'pchip');
            fexpNADH{i} = data.gapdhF.NADH{i};
            fsimGAPDHr{i} = interp1(simTime,simRateF,data.gapdhF.tempTime{i},'pchip');
            fexpGAPDHr{i} = data.gapdhF.Vprofs{i};
            
        end
        
        % create the cost function
        switch costfun
            case 1 % DF1
                    rerrorNADH = wDR * (rsimNADH{4} - rexpNADH{4});
                    ferrorNADH = wDF * (fsimNADH{4} - fexpNADH{4});
% % % %                         rerrorNADH1 = wDR * (rsimNADH{4} - rexpNADH{4});
% % % %                         rerrorNADH2 = wDR * (rsimNADH{3} - rexpNADH{3});
% % % %                         rerrorNADH = [rerrorNADH1; rerrorNADH2];
%                     ferrorNADH(end-10:end) = ferrorNADH(end-10:end) * 10;
% % % %                     ferrorNADH(round(length(ferrorNADH)/2)-5:round(length(ferrorNADH)/2)+5) = ferrorNADH(round(length(ferrorNADH)/2)-5:round(length(ferrorNADH)/2)+5) * 10;
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
                            errorHaldane = wH * (Keq - (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2) );
                        otherwise
                            errorHaldane = wH * (Keq - (vmf * kp1 * kp2) / (vmr * ks1 * ks2) );
                    end
                    errorReg = lambda * x_temp';
                        errorReg(1) = errorReg(1)*2; errorReg(6) = errorReg(6)*2;
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
            rerrorNADH;
            ferrorNADH;
            errorHaldane;
            errorReg;
            ];        
%         disp(lambda);
    otherwise
        disp('No enzyme has been selected in the cost function file');        
end

if plotEachSimCF == 1
    if setup.simAllProfiles == 0 % plots the fit of each data set
        figure
        subplot(2,2,1) % NADH reverse direction
        for j = DFstudy
            simRes = simResult;
            plot(simRes.t,simRes.y(:,8),'-')
            hold on
            plot(data.gapdhR.time{data.i,j},data.gapdhR.conc_mean{data.i,j},'k+')
        end
        title('rev NADH')
        subplot(2,2,2) % GAPDH reverse direction
        for j = DFstudy
            simRRs = vObs;
            plot(simRes.t,simRRs(:,1),'-')
            hold on
            plot(data.gapdhR.time{i,j},-data.gapdhR.RRs{i,j},'k+')
        end
        title('rev v_{GAPDH}')  
        subplot(2,2,3) % NADH forward direction
        for j = DFstudy
            simRes = simResult;
            plot(simRes.t,simRes.y(:,16),'-')
            hold on
            plot(data.gapdhF.time{data.i,j},data.gapdhF.conc_mean{data.i,j},'k+')
        end
        title('fwd NADH')
        subplot(2,2,4) % GAPDH forward direction
        for j = DFstudy
            simRRs = vObs;
            plot(simRes.t,simRRs(:,2),'-')
            hold on
            plot(data.gapdhF.time{i,j},-data.gapdhF.RRs{i,j},'k+')
        end
        title('fwd v_{GAPDH}')
        suptitle(erase(sprintf('pH = %d',pHarray(data.i)),"0000e+00"));
    elseif setup.simAllProfiles == 1 % 
        figure
        subplot(2,2,1) % NADH reverse direction
        for j = DFstudy
%             simRes = simResult;
            plot(simRes{j}.t,simRes{j}.y(:,8),'-')
            hold on
            plot(data.gapdhR.time{data.i,j},data.gapdhR.conc_mean{data.i,j},'k+')
        end
        xlim([0 300])
        title('rev NADH')
        subplot(2,2,2) % GAPDH reverse direction
        for j = DFstudy
%             simRRs = vObs;
            plot(simRes{j}.t,simRes{j}.v(:,1),'-')
            hold on
            plot(data.gapdhR.time{i,j},-data.gapdhR.RRs{i,j},'k+')
        end
        xlim([0 300])
        title('rev v_{GAPDH}')  
        subplot(2,2,3) % NADH forward direction
        for j = DFstudy
%             simRes = simResult;
            plot(simRes{j}.t,simRes{j}.y(:,16),'-')
            hold on
            plot(data.gapdhF.time{data.i,j},data.gapdhF.conc_mean{data.i,j},'k+')
        end
        xlim([0 300])
        title('fwd NADH')
        subplot(2,2,4) % GAPDH forward direction
        for j = DFstudy
%             simRRs = vObs;
            plot(simRes{j}.t,simRes{j}.v(:,2),'-')
            hold on
            plot(data.gapdhF.time{i,j},-data.gapdhF.RRs{i,j},'k+')
        end
        xlim([0 300])
        title('fwd v_{GAPDH}')
        suptitle(erase(sprintf('pH = %d',pHarray(data.i)),"0000e+00"));        
%         figure
%         subplot(1,2,1)
%         for j = DFstudy
% %             plot(data.tempTime{j}, simNADH{j},'-')
%             plot(simRes{j}.t, simRes{j}.y(:,8),'-')
%             hold on
%             plot(data.tempTime{j}, expNADH{j},'k+')
%             hold on
%         end
%         title('NADH')
%         subplot(1,2,2)
%         for j = DFstudy
% %             plot(data.tempTime{j}, simGAPDHr{j},'-')
%             plot(simRes{j}.t, simRes{j}.v, '-')
%             hold on
%             plot(data.tempTime{j}, expGAPDHr{j},'k+')
%             hold on
%         end
%         title('v_{apparent.GAPDHr}')
%         suptitle(erase(sprintf('pH = %d',pHarray(data.i)),"0000e+00"));
    end

end

end

                    
                    
                    