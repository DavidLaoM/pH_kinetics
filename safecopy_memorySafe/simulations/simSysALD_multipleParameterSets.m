function [error,FullSim] = simSysALD_multipleParameterSets(x_tempOri,data,setup)
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

[lenPset,~] = size(x_tempOri);
FullSim = cell(lenPset,1);
for o = 1:lenPset
    x_temp = x_tempOri(o,:);
    fprintf('Simulation %f of %f \n',o,lenPset);

    switch enzyme

        case 'ald'
            % zeros
            simNADH = cell(DFs,numpH);
            expNADH = cell(DFs,numpH);
            simGPD = cell(DFs,numpH);
            expGPD = cell(DFs,numpH);
            temp_simResult = cell(DFs,numpH);                   
            
            % simulations loop for each pH value
            for j = 1:numpH
    %         for j = 11:12
                % select required data
                data.Vmaxs = data.Vmax(j,:);
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
                    simResult.v = simRate;

                    simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                    simGPD{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                    expNADH{i,j} = data.NADH{i};
                    expGPD{i,j} = data.Vprofs{i};
                    
                    temp_simResult{i,j} = simResult;
                    
                end
            end

            for j = 1:numpH
                data.tempTime = data.time(j,:);
                if plotEachSimCF == 1
                    if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                        simulationVisualization;

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
                case 3 % selected DFs
                        errorHaldane = 0;
                        % DF1
                        errorNADH1_df1 = simNADH{4,1} - expNADH{4,1};
                        errorNADH2_df1 = simNADH{4,2} - expNADH{4,2};
                        errorNADH3_df1 = simNADH{4,3} - expNADH{4,3};
                        errorNADH4_df1 = simNADH{4,4} - expNADH{4,4};
                        errorNADH5_df1 = simNADH{4,5} - expNADH{4,5};
                        errorNADH6_df1 = simNADH{4,6} - expNADH{4,6};
                        errorNADH7_df1 = simNADH{4,7} - expNADH{4,7};
                        errorNADH8_df1 = simNADH{4,8} - expNADH{4,8};
                            % DF2
                            errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                            errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                            errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                            errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                            errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                            errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                            errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
                            errorNADH8_df2 = simNADH{3,8} - expNADH{3,8};
                                % DF4
                                errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                                errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                                errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                                errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                                errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                                errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                                errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
                                errorNADH8_df4 = simNADH{2,8} - expNADH{2,8};
                                    % DF8
                                    errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                    errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                    errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                    errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                    errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                    errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                    errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
                                    errorNADH8_df8 = simNADH{1,8} - expNADH{1,8};
                        wDesp_t = wDesp';
                        errorNADH = [...
                            wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                            wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                            wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                            wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                            wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                            wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                            wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1;
                            wDesp_t(1,8)*errorNADH8_df8; wDesp_t(2,8)*errorNADH8_df4; wDesp_t(3,8)*errorNADH8_df2; wDesp_t(4,8)*errorNADH8_df1];
                        errorReg = lambda * x_temp';
                otherwise
                    disp('No specific cost function has been appointed');
            end
            error = [
                wD * errorNADH;
                wH * errorHaldane;
                wL * errorReg;
                ];        

        otherwise
            disp('No enzyme has been selected in the cost function file');

    end
    
    FullSim{o,1} = temp_simResult;
        
end


end

% %% memoryDump
% eD1 = sum(abs(errorPEP1));
% eD2 = sum(abs(errorPEP2));
% eD3 = sum(abs(errorPEP3));
% eD4 = sum(abs(errorPEP4));
% eD5 = sum(abs(errorPEP5));
% eD6 = sum(abs(errorPEP6));
% eD7 = sum(abs(errorPEP7));
% eD8 = sum(abs(errorPEP8));
% eD9 = sum(abs(errorPEP9));
% eD10 = sum(abs(errorPEP10));
% eDarray = [eD1, eD2, eD3, eD4, eD5, eD6, eD7, eD8, eD9, eD10];
% 
% figure
% plot(setup.pH_vals, eDarray,'-o')
% xlabel('pH value')
% ylabel('errorData')

