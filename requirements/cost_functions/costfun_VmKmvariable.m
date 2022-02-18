function [error] = costfun_VmKmvariable(x_temp,data,setup)
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
                % prepare idxs
                idx1 = 0 + j; %idx1 = 1 + (j - 1) * 3;
                idx2 = 12 + j; %2 + (j - 1) * 3;
                idx3 = 24 + j; %3 + (j - 1) * 3;
                % set up the values inside
                xassay = zeros(1,3);
                xassay(1) = x_temp(idx1);
                xassay(2) = x_temp(idx2);
                xassay(3) = x_temp(idx3);
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
                    simulationVisualization;
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101)
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 1.2])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('PEP concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102)
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 1.5E-3])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('ENO reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end               
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
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';
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
    otherwise
        disp('No enzyme has been selected in the cost function file');
        
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

