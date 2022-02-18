function [error] = costfun_pH(x_temp,data,setup)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
enzyme = setup.enzymeName;
DFs = setup.DFactorsTotal;
DFstudy = setup.DFstudy; % default is the lowest dilution factor (usually DF1, location 4)
obsMet = setup.observableMetabolite;
costfun = setup.costfun; % default value is 1. No regularization
costfun2 = setup.costfun2; % default value is 0. Knadh reg. at 0.
lambda = setup.selectedLambda;
% lambda = 10;

% % % % % 
% x_temp = [x_temp, 0];
% % % % % 

switch enzyme
    
    case 'gapdhr'
        
        % simulations
        simNADH = cell(DFs,1);
        expNADH = cell(DFs,1);
        simGAPDHr = cell(DFs,1);
        expGAPDHr = cell(DFs,1);
        for i = DFstudy
            % recall vmax for the specific value and simulate
            data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,i); % vmax from the highest DF is taken and then divided
            data.chosenLink = data.DF(1,i);
            data.chosenNADini = data.NADH{i}(1);
            % check only one vmax
            if setup.constantVm == 1
                data.chosenVmax = 0.0021/data.DF(1,i); % fixing at the maximum
%                 data.chosenVmax = 0.0009/data.DF(1,i); % fixing at the minimum
%                 data.chosenVmax = 0.0015/data.DF(1,i); % middle value
            else
            end
            [simResult] = simSys(x_temp,data,setup);
            % calculation of reaction rates
            [vObs,~] = calcRates(x_temp,simResult,data,setup);   
            % cost function (NADHexp + vGAPDHr)
            simTime = simResult.t;
            simMet = simResult.y(:,obsMet);
            simRate = vObs;
            simNADH{i} = interp1(simTime,simMet,data.tempTime{i},'pchip');
            simGAPDHr{i} = interp1(simTime,simRate,data.tempTime{i},'pchip');
            expNADH{i} = data.NADH{i};
            expGAPDHr{i} = -data.Vprofs{i};
        end
        
        % create the cost function
        switch costfun
            case 1 % DF1
                    errorNADH = simNADH{4} - expNADH{4};
                    errorGAPDHr = 0;
                    errorKeq = 0;
                    errorVmax = 0;
                    errorLinking = 0;
                    errorKm = 0; 
                    errorReg = lambda * x_temp';
                    switch costfun2
                        case 1
%                             errorReg(5) = 1 * (x_temp(5) + 1);
%                             errorReg(5) = 1 * (x_temp(5) + 0.5);
%                             errorReg(5) = 1 * (x_temp(5) - 0);
%                             errorReg(5) = 1 * (x_temp(5) - 0.5);
                            errorReg(5) = 1 * (x_temp(5) - 1);
                        case 0
                        otherwise
                            disp('costfun2 active but not specified.');
                    end                        
            case 2 % DF1,2
                    errorNADH1 = simNADH{4} - expNADH{4};
                    errorNADH2 = simNADH{3} - expNADH{3};
                    errorNADH = [errorNADH1; errorNADH2];
                    errorGAPDHr = 0;
                    errorKeq = 0;
                    errorVmax = 0;
                    errorLinking = 0;
                    errorKm = 0; 
                    errorReg = lambda * x_temp';
                    switch costfun2
                        case 1
%                             errorReg(5) = 1 * (x_temp(5) + 1);
%                             errorReg(5) = 1 * (x_temp(5) + 0.5);
%                             errorReg(5) = 1 * (x_temp(5) - 0);
%                             errorReg(5) = 1 * (x_temp(5) - 0.5);
                            errorReg(5) = 1 * (x_temp(5) - 1);
                        case 0
                        otherwise
                            disp('costfun2 active but not specified.');
                    end
            case 3 % DF1,2,4,8 (all)
                    errorNADH1 = simNADH{4} - expNADH{4};
                    errorNADH2 = simNADH{3} - expNADH{3};
                    errorNADH3 = simNADH{2} - expNADH{2};
                    errorNADH4 = simNADH{1} - expNADH{1};
                    errorNADH = [errorNADH1; errorNADH2; errorNADH3; errorNADH4];
                    errorGAPDHr = 0;
                    errorKeq = 0;
                    errorVmax = 0;
                    errorLinking = 0;
                    errorKm = 0;  
                    errorReg = lambda * x_temp';
                    switch costfun2
                        case 1
%                             errorReg(5) = 1 * (x_temp(5) + 1);
%                             errorReg(5) = 1 * (x_temp(5) + 0.5);
%                             errorReg(5) = 1 * (x_temp(5) - 0);
%                             errorReg(5) = 1 * (x_temp(5) - 0.5);
                            errorReg(5) = 1 * (x_temp(5) - 1);
                        case 0
                        otherwise
                            disp('costfun2 active but not specified.');
                    end
            otherwise
                disp('No specific cost function has been appointed');
        end
        error = [
            errorNADH;
            errorGAPDHr;
            errorVmax;
            errorKeq;
            errorLinking;
            errorKm;
            errorReg;
            ];        
%         disp(lambda);
    otherwise
        disp('No enzyme has been selected in the cost function file');
        
end
% esize = sprintf('The size of the error array is %d',length(error));
% disp(esize);
% disp(lambda);
end