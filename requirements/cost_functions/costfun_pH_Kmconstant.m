function [error] = costfun_pH_Kmconstant(x_temp,data,setup)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%     x_temp(1) = Kgap
%     x_temp(2) = Kbpg
%     x_temp(3) = Knad
%     x_temp(4) = Knadh
%     x_temp([5:7]) = {Keq, Vm, pLink}, pH#1
%     x_temp([8:10]) = {Keq, Vm, pLink}, pH#2
%     x_temp([11:13]) = {Keq, Vm, pLink}, pH#3
%     x_temp([14:16]) = {Keq, Vm, pLink}, pH#4
%     x_temp([17:19]) = {Keq, Vm, pLink}, pH#5
%     x_temp([20:22]) = {Keq, Vm, pLink}, pH#6
%     x_temp([23:25]) = {Keq, Vm, pLink}, pH#7
%     x_temp([26:28]) = {Keq, Vm, pLink}, pH#8
%     x_temp([29:31]) = {Keq, Vm, pLink}, pH#9
%     x_temp([32:34]) = {Keq, Vm, pLink}, pH#10
enzyme = setup.enzymeName;
DFs = setup.DFactorsTotal;
DFstudy = setup.DFstudy; % default is the lowest dilution factor (usually DF1, location 4)
obsMet = setup.observableMetabolite;
costfun = setup.costfun; % default value is 1. No regularization
costfun2 = setup.costfun2; % default value is 0. Knadh reg. at 0.
lambda = setup.selectedLambda;
selectedLambda = setup.selectedLambda;
numpH = setup.numpHtested;



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

            % selecting the right parameters
            switch j
                case 1
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(5);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(7);
                case 2
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(8);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(9);
                    xassay(7) = x_temp(10);
                case 3
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(11);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(12);
                    xassay(7) = x_temp(13);
                case 4
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(14);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(15);
                    xassay(7) = x_temp(16);
                case 5
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(17);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(18);
                    xassay(7) = x_temp(19);
                case 6
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(20);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(21);
                    xassay(7) = x_temp(22);
                case 7
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(23);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(24);
                    xassay(7) = x_temp(25);
                case 8
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(26);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(27);
                    xassay(7) = x_temp(28);
                case 9
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(29);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(30);
                    xassay(7) = x_temp(31);
                case 10
                    xassay = zeros(1,7);
                    xassay(1) = x_temp(32);
                    xassay(2) = x_temp(1);
                    xassay(3) = x_temp(2);
                    xassay(4) = x_temp(3);
                    xassay(5) = x_temp(4);
                    xassay(6) = x_temp(33);
                    xassay(7) = x_temp(34);
                otherwise
                    disp('Something went wrong is selecting the pH value');
            end
            
            % simulations
            for i = DFstudy
                % recall vmax for the specific value and simulate
                data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,i); % vmax from the highest DF is taken and then divided
                data.chosenLink = data.DF(1,i);
                data.chosenNADini = data.NADH{i}(1);
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
                    errorGAPDHr = 0;
                    errorKeq = 0;
                    errorVmax = 0;
                    errorLinking = 0;
                    errorKm = 0; 
                    errorReg = lambda * x_temp';               
 
            otherwise
                disp('No specific cost function has been appointed');
        end
        error = [
            errorNADH1;
            errorNADH2;
            errorNADH3;
            errorNADH4;
            errorNADH5;
            errorNADH6;
            errorNADH7;
            errorNADH8;
            errorNADH9;
            errorNADH10;
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

end



