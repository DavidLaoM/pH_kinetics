function [error] = costfun_pH_FandR_fixed(x_temp,data,setup)
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

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
numpH = 10;
errorTemp = cell(numpH,1);
for o = 1:numpH 
    data.KeqGAPDH = setup.pH_Keq_gapdh_eQ(o); %setup.pH_Keq_gapdh(i);
    data.KeqPGK = setup.pH_Keq_pgk(o);
    setup.excessPGK = 1;
    data.gapdhR.NADH = data.gapdhR.conc_mean(o,:);
    data.gapdhR.Vprofs = data.gapdhR.RRs(o,:);
    data.gapdhR.tempTime = data.gapdhR.time(o,:);
    data.gapdhF.NADH = data.gapdhF.conc_mean(o,:);
    data.gapdhF.Vprofs = data.gapdhF.RRs(o,:);
    data.gapdhF.tempTime = data.gapdhF.time(o,:);
    data.i = o;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
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

    % error calculation
    rerrorNADH = wDR * (rsimNADH{4} - rexpNADH{4});
    ferrorNADH = wDF * (fsimNADH{4} - fexpNADH{4});
% % % %         rerrorNADH1 = wDR * (rsimNADH{4} - rexpNADH{4});
% % % %         rerrorNADH2 = wDR * (rsimNADH{3} - rexpNADH{3});
% % % %         rerrorNADH = [rerrorNADH1; rerrorNADH2];
                        
    Keq = data.chosenKeqGAPDH; %[]
    vmf = 10.^x_temp(1).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
    vmr = 10.^x_temp(6).*setup.exp_vmax_gapdhr(6); % mM s^{-1} %.*data.chosenVmax
    ks1 = 10 .^ x_temp(2) .* 2.48; % mM %k_gap
    ks2 = 10 .^ x_temp(4) .* 2.92; %mM %k_nad
    kp1 = 10 .^ x_temp(3) .* 1.18; % mM %k_bpg
    kp2 = 10 .^ x_temp(5) .* 0.022; % mM %k_nadh

    H_effect = 10^(setup.pH_vals(data.i) - setup.pH_vals(6));
    errorHaldane = wH * (Keq - (vmf * kp1 * kp2 * H_effect) / (vmr * ks1 * ks2) );

    errorReg = lambda * x_temp';
        errorReg(1) = errorReg(1)*2; errorReg(6) = errorReg(6)*2;
    
    % putting error together
    errorTemp{o} = [
        rerrorNADH;
        ferrorNADH;
        errorHaldane;
        errorReg;
        ];        
%         disp(lambda);


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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
end

% error = [errorTemp{1};
%         errorTemp{2};
%         errorTemp{3};
%         errorTemp{4};
%         errorTemp{5};
%         errorTemp{6};
%         errorTemp{7};
%         errorTemp{8};
%         errorTemp{9};
%         errorTemp{10}];

error = [sum(abs(errorTemp{1}));
        sum(abs(errorTemp{2}));
        sum(abs(errorTemp{3}));
        sum(abs(errorTemp{4}));
        sum(abs(errorTemp{5}));
        sum(abs(errorTemp{6}));
        sum(abs(errorTemp{7}));
        sum(abs(errorTemp{8}));
        sum(abs(errorTemp{9}));
        sum(abs(errorTemp{10}))];

end                   