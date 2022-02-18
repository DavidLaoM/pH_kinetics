function [error] = costfun_Kmfixed(x_temp,data,setup)
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
    
    case 'gapdhr'
        caseKm = setup.caseKm;
        
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
            
            % parameter setup
            switch caseKm
                case 'pH_independent'
                    % Mode 1. still keeping kms, but fixed
                    idx5 = 4 + j;
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(idx5);
                case 'pH_dependent'
                    idx1 = j + numpH * 0;
                    idx2 = j + numpH * 1;
                    idx3 = j + numpH * 2;
                    idx4 = j + numpH * 3;
                    idx5 = j + numpH * 4;
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(idx1);
                    xassay(2) = x_temp(idx2);
                    xassay(3) = x_temp(idx3);
                    xassay(4) = x_temp(idx4);
                    xassay(5) = x_temp(idx5);
                otherwise
                    disp('Warning: no specification has been made on Km being pH dependent or independent');
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
        end
        
% % % %         for j = 1:numpH
        for j = 4
            if plotEachSimCF == 1
                if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                    
% % % %                     simulationVisualization;

                     for ENOsim = 1
    %                     %% Simulation concentrations
                        % selection of the timeframe
                        data.tempTime = data.time(j,:);
                        % plottinh
                        if j == 4
% % % %                             h101 = figure(101);
                            h101 = figure(1000);
                            hold on
                        end
                        set(0,'CurrentFigure', h101);
% % % %                         p101 = subplot(3,4,j);
                        if setup.caseStudyENO == 1
                            p101 = subplot(3,3,2);
                            hold on
                        elseif setup.caseStudyGAPDHr == 1
                            p101 = subplot(3,3,5);
                            hold on
                        elseif setup.caseStudyPGM == 1
                            p101 = subplot(3,3,6);
                            hold on
                        end
                        
% % % %                         if setup.literatureValues == 0
% % % %                             p101 = subplot(3,3,2);
% % % %                         elseif setup.literatureValues == 1
% % % %                             p101 = subplot(3,3,5);
% % % %                         end
                        
                        for i = DFstudy
                            if setup.caseStudyENO == 1
% % % %                                 plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
% % % %                                 ylim([0 1.2])
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.2])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.2])
                                    hold on
                                end
                                
                                
                            elseif setup.caseStudyHXK == 1
                                plot(data.tempTime{i}, simNADPH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADPH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                            elseif setup.caseStudyGAPDH == 1
                        %         plot(data.time{j,i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                        %         hold on
                        %         plot(data.time{j,i}, expNADH{i,j},'k.','MarkerSize',4)
                        %         ylim([0 0.15])
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                        %         disp(sprintf('j=%d, i=%d',j,i));
                            elseif setup.caseStudyGAPDHr == 1
                                
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                end
                                
                                                                    
                            else
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                                if(setup.caseStudyPDC == 1)
                                    xlim([0 300])
                                elseif(setup.caseStudyPGM == 1)
                                    ylim([0 0.10])
                                    xlim([0 300])
                                elseif((setup.caseStudyPGI == 1)||(setup.caseStudyPFK == 1))
                                    ylim([0 0.15])
                                    xlim([0 1000])
                                elseif(setup.caseStudyTPI == 1)
                                    ylim([0 0.15])
                                    xlim([0 600])
                                end        
                            end
                        end
                        % ylabel left
                        if((j == 1)||(j == 5)||(j == 9))
                            if setup.caseStudyENO == 1
                                ylabel('PEP concentration [mM]')
                            elseif setup.caseStudyHXK == 1
                                ylabel('NAPDH concentration [mM]')
                            else 
                                ylabel('NADH concentration [mM]')        
                            end
                        end
                        % xlabel
                        if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
                            i_xlabel = 5:8;
                        else
                            i_xlabel = 9:12;
                        end
                        if((j == i_xlabel(1))||(j == i_xlabel(2))||(j == i_xlabel(3))||(j == i_xlabel(4)))
%                             xlabel('assay time [s]')
                        end
                        % pH text box
                        tempText = erase(sprintf('@pH %d', data.pH(j,1)),"0000e+00");
%                         text(30, p101.YLim(2)*0.9, tempText);
%                         text(300*0.95, p101.YLim(2)*0.9, tempText,...
%                             'HorizontalAlignment', 'right')
                        % suptitle
                        if j == numpH
% % % %                             textHere_pre = [setup.enzymeName,': concentration progression curve. Each box contains a different pH value.'];
% % % %                             textHere = {textHere_pre;...
% % % %                                 'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
% % % %                             suptitle(textHere)
% % % %                             set(101,'color','white'),
                        end
                        % title
% % % %                         title('train data: concentrations')
                        box on

    %                     %% Simulation reaction rates
                        if j == 4
% % % %                             h102 = figure(102);
                            h102 = figure(1000);
                            hold on
                        end
                        set(0,'CurrentFigure', h102);
% % % %                         p102 = subplot(3,4,j);
% % % %                         p102 = subplot(3,3,3);
% % % %                         hold on
                        
                        if setup.caseStudyENO == 1
                            p102 = subplot(3,3,3);
                            hold on
                        elseif setup.caseStudyGAPDHr == 1
                            p102 = subplot(3,3,6);
                            hold on
                        elseif setup.caseStudyPGM == 1
                            p102 = subplot(3,3,9);
                            hold on
                        end
                        
                        
% % % %                         if setup.literatureValues == 0
% % % %                             p102 = subplot(3,3,3);
% % % %                         elseif setup.literatureValues == 1
% % % %                             p102 = subplot(3,3,6);
% % % %                         end
                        
                        for i = DFstudy
                            if setup.caseStudyALD == 1
                                plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyENO == 1
% % % %                                 plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
% % % %                                 ylim([0 1.5E-3])
% % % %                                 hold on                               
                                
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.5E-3])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color', setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.5E-3])
                                    hold on
                                end
                                
                            elseif setup.caseStudyHXK == 1
                                plot(data.tempTime{i}, simG6PDH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expG6PDH{i,j},'k.','MarkerSize',4)
                                ylim([0 5E-4])
                            elseif setup.caseStudyPDC == 1
                                plot(data.tempTime{i}, simPDC{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expPDC{i,j},'k.','MarkerSize',4)
                                xlim([0 300])
                            elseif setup.caseStudyPYK == 1  
                                plot(data.tempTime{i}, simLDH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expLDH{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyPGM == 1  
                                plot(data.tempTime{i}, simPGM{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expPGM{i,j},'k.','MarkerSize',4)
                                xlim([0 300])
                            elseif setup.caseStudyPGI == 1
                                plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
                                xlim([0 1000])
                            elseif setup.caseStudyTPI == 1
                                plot(data.tempTime{i}, simTPI{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expTPI{i,j},'k.','MarkerSize',4)
                                xlim([0 600])
                            elseif setup.caseStudyPFK == 1
                                plot(data.tempTime{i}, simPFK{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expPFK{i,j},'k.','MarkerSize',4)
                                xlim([0 1000])
                            elseif setup.caseStudyGAPDH == 1
                                plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyGAPDHr == 1
% % % %                                 plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                
                                if setup.literatureValues == 0
                                    plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.003])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.003])
                                    hold on
                                end
                                
                            end
                        end
                        % ylabel left
                        if((j == 1)||(j == 5)||(j == 9))
                            if((setup.caseStudyALD == 1)||(setup.caseStudyPGI == 1)||(setup.caseStudyPFK == 1)||(setup.caseStudyTPI == 1))
                                ylabel('GPD reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyENO == 1
                                ylabel('ENO reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyHXK == 1
                                ylabel('G6PDH reaction rate [mM s^{-1}]')
                            elseif((setup.caseStudyPYK == 1)||(setup.caseStudyPDC == 1)||(setup.caseStudyPGM == 1))
                                ylabel('LDH reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyGAPDH== 1
                                ylabel('GAPDH_{fwd} reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyGAPDHr== 1
                                ylabel('GAPDH_{rev} reaction rate [mM s^{-1}]')
                            end
                        end
                        % xlabel
                        if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
                            i_xlabel = 5:8;
                        else
                            i_xlabel = 9:12;
                        end
                        if((j == i_xlabel(1))||(j == i_xlabel(2))||(j == i_xlabel(3))||(j == i_xlabel(4)))
%                             xlabel('assay time [s]')
                        end
                        % textbox
                        tempText = erase(sprintf('@pH %d', data.pH(j,1)),"0000e+00");
%                         text(30, p102.YLim(2)*0.9, tempText)
%                         text(300*0.95, p102.YLim(2)*0.9, tempText,...
%                             'HorizontalAlignment', 'right')
                        % suptitle
                        if j == numpH
% % % %                             textHere_pre = [setup.enzymeName,': reaction rate progression curve. Each box contains a different pH value.'];
% % % %                             textHere = {textHere_pre;...
% % % %                                 'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
% % % %                             suptitle(textHere);
% % % %                             set(102,'color','white'),
                        end
                        % title
% % % %                         title('test data: progression curve')
                        box on

                    end
                    

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
                    title('PEP')
                    subplot(1,2,2)
                    for i = DFstudy
                        plot(data.tempTime{i}, simNADH{i,j},'-')
%                         plot(simRes{j}.t, simRes{j}.v, '-')
                        hold on
                        plot(data.tempTime{i}, expNADH{i,j},'k+')
                        hold on
                    end
                    title('v_{apparent.GAPDHrev}')
                    % % % % ONLY FOR ENO
                    suptitle(erase(sprintf('pH = %d',setup.fullpHarray(j)),"0000e+00"));
                    % % % % ONLY FOR ENO
                end
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
                                H_effect = zeros(1:10);
                                for j = 1:numpH
                                    H_effect(j) = 10^(setup.pH_vals(j) - setup.pH_vals(6));
                                end
                                errorHaldane1 = (Keq(1) - (vmf1 * kp1 * kp2 * H_effect(1)) / (vmr1 * ks1 * ks2) );
                                errorHaldane2 = (Keq(2) - (vmf2 * kp1 * kp2 * H_effect(2)) / (vmr2 * ks1 * ks2) );
                                errorHaldane3 = (Keq(3) - (vmf3 * kp1 * kp2 * H_effect(3)) / (vmr3 * ks1 * ks2) );
                                errorHaldane4 = (Keq(4) - (vmf4 * kp1 * kp2 * H_effect(4)) / (vmr4 * ks1 * ks2) );
                                errorHaldane5 = (Keq(5) - (vmf5 * kp1 * kp2 * H_effect(5)) / (vmr5 * ks1 * ks2) );
                                errorHaldane6 = (Keq(6) - (vmf6 * kp1 * kp2 * H_effect(6)) / (vmr6 * ks1 * ks2) );
                                errorHaldane7 = (Keq(7) - (vmf7 * kp1 * kp2 * H_effect(7)) / (vmr7 * ks1 * ks2) );
                                errorHaldane8 = (Keq(8) - (vmf8 * kp1 * kp2 * H_effect(8)) / (vmr8 * ks1 * ks2) );
                                errorHaldane9 = (Keq(9) - (vmf9 * kp1 * kp2 * H_effect(9)) / (vmr9 * ks1 * ks2) );
                                errorHaldane10 =(Keq(10) - (vmf10 * kp1 * kp2 * H_effect(10)) / (vmr10 * ks1 * ks2) );
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
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:4)';
                    
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
                    errorNADH9_df1 = simNADH{4,9} - expNADH{4,9};
                    errorNADH10_df1 = simNADH{4,10} - expNADH{4,10};
                        % DF2
                        errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_df2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_df2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_df2 = simNADH{3,10} - expNADH{3,10};
                            % DF4
                            errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
                            errorNADH8_df4 = simNADH{2,8} - expNADH{2,8};
                            errorNADH9_df4 = simNADH{2,9} - expNADH{2,9};
                            errorNADH10_df4 = simNADH{2,10} - expNADH{2,10};
                                % DF8
                                errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
                                errorNADH8_df8 = simNADH{1,8} - expNADH{1,8};
                                errorNADH9_df8 = simNADH{1,9} - expNADH{1,9};
                                errorNADH10_df8 = simNADH{1,10} - expNADH{1,10};
                wDesp_t = wDesp';
                errorNADH = [...
                    wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                    wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                    wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                    wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                    wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                    wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                    wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1;
                    wDesp_t(1,8)*errorNADH8_df8; wDesp_t(2,8)*errorNADH8_df4; wDesp_t(3,8)*errorNADH8_df2; wDesp_t(4,8)*errorNADH8_df1;
                    wDesp_t(1,9)*errorNADH9_df8; wDesp_t(2,9)*errorNADH9_df4; wDesp_t(3,9)*errorNADH9_df2; wDesp_t(4,9)*errorNADH9_df1;
                    wDesp_t(1,10)*errorNADH10_df8; wDesp_t(2,10)*errorNADH10_df4; wDesp_t(3,10)*errorNADH10_df2; wDesp_t(4,10)*errorNADH10_df1];
% % % %                     errorReg = lambda * x_temp';
%                     errorReg = lambda * x_temp(1:4)'; 
                    caseKm = setup.caseKm;
                    switch caseKm
                        case 'pH_independent'
                            errorReg = lambda * x_temp(1:4)';
                        case 'pH_dependent'
                            errorReg = lambda * x_temp(1:40)';
                    end               
            otherwise
                disp('No specific cost function has been appointed');
        end
        error = [
            wD * errorNADH;
            wH * errorHaldane;
            wL * errorReg;
            ];        
%         disp(lambda);
    
    case 'gapdh'
        caseKm = setup.caseKm;

        simNADH = cell(DFs,numpH);
        expNADH = cell(DFs,numpH);
        simGAPDHr = cell(DFs,numpH);
        expGAPDHr = cell(DFs,numpH);
        temp_simResult = cell(DFs,numpH);
        
        % simulations loop for each pH value
        for j = 1:numpH
            % select required data
            data.Vmaxs = data.Vmax(j,:);
            data.NADH = data.conc_mean(j,:);
            data.Vprofs = data.RRs(j,:);
            data.tempTime = data.time(j,:);
            % inputs to be selected
%             data.KeqGAPDH = setup.pH_Keq_gapdh_eQ; %setup.pH_Keq_gapdh(i);
%             data.KeqPGK = setup.pH_Keq_pgk;
            data.KeqGAPDH = setup.pH_Keq_gapdh_eQ_fwd; %setup.pH_Keq_gapdh(i);
            data.KeqPGK = setup.pH_Keq_pgk_fwd;
            data.chosenKeqGAPDH = data.KeqGAPDH(j);
            data.chosenKeqPGK = data.KeqPGK(j);
            
            % selecting the right parameters
            switch caseKm
                case 'pH_independent'
                    % Mode 1. still keeping kms, but fixed
                    idx5 = 4 + j;
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(idx5);
                case 'pH_dependent'
                    idx1 = j + numpH * 0;
                    idx2 = j + numpH * 1;
                    idx3 = j + numpH * 2;
                    idx4 = j + numpH * 3;
                    idx5 = j + numpH * 4;
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(idx1);
                    xassay(2) = x_temp(idx2);
                    xassay(3) = x_temp(idx3);
                    xassay(4) = x_temp(idx4);
                    xassay(5) = x_temp(idx5);
                otherwise
                    disp('Warning: no specification has been made on Km being pH dependent or independent');
            end
            
            % simulations
            for i = DFstudy
                % recall vmax for the specific value and simulate
%                 data.chosenVmax = data.Vmaxs(1,4)/data.DF(1,i); % vmax from the highest DF is taken and then divided
                data.chosenVmax = data.Vmax(j,4)/data.DF(j,i); % vmax from the highest DF is taken and then divided
                
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
                expGAPDHr{i,j} = data.Vprofs{i};
                simResult.vObs = vObs;
                temp_simResult{i,j} = simResult;
            end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
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
% % % %                         plot(data.time{j,i}, simNADH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.time{j,i}, expNADH{i,j},'k.','MarkerSize',4)
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('NADH concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102)
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('GAPDHfwd reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end       
% % % %                     
                elseif((simAllProfiles == 0)&&(setup.plotEachSim == 1))
                    figure
                    subplot(1,2,1)
                    for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-')
                        plot(data.time{j,i}, simNADH{i,j},'-')
%                         plot(simRes{j}.t, simRes{j}.y(:,8),'-')
                        hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k+')
                        plot(data.time{j,i}, expNADH{i,j},'k+')
                        hold on
                    end
                    title('NADH')
                    subplot(1,2,2)
                    for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-')
                        plot(data.time{j,i}, simNADH{i,j},'-')
%                         plot(simRes{j}.t, simRes{j}.v, '-')
                        hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k+')
                        plot(data.time{j,i}, expGAPDHr{i,j},'k+')
                        hold on
                    end
                    title('v_{apparent.GAPDHfwd}')
                    % % % % ONLY FOR ENO
                    suptitle(erase(sprintf('pH = %d',setup.fullpHarray(j)),"0000e+00"));
                    % % % % ONLY FOR ENO
                end
            end
        end

%         %%
%         tempTitle = {'P3G','ATP','BPG','ADP','NAD','GAP','PHOS','NADH'};
%         
%         figure,
%         for i = 1:8
%             subplot(3,3,i)
%             plot(temp_simResult{2,5}.t, temp_simResult{2,5}.y(:,i),'b.-')
%             hold on
%             plot(temp_simResult{2,6}.t, temp_simResult{2,6}.y(:,i),'r.-')
%             title(tempTitle{i})
%         end
%         subplot(3,3,9)
%         plot(temp_simResult{2,5}.t, temp_simResult{2,5}.vObs,'b.-')
%         hold on
%         plot(temp_simResult{2,6}.t, temp_simResult{2,6}.vObs,'r.-')
%         title('v_{obs}')
%         %%      
        
        % calculation cost function
        switch costfun
            case 1 % DF1
%                     errorNADH1 = simNADH{4,1} - expNADH{4,1};
%                     errorNADH2 = simNADH{4,2} - expNADH{4,2};
%                     errorNADH3 = simNADH{4,3} - expNADH{4,3};
%                     errorNADH4 = simNADH{4,4} - expNADH{4,4};
%                     errorNADH5 = simNADH{4,5} - expNADH{4,5};
%                     errorNADH6 = simNADH{4,6} - expNADH{4,6};
%                     errorNADH7 = simNADH{4,7} - expNADH{4,7};
%                     errorNADH8 = simNADH{4,8} - expNADH{4,8};
%                     errorNADH9 = simNADH{4,9} - expNADH{4,9};
%                     errorNADH10 = simNADH{4,10} - expNADH{4,10};
%                     errorNADH = [...
%                         wDesp(1)*errorNADH1;
%                         wDesp(2)*errorNADH2;
%                         wDesp(3)*errorNADH3;
%                         wDesp(4)*errorNADH4;
%                         wDesp(5)*errorNADH5;
%                         wDesp(6)*errorNADH6;
%                         wDesp(7)*errorNADH7;
%                         wDesp(8)*errorNADH8;
%                         wDesp(9)*errorNADH9;
%                         wDesp(10)*errorNADH10];
%                     for temp1 = 1                    
%                         Keq = setup.pH_Keq_gapdh_eQ; %[]
% %                         Keq = data.chosenKeqGAPDH; %[]
%                         switch sourceVm
%                             case 'experimentalSlopesFixed'
%                                 vmf1 = 10.^x_temp(5).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr1 = 10.^x_temp(6).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf2 = 10.^x_temp(7).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr2 = 10.^x_temp(8).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf3 = 10.^x_temp(9).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr3 = 10.^x_temp(10).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf4 = 10.^x_temp(11).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr4 = 10.^x_temp(12).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf5 = 10.^x_temp(13).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr5 = 10.^x_temp(14).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf6 = 10.^x_temp(15).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr6 = 10.^x_temp(16).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf7 = 10.^x_temp(17).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr7 = 10.^x_temp(18).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf8 = 10.^x_temp(19).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr8 = 10.^x_temp(20).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf9 = 10.^x_temp(21).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr9 = 10.^x_temp(22).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                                 vmf10 = 10.^x_temp(23).*setup.exp_vmax_gapdhf(6);% mM s^{-1}
%                                 vmr10 = 10.^x_temp(24).*setup.exp_vmax_gapdhr(6); % mM s^{-1}
%                             otherwise
%                                 disp('No source for vmax has been selected');
%                         end
%                         ks1 = 10 .^ x_temp(1) .* 2.48; % mM %k_gap
%                         ks2 = 10 .^ x_temp(3) .* 2.92; %mM %k_nad
%                         kp1 = 10 .^ x_temp(2) .* 1.18; % mM %k_bpg
%                         kp2 = 10 .^ x_temp(4) .* 0.022; % mM %k_nadh
%                         switch ode_pH
%                             case 'on'
%                                 H_effect = zeros(1:10);
%                                 for j = 1:numpH
%                                     H_effect(j) = 10^(setup.pH_vals(j) - setup.pH_vals(6));
%                                 end
%                                 errorHaldane1 = (Keq(1) - (vmf1 * kp1 * kp2 * H_effect(1)) / (vmr1 * ks1 * ks2) );
%                                 errorHaldane2 = (Keq(2) - (vmf2 * kp1 * kp2 * H_effect(2)) / (vmr2 * ks1 * ks2) );
%                                 errorHaldane3 = (Keq(3) - (vmf3 * kp1 * kp2 * H_effect(3)) / (vmr3 * ks1 * ks2) );
%                                 errorHaldane4 = (Keq(4) - (vmf4 * kp1 * kp2 * H_effect(4)) / (vmr4 * ks1 * ks2) );
%                                 errorHaldane5 = (Keq(5) - (vmf5 * kp1 * kp2 * H_effect(5)) / (vmr5 * ks1 * ks2) );
%                                 errorHaldane6 = (Keq(6) - (vmf6 * kp1 * kp2 * H_effect(6)) / (vmr6 * ks1 * ks2) );
%                                 errorHaldane7 = (Keq(7) - (vmf7 * kp1 * kp2 * H_effect(7)) / (vmr7 * ks1 * ks2) );
%                                 errorHaldane8 = (Keq(8) - (vmf8 * kp1 * kp2 * H_effect(8)) / (vmr8 * ks1 * ks2) );
%                                 errorHaldane9 = (Keq(9) - (vmf9 * kp1 * kp2 * H_effect(9)) / (vmr9 * ks1 * ks2) );
%                                 errorHaldane10 =(Keq(10) - (vmf10 * kp1 * kp2 * H_effect(10)) / (vmr10 * ks1 * ks2) );
%                                 errorHaldane = [...
%                                     errorHaldane1;
%                                     errorHaldane2;
%                                     errorHaldane3;
%                                     errorHaldane4;
%                                     errorHaldane5;
%                                     errorHaldane6;
%                                     errorHaldane7;
%                                     errorHaldane8;
%                                     errorHaldane9;
%                                     errorHaldane10];
%                         %         errorHaldane = sum(abs(errorHaldane));
%                             otherwise
%                         % % % %         errorHaldane = wH * (Keq - (vmf * kp1 * kp2) / (vmr * ks1 * ks2) );
%                                 disp('No source for vmax has been selected');
%                         end
%                     end
%                     errorReg = lambda * x_temp';
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
                    errorNADH9_df1 = simNADH{4,9} - expNADH{4,9};
                    errorNADH10_df1 = simNADH{4,10} - expNADH{4,10};
                    errorNADH11_df1 = simNADH{4,11} - expNADH{4,11};
                    errorNADH12_df1 = simNADH{4,12} - expNADH{4,12};
                        % DF2
                        errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_df2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_df2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_df2 = simNADH{3,10} - expNADH{3,10};
                        errorNADH11_df2 = simNADH{3,11} - expNADH{3,11};
                        errorNADH12_df2 = simNADH{3,12} - expNADH{3,12};
                            % DF4
                            errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
                            errorNADH8_df4 = simNADH{2,8} - expNADH{2,8};
                            errorNADH9_df4 = simNADH{2,9} - expNADH{2,9};
                            errorNADH10_df4 = simNADH{2,10} - expNADH{2,10};
                            errorNADH11_df4 = simNADH{2,11} - expNADH{2,11};
                            errorNADH12_df4 = simNADH{2,12} - expNADH{2,12};
                                % DF8
                                errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
                                errorNADH8_df8 = simNADH{1,8} - expNADH{1,8};
                                errorNADH9_df8 = simNADH{1,9} - expNADH{1,9};
                                errorNADH10_df8 = simNADH{1,10} - expNADH{1,10};
                                errorNADH11_df8 = simNADH{1,11} - expNADH{1,11};
                                errorNADH12_df8 = simNADH{1,12} - expNADH{1,12};
                wDesp_t = wDesp';
                errorNADH = [...
                    wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                    wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                    wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                    wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                    wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                    wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                    wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1;
                    wDesp_t(1,8)*errorNADH8_df8; wDesp_t(2,8)*errorNADH8_df4; wDesp_t(3,8)*errorNADH8_df2; wDesp_t(4,8)*errorNADH8_df1;
                    wDesp_t(1,9)*errorNADH9_df8; wDesp_t(2,9)*errorNADH9_df4; wDesp_t(3,9)*errorNADH9_df2; wDesp_t(4,9)*errorNADH9_df1;
                    wDesp_t(1,10)*errorNADH10_df8; wDesp_t(2,10)*errorNADH10_df4; wDesp_t(3,10)*errorNADH10_df2; wDesp_t(4,10)*errorNADH10_df1;
                    wDesp_t(1,11)*errorNADH11_df8; wDesp_t(2,11)*errorNADH11_df4; wDesp_t(3,11)*errorNADH11_df2; wDesp_t(4,11)*errorNADH11_df1;
                    wDesp_t(1,12)*errorNADH12_df8; wDesp_t(2,12)*errorNADH12_df4; wDesp_t(3,12)*errorNADH12_df2; wDesp_t(4,12)*errorNADH12_df1];
%                 errorNADH = [...
%                     wDesp_t(1,1)*errorNADH1_df8(1:21); wDesp_t(2,1)*errorNADH1_df4(1:21); wDesp_t(3,1)*errorNADH1_df2(1:21); wDesp_t(4,1)*errorNADH1_df1(1:21); 
%                     wDesp_t(1,2)*errorNADH2_df8(1:21); wDesp_t(2,2)*errorNADH2_df4(1:21); wDesp_t(3,2)*errorNADH2_df2(1:21); wDesp_t(4,2)*errorNADH2_df1(1:21);
%                     wDesp_t(1,3)*errorNADH3_df8(1:21); wDesp_t(2,3)*errorNADH3_df4(1:21); wDesp_t(3,3)*errorNADH3_df2(1:21); wDesp_t(4,3)*errorNADH3_df1(1:21);
%                     wDesp_t(1,4)*errorNADH4_df8(1:21); wDesp_t(2,4)*errorNADH4_df4(1:21); wDesp_t(3,4)*errorNADH4_df2(1:21); wDesp_t(4,4)*errorNADH4_df1(1:21);
%                     wDesp_t(1,5)*errorNADH5_df8(1:21); wDesp_t(2,5)*errorNADH5_df4(1:21); wDesp_t(3,5)*errorNADH5_df2(1:21); wDesp_t(4,5)*errorNADH5_df1(1:21);
%                     wDesp_t(1,6)*errorNADH6_df8(1:21); wDesp_t(2,6)*errorNADH6_df4(1:21); wDesp_t(3,6)*errorNADH6_df2(1:21); wDesp_t(4,6)*errorNADH6_df1(1:21);
%                     wDesp_t(1,7)*errorNADH7_df8(1:21); wDesp_t(2,7)*errorNADH7_df4(1:21); wDesp_t(3,7)*errorNADH7_df2(1:21); wDesp_t(4,7)*errorNADH7_df1(1:21);
%                     wDesp_t(1,8)*errorNADH8_df8(1:21); wDesp_t(2,8)*errorNADH8_df4(1:21); wDesp_t(3,8)*errorNADH8_df2(1:21); wDesp_t(4,8)*errorNADH8_df1(1:21);
%                     wDesp_t(1,9)*errorNADH9_df8(1:21); wDesp_t(2,9)*errorNADH9_df4(1:21); wDesp_t(3,9)*errorNADH9_df2(1:21); wDesp_t(4,9)*errorNADH9_df1(1:21);
%                     wDesp_t(1,10)*errorNADH10_df8(1:21); wDesp_t(2,10)*errorNADH10_df4(1:21); wDesp_t(3,10)*errorNADH10_df2(1:21); wDesp_t(4,10)*errorNADH10_df1(1:21);
%                     wDesp_t(1,11)*errorNADH11_df8(1:21); wDesp_t(2,11)*errorNADH11_df4(1:21); wDesp_t(3,11)*errorNADH11_df2(1:21); wDesp_t(4,11)*errorNADH11_df1(1:21);
%                     wDesp_t(1,12)*errorNADH12_df8(1:21); wDesp_t(2,12)*errorNADH12_df4(1:21); wDesp_t(3,12)*errorNADH12_df2(1:21); wDesp_t(4,12)*errorNADH12_df1(1:21)];
% % % %                 errorReg = lambda * x_temp';
%                 errorReg = lambda * x_temp(1:4)';
                caseKm = setup.caseKm;
                switch caseKm
                    case 'pH_independent'
                        errorReg = lambda * x_temp(1:4)';
                    case 'pH_dependent'
                        errorReg = lambda * x_temp(1:48)';
                end
                
            otherwise
                disp('No specific cost function has been appointed');
        end
        error = [
            wD * errorNADH;
            wH * errorHaldane;
            wL * errorReg;
            ];        
%         disp(lambda);

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
            switch j
                case 1
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                case 2
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(4);
                case 3
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(5);
                case 4
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(6);
                case 5
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(7);
                case 6
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(8);
                case 7
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(9);
                case 8
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(10);
                case 9
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(11);
                case 10
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(12);
                case 11
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(13);
                case 12
                    xassay = zeros(1,3);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(14);
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
        
%         for j = 1:numpH
        for j = 12
            if plotEachSimCF == 1
                if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
% % % %                     simulationVisualization;
                    for ENOsim = 1
    %                     %% Simulation concentrations
                        % selection of the timeframe
                        data.tempTime = data.time(j,:);
                        % plottinh
                        if j == 12
% % % %                             h101 = figure(101);
                            h101 = figure(1000);
                            hold on
                        end
                        set(0,'CurrentFigure', h101);
% % % %                         p101 = subplot(3,4,j);
                        if setup.caseStudyENO == 1
                            p101 = subplot(3,3,2);
                            hold on
                        elseif setup.caseStudyGAPDHr == 1
                            p101 = subplot(3,3,5);
                            hold on
                        elseif setup.caseStudyPGM == 1
                            p101 = subplot(3,3,8);
                            hold on
                        end
                        
% % % %                         if setup.literatureValues == 0
% % % %                             p101 = subplot(3,3,2);
% % % %                         elseif setup.literatureValues == 1
% % % %                             p101 = subplot(3,3,5);
% % % %                         end
                        
                        for i = DFstudy
                            if setup.caseStudyENO == 1
% % % %                                 plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
% % % %                                 ylim([0 1.2])
                                if setup.literatureValues == 0
%                                     plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',setup.colLightSkyBlue)
                                    plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}(linspace(1,121,16)), expPEP{i,j}(linspace(1,121,16)),'k.','MarkerSize',4)
                                    ylim([0 1.2])
                                    hold on
                                elseif setup.literatureValues == 1
%                                     plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',setup.col4)
                                    plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}(linspace(1,121,16)), expPEP{i,j}(linspace(1,121,16)),'k.','MarkerSize',4)
                                    ylim([0 1.2])
                                    hold on
                                end
                                
                                
                            elseif setup.caseStudyHXK == 1
                                plot(data.tempTime{i}, simNADPH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADPH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                            elseif setup.caseStudyGAPDH == 1
                        %         plot(data.time{j,i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                        %         hold on
                        %         plot(data.time{j,i}, expNADH{i,j},'k.','MarkerSize',4)
                        %         ylim([0 0.15])
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                        %         disp(sprintf('j=%d, i=%d',j,i));
                            elseif setup.caseStudyGAPDHr == 1
                                
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                end
                                
                                                                    
                            else
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                                if(setup.caseStudyPDC == 1)
                                    xlim([0 300])
                                elseif(setup.caseStudyPGM == 1)
                                    ylim([0 0.10])
                                    xlim([0 300])
                                elseif((setup.caseStudyPGI == 1)||(setup.caseStudyPFK == 1))
                                    ylim([0 0.15])
                                    xlim([0 1000])
                                elseif(setup.caseStudyTPI == 1)
                                    ylim([0 0.15])
                                    xlim([0 600])
                                end        
                            end
                        end
                        % ylabel left
                        if((j == 1)||(j == 5)||(j == 9))
                            if setup.caseStudyENO == 1
                                ylabel('PEP concentration [mM]')
                            elseif setup.caseStudyHXK == 1
                                ylabel('NAPDH concentration [mM]')
                            else 
                                ylabel('NADH concentration [mM]')        
                            end
                        end
                        % xlabel
                        if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
                            i_xlabel = 5:8;
                        else
                            i_xlabel = 9:12;
                        end
                        if((j == i_xlabel(1))||(j == i_xlabel(2))||(j == i_xlabel(3))||(j == i_xlabel(4)))
%                             xlabel('assay time [s]')
                        end
                        % pH text box
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         
                        % Only the first time that the function is called
% % % %                         temp = get(p101,'Children');
% % % %                         temp2 = findobj(temp, 'Type','text');
% % % %                         if sum(abs(size(temp2))) ~= 0%exist('temp2','var')
    %                             temp2 = temp(u);
    %                         end

                            tempText = erase(sprintf('@pH %d', data.pH(j,1)),"0000e+00");
    %                         text(30, p101.YLim(2)*0.9, tempText);
    %                         if i == 1 
%                                 text(600*0.95, p101.YLim(2)*0.9, tempText,...
%                                     'HorizontalAlignment', 'right')
    %                         end
                            % suptitle
                            if j == numpH
    % % % %                             textHere_pre = [setup.enzymeName,': concentration progression curve. Each box contains a different pH value.'];
    % % % %                             textHere = {textHere_pre;...
    % % % %                                 'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
    % % % %                             suptitle(textHere)
    % % % %                             set(101,'color','white'),
                            end
                            % title
    %                         title({'training data:';'concentrations vs time';'(mM, s)'}, ...%)
    %                             'FontSize',10)
    %                         if i == 1
%                                 text(300, p101.YLim(2) + (p101.YLim(2) - p101.YLim(1)) * 0.3, ...
%                                     {'training data:';'concentrations vs time';'(mM, s)'}, ...
%                                     'FontSize', 10, ...
%                                     'HorizontalAlignment', 'center')
    %                         end
% % % %                         end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                             

                        box on

    %                     %% Simulation reaction rates
                        if j == 12
% % % %                             h102 = figure(102);
                            h102 = figure(1000);
                            hold on
                        end
                        set(0,'CurrentFigure', h102);
% % % %                         p102 = subplot(3,4,j);
% % % %                         p102 = subplot(3,3,3);
% % % %                         hold on
                        
                        if setup.caseStudyENO == 1
                            p102 = subplot(3,3,3);
                            hold on
                        elseif setup.caseStudyGAPDHr == 1
                            p102 = subplot(3,3,6);
                            hold on
                        elseif setup.caseStudyPGM == 1
                            p102 = subplot(3,3,9);
                            hold on
                        end
                        
                        
% % % %                         if setup.literatureValues == 0
% % % %                             p102 = subplot(3,3,3);
% % % %                         elseif setup.literatureValues == 1
% % % %                             p102 = subplot(3,3,6);
% % % %                         end
                        
                        for i = DFstudy
                            if setup.caseStudyALD == 1
                                plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyENO == 1
% % % %                                 plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
% % % %                                 ylim([0 1.5E-3])
% % % %                                 hold on                               
                                
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.5E-3])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color', setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.5E-3])
                                    hold on
                                end
                                
                            elseif setup.caseStudyHXK == 1
                                plot(data.tempTime{i}, simG6PDH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expG6PDH{i,j},'k.','MarkerSize',4)
                                ylim([0 5E-4])
                            elseif setup.caseStudyPDC == 1
                                plot(data.tempTime{i}, simPDC{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expPDC{i,j},'k.','MarkerSize',4)
                                xlim([0 300])
                            elseif setup.caseStudyPYK == 1  
                                plot(data.tempTime{i}, simLDH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expLDH{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyPGM == 1  
                                plot(data.tempTime{i}, simPGM{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expPGM{i,j},'k.','MarkerSize',4)
                                xlim([0 300])
                            elseif setup.caseStudyPGI == 1
                                plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
                                xlim([0 1000])
                            elseif setup.caseStudyTPI == 1
                                plot(data.tempTime{i}, simTPI{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expTPI{i,j},'k.','MarkerSize',4)
                                xlim([0 600])
                            elseif setup.caseStudyPFK == 1
                                plot(data.tempTime{i}, simPFK{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expPFK{i,j},'k.','MarkerSize',4)
                                xlim([0 1000])
                            elseif setup.caseStudyGAPDH == 1
                                plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyGAPDHr == 1
% % % %                                 plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                
                                if setup.literatureValues == 0
                                    plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                end
                                
                            end
                        end
                        % ylabel left
                        if((j == 1)||(j == 5)||(j == 9))
                            if((setup.caseStudyALD == 1)||(setup.caseStudyPGI == 1)||(setup.caseStudyPFK == 1)||(setup.caseStudyTPI == 1))
                                ylabel('GPD reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyENO == 1
                                ylabel('ENO reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyHXK == 1
                                ylabel('G6PDH reaction rate [mM s^{-1}]')
                            elseif((setup.caseStudyPYK == 1)||(setup.caseStudyPDC == 1)||(setup.caseStudyPGM == 1))
                                ylabel('LDH reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyGAPDH== 1
                                ylabel('GAPDH_{fwd} reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyGAPDHr== 1
                                ylabel('GAPDH_{rev} reaction rate [mM s^{-1}]')
                            end
                        end
                        % xlabel
                        if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
                            i_xlabel = 5:8;
                        else
                            i_xlabel = 9:12;
                        end
                        if((j == i_xlabel(1))||(j == i_xlabel(2))||(j == i_xlabel(3))||(j == i_xlabel(4)))
%                             xlabel('assay time [s]')
                        end
                        % textbox
                        tempText = erase(sprintf('@pH %d', data.pH(j,1)),"0000e+00");
%                         if i == 1
%                             text(600*0.95, p102.YLim(2)*0.9, tempText,...
%                                 'HorizontalAlignment', 'right')
%                         end
                        % suptitle
                        if j == numpH
% % % %                             textHere_pre = [setup.enzymeName,': reaction rate progression curve. Each box contains a different pH value.'];
% % % %                             textHere = {textHere_pre;...
% % % %                                 'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
% % % %                             suptitle(textHere);
% % % %                             set(102,'color','white'),
                        end
                        % title
%                         title({'test data:';'reaction rates vs time';'(mM s^{-1}, s)'},...%)
%                             'FontSize',10)
%                         if i == 1
%                             text(300, p102.YLim(2) + (p102.YLim(2) - p102.YLim(1)) * 0.3, ...
%                                 {'test data:';'reaction rates vs time'; '(mM s^{-1}, s)'}, ...
%                                 'FontSize', 10, ...
%                                 'HorizontalAlignment', 'center')
%                         end
%                         title('test data: progression curve')
                        box on

                    end
                    
                    
                    
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
    case 'hxk'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = data.Vmax(j,:);
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
            switch j
                case 1
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                case 2
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(6);
                case 3
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(7);
                case 4
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(8);
                case 5
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(9);
                case 6
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(10);
                case 7
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(11);
                case 8
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(12);
                case 9
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(13);
                case 10
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(14);
                case 11
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(15);
                case 12
                    xassay = zeros(1,5);
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(16);
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
                    simulationVisualization;
% % % %                     
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADPH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expNADPH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 0.15])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('NADPH concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simG6PDH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expG6PDH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 5E-4])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('G6PDH reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end       
% % % %                     
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
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:4)';
            case 3 % selected DFs
                    errorHaldane = 0;
                    % DF1
                    errorNADPH1_df1 = simNADPH{4,1} - expNADPH{4,1};
                    errorNADPH2_df1 = simNADPH{4,2} - expNADPH{4,2};
                    errorNADPH3_df1 = simNADPH{4,3} - expNADPH{4,3};
                    errorNADPH4_df1 = simNADPH{4,4} - expNADPH{4,4};
                    errorNADPH5_df1 = simNADPH{4,5} - expNADPH{4,5};
                    errorNADPH6_df1 = simNADPH{4,6} - expNADPH{4,6};
                    errorNADPH7_df1 = simNADPH{4,7} - expNADPH{4,7};
                    errorNADPH8_df1 = simNADPH{4,8} - expNADPH{4,8};
                    errorNADPH9_df1 = simNADPH{4,9} - expNADPH{4,9};
                    errorNADPH10_df1 = simNADPH{4,10} - expNADPH{4,10};
                    errorNADPH11_df1 = simNADPH{4,11} - expNADPH{4,11};
                    errorNADPH12_df1 = simNADPH{4,12} - expNADPH{4,12};
                        % DF2
                        errorNADPH1_df2 = simNADPH{3,1} - expNADPH{3,1};
                        errorNADPH2_df2 = simNADPH{3,2} - expNADPH{3,2};
                        errorNADPH3_df2 = simNADPH{3,3} - expNADPH{3,3};
                        errorNADPH4_df2 = simNADPH{3,4} - expNADPH{3,4};
                        errorNADPH5_df2 = simNADPH{3,5} - expNADPH{3,5};
                        errorNADPH6_df2 = simNADPH{3,6} - expNADPH{3,6};
                        errorNADPH7_df2 = simNADPH{3,7} - expNADPH{3,7};
                        errorNADPH8_df2 = simNADPH{3,8} - expNADPH{3,8};
                        errorNADPH9_df2 = simNADPH{3,9} - expNADPH{3,9};
                        errorNADPH10_df2 = simNADPH{3,10} - expNADPH{3,10};
                        errorNADPH11_df2 = simNADPH{3,11} - expNADPH{3,11};
                        errorNADPH12_df2 = simNADPH{3,12} - expNADPH{3,12};
                            % DF4
                            errorNADPH1_df4 = simNADPH{2,1} - expNADPH{2,1};
                            errorNADPH2_df4 = simNADPH{2,2} - expNADPH{2,2};
                            errorNADPH3_df4 = simNADPH{2,3} - expNADPH{2,3};
                            errorNADPH4_df4 = simNADPH{2,4} - expNADPH{2,4};
                            errorNADPH5_df4 = simNADPH{2,5} - expNADPH{2,5};
                            errorNADPH6_df4 = simNADPH{2,6} - expNADPH{2,6};
                            errorNADPH7_df4 = simNADPH{2,7} - expNADPH{2,7};
                            errorNADPH8_df4 = simNADPH{2,8} - expNADPH{2,8};
                            errorNADPH9_df4 = simNADPH{2,9} - expNADPH{2,9};
                            errorNADPH10_df4 = simNADPH{2,10} - expNADPH{2,10};
                            errorNADPH11_df4 = simNADPH{2,11} - expNADPH{2,11};
                            errorNADPH12_df4 = simNADPH{2,12} - expNADPH{2,12};
                                % DF8
                                errorNADPH1_df8 = simNADPH{1,1} - expNADPH{1,1};
                                errorNADPH2_df8 = simNADPH{1,2} - expNADPH{1,2};
                                errorNADPH3_df8 = simNADPH{1,3} - expNADPH{1,3};
                                errorNADPH4_df8 = simNADPH{1,4} - expNADPH{1,4};
                                errorNADPH5_df8 = simNADPH{1,5} - expNADPH{1,5};
                                errorNADPH6_df8 = simNADPH{1,6} - expNADPH{1,6};
                                errorNADPH7_df8 = simNADPH{1,7} - expNADPH{1,7};
                                errorNADPH8_df8 = simNADPH{1,8} - expNADPH{1,8};
                                errorNADPH9_df8 = simNADPH{1,9} - expNADPH{1,9};
                                errorNADPH10_df8 = simNADPH{1,10} - expNADPH{1,10};
                                errorNADPH11_df8 = simNADPH{1,11} - expNADPH{1,11};
                                errorNADPH12_df8 = simNADPH{1,12} - expNADPH{1,12};
                    wDesp_t = wDesp';
                    errorNADPH = [...
                        wDesp_t(1,1)*errorNADPH1_df8; wDesp_t(2,1)*errorNADPH1_df4; wDesp_t(3,1)*errorNADPH1_df2; wDesp_t(4,1)*errorNADPH1_df1; 
                        wDesp_t(1,2)*errorNADPH2_df8; wDesp_t(2,2)*errorNADPH2_df4; wDesp_t(3,2)*errorNADPH2_df2; wDesp_t(4,2)*errorNADPH2_df1;
                        wDesp_t(1,3)*errorNADPH3_df8; wDesp_t(2,3)*errorNADPH3_df4; wDesp_t(3,3)*errorNADPH3_df2; wDesp_t(4,3)*errorNADPH3_df1;
                        wDesp_t(1,4)*errorNADPH4_df8; wDesp_t(2,4)*errorNADPH4_df4; wDesp_t(3,4)*errorNADPH4_df2; wDesp_t(4,4)*errorNADPH4_df1;
                        wDesp_t(1,5)*errorNADPH5_df8; wDesp_t(2,5)*errorNADPH5_df4; wDesp_t(3,5)*errorNADPH5_df2; wDesp_t(4,5)*errorNADPH5_df1;
                        wDesp_t(1,6)*errorNADPH6_df8; wDesp_t(2,6)*errorNADPH6_df4; wDesp_t(3,6)*errorNADPH6_df2; wDesp_t(4,6)*errorNADPH6_df1;
                        wDesp_t(1,7)*errorNADPH7_df8; wDesp_t(2,7)*errorNADPH7_df4; wDesp_t(3,7)*errorNADPH7_df2; wDesp_t(4,7)*errorNADPH7_df1;
                        wDesp_t(1,8)*errorNADPH8_df8; wDesp_t(2,8)*errorNADPH8_df4; wDesp_t(3,8)*errorNADPH8_df2; wDesp_t(4,8)*errorNADPH8_df1;
                        wDesp_t(1,9)*errorNADPH9_df8; wDesp_t(2,9)*errorNADPH9_df4; wDesp_t(3,9)*errorNADPH9_df2; wDesp_t(4,9)*errorNADPH9_df1;
                        wDesp_t(1,10)*errorNADPH10_df8; wDesp_t(2,10)*errorNADPH10_df4; wDesp_t(3,10)*errorNADPH10_df2; wDesp_t(4,10)*errorNADPH10_df1;
                        wDesp_t(1,11)*errorNADPH11_df8; wDesp_t(2,11)*errorNADPH11_df4; wDesp_t(3,11)*errorNADPH11_df2; wDesp_t(4,11)*errorNADPH11_df1;
                        wDesp_t(1,12)*errorNADPH12_df8; wDesp_t(2,12)*errorNADPH12_df4; wDesp_t(3,12)*errorNADPH12_df2; wDesp_t(4,12)*errorNADPH12_df1];
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:4)';
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
% % % %                     xassay(5) = x_temp(12);
% % % %                     xassay(6) = x_temp(13);
                case 2
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(5);
% % % %                     xassay(5) = x_temp(12);
% % % %                     xassay(6) = x_temp(13);
                case 3
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(6);
% % % %                     xassay(5) = x_temp(12);
% % % %                     xassay(6) = x_temp(13);
                case 4
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(7);
% % % %                     xassay(5) = x_temp(12);
% % % %                     xassay(6) = x_temp(13);
                case 5
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(8);
% % % %                     xassay(5) = x_temp(12);
% % % %                     xassay(6) = x_temp(13);
                case 6
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(9);
% % % %                     xassay(5) = x_temp(12);
% % % %                     xassay(6) = x_temp(13);
                case 7
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(10);
% % % %                     xassay(5) = x_temp(12);
% % % %                     xassay(6) = x_temp(13);
                case 8
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(11);
% % % %                     xassay(5) = x_temp(12);
% % % %                     xassay(6) = x_temp(13);
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
                    simulationVisualization;
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101);
% % % %                     p101 = subplot(3,4,j);
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 0.15])
% % % %                     end
% % % %                     % ylabel left
% % % %                     if((j == 1)||(j == 5))
% % % %                         ylabel('NADH concentration [mM]')
% % % %                     end
% % % %                     % xlabel
% % % %                     if((j == 5)||(j == 6)||(j == 7)||(j == 8))
% % % %                         xlabel('assay time [s]')
% % % %                     end
% % % %                     % pH text box
% % % %                     tempText = erase(sprintf('pH %d', data.pH(j,1)),"0000e+00");
% % % %                     text(30, p101.YLim(2)*0.9, tempText);
% % % %                     % suptitle
% % % %                     if j == numpH
% % % %                         textHere = {'Aldolase: concentration progression curve. Each box contains a different pH value.';...
% % % %                             'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
% % % %                         suptitle(textHere)
% % % %                         set(101,'color','white'),
% % % %                     end
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102);
% % % %                     p102 = subplot(3,4,j);
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
% % % %                     end
% % % %                     % ylabel left
% % % %                     if((j == 1)||(j == 5))
% % % %                         ylabel('GPD reaction rate [mM s^{-1}]')
% % % %                     end
% % % %                     % xlabel
% % % %                     if((j == 5)||(j == 6)||(j == 7)||(j == 8))
% % % %                         xlabel('assay time [s]')
% % % %                     end
% % % %                     % textbox
% % % %                     tempText = erase(sprintf('pH %d', data.pH(j,1)),"0000e+00");
% % % %                     text(30, p102.YLim(2)*0.9, tempText)
% % % %                     % % % % Weird layout: stop here.
% % % %                     % suptitle
% % % %                     if j == numpH
% % % %                         textHere = {'Aldolase: reaction rate progression curve. Each box contains a different pH value.';...
% % % %                             'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
% % % %                         suptitle(textHere);
% % % %                         set(102,'color','white'),
% % % %                     end                    
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
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:3)';
                     
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
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:3)';
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
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:3)';
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
    case 'pyk'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = data.Vmax(j,:);
%             data.Vmaxs = data.Vmax(8,:);
            data.NADH = data.conc_mean(j,:);
            data.Vprofs = data.RRs(j,:);
            data.tempTime = data.time(j,:);
            % inputs to be selected
%             data.chosenKeq = setup.keq(j); 
            data.chosenKeq_PYK = setup.Keq_PYK(j);% = [1/(3.1E-6) 1/(3.5E-6) 1/(3.9E-6) 1/(4.5E-6) 1/(6.5E-6) 1/(9.7E-6) 1/(1.6E-5) 1/(2.5E-5) 1/(4.1E-5) 1/(5.9E-5) 1/(7.9E-5) 1/(9.6E-5)];  %dir-
            data.chosenKeq_LDH = setup.Keq_LDH(j);%
            data.i = j;
            % selecting the right parameters
            for temp11 = 1
                xassay = zeros(1,7);
            switch j
                case 1
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(7);
                case 2
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(8);
                case 3
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(9);
                case 4
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(10);
                case 5
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(11);
                case 6
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(12);
                case 7
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(13);
                case 8
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(14);
                case 9
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(15);
                case 10
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(16);
                case 11
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(17);
                case 12
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    xassay(4) = x_temp(4);
                    xassay(5) = x_temp(5);
                    xassay(6) = x_temp(6);
                    xassay(7) = x_temp(18);
                otherwise
                    disp('Something went wrong is selecting the pH value');
            end
            end
            
            % simulations
%             simRes = cell(numpH,DFstudy);
%             j, data.Vmaxs,
            for i = DFstudy
                % recall vmax for the specific value and simulate
                data.chosenDF = data.DF(j,i);
                data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
                data.chosenNADHini = data.NADH{i}(1);                
                % simulate metabolites
% % % %                 disp(j);
                [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                % calculation of reaction rates
                [vObs,~] = calcRates(xassay,simResult,data,setup);   
                % cost function (NADHexp + vGAPDHr)
                simTime = simResult.t;
                simMet = simResult.y(:,obsMet);
                simRate = vObs(:,2);

                simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                simLDH{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                expNADH{i,j} = data.NADH{i};
                expLDH{i,j} = data.Vprofs{i};
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
                    simulationVisualization;
% % % %                     
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 0.15])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('NADH concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simLDH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, -expLDH{i,j},'k.','MarkerSize',4)
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('LDH reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end
% % % %                     
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
                        plot(data.tempTime{i}, simLDH{i,j},'-')
%                         plot(simRes{j}.t, simRes{j}.v, '-')
                        hold on
                        plot(data.tempTime{i}, -expLDH{i,j},'k+')
                        hold on
                    end
                    title('v_{apparent.LDH}')
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
                        errorNADH1_2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_2 = simNADH{3,10} - expNADH{3,10};
                        errorNADH11_2 = simNADH{3,11} - expNADH{3,11};
                        errorNADH12_2 = simNADH{3,12} - expNADH{3,12};
                    errorNADH = [...
                        wDesp(1)*errorNADH1; wDesp(1)*errorNADH1_2; 
                        wDesp(2)*errorNADH2; wDesp(2)*errorNADH2_2;
                        wDesp(3)*errorNADH3; wDesp(3)*errorNADH3_2;
                        wDesp(4)*errorNADH4; wDesp(4)*errorNADH4_2;
                        wDesp(5)*errorNADH5; wDesp(5)*errorNADH5_2;
                        wDesp(6)*errorNADH6; wDesp(6)*errorNADH6_2;
                        wDesp(7)*errorNADH7; wDesp(7)*errorNADH7_2;
                        wDesp(8)*errorNADH8; wDesp(8)*errorNADH8_2;
                        wDesp(9)*errorNADH9; wDesp(9)*errorNADH9_2;
                        wDesp(10)*errorNADH10; wDesp(10)*errorNADH10_2;
                        wDesp(11)*errorNADH11; wDesp(11)*errorNADH11_2;
                        wDesp(12)*errorNADH12; wDesp(12)*errorNADH12_2];
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:6)';
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
                    errorNADH9_df1 = simNADH{4,9} - expNADH{4,9};
                    errorNADH10_df1 = simNADH{4,10} - expNADH{4,10};
                    errorNADH11_df1 = simNADH{4,11} - expNADH{4,11};
                    errorNADH12_df1 = simNADH{4,12} - expNADH{4,12};
                        % DF2
                        errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_df2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_df2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_df2 = simNADH{3,10} - expNADH{3,10};
                        errorNADH11_df2 = simNADH{3,11} - expNADH{3,11};
                        errorNADH12_df2 = simNADH{3,12} - expNADH{3,12};
                            % DF4
                            errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
                            errorNADH8_df4 = simNADH{2,8} - expNADH{2,8};
                            errorNADH9_df4 = simNADH{2,9} - expNADH{2,9};
                            errorNADH10_df4 = simNADH{2,10} - expNADH{2,10};
                            errorNADH11_df4 = simNADH{2,11} - expNADH{2,11};
                            errorNADH12_df4 = simNADH{2,12} - expNADH{2,12};
                                % DF8
                                errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
                                errorNADH8_df8 = simNADH{1,8} - expNADH{1,8};
                                errorNADH9_df8 = simNADH{1,9} - expNADH{1,9};
                                errorNADH10_df8 = simNADH{1,10} - expNADH{1,10};
                                errorNADH11_df8 = simNADH{1,11} - expNADH{1,11};
                                errorNADH12_df8 = simNADH{1,12} - expNADH{1,12};
                wDesp_t = wDesp';
                errorNADH = [...
                    wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                    wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                    wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                    wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                    wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                    wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                    wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1;
                    wDesp_t(1,8)*errorNADH8_df8; wDesp_t(2,8)*errorNADH8_df4; wDesp_t(3,8)*errorNADH8_df2; wDesp_t(4,8)*errorNADH8_df1;
                    wDesp_t(1,9)*errorNADH9_df8; wDesp_t(2,9)*errorNADH9_df4; wDesp_t(3,9)*errorNADH9_df2; wDesp_t(4,9)*errorNADH9_df1;
                    wDesp_t(1,10)*errorNADH10_df8; wDesp_t(2,10)*errorNADH10_df4; wDesp_t(3,10)*errorNADH10_df2; wDesp_t(4,10)*errorNADH10_df1;
                    wDesp_t(1,11)*errorNADH11_df8; wDesp_t(2,11)*errorNADH11_df4; wDesp_t(3,11)*errorNADH11_df2; wDesp_t(4,11)*errorNADH11_df1;
                    wDesp_t(1,12)*errorNADH12_df8; wDesp_t(2,12)*errorNADH12_df4; wDesp_t(3,12)*errorNADH12_df2; wDesp_t(4,12)*errorNADH12_df1];
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:6)';
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
    case 'pgi'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = data.Vmax(j,:);
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
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                    % % % % xassay(4) = x_temp(15);
                case 2
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(4);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(16);
                case 3
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(5);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(17);
                case 4
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(6);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(18);
                case 5
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(7);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(19);
                case 6
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(8);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(20);
                case 7
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(9);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(21);
                case 8
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(10);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(22);
                case 9
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(11);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(23);
                case 10
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(12);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(24);
                case 11
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(13);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(25);
                case 12
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(14);
% % % %                     % % % % xassay(4) = x_temp(15);
                    % % % % xassay(4) = x_temp(26);
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
                    simulationVisualization;
% % % %                     
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 0.15])
% % % %                         xlim([0 1000])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('NADH concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
% % % %                     end
% % % %                     xlim([0 1000])
% % % %                     if j == numpH
% % % %                         suptitle('GPD reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end                    
% % % % 
% % % %                     if j == 1
% % % %                         h103 = figure(103);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h103);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simG6P{i,j},'-','LineWidth',2)
% % % %                         ylim([4 5])
% % % %                         xlim([0 1000])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('G6P concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
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
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';
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
                    errorNADH9_df1 = simNADH{4,9} - expNADH{4,9};
                    errorNADH10_df1 = simNADH{4,10} - expNADH{4,10};
                    errorNADH11_df1 = simNADH{4,11} - expNADH{4,11};
                    errorNADH12_df1 = simNADH{4,12} - expNADH{4,12};
                        % DF2
                        errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_df2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_df2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_df2 = simNADH{3,10} - expNADH{3,10};
                        errorNADH11_df2 = simNADH{3,11} - expNADH{3,11};
                        errorNADH12_df2 = simNADH{3,12} - expNADH{3,12};
                            % DF4
                            errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
                            errorNADH8_df4 = simNADH{2,8} - expNADH{2,8};
                            errorNADH9_df4 = simNADH{2,9} - expNADH{2,9};
                            errorNADH10_df4 = simNADH{2,10} - expNADH{2,10};
                            errorNADH11_df4 = simNADH{2,11} - expNADH{2,11};
                            errorNADH12_df4 = simNADH{2,12} - expNADH{2,12};
                                % DF8
                                errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
                                errorNADH8_df8 = simNADH{1,8} - expNADH{1,8};
                                errorNADH9_df8 = simNADH{1,9} - expNADH{1,9};
                                errorNADH10_df8 = simNADH{1,10} - expNADH{1,10};
                                errorNADH11_df8 = simNADH{1,11} - expNADH{1,11};
                                errorNADH12_df8 = simNADH{1,12} - expNADH{1,12};
                wDesp_t = wDesp';
                errorNADH = [...
                    wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                    wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                    wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                    wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                    wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                    wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                    wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1;
                    wDesp_t(1,8)*errorNADH8_df8; wDesp_t(2,8)*errorNADH8_df4; wDesp_t(3,8)*errorNADH8_df2; wDesp_t(4,8)*errorNADH8_df1;
                    wDesp_t(1,9)*errorNADH9_df8; wDesp_t(2,9)*errorNADH9_df4; wDesp_t(3,9)*errorNADH9_df2; wDesp_t(4,9)*errorNADH9_df1;
                    wDesp_t(1,10)*errorNADH10_df8; wDesp_t(2,10)*errorNADH10_df4; wDesp_t(3,10)*errorNADH10_df2; wDesp_t(4,10)*errorNADH10_df1;
                    wDesp_t(1,11)*errorNADH11_df8; wDesp_t(2,11)*errorNADH11_df4; wDesp_t(3,11)*errorNADH11_df2; wDesp_t(4,11)*errorNADH11_df1;
                    wDesp_t(1,12)*errorNADH12_df8; wDesp_t(2,12)*errorNADH12_df4; wDesp_t(3,12)*errorNADH12_df2; wDesp_t(4,12)*errorNADH12_df1];
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';
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
    case 'pdc'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = data.Vmax(j,:);
            data.NADH = data.conc_mean(j,:);
            data.Vprofs = data.RRs(j,:);
            data.tempTime = data.time(j,:);
            % inputs to be selected
        	data.chosenKeq_ADH = setup.Keq_ADH(j);%
            data.i = j;
            % selecting the right parameters
            for temp11 = 1
                xassay = zeros(1,3);
            switch j
                case 1
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                case 2
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(4);
                case 3
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(5);
                case 4
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(6);
                case 5
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(7);
                case 6
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(8);
                case 7
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(9);
                case 8
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(10);
                case 9
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(11);
                case 10
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(12);
                case 11
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(13);
                case 12
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(14);
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
                data.chosenNADHini = data.NADH{i}(1);                
                % simulate metabolites
                [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                % calculation of reaction rates
                [vObs,~] = calcRates(xassay,simResult,data,setup);   
                % cost function (NADHexp + vGAPDHr)
                simTime = simResult.t;
                simMet = simResult.y(:,obsMet);
                simRate = vObs(:,2);
                % locate values in structure
                simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                simPDC{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                expNADH{i,j} = data.NADH{i};
                expPDC{i,j} = data.Vprofs{i};
            end
        end
        
        for j = 1:numpH
            data.tempTime = data.time(j,:);
            if plotEachSimCF == 1
                if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                    simulationVisualization;
% % % %                     
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 0.15])
% % % %                         xlim([0 300])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('NADH concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simPDC{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, -expPDC{i,j},'k.','MarkerSize',4)
% % % % %                         ylim([0 5E-4])
% % % %                     end
% % % %                     xlim([0 300])
% % % %                     if j == numpH
% % % %                         suptitle('PDC reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end
% % % %                     
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
                    xlim([0 300])
                    ylim([0 0.15])
                    subplot(1,2,2)
                    for i = DFstudy
                        plot(data.tempTime{i}, simPDC{i,j},'-')
                        hold on
                        plot(data.tempTime{i}, -expPDC{i,j},'k+')
                        hold on
                    end
                    title('v_{apparent.PDC}')
                    xlim([0 300])
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
                        errorNADH1_2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_2 = simNADH{3,10} - expNADH{3,10};
                        errorNADH11_2 = simNADH{3,11} - expNADH{3,11};
                        errorNADH12_2 = simNADH{3,12} - expNADH{3,12};
                    errorNADH = [...
                        wDesp(1)*errorNADH1; (1/2)*wDesp(1)*errorNADH1_2; 
                        wDesp(2)*errorNADH2; (1/2)*wDesp(2)*errorNADH2_2;
                        wDesp(3)*errorNADH3; (1/2)*wDesp(3)*errorNADH3_2;
                        wDesp(4)*errorNADH4; (1/2)*wDesp(4)*errorNADH4_2;
                        wDesp(5)*errorNADH5; (1/2)*wDesp(5)*errorNADH5_2;
                        wDesp(6)*errorNADH6; (1/2)*wDesp(6)*errorNADH6_2;
                        wDesp(7)*errorNADH7; (1/2)*wDesp(7)*errorNADH7_2;
                        wDesp(8)*errorNADH8; (1/2)*wDesp(8)*errorNADH8_2;
                        wDesp(9)*errorNADH9; (1/2)*wDesp(9)*errorNADH9_2;
                        wDesp(10)*errorNADH10; (1/2)*wDesp(10)*errorNADH10_2;
                        wDesp(11)*errorNADH11; (1/2)*wDesp(11)*errorNADH11_2;
                        wDesp(12)*errorNADH12; (1/2)*wDesp(12)*errorNADH12_2];
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';
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
                    errorNADH9_df1 = simNADH{4,9} - expNADH{4,9};
                    errorNADH10_df1 = simNADH{4,10} - expNADH{4,10};
                    errorNADH11_df1 = simNADH{4,11} - expNADH{4,11};
                    errorNADH12_df1 = simNADH{4,12} - expNADH{4,12};
                        % DF2
                        errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_df2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_df2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_df2 = simNADH{3,10} - expNADH{3,10};
                        errorNADH11_df2 = simNADH{3,11} - expNADH{3,11};
                        errorNADH12_df2 = simNADH{3,12} - expNADH{3,12};
                            % DF4
                            errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
                            errorNADH8_df4 = simNADH{2,8} - expNADH{2,8};
                            errorNADH9_df4 = simNADH{2,9} - expNADH{2,9};
                            errorNADH10_df4 = simNADH{2,10} - expNADH{2,10};
                            errorNADH11_df4 = simNADH{2,11} - expNADH{2,11};
                            errorNADH12_df4 = simNADH{2,12} - expNADH{2,12};
                                % DF8
                                errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
                                errorNADH8_df8 = simNADH{1,8} - expNADH{1,8};
                                errorNADH9_df8 = simNADH{1,9} - expNADH{1,9};
                                errorNADH10_df8 = simNADH{1,10} - expNADH{1,10};
                                errorNADH11_df8 = simNADH{1,11} - expNADH{1,11};
                                errorNADH12_df8 = simNADH{1,12} - expNADH{1,12};
                wDesp_t = wDesp';
                errorNADH = [...
                    wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                    wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                    wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                    wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                    wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                    wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                    wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1;
                    wDesp_t(1,8)*errorNADH8_df8; wDesp_t(2,8)*errorNADH8_df4; wDesp_t(3,8)*errorNADH8_df2; wDesp_t(4,8)*errorNADH8_df1;
                    wDesp_t(1,9)*errorNADH9_df8; wDesp_t(2,9)*errorNADH9_df4; wDesp_t(3,9)*errorNADH9_df2; wDesp_t(4,9)*errorNADH9_df1;
                    wDesp_t(1,10)*errorNADH10_df8; wDesp_t(2,10)*errorNADH10_df4; wDesp_t(3,10)*errorNADH10_df2; wDesp_t(4,10)*errorNADH10_df1;
                    wDesp_t(1,11)*errorNADH11_df8; wDesp_t(2,11)*errorNADH11_df4; wDesp_t(3,11)*errorNADH11_df2; wDesp_t(4,11)*errorNADH11_df1;
                    wDesp_t(1,12)*errorNADH12_df8; wDesp_t(2,12)*errorNADH12_df4; wDesp_t(3,12)*errorNADH12_df2; wDesp_t(4,12)*errorNADH12_df1];
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';
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
    case 'pgm'
        % simulations loop for each pH value
% % % %         for j = 1:numpH
        for j = 7
%         for j = 11:12
            % select required data
            data.Vmaxs = data.Vmax(j,:);
            data.NADH = data.conc_mean(j,:);
            data.Vprofs = data.RRs(j,:);
            data.tempTime = data.time(j,:);
            % inputs to be selected
            data.chosenKeq_PGM = setup.Keq_PGM(j);
            data.chosenKeq_ENO = setup.Keq_ENO(j);
            data.chosenKeq_PYK = setup.Keq_PYK(j);
            data.chosenKeq_LDH = setup.Keq_LDH(j);
            data.i = j;
            % selecting the right parameters
            for temp11 = 1
                xassay = zeros(1,3);
            switch j
                case 1
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                case 2
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(4);
                case 3
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(5);
                case 4
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(6);
                case 5
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(7);
                case 6
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(8);
                case 7
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(9);
                case 8
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(10);
                case 9
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(11);
                case 10
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(12);
                case 11
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(13);
                case 12
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(14);
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
                data.chosenNADHini = data.NADH{i}(1);                
                % simulate metabolites
                [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                % calculation of reaction rates
                [vObs,~] = calcRates(xassay,simResult,data,setup);   
                % cost function (NADHexp + vGAPDHr)
                simTime = simResult.t;
                simMet = simResult.y(:,obsMet);
                simRate = vObs(:,4);
                % locate values in structure
                simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                simPGM{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                expNADH{i,j} = data.NADH{i};
                expPGM{i,j} = data.Vprofs{i};
            end
        end
        
% % % %         for j = 1:numpH
        for j = 7
            data.tempTime = data.time(j,:);
            if plotEachSimCF == 1
                if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                    
% % % %                     simulationVisualization;
                    for ENOsim = 1
    %                     %% Simulation concentrations
                        % selection of the timeframe
                        data.tempTime = data.time(j,:);
                        % plottinh
                        if j == 7
% % % %                             h101 = figure(101);
                            h101 = figure(1000);
                            hold on
                        end
                        set(0,'CurrentFigure', h101);
% % % %                         p101 = subplot(3,4,j);
                        if setup.caseStudyENO == 1
                            p101 = subplot(3,3,2);
                            hold on
                        elseif setup.caseStudyGAPDHr == 1
                            p101 = subplot(3,3,5);
                            hold on
                        elseif setup.caseStudyPGM == 1
                            p101 = subplot(3,3,8);
                            hold on
                            if((isfield(setup, 'only2DFs'))&&(setup.only2DFs == 1))
                                DFstudy([2 3]) = [];
                            end
                        end
                        
% % % %                         if setup.literatureValues == 0
% % % %                             p101 = subplot(3,3,2);
% % % %                         elseif setup.literatureValues == 1
% % % %                             p101 = subplot(3,3,5);
% % % %                         end
                        
                        for i = DFstudy
                            if setup.caseStudyENO == 1
% % % %                                 plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
% % % %                                 ylim([0 1.2])
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.2])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.2])
                                    hold on
                                end
                                
                                
                            elseif setup.caseStudyHXK == 1
                                plot(data.tempTime{i}, simNADPH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADPH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                            elseif setup.caseStudyGAPDH == 1
                        %         plot(data.time{j,i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                        %         hold on
                        %         plot(data.time{j,i}, expNADH{i,j},'k.','MarkerSize',4)
                        %         ylim([0 0.15])
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                        %         disp(sprintf('j=%d, i=%d',j,i));
                            elseif setup.caseStudyGAPDHr == 1
                                
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                end
                                
                            elseif setup.caseStudyPGM == 1
                                
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}(linspace(1,61,16)), expNADH{i,j}(linspace(1,61,16)),'k.','MarkerSize',4)
%                                     ylim([0 0.15])
                                    ylim([0 0.10])
                                    xlim([0 300])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}(linspace(1,61,16)), expNADH{i,j}(linspace(1,61,16)),'k.','MarkerSize',4)
%                                     ylim([0 0.15])
                                    ylim([0 0.10])
                                    xlim([0 300])
                                    hold on
                                end
                                                                    
                            else
                                plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                                ylim([0 0.15])
                                if(setup.caseStudyPDC == 1)
                                    xlim([0 300])
                                elseif(setup.caseStudyPGM == 1)
                                    ylim([0 0.10])
                                    xlim([0 300])
                                elseif((setup.caseStudyPGI == 1)||(setup.caseStudyPFK == 1))
                                    ylim([0 0.15])
                                    xlim([0 1000])
                                elseif(setup.caseStudyTPI == 1)
                                    ylim([0 0.15])
                                    xlim([0 600])
                                end        
                            end
                        end
                        % ylabel left
                        if((j == 1)||(j == 5)||(j == 9))
                            if setup.caseStudyENO == 1
                                ylabel('PEP concentration [mM]')
                            elseif setup.caseStudyHXK == 1
                                ylabel('NAPDH concentration [mM]')
                            else 
                                ylabel('NADH concentration [mM]')        
                            end
                        end
                        % xlabel
                        if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
                            i_xlabel = 5:8;
                        else
                            i_xlabel = 9:12;
                        end
                        if((j == i_xlabel(1))||(j == i_xlabel(2))||(j == i_xlabel(3))||(j == i_xlabel(4)))
%                             xlabel('assay time [s]')
                        end
                        % pH text box
                        tempText = erase(sprintf('@pH %d', data.pH(j,1)),"0000e+00");
%                         text(30, p101.YLim(2)*0.9, tempText);
%                         text(300*0.95, p101.YLim(2)*0.9, tempText,...
%                             'HorizontalAlignment', 'right')
                        % suptitle
                        if j == numpH
% % % %                             textHere_pre = [setup.enzymeName,': concentration progression curve. Each box contains a different pH value.'];
% % % %                             textHere = {textHere_pre;...
% % % %                                 'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
% % % %                             suptitle(textHere)
% % % %                             set(101,'color','white'),
                        end
                        % title
% % % %                         title('train data: concentrations')
                        box on

    %                     %% Simulation reaction rates
                        if j == 7
% % % %                             h102 = figure(102);
                            h102 = figure(1000);
                            hold on
                        end
                        set(0,'CurrentFigure', h102);
% % % %                         p102 = subplot(3,4,j);
% % % %                         p102 = subplot(3,3,3);
% % % %                         hold on
                        
                        if setup.caseStudyENO == 1
                            p102 = subplot(3,3,3);
                            hold on
                        elseif setup.caseStudyGAPDHr == 1
                            p102 = subplot(3,3,6);
                            hold on
                        elseif setup.caseStudyPGM == 1
                            p102 = subplot(3,3,9);
                            hold on
                        end
                        
                        
% % % %                         if setup.literatureValues == 0
% % % %                             p102 = subplot(3,3,3);
% % % %                         elseif setup.literatureValues == 1
% % % %                             p102 = subplot(3,3,6);
% % % %                         end
                        
                        for i = DFstudy
                            if setup.caseStudyALD == 1
                                plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyENO == 1
% % % %                                 plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
% % % %                                 ylim([0 1.5E-3])
% % % %                                 hold on                               
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.5E-3])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color', setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
                                    ylim([0 1.5E-3])
                                    hold on
                                end
                                
                            elseif setup.caseStudyHXK == 1
                                plot(data.tempTime{i}, simG6PDH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, expG6PDH{i,j},'k.','MarkerSize',4)
                                ylim([0 5E-4])
                            elseif setup.caseStudyPDC == 1
                                plot(data.tempTime{i}, simPDC{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expPDC{i,j},'k.','MarkerSize',4)
                                xlim([0 300])
                            elseif setup.caseStudyPYK == 1  
                                plot(data.tempTime{i}, simLDH{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expLDH{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyPGM == 1  
% % % %                                 plot(data.tempTime{i}, simPGM{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.tempTime{i}, -expPGM{i,j},'k.','MarkerSize',4)
% % % %                                 xlim([0 300])
                                if setup.literatureValues == 0
                                    plot(data.tempTime{i}, simPGM{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.tempTime{i}, -expPGM{i,j},'k.','MarkerSize',4)
                                    xlim([0 300])
                                    ylim([0 0.0015])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.tempTime{i}, simPGM{i,j},'-','LineWidth',2,'color', setup.c_royalBlue)
                                    hold on
                                    plot(data.tempTime{i}, -expPGM{i,j},'k.','MarkerSize',4)
                                    xlim([0 300])
                                    ylim([0 0.0015])
                                    hold on
                                end
                                
                                
                            elseif setup.caseStudyPGI == 1
                                plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
                                xlim([0 1000])
                            elseif setup.caseStudyTPI == 1
                                plot(data.tempTime{i}, simTPI{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expTPI{i,j},'k.','MarkerSize',4)
                                xlim([0 600])
                            elseif setup.caseStudyPFK == 1
                                plot(data.tempTime{i}, simPFK{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.tempTime{i}, -expPFK{i,j},'k.','MarkerSize',4)
                                xlim([0 1000])
                            elseif setup.caseStudyGAPDH == 1
                                plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
                                hold on
                                plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                            elseif setup.caseStudyGAPDHr == 1
% % % %                                 plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
% % % %                                 hold on
% % % %                                 plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                
                                if setup.literatureValues == 0
                                    plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',setup.c_midnightblue)
                                    hold on
                                    plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                elseif setup.literatureValues == 1
                                    plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',setup.c_royalBlue)
                                    hold on
                                    plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                                    ylim([0 0.15])
                                    hold on
                                end
                                
                            end
                        end
                        % ylabel left
                        if((j == 1)||(j == 5)||(j == 9))
                            if((setup.caseStudyALD == 1)||(setup.caseStudyPGI == 1)||(setup.caseStudyPFK == 1)||(setup.caseStudyTPI == 1))
                                ylabel('GPD reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyENO == 1
                                ylabel('ENO reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyHXK == 1
                                ylabel('G6PDH reaction rate [mM s^{-1}]')
                            elseif((setup.caseStudyPYK == 1)||(setup.caseStudyPDC == 1)||(setup.caseStudyPGM == 1))
                                ylabel('LDH reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyGAPDH== 1
                                ylabel('GAPDH_{fwd} reaction rate [mM s^{-1}]')
                            elseif setup.caseStudyGAPDHr== 1
                                ylabel('GAPDH_{rev} reaction rate [mM s^{-1}]')
                            end
                        end
                        % xlabel
                        if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
                            i_xlabel = 5:8;
                        else
                            i_xlabel = 9:12;
                        end
                        if((j == i_xlabel(1))||(j == i_xlabel(2))||(j == i_xlabel(3))||(j == i_xlabel(4)))
%                             xlabel('assay time [s]')
                        end
                        % textbox
                        tempText = erase(sprintf('@pH %d', data.pH(j,1)),"0000e+00");
%                         text(30, p102.YLim(2)*0.9, tempText)
%                         text(300*0.95, p102.YLim(2)*0.9, tempText,...
%                             'HorizontalAlignment', 'right')
                        % suptitle
                        if j == numpH
% % % %                             textHere_pre = [setup.enzymeName,': reaction rate progression curve. Each box contains a different pH value.'];
% % % %                             textHere = {textHere_pre;...
% % % %                                 'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
% % % %                             suptitle(textHere);
% % % %                             set(102,'color','white'),
                        end
                        % title
% % % %                         title('test data: progression curve')
                        box on

                    end
                    
                    
                    
% % % %                     
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 0.10])
% % % %                         xlim([0 300])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('NADH concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simPGM{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, -expPGM{i,j},'k.','MarkerSize',4)
% % % %                     end
% % % %                     xlim([0 300])
% % % %                     if j == numpH
% % % %                         suptitle('PGM reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end                    
% % % % 
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
                    xlim([0 300])
                    ylim([0 0.10])
                    subplot(1,2,2)
                    for i = DFstudy
                        plot(data.tempTime{i}, simPGM{i,j},'-')
                        hold on
                        plot(data.tempTime{i}, -expPGM{i,j},'k+')
                        hold on
                    end
                    title('v_{apparent.PGM}')
                    xlim([0 300])
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
% % % %                     errorNADH8 = simNADH{4,8} - expNADH{4,8};
% % % %                     errorNADH9 = simNADH{4,9} - expNADH{4,9};
% % % %                     errorNADH10 = simNADH{4,10} - expNADH{4,10};
% % % %                     errorNADH11 = simNADH{4,11} - expNADH{4,11};
% % % %                     errorNADH12 = simNADH{4,12} - expNADH{4,12};
                        errorNADH1_2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_2 = simNADH{3,7} - expNADH{3,7};
% % % %                         errorNADH8_2 = simNADH{3,8} - expNADH{3,8};
% % % %                         errorNADH9_2 = simNADH{3,9} - expNADH{3,9};
% % % %                         errorNADH10_2 = simNADH{3,10} - expNADH{3,10};
% % % %                         errorNADH11_2 = simNADH{3,11} - expNADH{3,11};
% % % %                         errorNADH12_2 = simNADH{3,12} - expNADH{3,12};
                            errorNADH1_3 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_3 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_3 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_3 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_3 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_3 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_3 = simNADH{2,7} - expNADH{2,7};
% % % %                             errorNADH8_3 = simNADH{2,8} - expNADH{2,8};
% % % %                             errorNADH9_3 = simNADH{2,9} - expNADH{2,9};
% % % %                             errorNADH10_3 = simNADH{2,10} - expNADH{2,10};
% % % %                             errorNADH11_3 = simNADH{2,11} - expNADH{2,11};
% % % %                             errorNADH12_3 = simNADH{2,12} - expNADH{2,12};
                    errorNADH = [...
                        wDesp(1)*errorNADH1; (1/2)*wDesp(1)*errorNADH1_2; (1/4)*wDesp(1)*errorNADH1_3; % 
                        wDesp(2)*errorNADH2; (1/2)*wDesp(2)*errorNADH2_2; (1/4)*wDesp(2)*errorNADH2_3; %
                        wDesp(3)*errorNADH3; (1/2)*wDesp(3)*errorNADH3_2; (1/4)*wDesp(3)*errorNADH3_3; %
                        wDesp(4)*errorNADH4; (1/2)*wDesp(4)*errorNADH4_2; (1/4)*wDesp(4)*errorNADH4_3; %
                        wDesp(5)*errorNADH5; (1/2)*wDesp(5)*errorNADH5_2; (1/4)*wDesp(5)*errorNADH5_3; %
                        wDesp(6)*errorNADH6; (1/2)*wDesp(6)*errorNADH6_2; (1/4)*wDesp(6)*errorNADH6_3; %
                        wDesp(7)*errorNADH7; (1/2)*wDesp(7)*errorNADH7_2; (1/4)*wDesp(7)*errorNADH7_3]; %
% % % %                         wDesp(8)*errorNADH8; (1/2)*wDesp(8)*errorNADH8_2; (1/4)*wDesp(8)*errorNADH8_3; %
% % % %                         wDesp(9)*errorNADH9; (1/2)*wDesp(9)*errorNADH9_2; (1/4)*wDesp(9)*errorNADH9_3; %
% % % %                         wDesp(10)*errorNADH10; (1/2)*wDesp(10)*errorNADH10_2; (1/4)*wDesp(10)*errorNADH10_3;
% % % %                         wDesp(11)*errorNADH11; (1/2)*wDesp(11)*errorNADH11_2; (1/4)*wDesp(11)*errorNADH11_3;
% % % %                         wDesp(12)*errorNADH12; (1/2)*wDesp(12)*errorNADH12_2; (1/4)*wDesp(12)*errorNADH12_3];
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';
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
% % % %                     errorNADH8_df1 = simNADH{4,8} - expNADH{4,8};
% % % %                     errorNADH9_df1 = simNADH{4,9} - expNADH{4,9};
% % % %                     errorNADH10_df1 = simNADH{4,10} - expNADH{4,10};
% % % %                     errorNADH11_df1 = simNADH{4,11} - expNADH{4,11};
% % % %                     errorNADH12_df1 = simNADH{4,12} - expNADH{4,12};
                        % DF2
                        errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
% % % %                         errorNADH8_df2 = simNADH{3,8} - expNADH{3,8};
% % % %                         errorNADH9_df2 = simNADH{3,9} - expNADH{3,9};
% % % %                         errorNADH10_df2 = simNADH{3,10} - expNADH{3,10};
% % % %                         errorNADH11_df2 = simNADH{3,11} - expNADH{3,11};
% % % %                         errorNADH12_df2 = simNADH{3,12} - expNADH{3,12};
                            % DF4
                            errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
% % % %                             errorNADH8_df4 = simNADH{2,8} - expNADH{2,8};
% % % %                             errorNADH9_df4 = simNADH{2,9} - expNADH{2,9};
% % % %                             errorNADH10_df4 = simNADH{2,10} - expNADH{2,10};
% % % %                             errorNADH11_df4 = simNADH{2,11} - expNADH{2,11};
% % % %                             errorNADH12_df4 = simNADH{2,12} - expNADH{2,12};
                                % DF8
                                errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
% % % %                                 errorNADH8_df8 = simNADH{1,8} - expNADH{1,8};
% % % %                                 errorNADH9_df8 = simNADH{1,9} - expNADH{1,9};
% % % %                                 errorNADH10_df8 = simNADH{1,10} - expNADH{1,10};
% % % %                                 errorNADH11_df8 = simNADH{1,11} - expNADH{1,11};
% % % %                                 errorNADH12_df8 = simNADH{1,12} - expNADH{1,12};
                wDesp_t = wDesp';
                errorNADH = [...
                    wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                    wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                    wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                    wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                    wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                    wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                    wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1];
% % % %                     wDesp_t(1,8)*errorNADH8_df8; wDesp_t(2,8)*errorNADH8_df4; wDesp_t(3,8)*errorNADH8_df2; wDesp_t(4,8)*errorNADH8_df1;
% % % %                     wDesp_t(1,9)*errorNADH9_df8; wDesp_t(2,9)*errorNADH9_df4; wDesp_t(3,9)*errorNADH9_df2; wDesp_t(4,9)*errorNADH9_df1;
% % % %                     wDesp_t(1,10)*errorNADH10_df8; wDesp_t(2,10)*errorNADH10_df4; wDesp_t(3,10)*errorNADH10_df2; wDesp_t(4,10)*errorNADH10_df1;
% % % %                     wDesp_t(1,11)*errorNADH11_df8; wDesp_t(2,11)*errorNADH11_df4; wDesp_t(3,11)*errorNADH11_df2; wDesp_t(4,11)*errorNADH11_df1;
% % % %                     wDesp_t(1,12)*errorNADH12_df8; wDesp_t(2,12)*errorNADH12_df4; wDesp_t(3,12)*errorNADH12_df2; wDesp_t(4,12)*errorNADH12_df1];
% % % %                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';                
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
    case 'tpi'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = data.Vmax(j,:);
            data.NADH = data.conc_mean(j,:);
            data.Vprofs = data.RRs(j,:);
            data.tempTime = data.time(j,:);
            % inputs to be selected
            data.chosenKeq_TPI = setup.Keq_TPI(j);
            data.chosenKeq_GPD = setup.Keq_GPD(j);
            data.i = j;
            % selecting the right parameters
            for temp11 = 1
                xassay = zeros(1,3);
            switch j
                case 1
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(3);
                case 2
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(4);
                case 3
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(5);
                case 4
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(6);
                case 5
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(7);
                case 6
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(8);
                case 7
                    xassay(1) = x_temp(1);
                    xassay(2) = x_temp(2);
                    xassay(3) = x_temp(9);
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
                data.chosenNADHini = data.NADH{i}(1);                
                % simulate metabolites
                [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                % calculation of reaction rates
                [vObs,~] = calcRates(xassay,simResult,data,setup);   
                % cost function (NADHexp + vGAPDHr)
                simTime = simResult.t;
                simMet = simResult.y(:,obsMet);
                simRate = vObs(:,2);
                % locate values in structure
                simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                simTPI{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                expNADH{i,j} = data.NADH{i};
                expTPI{i,j} = data.Vprofs{i};
            end
        end
        
        for j = 1:numpH
            data.tempTime = data.time(j,:);
            if plotEachSimCF == 1
                if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                    simulationVisualization;
% % % %                     
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 0.15])
% % % %                         xlim([0 600])
% % % %                     end
% % % %                     title(erase(sprintf('pH = %d',data.pH(j,1)),"0000e+00"));
% % % %                     if j == numpH
% % % %                         suptitle('NADH concentration [mM] vs assay time [s]')
% % % %                     end
% % % % 
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simTPI{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, -expTPI{i,j},'k.','MarkerSize',4)
% % % %                     end
% % % %                     xlim([0 600])
% % % %                     title(erase(sprintf('pH = %d',data.pH(j,1)),"0000e+00"));
% % % %                     if j == numpH
% % % %                         suptitle('TPI reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end                    
% % % % 
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
                    xlim([0 600])
                    ylim([0 0.15])
                    subplot(1,2,2)
                    for i = DFstudy
                        plot(data.tempTime{i}, simTPI{i,j},'-')
                        hold on
                        plot(data.tempTime{i}, -expTPI{i,j},'k+')
                        hold on
                    end
                    title('v_{apparent.TPI}')
                    xlim([0 600])
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
                        errorNADH1_2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_2 = simNADH{3,7} - expNADH{3,7};
                            errorNADH1_3 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_3 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_3 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_3 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_3 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_3 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_3 = simNADH{2,7} - expNADH{2,7};
                                errorNADH1_4 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_4 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_4 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_4 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_4 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_4 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_4 = simNADH{1,7} - expNADH{1,7};
                    errorNADH = [...
                        wDesp(1)*errorNADH1; (1/2)*wDesp(1)*errorNADH1_2; (1/4)*wDesp(1)*errorNADH1_3; (1/8)*wDesp(1)*errorNADH1_4; 
                        wDesp(2)*errorNADH2; (1/2)*wDesp(2)*errorNADH2_2; (1/4)*wDesp(2)*errorNADH2_3; (1/8)*wDesp(2)*errorNADH2_4;
                        wDesp(3)*errorNADH3; (1/2)*wDesp(3)*errorNADH3_2; (1/4)*wDesp(3)*errorNADH3_3; (1/8)*wDesp(3)*errorNADH3_4;
                        wDesp(4)*errorNADH4; (1/2)*wDesp(4)*errorNADH4_2; (1/4)*wDesp(4)*errorNADH4_3; (1/8)*wDesp(4)*errorNADH4_4;
                        wDesp(5)*errorNADH5; (1/2)*wDesp(5)*errorNADH5_2; (1/4)*wDesp(5)*errorNADH5_3; (1/8)*wDesp(5)*errorNADH5_4;
                        wDesp(6)*errorNADH6; (1/2)*wDesp(6)*errorNADH6_2; (1/4)*wDesp(6)*errorNADH6_3; (1/8)*wDesp(6)*errorNADH6_4;
                        wDesp(7)*errorNADH7; (1/2)*wDesp(7)*errorNADH7_2; (1/4)*wDesp(7)*errorNADH7_3; (1/8)*wDesp(7)*errorNADH7_4];
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';
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
                        % DF2
                        errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
                            % DF4
                            errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
                                % DF8
                                errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
                wDesp_t = wDesp';
                errorNADH = [...
                    wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                    wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                    wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                    wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                    wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                    wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                    wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1];
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:2)';
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
    case 'pfk'
        % simulations loop for each pH value
        for j = 1:numpH
%         for j = 11:12
            % select required data
            data.Vmaxs = data.Vmax(j,:);
            data.NADH = data.conc_mean(j,:);
            data.Vprofs = data.RRs(j,:);
            data.tempTime = data.time(j,:);
            % inputs to be selected
            data.chosenKeq_FBA = setup.Keq_FBA(j);
            data.chosenKeq_GPD = setup.Keq_GPD(j);
            data.chosenKeq_TPI = setup.Keq_TPI(j);
            data.chosenKeq_PFK = setup.Keq_PFK(j);
            data.i = j;
            % selecting the right parameters
            problemStudy = setup.problemStudy;
            switch problemStudy
                
                case 'onlyVmax' % case 'onlyVmax'
                for temp11 = 1
    %                 xassay = zeros(1,14);
                    xassay = zeros(1,1);
                    switch j
                        case 1
                            xassay(1) = x_temp(1);
                        case 2
                            xassay(1) = x_temp(2);
                        case 3
                            xassay(1) = x_temp(3);
                        case 4
                            xassay(1) = x_temp(4);
                        case 5
                            xassay(1) = x_temp(5);
                        case 6
                            xassay(1) = x_temp(6);
                        case 7
                            xassay(1) = x_temp(7);
                        case 8
                            xassay(1) = x_temp(8);
                        case 9
                            xassay(1) = x_temp(9);
                        case 10
                            xassay(1) = x_temp(10);
                        case 11
                            xassay(1) = x_temp(11);
                        case 12
                            xassay(1) = x_temp(12);
                        otherwise
                            disp('Something went wrong is selecting the pH value');
                    end
                end

                case 'fullKinetics_paramsPartFixed'
                for temp11 = 2
                    xassay = zeros(1,14);
                    switch j
                        case 1
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(14);
                        case 2
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(15);
                        case 3
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(16);
                        case 4
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(17);
                        case 5
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(18);
                        case 6
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(19);
                        case 7
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(20);
                        case 8
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(21);
                        case 9
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(22);
                        case 10
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(23);
                        case 11
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(24);
                        case 12
                            xassay(1) = x_temp(1);
                            xassay(2) = x_temp(2);
                            xassay(3) = x_temp(3);
                            xassay(4) = x_temp(4);
                            xassay(5) = x_temp(5);
                            xassay(6) = x_temp(6);
                            xassay(7) = x_temp(7);
                            xassay(8) = x_temp(8);
                            xassay(9) = x_temp(9);
                            xassay(10) = x_temp(10);
                            xassay(11) = x_temp(11);
                            xassay(12) = x_temp(12);
                            xassay(13) = x_temp(13);
                            xassay(14) = x_temp(25);
                        otherwise
                            disp('Something went wrong is selecting the pH value');
                    end
                end

                case 'fullKinetics_paramsAllFlexible'
                for temp11 = 3
                    xassay = zeros(1,14);
                    % automatically locating idxs
                    idx1 = j;
                    idx2 = j + 12 * 1;
                    idx3 = j + 12 * 2;
                    idx4 = j + 12 * 3;
                    idx5 = j + 12 * 4;
                    idx6 = j + 12 * 5;
                    idx7 = j + 12 * 6;
                    idx8 = j + 12 * 7;
                    idx9 = j + 12 * 8;
                    idx10 = j + 12 * 9;
                    idx11 = j + 12 * 10;
                    idx12 = j + 12 * 11;
                    idx13 = j + 12 * 12;
                    idx14 = j + 12 * 13;
                    % selecting xassay based on the idxs
                    xassay(1) = x_temp(idx1);
                    xassay(2) = x_temp(idx2);
                    xassay(3) = x_temp(idx3);
                    xassay(4) = x_temp(idx4);
                    xassay(5) = x_temp(idx5);
                    xassay(6) = x_temp(idx6);
                    xassay(7) = x_temp(idx7);
                    xassay(8) = x_temp(idx8);
                    xassay(9) = x_temp(idx9);
                    xassay(10) = x_temp(idx10);
                    xassay(11) = x_temp(idx11);
                    xassay(12) = x_temp(idx12);
                    xassay(13) = x_temp(idx13);
                    xassay(14) = x_temp(idx14);
                end
                
                otherwise
                    disp('Warning: problem study (reaction kinetics) have not been selected');
            end
            
            % simulations
%             simRes = cell(numpH,DFstudy);
            for i = DFstudy
                % recall vmax for the specific value and simulate
                data.chosenDF = data.DF(j,i);
                data.chosenVmax = data.Vmaxs(1,4)/data.chosenDF; % vmax from the highest DF is taken and then divided
                data.chosenNADHini = data.NADH{i}(1);                
                % simulate metabolites
                [simResult] = simSys(xassay,data,setup); % simRes{i,j} = simResult;
                % calculation of reaction rates
                [vObs,~] = calcRates(xassay,simResult,data,setup);   
                % cost function (NADHexp + vGAPDHr)
                simTime = simResult.t;
                simMet = simResult.y(:,obsMet);
                simRate = vObs(:,3);
                % locate values in structure
                simNADH{i,j} = interp1(simTime,simMet,data.tempTime{i},'pchip');
                simPFK{i,j} = interp1(simTime,simRate,data.tempTime{i},'pchip');
                expNADH{i,j} = data.NADH{i};
                expPFK{i,j} = data.Vprofs{i};
            end
        end
        
        for j = 1:numpH
            data.tempTime = data.time(j,:);
            if plotEachSimCF == 1
                if((simAllProfiles == 1)&&(setup.plotEachSim == 0))
                    simulationVisualization;
% % % %                     
% % % %                     if j == 1
% % % %                         h101 = figure(101);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h101);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
% % % %                         ylim([0 0.15])
% % % %                         xlim([0 1000])
% % % %                     end
% % % %                     if j == numpH
% % % %                         suptitle('NADH concentration [mM] vs assay time [s]')
% % % %                     end
% % % %                     
% % % %                     if j == 1
% % % %                         h102 = figure(102);
% % % %                     end
% % % %                     set(0,'CurrentFigure', h102);
% % % %                     subplot(3,4,j)
% % % %                     for i = DFstudy
% % % %                         plot(data.tempTime{i}, simPFK{i,j},'-','LineWidth',2)
% % % %                         hold on
% % % %                         plot(data.tempTime{i}, -expPFK{i,j},'k.','MarkerSize',4)
% % % %                     end
% % % %                     xlim([0 1000])
% % % %                     if j == numpH
% % % %                         suptitle('PFK reaction rate [mM s^{-1}] vs assay time [s]')
% % % %                     end                    
% % % % 
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
                        plot(data.tempTime{i}, simPFK{i,j},'-')
                        hold on
                        plot(data.tempTime{i}, -expPFK{i,j},'k+')
                        hold on
                    end
                    title('v_{apparent.PFK}')
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
%                     errorNADH1 = simNADH{1,1} - expNADH{1,1};
%                     errorNADH2 = simNADH{1,2} - expNADH{1,2};
%                     errorNADH3 = simNADH{1,3} - expNADH{1,3};
%                     errorNADH4 = simNADH{1,4} - expNADH{1,4};
%                     errorNADH5 = simNADH{1,5} - expNADH{1,5};
%                     errorNADH6 = simNADH{1,6} - expNADH{1,6};
%                     errorNADH7 = simNADH{1,7} - expNADH{1,7};
%                     errorNADH8 = simNADH{1,8} - expNADH{1,8};
%                     errorNADH9 = simNADH{1,9} - expNADH{1,9};
%                     errorNADH10 = simNADH{1,10} - expNADH{1,10};
%                     errorNADH11 = simNADH{1,11} - expNADH{1,11};
%                     errorNADH12 = simNADH{1,12} - expNADH{1,12};
                        errorNADH1_2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_2 = simNADH{3,10} - expNADH{3,10};
                        errorNADH11_2 = simNADH{3,11} - expNADH{3,11};
                        errorNADH12_2 = simNADH{3,12} - expNADH{3,12};
                            errorNADH1_3 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_3 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_3 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_3 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_3 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_3 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_3 = simNADH{2,7} - expNADH{2,7};
                            errorNADH8_3 = simNADH{2,8} - expNADH{2,8};
                            errorNADH9_3 = simNADH{2,9} - expNADH{2,9};
                            errorNADH10_3 = simNADH{2,10} - expNADH{2,10};
                            errorNADH11_3 = simNADH{2,11} - expNADH{2,11};
                            errorNADH12_3 = simNADH{2,12} - expNADH{2,12};
                    errorNADH = [...
                        wDesp(1)*errorNADH1; %(1/2)*wDesp(1)*errorNADH1_2; (1/4)*wDesp(1)*errorNADH1_3; 
                        wDesp(2)*errorNADH2; %(1/2)*wDesp(2)*errorNADH2_2; (1/4)*wDesp(2)*errorNADH2_3;
                        wDesp(3)*errorNADH3; %(1/2)*wDesp(3)*errorNADH3_2; (1/4)*wDesp(3)*errorNADH3_3;
                        wDesp(4)*errorNADH4; %(1/2)*wDesp(4)*errorNADH4_2; (1/4)*wDesp(4)*errorNADH4_3;
                        wDesp(5)*errorNADH5; %(1/2)*wDesp(5)*errorNADH5_2; (1/4)*wDesp(5)*errorNADH5_3;
                        wDesp(6)*errorNADH6; %(1/2)*wDesp(6)*errorNADH6_2; (1/4)*wDesp(6)*errorNADH6_3;
                        wDesp(7)*errorNADH7; %(1/2)*wDesp(7)*errorNADH7_2; (1/4)*wDesp(7)*errorNADH7_3;
                        wDesp(8)*errorNADH8; %(1/2)*wDesp(8)*errorNADH8_2; (1/4)*wDesp(8)*errorNADH8_3;
                        wDesp(9)*errorNADH9; %(1/2)*wDesp(9)*errorNADH9_2; (1/4)*wDesp(9)*errorNADH9_3;
                        wDesp(10)*errorNADH10; %(1/2)*wDesp(10)*errorNADH10_2; (1/4)*wDesp(10)*errorNADH10_3;
                        wDesp(11)*errorNADH11; %(1/2)*wDesp(11)*errorNADH11_2; (1/4)*wDesp(11)*errorNADH11_3;
                        wDesp(12)*errorNADH12]; %(1/2)*wDesp(12)*errorNADH12_2; (1/4)*wDesp(12)*errorNADH12_3];
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:13)';
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
                    errorNADH9_df1 = simNADH{4,9} - expNADH{4,9};
                    errorNADH10_df1 = simNADH{4,10} - expNADH{4,10};
                    errorNADH11_df1 = simNADH{4,11} - expNADH{4,11};
                    errorNADH12_df1 = simNADH{4,12} - expNADH{4,12};
                        % DF2
                        errorNADH1_df2 = simNADH{3,1} - expNADH{3,1};
                        errorNADH2_df2 = simNADH{3,2} - expNADH{3,2};
                        errorNADH3_df2 = simNADH{3,3} - expNADH{3,3};
                        errorNADH4_df2 = simNADH{3,4} - expNADH{3,4};
                        errorNADH5_df2 = simNADH{3,5} - expNADH{3,5};
                        errorNADH6_df2 = simNADH{3,6} - expNADH{3,6};
                        errorNADH7_df2 = simNADH{3,7} - expNADH{3,7};
                        errorNADH8_df2 = simNADH{3,8} - expNADH{3,8};
                        errorNADH9_df2 = simNADH{3,9} - expNADH{3,9};
                        errorNADH10_df2 = simNADH{3,10} - expNADH{3,10};
                        errorNADH11_df2 = simNADH{3,11} - expNADH{3,11};
                        errorNADH12_df2 = simNADH{3,12} - expNADH{3,12};
                            % DF4
                            errorNADH1_df4 = simNADH{2,1} - expNADH{2,1};
                            errorNADH2_df4 = simNADH{2,2} - expNADH{2,2};
                            errorNADH3_df4 = simNADH{2,3} - expNADH{2,3};
                            errorNADH4_df4 = simNADH{2,4} - expNADH{2,4};
                            errorNADH5_df4 = simNADH{2,5} - expNADH{2,5};
                            errorNADH6_df4 = simNADH{2,6} - expNADH{2,6};
                            errorNADH7_df4 = simNADH{2,7} - expNADH{2,7};
                            errorNADH8_df4 = simNADH{2,8} - expNADH{2,8};
                            errorNADH9_df4 = simNADH{2,9} - expNADH{2,9};
                            errorNADH10_df4 = simNADH{2,10} - expNADH{2,10};
                            errorNADH11_df4 = simNADH{2,11} - expNADH{2,11};
                            errorNADH12_df4 = simNADH{2,12} - expNADH{2,12};
                                % DF8
                                errorNADH1_df8 = simNADH{1,1} - expNADH{1,1};
                                errorNADH2_df8 = simNADH{1,2} - expNADH{1,2};
                                errorNADH3_df8 = simNADH{1,3} - expNADH{1,3};
                                errorNADH4_df8 = simNADH{1,4} - expNADH{1,4};
                                errorNADH5_df8 = simNADH{1,5} - expNADH{1,5};
                                errorNADH6_df8 = simNADH{1,6} - expNADH{1,6};
                                errorNADH7_df8 = simNADH{1,7} - expNADH{1,7};
                                errorNADH8_df8 = simNADH{1,8} - expNADH{1,8};
                                errorNADH9_df8 = simNADH{1,9} - expNADH{1,9};
                                errorNADH10_df8 = simNADH{1,10} - expNADH{1,10};
                                errorNADH11_df8 = simNADH{1,11} - expNADH{1,11};
                                errorNADH12_df8 = simNADH{1,12} - expNADH{1,12};
                wDesp_t = wDesp';
                errorNADH = [...
                    wDesp_t(1,1)*errorNADH1_df8; wDesp_t(2,1)*errorNADH1_df4; wDesp_t(3,1)*errorNADH1_df2; wDesp_t(4,1)*errorNADH1_df1; 
                    wDesp_t(1,2)*errorNADH2_df8; wDesp_t(2,2)*errorNADH2_df4; wDesp_t(3,2)*errorNADH2_df2; wDesp_t(4,2)*errorNADH2_df1;
                    wDesp_t(1,3)*errorNADH3_df8; wDesp_t(2,3)*errorNADH3_df4; wDesp_t(3,3)*errorNADH3_df2; wDesp_t(4,3)*errorNADH3_df1;
                    wDesp_t(1,4)*errorNADH4_df8; wDesp_t(2,4)*errorNADH4_df4; wDesp_t(3,4)*errorNADH4_df2; wDesp_t(4,4)*errorNADH4_df1;
                    wDesp_t(1,5)*errorNADH5_df8; wDesp_t(2,5)*errorNADH5_df4; wDesp_t(3,5)*errorNADH5_df2; wDesp_t(4,5)*errorNADH5_df1;
                    wDesp_t(1,6)*errorNADH6_df8; wDesp_t(2,6)*errorNADH6_df4; wDesp_t(3,6)*errorNADH6_df2; wDesp_t(4,6)*errorNADH6_df1;
                    wDesp_t(1,7)*errorNADH7_df8; wDesp_t(2,7)*errorNADH7_df4; wDesp_t(3,7)*errorNADH7_df2; wDesp_t(4,7)*errorNADH7_df1;
                    wDesp_t(1,8)*errorNADH8_df8; wDesp_t(2,8)*errorNADH8_df4; wDesp_t(3,8)*errorNADH8_df2; wDesp_t(4,8)*errorNADH8_df1;
                    wDesp_t(1,9)*errorNADH9_df8; wDesp_t(2,9)*errorNADH9_df4; wDesp_t(3,9)*errorNADH9_df2; wDesp_t(4,9)*errorNADH9_df1;
                    wDesp_t(1,10)*errorNADH10_df8; wDesp_t(2,10)*errorNADH10_df4; wDesp_t(3,10)*errorNADH10_df2; wDesp_t(4,10)*errorNADH10_df1;
                    wDesp_t(1,11)*errorNADH11_df8; wDesp_t(2,11)*errorNADH11_df4; wDesp_t(3,11)*errorNADH11_df2; wDesp_t(4,11)*errorNADH11_df1;
                    wDesp_t(1,12)*errorNADH12_df8; wDesp_t(2,12)*errorNADH12_df4; wDesp_t(3,12)*errorNADH12_df2; wDesp_t(4,12)*errorNADH12_df1];
%                     errorReg = lambda * x_temp';
                    errorReg = lambda * x_temp(1:13)';
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

