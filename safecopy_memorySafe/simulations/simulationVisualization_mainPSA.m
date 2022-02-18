% % SimulationVisualization


%% Simulation concentrations
% selection of the timeframe
data.tempTime = data.time(j,:);
% plottinh
if j == 1
    h101 = figure(101);
end
set(0,'CurrentFigure', h101);
p101 = subplot(3,4,j);
% % % % % Added for the Keq
% colArray = jet(length(simResults));
colArray = cool(length(simResults));
% colArray = parula(length(simResults));
% colArray = winter(length(simResults));
% colArray = copper(length(simResults));
% colArray = summer(length(simResults));
for k = 1:length(simResults)
    % load back simResult
    if setup.caseStudyENO == 1
        simPEP = simResults{k}.simPEP; %: {4×12 cell}
%         simENO = simResults{k}.simENO; %: {4×12 cell}
        expPEP = simResults{k}.expPEP; %: {4×12 cell}
%         expENO = simResults{k}.expENO; %: {4×12 cell}
    else
        simNADH = simResults{k}.simNADH;
        expNADH = simResults{k}.expNADH;
    end
    DFstudy = simResult.DFstudy;
    data.tempTime = data.time(j,:);
    
    for i = DFstudy
        if setup.caseStudyENO == 1
%             plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',colArray(k,:))
            timeInterp = [data.tempTime{i}; setup.addedTime];
            plot(timeInterp, simPEP{i,j},'-','LineWidth',2,'color',colArray(k,:))
            hold on
            if k == length(simResults)
                % % addition for DoE
                % DLM 2020-11-19
                if isfield(setup,'PSAstudy') % PSA?
                    % which PSA?
                    if isfield(setup,'PSAstudy_ENO_kp2g')
                        if setup.PSAstudy_ENO_k2pg == 1
                            tempVal = expPEP{i,j} - expPEP{i,j}(1);
                            plot(data.tempTime{i}, tempVal,'k.','MarkerSize',4)
                        end
                    end
                    if isfield(setup,'PSAstudy_ENO_kpep')
                        if setup.PSAstudy_ENO_kpep == 1
                            tempVal = expPEP{i,j} - expPEP{i,j}(1);
                            plot(data.tempTime{i}, tempVal,'k.','MarkerSize',4)
                        end
                    end
                else
                    plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
                end
                % % 
% % % % % % % % % %                 plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
                if j == numpHtested
                    legend off
% % % %                     set(gcf, 'Colormap', colArray);
% % % %                     cb = colorbar;
% % % %                     cb.TickLabels = tempTickLabels;
                end
            end
            
            % % addition for DoE
            % DLM 2020-11-19
            if isfield(setup,'PSAstudy') % PSA?
                % which PSA?
                if isfield(setup,'PSAstudy_ENO_k2pg')
                    if setup.PSAstudy_ENO_k2pg == 1
                        ylim([0 0.5])
                    end
                end
                if isfield(setup,'PSAstudy_ENO_pep')
                    if setup.PSAstudy_ENO_kpep == 1
                        ylim([0 0.5])
                    end
                end
            else
                ylim([0 1.2])
            end
            % % 
% % % % % % % % % %         elseif setup.caseStudyHXK == 1
% % % % % % % % % %             plot(data.tempTime{i}, simNADPH{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             if k == length(simResults)
% % % % % % % % % %                 plot(data.tempTime{i}, expNADPH{i,j},'k.','MarkerSize',4)
% % % % % % % % % %             end
% % % % % % % % % %             ylim([0 0.15])
        elseif setup.caseStudyGAPDH == 1
%             plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',colArray(k,:))
            timeInterp = [data.tempTime{i}; setup.addedTime];
            plot(timeInterp, simNADH{i,j},'-','LineWidth',2,'color',colArray(k,:))
            hold on
% %             
            if k == length(simResults)
                % % addition for DoE
                % DLM 2020-11-19
                if isfield(setup,'PSAstudy') % PSA?
                    % which PSA?
                    if isfield(setup,'PSAstudy_GAPDH_kgap')
                        if setup.PSAstudy_GAPDH_kgap == 1
                            tempVal = expNADH{i,j} - 0;%expPEP{i,j}(1);
                            plot(data.tempTime{i}, tempVal,'k.','MarkerSize',4)
                        end
                    end
%                     if isfield(setup,'PSAstudy_ENO_kpep')
%                         if setup.PSAstudy_ENO_kpep == 1
%                             tempVal = expPEP{i,j} - expPEP{i,j}(1);
%                             plot(data.tempTime{i}, tempVal,'k.','MarkerSize',4)
%                         end
%                     end
                else
                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                end
                if j == numpHtested
                    legend off
                end
            end
            ylim([0 0.15])  
        elseif setup.caseStudyPFK == 1
%             plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',colArray(k,:))
            timeInterp = [data.tempTime{i}; setup.addedTime];
            plot(timeInterp, simNADH{i,j},'-','LineWidth',2,'color',colArray(k,:))
            hold on
            if k == length(simResults)
                % % addition for DoE
                % DLM 2020-11-19
                if isfield(setup,'PSAstudy') % PSA?
                    % which PSA?
                    if isfield(setup,'PSAstudy_PFK_kf6p')
                        if setup.PSAstudy_PFK_kf6p == 1
% % % %                             if j <= 6
                                tempVal = expNADH{i,j} - 0;%expPEP{i,j}(1);
                                tempTime = data.tempTime{i};
% % % %                             else
% % % %                                 tempVal = expNADH{i,j} - 0;%expPEP{i,j}(1);
% % % %                                 tempTime = data.tempTime{i};
% % % %                             end
                            plot(tempTime, tempVal,'k.','MarkerSize',4)
                        end
                    end
                else
                    plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
                end
                if j == numpHtested
                    legend off
                end
            end
            ylim([0 0.15])
            xlim([0 1000])
% % % % % % % % % %         else
% % % % % % % % % %             plot(data.tempTime{i}, simNADH{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             plot(data.tempTime{i}, expNADH{i,j},'k.','MarkerSize',4)
% % % % % % % % % %             ylim([0 0.15])
% % % % % % % % % %             if(setup.caseStudyPDC == 1)
% % % % % % % % % %                 xlim([0 300])
% % % % % % % % % %             elseif(setup.caseStudyPGM == 1)
% % % % % % % % % %                 ylim([0 0.10])
% % % % % % % % % %                 xlim([0 300])
% % % % % % % % % %             elseif((setup.caseStudyPGI == 1)||(setup.caseStudyPFK == 1))
% % % % % % % % % %                 ylim([0 0.15])
% % % % % % % % % %                 xlim([0 1000])
% % % % % % % % % %             elseif(setup.caseStudyTPI == 1)
% % % % % % % % % %                 ylim([0 0.15])
% % % % % % % % % %                 xlim([0 600])
% % % % % % % % % %             end        
        end
    end
    % ylabel left
    if((j == 1)||(j == 5)||(j == 9))
        if setup.caseStudyENO == 1
            ylabel('PEP concentration [mM]')
% % % % % % % % % %         elseif setup.caseStudyHXK == 1
% % % % % % % % % %             ylabel('NAPDH concentration [mM]')
% % % % % % % % % %         else 
% % % % % % % % % %             ylabel('NADH concentration [mM]')        
        end
    end
    % xlabel
    if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
        i_xlabel = 5:8;
    else
        i_xlabel = 9:12;
    end
    if((j == i_xlabel(1))||(j == i_xlabel(2))||(j == i_xlabel(3))||(j == i_xlabel(4)))
        xlabel('assay time [s]')
    end
    % pH text box
    tempText = erase(sprintf('pH %d', data.pH(j,1)),"0000e+00");
    text(30, p101.YLim(2)*0.9, tempText);
    % suptitle
    if j == numpH
% % % % % % % % % %         textHere_pre = [setup.enzymeName,': concentration progression curve. Each box contains a different pH value.'];
% % % % % % % % % %         textHere = {textHere_pre;...
% % % % % % % % % %             'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
%         get(101);
%         gcf(101);
        figure(101);
%         fig = get(101,'CurrentFigure');
        textHere = setup.supTitleText;
        hST = suptitle(textHere);
        hST.Position = [0.1500 -0.0500 0];
        set(101,'color','white'),
    end
    %
    %
    if i == 4
        % setup.plotBackColor = [.5 .5 1]; % blue
        % setup.plotBackColor = [.5 1 .5]; % green
        % setup.plotBackColor = [1 .5 .5]; % red
        % setup.plotBackColor = [.94 .94 .94]; % standard grey
        val = 0.4;
        if((j == 1)||(j == 6))
            set(gca,'Color',[.5 1 .5] + [val 0 val]);
        elseif((j == 7)||(j == 8)||(j == 12))
            set(gca,'Color',[.5 .5 1] + [val val 0]);
        end
    end
    tempNum = numpH - 1;
    if((j == tempNum)&&(k == length(simResults)))
        figure(101);
%         legNames = {'uno','dos','tres'};
%         hL = legend(legNames);
        hL = legend(setup.legNames);
        set(hL,'orientation','horizontal');
        set(hL,'position',setup.legendLocation);
    end
end

%% Simulation reaction rates
if j == 1
    h102 = figure(102);
end
set(0,'CurrentFigure', h102);
p102 = subplot(3,4,j);
for k = 1:length(simResults)
    % load back simResult
    if setup.caseStudyENO == 1
%         simPEP = simResults{k}.simPEP; %: {4×12 cell}
        simENO = simResults{k}.simENO; %: {4×12 cell}
%         expPEP = simResults{k}.expPEP; %: {4×12 cell}
        expENO = simResults{k}.expENO; %: {4×12 cell}
    elseif setup.caseStudyPFK == 1
        simPFK = simResults{k}.simPFK; %: {4×12 cell}
        expPFK = simResults{k}.expPFK; %: {4×12 cell}
    else
        simGAPDHr = simResults{k}.simGAPDHr;
        expGAPDHr = simResults{k}.expGAPDHr;
    end
%     DFstudy = simResult.DFstudy;
%     data.tempTime = data.time(j,:);
    
    for i = DFstudy
        if setup.caseStudyALD == 1
            plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',colArray(k,:))
            hold on
            plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
        elseif setup.caseStudyENO == 1
%             plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color',colArray(k,:))
            timeInterp = [data.tempTime{i}; setup.addedTime];
            plot(timeInterp, simENO{i,j},'-','LineWidth',2,'color',colArray(k,:))
            hold on
            if k == length(simResults)
                plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
                if j == numpHtested
                    legend off
% % % %                     set(gcf, 'Colormap', colArray);
% % % %                     cb = colorbar;
                end
            end
            ylim([0 1.5E-3])
        elseif setup.caseStudyGAPDH == 1
            timeInterp = [data.tempTime{i}; setup.addedTime];
            plot(timeInterp, simGAPDHr{i,j},'-','LineWidth',2,'color',colArray(k,:))
            hold on
            if k == length(simResults)
                plot(data.tempTime{i}, expGAPDHr{i,j},'k.','MarkerSize',4)
                if j == numpHtested
                    legend off
% % % %                     set(gcf, 'Colormap', colArray);
% % % %                     cb = colorbar;
                end
            end
            ylim([0 1.5E-3])
        elseif setup.caseStudyPFK == 1
            timeInterp = [data.tempTime{i}; setup.addedTime];
            plot(timeInterp, simPFK{i,j},'-','LineWidth',2,'color',colArray(k,:))
            hold on
            if k == length(simResults)
                plot(data.tempTime{i}, -expPFK{i,j},'k.','MarkerSize',4)
                if j == numpHtested
                    legend off
                end
            end
            xlim([0 1000])
            ylim([0 1.5E-4])
% % % % % % % % % %         elseif setup.caseStudyPYK == 1  
% % % % % % % % % %             plot(data.tempTime{i}, simLDH{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             plot(data.tempTime{i}, -expLDH{i,j},'k.','MarkerSize',4)
% % % % % % % % % %         elseif setup.caseStudyHXK == 1
% % % % % % % % % %             plot(data.tempTime{i}, simG6PDH{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             plot(data.tempTime{i}, expG6PDH{i,j},'k.','MarkerSize',4)
% % % % % % % % % %             ylim([0 5E-4])
% % % % % % % % % %         elseif setup.caseStudyPDC == 1
% % % % % % % % % %             plot(data.tempTime{i}, simPDC{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             plot(data.tempTime{i}, -expPDC{i,j},'k.','MarkerSize',4)
% % % % % % % % % %             xlim([0 300])
% % % % % % % % % %         elseif setup.caseStudyPGM == 1  
% % % % % % % % % %             plot(data.tempTime{i}, simPGM{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             plot(data.tempTime{i}, -expPGM{i,j},'k.','MarkerSize',4)
% % % % % % % % % %             xlim([0 300])
% % % % % % % % % %         elseif setup.caseStudyPGI == 1
% % % % % % % % % %             plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
% % % % % % % % % %             xlim([0 1000])
% % % % % % % % % %         elseif setup.caseStudyTPI == 1
% % % % % % % % % %             plot(data.tempTime{i}, simTPI{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             plot(data.tempTime{i}, -expTPI{i,j},'k.','MarkerSize',4)
% % % % % % % % % %             xlim([0 600])
% % % % % % % % % %         elseif setup.caseStudyGAPDHr == 1
% % % % % % % % % %             plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',colArray(k,:))
% % % % % % % % % %             hold on
% % % % % % % % % %             plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
        end
    end
    % ylabel left
    if((j == 1)||(j == 5)||(j == 9))
        if((setup.caseStudyALD == 1)||(setup.caseStudyPGI == 1)||(setup.caseStudyPFK == 1)||(setup.caseStudyTPI == 1))
            ylabel('GPD reaction rate [mM s^{-1}]')
        elseif setup.caseStudyENO == 1
            ylabel('ENO reaction rate [mM s^{-1}]')
% % % % % % % % % %         elseif setup.caseStudyHXK == 1
% % % % % % % % % %             ylabel('G6PDH reaction rate [mM s^{-1}]')
% % % % % % % % % %         elseif((setup.caseStudyPYK == 1)||(setup.caseStudyPDC == 1)||(setup.caseStudyPGM == 1))
% % % % % % % % % %             ylabel('LDH reaction rate [mM s^{-1}]')
% % % % % % % % % %         elseif setup.caseStudyGAPDH== 1
% % % % % % % % % %             ylabel('GAPDH_{fwd} reaction rate [mM s^{-1}]')
% % % % % % % % % %         elseif setup.caseStudyGAPDHr== 1
% % % % % % % % % %             ylabel('GAPDH_{rev} reaction rate [mM s^{-1}]')
        end
    end
    % xlabel
    if((setup.caseStudyALD == 1)||(setup.caseStudyTPI == 1))
        i_xlabel = 5:8;
    else
        i_xlabel = 9:12;
    end
    if((j == i_xlabel(1))||(j == i_xlabel(2))||(j == i_xlabel(3))||(j == i_xlabel(4)))
        xlabel('assay time [s]')
    end
    % textbox
    tempText = erase(sprintf('pH %d', data.pH(j,1)),"0000e+00");
    text(30, p102.YLim(2)*0.9, tempText)
    % suptitle
    if j == numpH
        textHere_pre = [setup.enzymeName,': reaction rate progression curve. Each box contains a different pH value.'];
        textHere = {textHere_pre;...
            'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
        suptitle(textHere);
        set(102,'color','white'),
    end
    % 
    % 
    if i == 4
        % setup.plotBackColor = [.5 .5 1]; % blue
        % setup.plotBackColor = [.5 1 .5]; % green
        % setup.plotBackColor = [1 .5 .5]; % red
        % setup.plotBackColor = [.94 .94 .94]; % standard grey
        val = 0.4;
        if((j == 1)||(j == 6))
            set(gca,'Color',[.5 1 .5] + [val 0 val]);
        elseif((j == 7)||(j == 8)||(j == 12))
            set(gca,'Color',[.5 .5 1] + [val val 0]);
        end
    end
    
end

