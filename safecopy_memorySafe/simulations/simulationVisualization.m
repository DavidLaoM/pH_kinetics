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
for i = DFstudy
    if setup.caseStudyENO == 1
        plot(data.tempTime{i}, simPEP{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
        hold on
        plot(data.tempTime{i}, expPEP{i,j},'k.','MarkerSize',4)
        ylim([0 1.2])
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
    xlabel('assay time [s]')
end
% pH text box
tempText = erase(sprintf('pH %d', data.pH(j,1)),"0000e+00");
text(30, p101.YLim(2)*0.9, tempText);
% suptitle
if j == numpH
    textHere_pre = [setup.enzymeName,': concentration progression curve. Each box contains a different pH value.'];
    textHere = {textHere_pre;...
        'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
    suptitle(textHere)
    set(101,'color','white'),
end


%% Simulation reaction rates
if j == 1
    h102 = figure(102);
end
set(0,'CurrentFigure', h102);
p102 = subplot(3,4,j);
for i = DFstudy
    if setup.caseStudyALD == 1
        plot(data.tempTime{i}, simGPD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
        hold on
        plot(data.tempTime{i}, -expGPD{i,j},'k.','MarkerSize',4)
    elseif setup.caseStudyENO == 1
        plot(data.tempTime{i}, simENO{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
        hold on
        plot(data.tempTime{i}, expENO{i,j},'k.','MarkerSize',4)
        ylim([0 1.5E-3])
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
        plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
        hold on
        plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
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

