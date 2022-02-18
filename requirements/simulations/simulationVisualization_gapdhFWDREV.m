% % SimulationVisualization


%% Simulation concentrations
% selection of the timeframe
data.FWD.tempTime = data.FWD.time(j,:);
% data.REV.tempTime = data.REV.time(j,:);
% plottinh
if j == 1
    h101 = figure(101);
end

% fowward reaction
set(0,'CurrentFigure', h101);
idx = j*2-1;
p101 = subplot(3,8,idx);
for i = DFstudy
        plot(data.FWD.tempTime{i}, simNADH_FWD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
        hold on
        plot(data.FWD.tempTime{i}, expNADH_FWD{i,j},'k.','MarkerSize',4)
        ylim([0 0.15])
        xlim([0 300])
        if i == DFstudy(end)
            % pH text box
            tempText = erase(sprintf('pH %d', data.FWD.pH(j,1)),"0000e+00");
            text(30, p101.YLim(2)*0.9, tempText);
        end
end


% reverse reaction
tempIdx = find(idxsTranspose == j);
if tempIdx >= 1
    k = j;%idxsTranspose(j);
    k2 = find(idxsTranspose == k);
    data.REV.tempTime = data.REV.time(k2,:);
    idx2 = k*2;
    p102 = subplot(3,8,idx2);
    for i = DFstudy
        plot(data.REV.tempTime{i}, simNADH_REV{i,k2},'-','LineWidth',2,'color',[0.5 0.5 1])
        hold on
        plot(data.REV.tempTime{i}, expNADH_REV{i,k2},'k.','MarkerSize',4)
        ylim([0 0.15])
        xlim([0 300])
    end
end

% suptitle
if j == numpH
%     textHere_pre = [setup.enzymeName,': concentration progression curve. Each box contains a different pH value.'];
%     textHere = {textHere_pre;...
%         'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
%     suptitle(textHere)
    set(101,'color','white'),
end







%% Simulation reaction rates
if j == 1
    h102 = figure(102);
end
set(0,'CurrentFigure', h102);
idx = j*2-1;
p102= subplot(3,8,idx);
for i = DFstudy
        plot(data.FWD.time{j,i}, simGAPDHr_FWD{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
        hold on
        plot(data.FWD.time{j,i}, expGAPDHr_FWD{i,j},'k.','MarkerSize',4)
        xlim([0 300])
%     elseif setup.caseStudyGAPDHr == 1
%         plot(data.time{j,i}, simGAPDHr{i,j},'-','LineWidth',2,'color',[0.5 0.5 1])
%         hold on
%         plot(data.time{j,i}, expGAPDHr{i,j},'k.','MarkerSize',4)
end
% textbox
tempText = erase(sprintf('pH %d', data.FWD.pH(j,1)),"0000e+00");
text(30, p102.YLim(2)*0.9, tempText)
% suptitle
if j == numpH
%     textHere_pre = [setup.enzymeName,': reaction rate progression curve. Each box contains a different pH value.'];
%     textHere = {textHere_pre;...
%         'Experimental data is shown in black dots, while simulations in blue continuous lines.'};
%     suptitle(textHere);
    set(102,'color','white'),
end

