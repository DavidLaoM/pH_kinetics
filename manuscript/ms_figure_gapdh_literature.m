% % % % function x=dev_pEst_gapdh_rev
% % PEST_ALD.m
% Parameter estimation for the data in the Enolase assay.
% Vm are estimated changing with pH
% Kms are assumed constant
% Keq from the eQuilibrator


%% (0) Setup and data load
% clear, close all
% set_paths_pHstudy;
% dbstop if error
for step0 = 1
    % select specific case and recall data
    setup.caseStudyALD = 0;
    setup.caseStudyENO = 0;
    setup.caseStudyGAPDH = 0;
    setup.caseStudyGAPDHr = 1;
    setup.caseStudyHXK = 0;
    setup.caseStudyPDC = 0;
    setup.caseStudyPFK = 0;
    setup.caseStudyPGI = 0;
    setup.caseStudyPGM = 0;
    setup.caseStudyPYK = 0;
    setup.caseStudyTPI = 0;
    setup.caseStudyENO = 0;
    selectSetup_pH;
    % added
    setup.saveOutput = 0;
    
    load('expData.mat','expData');
    import_gapdhr = expData.gapdhr;
    
    DFs = setup.DFactorsTotal;
    pHtested = setup.pHtested;
    numpHtested = nnz(pHtested);
    pHs = numpHtested;
    blank = zeros(pHs,DFs);
    blankCell = cell(pHs,DFs);
    
    % data reorganization
    pHTemp = blank';
    DFTemp = blank';
    abs_meanTemp = blankCell';
    abs_stdTemp = blankCell';
    conc_meanTemp = blankCell';
    conc_stdTemp = blankCell';
    timeTemp = blankCell';
    RRsTemp = blankCell';
    
    pHarray = unique(import_gapdhr.treatedData.pH_corrected);
    for i = 1:numpHtested
        pHval = pHarray(i);
        tempID = find(import_gapdhr.treatedData.pH_corrected==pHval);
        pHTemp(:,i) = import_gapdhr.treatedData.pH_corrected(tempID);
        DFTemp(:,i) = import_gapdhr.treatedData.dilution_corrected(tempID);
        for j = 1:4
            abs_meanTemp{j,i} = import_gapdhr.treatedData.absorbance_mean{tempID(j)};
            abs_stdTemp{j,i} = import_gapdhr.treatedData.absorbance_std{tempID(j)};
            conc_meanTemp{j,i} = import_gapdhr.treatedData.concentration_mean{tempID(j)};
            conc_stdTemp{j,i} = import_gapdhr.treatedData.concentration_std{tempID(j)};
            timeTemp{j,i} = import_gapdhr.treatedData.time{tempID(j)};
            RRsTemp{j,i} = import_gapdhr.treatedData.reaction_rate{tempID(j)};
        end
    end

    pH = pHTemp';
    DF = DFTemp';
    abs_mean = abs_meanTemp';
    abs_std = abs_stdTemp';

    conc_mean = conc_meanTemp';
    conc_std = conc_stdTemp';
    time = timeTemp';
    RRs = RRsTemp';
    Vmax = blank';

    NADH = blankCell;
    Vmax = blank;
    for i = 1:(DFs*numpHtested)
        tempDiff = conc_mean{i} - conc_mean{i}(1); % all stoichiometries are 1-to-1.
        NADH{i} = conc_mean{i};
        %     RRs2 = RRs';
        % Option 1. Vmax from the values obtained
        Vmax(i) = max(abs(RRs{i}));
%         % Option 2. Vmax naive approach. First datapoints
%         Vmax(i) = (conc_mean{i}(end) - conc_mean{i}(1)) ./ (time{i}(end) - time{i}(1)); 
    end

    pH = pHTemp';
    DF = DFTemp';
    abs_mean = abs_meanTemp';
    abs_std = abs_stdTemp';
    conc_mean = conc_meanTemp';
    conc_std = conc_stdTemp';
    time = timeTemp';
    RRs = RRsTemp';
    clear pHTemp DFTemp abs_meanTemp abs_stdTemp conc_meanTemp conc_stdTemp timeTemp RRsTemp

    % save in data
    data.pH = pH;
    data.DF = DF;
    data.abs_mean = abs_mean;
    data.abs_std = abs_std;
    data.conc_mean = conc_mean;
    data.conc_std = conc_std;
    data.time = time;
    data.RRs = RRs;
    data.Vmax = Vmax;
    data.chosenVmax = max(max(Vmax));
    data.chosenNADHini = 0.15;
    temp1 = import_gapdhr.rawData.absorbance_corrected{4,4};
    temp2 = import_gapdhr.rawData.absorbance_corrected{5,4};
    temp3 = import_gapdhr.rawData.absorbance_corrected{6,4};
    data.raw.conc = [temp1, temp2, temp3]*setup.extinction_coefficient;
    data.raw.time = import_gapdhr.rawData.time{1};
    
        % (1) Correct for minimum value
        % (2) Bring the minimum to zero (apply to all)
        % (3) In principle, use the 3 first dilution rates
        % (4) Watch out with the dilution factors (first 2 cases are
        % reversed)
        
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Directly changing the concentration here, sicne the extinction
    % coefficient did not change.
    dps = length(NADH{1,1});
    for i = 1:DFs
        for j = 1:numpHtested
            for k = 1:dps
                switch j
                    case {1,2,3,4,5,6,7}
                        NADH{j,i}(k) = NADH{j,i}(k) - NADH{j,DFs}(dps);
                    case {8,9,10,11,12}
                        NADH{j,i}(k) = NADH{j,i}(k) - 0.0453;
    %                     NADH{j,i}(k) = NADH{j,i}(k) - NADH{6,DFs}(dps);
                end
            end
        end
    end
    data.conc_mean = NADH;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    pHvals = unique(import_gapdhr.treatedData.pH_corrected);
%     % visualize: check calculations made
%     figure('units','normalized','outerposition',[0 0 1 1])
%     for i = 1:numpHtested
%         subplot(3,4,i)
%         for j = 1:DFs
%             plot(time{i,j},NADH{i,j},'.-')
%             hold on
%         end
%         title(erase(sprintf('pH = %d', pHvals(i)),"0000e+00"))
%         if i == numpHtested
%             if setup.caseStudyGAPDHr == 1
%                 legend('DF 8','DF 4','DF 2','DF 1')
%             end
%         end
%         if setup.caseStudyENO == 1
%             ylim([0 1.5])
%         end
%         if setup.caseStudyHXK == 1
%             ylim([0 0.15])
%         end
%         if setup.caseStudyALD == 1
%             ylim([0 0.15])
%         end
%         if setup.caseStudyGAPDHr == 1
%             ylim([0 0.15])
%         end
%     end
%     suptitleName = ['Enzyme ', setup.enzymeName, ': NADH concentration profile'];
%     suptitle(suptitleName);
    
%     figure
%     plot(pHvals, Vmax(:,4),'.-')
%     title('Starting estimate: Vmax [mM s-1] vs pH')    
end


%% (0.1) Experimental Vmax determination
setup.plotOutput = 0;
% %% (0.1) Calculation of rates: moving window
    % intial things that could be in the setup
    minwindow = 6; % minimum size of the window
    limRates = [0 2E-3]; %Ylims plot vmaxs
    limR2 = [0 1]; %Ylims plot R2
    limcConc = [0 0.15];  %Ylims plot conc
% select start point (this needs manual selection deppending on previous plots)
dp_start = ones(size(data.conc_mean));
% blank cell total length
total_len = zeros(size(dp_start));
% DFs considered
DFarray = [1/8 1/4 1/2 1/1];
% idxs2consider
idxs2consider = ones(size(DF));

% Experimental rates determination and plotting
expRatesDetermination;

% % %% (0.3) saving
% if setup.plotOutput == 1
%     save(['results/',setup.enzymeName,'/',setup.enzymeName, '_initial_variables.mat'],'Vmax_mw_opt_corr','idxs2consider','DF','pH');
%     set(1,'color','white'), savefig(1,['results/',setup.enzymeName,'/',setup.enzymeName, '_concentrations_basezero.fig']);
%     set(2,'color','white'), savefig(2,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_df.fig']);
%     set(3,'color','white'), savefig(3,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_movingWindow.fig']);
%     set(4,'color','white'), savefig(4,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_R2_vs_movingWindow.fig']);
%     set(5,'color','white'), savefig(5,['results/',setup.enzymeName,'/',setup.enzymeName, '_mw_iniPoint_vs_movingWindow.fig']);
%     set(6,'color','white'), savefig(6,['results/',setup.enzymeName,'/',setup.enzymeName, '_experimental_vmax_vs_pH.fig']);
% end


%% (1.1) Simple parameter fit. Parameter estimation
% % % % setup.ode = 'vanHeerden2014';
setup.ode = 'gapdh_rev_simplified';
setup.sourceVm = 'experimentalSlopes';
setup.ode_pH = 'on';
setup.caseKm = 'pH_independent';
% setup.caseKm = 'pH_dependent';

setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

setup.numpHtested = numpHtested;
setup.DFstudy = 1:4;
setup.costfun = 3;

setup.weightData = 1;
% setup.weightDataEsp = ones(1,numpHtested);
setup.weightDataEsp = idxs2consider;
setup.weightHaldane = 0; % off for this case (Keq fixed)
setup.selectedLambda = 0; % by now just testing

% Km fixed
% caseKm = setup.caseKm;
% switch caseKm
%     case 'pH_independent'
        % plength = 10; % Kms (4*0) + Vms (1) * numpH (10) % overly simplified 2020-09-30
plength = 14; % Kms (4*1) + Vms (1) * numpH (10) % still keeping kms
%     case 'pH_dependent'
%         plength = 50; % ( Kms (4) + Vms (1) ) * numpH (10) % still keeping kms
%     otherwise
%         disp('Warning: no specification has been made on Km being pH dependent or independent');
% end
optfun = @costfun_Kmfixed;
x_temp = zeros(1,plength);
ub = 3*ones(1,plength);
lb = -3*ones(1,plength);
options = optimset('Display','iter');

%% test simulation
% 
setup.plotResults = 1;
setup.plotEachSimCF = 1;
setup.simAllProfiles = 1;
setup.plotEachSim = 0;
% 
setup.selectedLambda = 0;
[raw_error] = optfun(x_temp,data,setup);
% 
setup.plotResults = 0;
setup.plotEachSimCF = 0;
setup.simAllProfiles = 0;
setup.plotEachSim = 1;

%% adjusting looks
close(102)
c_classicBlue = [0 0.4470 0.7410];

for i = 7
% for i = 2 %[1 2] %[3 4] 
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % load the figure + rescale (train/concentrations)
    set(101, 'units','normalized','outerposition',[0 0 0.5 1]);
    fh1 = figure(101);
    fh1_children = get(fh1,'Children');
    for j = 1:length(fh1_children)
        
        % delete suptitle and some reordering
        if fh1_children(j).Tag == "suptitle"
            delete(fh1_children(j));
        else
            fh1_children(j).Children(1).Position(2) = fh1_children(j).YLim(2)*0.85;
            fh1_children(j).Children(1).Position(1) = fh1_children(j).XLim(2)*0.25;
        end
        
        % new changes in each subplot
        if IsAxes(fh1_children(j))
            % adding grid
            fh1_children(j).XGrid = 'on';
            fh1_children(j).YGrid = 'on';
%             % change the color
            fh1_children(j).Children(3).Color = c_classicBlue;
            fh1_children(j).Children(5).Color = c_classicBlue;
            fh1_children(j).Children(7).Color = c_classicBlue;
            fh1_children(j).Children(9).Color = c_classicBlue;
            % assay time -> time
            if contains(fh1_children(j).XLabel.String,'assay time [s]')
                fh1_children(j).XLabel.String = 'Time (s)';
            end
            % [mM] -> (mM) in NADPH
            if contains(fh1_children(j).YLabel.String,'NAPDH concentration [mM]')
                fh1_children(j).YLabel.String = 'NADPH concentration (mM)';
            end
            % [mM] -> (mM) in NADH
            if contains(fh1_children(j).YLabel.String,'NADH concentration [mM]')
                fh1_children(j).YLabel.String = 'NADH concentration (mM)';
            end
            % enzyme assay specific
            % axes
            if i == 1 % HXK
                fh1_children(j).XLim = [0 300];
                fh1_children(j).YLim = [0 0.20];
                fh1_children(j).XTick = [0 150 300];
            elseif i == 2 % PGI
                fh1_children(j).XLim = [0 500];
                fh1_children(j).YLim = [0 0.15];
                fh1_children(j).XTick = [0 250 500];
                fh1_children(j).YTickLabelRotation = 45;
            elseif i == 3 % PFK
                fh1_children(j).XLim = [0 1000];
                fh1_children(j).YLim = [0 0.15];
                fh1_children(j).XTick = [0 500 1000];
                fh1_children(j).YTickLabelRotation = 45;
            elseif i == 4 % ALD
                fh1_children(j).XLim = [0 400];
                fh1_children(j).YLim = [0 0.15];
                fh1_children(j).XTick = [0 200 400];
                fh1_children(j).YTickLabelRotation = 45;
            elseif i == 5 % TPI
                fh1_children(j).XLim = [0 400];
                fh1_children(j).YLim = [0 0.10];
                fh1_children(j).XTick = [0 200 400];
            elseif i == 6 % GAPDH
                fh1_children(j).XLim = [0 300];
                fh1_children(j).YLim = [0 0.15];
                fh1_children(j).XTick = [0 150 300];
                fh1_children(j).YTickLabelRotation = 45;
            elseif i == 7 % GAPDHR
                fh1_children(j).XLim = [0 300];
                fh1_children(j).YLim = [0 0.12];
                fh1_children(j).XTick = [0 150 300];
            elseif i == 8 % PGM
                fh1_children(j).XLim = [0 300];
                fh1_children(j).YLim = [0 0.1];
                fh1_children(j).XTick = [0 150 300];
            elseif i == 9 % ENO
                fh1_children(j).XLim = [0 600];
                fh1_children(j).YLim = [0 1.5];
                fh1_children(j).XTick = [0 300 600];
            elseif i == 10 % PYK
                fh1_children(j).XLim = [0 600];
                fh1_children(j).YLim = [0 0.12];
                fh1_children(j).XTick = [0 300 600];
            elseif i == 11 % PDC
                fh1_children(j).XLim = [0 300];
                fh1_children(j).YLim = [0 0.12];
                fh1_children(j).XTick = [0 150 300];
            end
            % common edits
            fh1_children(j).YTick = [fh1_children(j).YLim(1) 1/4*fh1_children(j).YLim(2) 2/4*fh1_children(j).YLim(2) 3/4*fh1_children(j).YLim(2) fh1_children(j).YLim(2)];
            fh1_children(j).XTickLabel = cellstr(num2str(fh1_children(j).XTick'));
%             if i == 6 % specific for GAPDH
%                 fh1_children(j).YTickLabel = cellstr(num2str(fh1_children(j).YTick','%7f'));
%             else
                fh1_children(j).YTickLabel = cellstr(num2str(fh1_children(j).YTick'));
%             end
            % text label
%             fh1_children(j).Children(1).Position = [0.75*fh1_children(j).XLim(2) 0.90*fh1_children(j).YLim(2) 0];
            fh1_children(j).Children(1).Position = [0.75*fh1_children(j).XLim(2) 1.10*fh1_children(j).YLim(2) 0];
            fh1_children(j).Children(1).HorizontalAlignment = 'center';
            
        end
        
    end
%     fh1.Position = [fh1.Position(1) fh1.Position(2) 0.9*fh1.Position(3) 0.9*fh1.Position(4) ];
    fh1.Position = [0.025    0.075    0.4425    0.8225];
    
    
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%     % Saving: save them back to .png in location
%     saveLoc = 'D:\OneDrive - TU Eindhoven\Documents\ch3_pHkinetics\results\manuscriptFigures\';
%     saveName_train_pdf = [saveLoc,'1appendix_trainData_gapdhr_LiteratureVals_fit_metabolites_reg.pdf'];
%     saveName_train_png = [saveLoc,'1appendix_trainData_gapdhr_LiteratureVals_fit_metabolites_reg.png'];
%     saveas(101,saveName_train_pdf);
%     saveas(101,saveName_train_png);

    % 
    if setup.saveOutput == 1
        % 
        savefig(101,'1appendix_trainData_gapdhr_LiteratureVals_fit_metabolites_reg')
        figure(101)
        % 
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        print -dpdf -painters 1appendix_trainData_gapdhr_LiteratureVals_fit_metabolites_reg
        print -dpng -painters 1appendix_trainData_gapdhr_LiteratureVals_fit_metabolites_reg
    end


end

