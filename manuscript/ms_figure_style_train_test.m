clear
% % recalling figures
% c_royalBlue = [65	105	225]/255; % royalblue
% c_midnightblue = [25	25	112]/255; % midnightblue
c_classicBlue = [0 0.4470 0.7410];

enzymeName = {'hxk';... %1 hxk
    'pgi';... %2 pgi
    'pfk';... %3 pfk
    'ald';... %4 ald
    'tpi';... %5 tpi
    'gapdh';... %6 gapdh_fwd
    'gapdhr';... %7 gapdh_rev
    'pgm';... %8 pgm
    'eno';... %9 eno_kmfixed
    'pyk';... %10 pyk
    'pdc'}; %11 pdc

%%
for i = 1:length(enzymeName)
% for i = 2 %[1 2] %[3 4] 
    
    % close to call by the same numbers
    close all
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % load the figure + rescale (train/concentrations)
    loadName_train = [enzymeName{i}, '_trainData_fit_metabolites_reg.fig'];
    fh1 = openfig(loadName_train);
    set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);
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
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % load the figure + rescale (test/reaction rate)
    loadName_test = [enzymeName{i}, '_testData_fit_fluxes_reg.fig'];
    fh2 = openfig(loadName_test);
    set(gcf, 'units','normalized','outerposition',[0.5 0 0.5 1]);
    fh2_children = get(fh2,'Children');
    for j = 1:length(fh2_children)
        
        % delete suptitle and some reordering
        if fh2_children(j).Tag == "suptitle"
            delete(fh2_children(j));
        else
            fh2_children(j).Children(1).Position(2) = fh2_children(j).YLim(2)*0.85;
            fh2_children(j).Children(1).Position(1) = fh2_children(j).XLim(2)*0.25;
        end
        
        % new changes in each subplot
        if IsAxes(fh2_children(j))
            % adding grid
            fh2_children(j).XGrid = 'on';
            fh2_children(j).YGrid = 'on';
%             % change the color
            fh2_children(j).Children(3).Color = c_classicBlue;
            fh2_children(j).Children(5).Color = c_classicBlue;
            fh2_children(j).Children(7).Color = c_classicBlue;
            fh2_children(j).Children(9).Color = c_classicBlue;
            % assay time -> time
            if contains(fh2_children(j).XLabel.String,'assay time [s]')
                fh2_children(j).XLabel.String = 'Time (s)';
            end
            % For reaction rate: name + [mM] -> (mM) 
            if contains(fh2_children(j).YLabel.String,'G6PDH reaction rate [mM s^{-1}]')
                fh2_children(j).YLabel.String = 'G6PDH reaction rate (mM s^{-1})';
            end
            % 
            if contains(fh2_children(j).YLabel.String,'GPD reaction rate [mM s^{-1}]')
                fh2_children(j).YLabel.String = 'GPD reaction rate (mM s^{-1})';
            end
            % 
            if contains(fh2_children(j).YLabel.String,'GAPDH_{fwd} reaction rate [mM s^{-1}]')
                fh2_children(j).YLabel.String = 'GAPDH reaction rate (mM s^{-1})';
            end
            % 
            if contains(fh2_children(j).YLabel.String,'GAPDH_{rev} reaction rate [mM s^{-1}]')
                fh2_children(j).YLabel.String = 'GAPDHR reaction rate (mM s^{-1})';
            end
            % 
            if contains(fh2_children(j).YLabel.String,'LDH reaction rate [mM s^{-1}]')
                fh2_children(j).YLabel.String = 'LDH reaction rate (mM s^{-1})';
            end
            % 
            if contains(fh2_children(j).YLabel.String,'ENO reaction rate [mM s^{-1}]')
                fh2_children(j).YLabel.String = 'ENO reaction rate (mM s^{-1})';
            end
            % enzyme assay specific
            % axes
            if i == 1 % HXK
                fh2_children(j).XLim = [0 300];
                fh2_children(j).YLim = [0 0.0005];
                fh2_children(j).XTick = [0 150 300];
            elseif i == 2 % PGI
                fh2_children(j).XLim = [0 500];
                fh2_children(j).YLim = [0 0.0015];
                fh2_children(j).XTick = [0 250 500];
            elseif i == 3 % PFK
                fh2_children(j).XLim = [0 1000];
                fh2_children(j).YLim = [0 0.00015];
                fh2_children(j).XTick = [0 500 1000];
            elseif i == 4 % PFK
                fh2_children(j).XLim = [0 400];
                fh2_children(j).YLim = [0 0.0020];
                fh2_children(j).XTick = [0 200 400];
            elseif i == 5 % TPI
                fh2_children(j).XLim = [0 400];
                fh2_children(j).YLim = [0 0.004];
                fh2_children(j).XTick = [0 200 400];
            elseif i == 6 % GAPDH
                fh2_children(j).XLim = [0 300];
                fh2_children(j).YLim = [0 0.0005];
                fh2_children(j).XTick = [0 150 300];
            elseif i == 7 % GAPDHR
                fh2_children(j).XLim = [0 300];
                fh2_children(j).YLim = [0 0.0015];
                fh2_children(j).XTick = [0 150 300];
            elseif i == 8 % PGM
                fh2_children(j).XLim = [0 300];
                fh2_children(j).YLim = [0 0.002];
                fh2_children(j).XTick = [0 150 300];
            elseif i == 9 % ENO
                fh2_children(j).XLim = [0 600];
                fh2_children(j).YLim = [0 0.0015];
                fh2_children(j).XTick = [0 300 600];
            elseif i == 10 % PYK
                fh2_children(j).XLim = [0 600];
                fh2_children(j).YLim = [0 0.002];
                fh2_children(j).XTick = [0 300 600];
            elseif i == 11 % PDC
                fh2_children(j).XLim = [0 300];
                fh2_children(j).YLim = [0 0.0015];
                fh2_children(j).XTick = [0 150 300];
            end
            % common edits
            fh2_children(j).XTickLabel = cellstr(num2str(fh2_children(j).XTick'));
            fh2_children(j).YTick = [fh2_children(j).YLim(1) 1/4*fh2_children(j).YLim(2) 2/4*fh2_children(j).YLim(2) 3/4*fh2_children(j).YLim(2) fh2_children(j).YLim(2)];
%             fh2_children(j).YTickLabel = cellstr(num2str(fh2_children(j).YTick'));
%             fh2_children(j).YTickLabel = cellstr(num2str(fh2_children(j).YTick','%.E'));
%             for k = 1:length(fh2_children(j).YTickLabel)
%                 fh2_children(j).YTickLabel{k} = fh2_children(j).YTickLabel{k}(1);
%             end
%             text(0.5,0.5,'text')
            
            % text label
%             fh2_children(j).Children(1).Position = [0.75*fh2_children(j).XLim(2) 0.90*fh2_children(j).YLim(2) 0];
            fh2_children(j).Children(1).Position = [0.75*fh2_children(j).XLim(2) 1.10*fh2_children(j).YLim(2) 0];
            fh2_children(j).Children(1).HorizontalAlignment = 'center';
            % correcting positions
            fh2_children(j).Position = fh1_children(j).Position;
        end
        
    end
    %linkprop(fh2_children, {'GridColor','GridLineStyle','GridAlpha'});
%     fh2.Position = [fh2.Position(1) fh2.Position(2) 0.9*fh2.Position(3) 0.9*fh2.Position(4) ];
    fh2.Position = [0.500    0.075    0.4425    0.8225];

    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % Saving: save them back to .png in location
    if setup.saveOutput == 1
    %     saveLoc = 'D:\OneDrive - TU Eindhoven\Documents\ch3_pHkinetics\results\manuscriptFigures\';
        saveName_train_png = ['1appendix_trainData_', enzymeName{i}, '_fit_metabolites_reg.png'];
        saveName_train_pdf = ['1appendix_trainData_', enzymeName{i}, '_fit_metabolites_reg.pdf'];
    %     saveName_train_pdf = [saveLoc,'1appendix_trainData_', enzymeName{i}, '_fit_metabolites_reg.pdf'];
    %     saveName_train_png = [saveLoc,'1appendix_trainData_', enzymeName{i}, '_fit_metabolites_reg.png'];
        saveas(1,saveName_train_pdf);
        saveas(1,saveName_train_png);
        % 
        saveName_test_png = ['1appendix_testData_', enzymeName{i}, '_fit_fluxes_reg.png'];
        saveName_test_pdf = ['1appendix_testData_', enzymeName{i}, '_fit_fluxes_reg.pdf'];
    %     saveName_test_pdf = [saveLoc,'1appendix_testData_', enzymeName{i}, '_fit_fluxes_reg.pdf'];
    %     saveName_test_png = [saveLoc,'1appendix_testData_', enzymeName{i}, '_fit_fluxes_reg.png'];
        saveas(2,saveName_test_pdf);
        saveas(2,saveName_test_png);
    end    
end




% %%
% % set(gcf, 'Position',  [0 0 1 1])
% % set(gcf, 'units','normalized','outerposition',[0 0 1 1])
% % figure('units','normalized','outerposition',[0 0 1 1])
% set(gcf, 'units','normalized','outerposition',[0.5 0 0.5 1])