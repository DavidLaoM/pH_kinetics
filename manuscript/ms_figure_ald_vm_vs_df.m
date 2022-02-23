% 
clear
close all
% manually made, the 2-y axis thing was taking too long

%%
% for i = 1:length(enzymeName)
for i = 1 %[1 2] %[3 4] 
    
    % close to call by the same numbers
    close all
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % load the figure + rescale
    fh1 = openfig('ald_mw_vmax_vs_df.fig');
    set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);
    fh1_children = get(fh1,'Children');
    for j = 1:length(fh1_children)
        
        % shift position from (2,4,#) to (4,2,#)... noe trying (3,3,#) rather.
        if j == 1, fh1_children(j).Position = [0.4108    0.1100    0.2134    0.2157]; %[0.5703    0.1100    0.3347    0.1577]; 
        elseif j == 3, fh1_children(j).Position = [0.1300    0.1100    0.2134    0.2157]; %[0.1300    0.1100    0.3347    0.1577]; 
        elseif j == 4, fh1_children(j).Position = [0.6916    0.4096    0.2134    0.2157]; %[0.5703    0.3291    0.3347    0.1577]; 
        elseif j == 5, fh1_children(j).Position = [0.4108    0.4096    0.2134    0.2157]; %[0.1300    0.3291    0.3347    0.1577]; 
        elseif j == 6, fh1_children(j).Position = [0.1300    0.4096    0.2134    0.2157]; %[0.5703    0.5482    0.3347    0.1577]; 
        elseif j == 7, fh1_children(j).Position = [0.6916    0.7093    0.2134    0.2157]; %[0.1300    0.5482    0.3347    0.1577]; 
        elseif j == 8, fh1_children(j).Position = [0.4108    0.7093    0.2134    0.2157]; %[0.5703    0.7673    0.3347    0.1577]; 
        elseif j == 9, fh1_children(j).Position = [0.1300    0.7093    0.2134    0.2157]; %[0.1300    0.7673    0.3347    0.1577]; 
        end
        
        % new changes in each subplot
        if IsAxes(fh1_children(j))
            % normalize content inside plots
%             % adding grid
%             fh1_children(j).XGrid = 'on';
%             fh1_children(j).YGrid = 'on';
            % x axis label
            if((j == 1)||(j == 3))
                fh1_children(j).XLabel.String = 'dilution factor';
            else
                fh1_children(j).XLabel.String = '';
            end
            % y axis label
            if((j == 1)||(j == 4)||(j == 7))
                fh1_children(j).YAxis(2).Label.String = {'Vmax_{corrected}','(umol mg_{P}^{-1} min^{-1})'};
                fh1_children(j).YAxis(1).Label.String = '';
            elseif((j == 3)||(j == 6)||(j == 9))
                fh1_children(j).YAxis(2).Label.String = {''};
                fh1_children(j).YAxis(1).Label.String ={'Vmax_{uncorrected}','(umol mg_{P}^{-1} min^{-1})'};
            elseif(j == 5)
                fh1_children(j).YAxis(2).Label.String = {''};
                fh1_children(j).YAxis(1).Label.String = {''};
            end
            fh1_children(j).YAxis(2).Label.FontSize = 13;
            fh1_children(j).YAxis(1).Label.FontSize = 13;
%             % axes
%             fh1_children(j).XLim = [0 1];
%             fh1_children(j).YLim = [0 1];
%             fh1_children(j).XTick = [0 0.125 0.25 0.5 1];
%             % common edits
%             fh1_children(j).YTick = [fh1_children(j).YLim(1) 1/4*fh1_children(j).YLim(2) 2/4*fh1_children(j).YLim(2) 3/4*fh1_children(j).YLim(2) fh1_children(j).YLim(2)];
%             fh1_children(j).XTickLabel = cellstr(num2str(fh1_children(j).XTick'));
%             fh1_children(j).YTickLabel = cellstr(num2str(fh1_children(j).YTick'));
%             % text 'pH' <- title label
%             fh1_children(j).Title.FontWeight = 'normal';
%             fh1_children(j).Title.FontSize = 11;
%             fh1_children(j).Title.Position = [0.25*fh1_children(j).XLim(2) 1.10*fh1_children(j).YLim(2) 0];
        end
        
% % % %         % remove extra DF-labels
% % % %         fh1_children(1).Children(10).String = '';
% % % %         fh1_children(3).Children(10).String = '';
% % % %         fh1_children(4).Children(10).String = '';
% % % %         fh1_children(8).Children(10).String = '';
% % % %         fh1_children(9).Children(10).String = '';
% % % %         fh1_children(9).Children(7).String = '';
        
        % delete suptitle and some reordering
        delete(fh1_children(2));
        
    end
%     fh1.Position = [fh1.Position(1) fh1.Position(2) 0.9*fh1.Position(3) 0.9*fh1.Position(4) ];
    fh1.Position = [0.025    0.075    0.6    0.8225]; % = [0.025    0.075    0.3    0.8225];
    
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%     % Saving: save them back to .png in location
%     savefig(1,'1appendix_ald_vm_dil');   
%     saveLoc = 'D:\OneDrive - TU Eindhoven\Documents\ch3_pHkinetics\results\manuscriptFigures\';
%     saveName_train_pdf = [saveLoc,'1appendix_ald_vm_dil.pdf'];
%     saveName_train_png = [saveLoc,'1appendix_ald_vm_dil.png'];
%     saveas(1,saveName_train_pdf);
%     saveas(1,saveName_train_png);

    % 
    if setup.saveOutput == 1
        % save
        savefig(1,'1appendix_ald_vm_dil'); 
        % specs printing (method 3)
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        print -dpdf -painters 1appendix_ald_vm_dil
        print -dpng -painters 1appendix_ald_vm_dil
    end
    
end

% %%
% figure(2)
% for j = 1:8
%     sp_temp = subplot(3,3,j);
%     disp(sp_temp.Position)
% end

% %%
% temp = fh1_children(1);
% for j = 1:6
%     figure
%     plot(temp.Children(j).XData,...
%     temp.Children(j).YData,'o-')
% end
% % fh1_children(1).YAxis(2).Children

