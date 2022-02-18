% 
clear
close all

%%
% for i = 1:length(enzymeName)
for i = 1 %[1 2] %[3 4] 
    
    % close to call by the same numbers
    close all
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % load the figure + rescale
    fh1 = openfig('ald_mw_R2_vs_movingWindow.fig');
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
            % x axis label
            if((j == 1)||(j == 3))
                fh1_children(j).XLabel.String = 'moving window size';
            else
                fh1_children(j).XLabel.String = '';
            end
            % y axis label
            if((j == 3)||(j == 6)||(j == 9))
                fh1_children(j).YLabel.String = 'R2';
            else
                fh1_children(j).YLabel.String = '';
            end
            % axes
            fh1_children(j).XLim = [0 100];
            fh1_children(j).YLim = [0 1];
            fh1_children(j).XTick = [0 50 100];
            % common edits
            fh1_children(j).YTick = [fh1_children(j).YLim(1) 2/4*fh1_children(j).YLim(2) fh1_children(j).YLim(2)];
            fh1_children(j).XTickLabel = cellstr(num2str(fh1_children(j).XTick'));
            fh1_children(j).YTickLabel = cellstr(num2str(fh1_children(j).YTick'));
            % text 'DF' <- text label (1 4 7 10)
            fh1_children(j).Children(1).Position(1) = 80;
            fh1_children(j).Children(4).Position(1) = 80;
            fh1_children(j).Children(7).Position(1) = 80;
            fh1_children(j).Children(10).Position(1) = 80;
            % text 'pH' <- title label
            fh1_children(j).Title.FontWeight = 'normal';
            fh1_children(j).Title.FontSize = 11;
            fh1_children(j).Title.Position = [0.25*fh1_children(j).XLim(2) 1.10*fh1_children(j).YLim(2) 0];
            
%             fh1_children(j).Children(1).Position = [0.75*fh1_children(j).XLim(2) 1.10*fh1_children(j).YLim(2) 0];
%             fh1_children(j).Children(1).HorizontalAlignment = 'center';
        end
        
        % remove extra DF-labels
        fh1_children(1).Children(10).String = '';
        fh1_children(3).Children(10).String = '';
        fh1_children(4).Children(10).String = '';
        fh1_children(8).Children(10).String = '';
        fh1_children(9).Children(10).String = '';
        fh1_children(9).Children(7).String = '';
        
        % delete suptitle and some reordering
        delete(fh1_children(2));
        
    end
%     fh1.Position = [fh1.Position(1) fh1.Position(2) 0.9*fh1.Position(3) 0.9*fh1.Position(4) ];
    fh1.Position = [0.025    0.075    0.6    0.8225]; % = [0.025    0.075    0.3    0.8225];
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % Saving: save them back to .png in location
    savefig(1,'1appendix_ald_mw_r2');   
    saveLoc = 'D:\OneDrive - TU Eindhoven\Documents\ch3_pHkinetics\results\manuscriptFigures\';
    saveName_train_pdf = [saveLoc,'1appendix_ald_mw_r2.pdf'];
    saveName_train_png = [saveLoc,'1appendix_ald_mw_r2.png'];
    saveas(1,saveName_train_pdf);
    saveas(1,saveName_train_png);
    
%     % save %% <- to use the PDF, fonts will be have to be adjusted, or change this save command.
%     savefig(1,'1appendix_ald_mw_r2'); 
%     % specs printing (method 3)
%     set(gcf,'Units','inches');
%     screenposition = get(gcf,'Position');
%     set(gcf,...
%         'PaperPosition',[0 0 screenposition(3:4)],...
%         'PaperSize',[screenposition(3:4)]);
%     print -dpdf -painters 1appendix_ald_mw_r2
%     print -dpng -painters 1appendix_ald_mw_r2
    

end

% %%
% figure(2)
% for j = 1:8
%     sp_temp = subplot(3,3,j);
%     disp(sp_temp.Position)
% end



