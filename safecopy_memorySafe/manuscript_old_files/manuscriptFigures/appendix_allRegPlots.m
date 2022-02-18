% % APPENDIX_ALLREGPLOTS
% 
enzymeNames = {'hxk';
                'pgi';
                'pfk';
                'ald';
                'tpi';
                'gapdh';
                'gapdhr';
                'pgm';
                'eno';
                'pyk';
                'pdc'};
selLamPosArray = [14;
                  14;
                  8;
                  17;
                  16;
                  15;
                  15;
                  13;
                  18;
                  17;
                  17];

figure(201)
for j = 1:length(enzymeNames)
    tempName = [enzymeNames{j},'_regularizationResults.mat']; % get the name
    load(tempName);% load the regularization results
    p1 = subplot(4,3,j);% subplot
    selLambdaPos = selLamPosArray(j);% selLamPos
    regularizationSimple;% run regularization simple
    %p1.YLim = [0 p1.YLim(2)];
    tempText = [enzymeNames{j}, ', lam=', sprintf('%d',lambdalist(selLambdaPos))];
    tempText2 = erase(tempText,".000000");
%     tempText = [enzymeNames{j}, ', lam=', erase(sprintf('%d',lambdalist(selLambdaPos),".000000"))];
    text(1E0,p1.YLim(2)*1.1,tempText2);
end
% suptitle('Regularization plots for all the enzymes')
set(201,'color','white')
savefig(201,'results/manuscriptFigures/appendixes_allRegPlots.fig');

%%
% ylabel('error_{Parameters} []') %left
% 
% ylabel('error_{Data} []') %right

%% memoryDump
% Cannot use copyobj since it cannot work with multiple axes (yyaxis left,
% right) and it seems it neither does for log scale.

% c=hgload('ald_regularization.fig');
% k=hgload('eno_regularization.fig');
% % Prepare subplots
% figure
% h(1)=subplot(1,2,1);
% h(2)=subplot(1,2,2);
% % Paste figures on the subplots
% copyobj(allchild(get(c,'CurrentAxes')),h(1));
% copyobj(allchild(get(k,'CurrentAxes')),h(2));
% % Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')
% 
% %%
% plot(peaks)
% a1 = gca
% f2 = figure
% a2 = copyobj(a1,f2) 
% 
% %%
% %%Create double yaxis
% ax1=axes('yaxislocation','left');hold on
% ax2=axes('yaxislocation','right','xcolor','none');hold on
% set([ax1 ax2],'color','none')
% linkaxes([ax1 ax2],'x')
% %%Call axes and plot like this
% axes(ax1)
% h1=semilogy([0 1],[0 1])
% axes(ax2)
% h2=semilogy([0 1],[1 0])
% %%Copy axes to new figure
% figure;
% copyobj(ax1,gcf)
% copyobj(ax2,gcf)