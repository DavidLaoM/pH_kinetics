% to locate better here
% dev_finalFigure_hxk_kEq_effect;
load('new_tempRes_figure_Keq_hxk.mat')
% dev_finalFigure_pyk_kEq_effect;
load('new_tempRes_figure_Keq_pyk.mat');


%% final plot
blueTriplet = [0, 0.4470, 0.7410];
orangeTriplet = [0.8500, 0.3250, 0.0980];
greyTriplet = [0.3 0.3 0.3];
rotationValue = 40 * ones(1,18);
    rotationValue(1) = 15;
    rotationValue(2) = 20;
    rotationValue(3) = 25;
% clf(100), 
figure(100)

% subplot(2,2,1); % vmHXK
subplot(2,2,1); % vmHXK
    load('new_tempRes_figure_Keq_hxk.mat');
% % % %     load('new_tempRes_figure_Keq_hxk_mod2.mat');
    legNames = cell(length(tempResult),1);
    c = cool(length(tempResult));
    c(3,:) = [0 0 0];
        
    % highlight in grey the meaningful area
%     x2 = [tempResult{end}.t', fliplr(tempResult{end}.t')];
    x2 = [tempResult{8}.t', fliplr(tempResult{2}.t')];
%     inBetween = [tempResult{2}.y(:,2)'/2, tempResult{6}.y(:,2)'/2];
    inBetween = [tempResult{8}.y(:,7)', fliplr(tempResult{2}.y(:,7)')];
    f1 = fill(x2, inBetween, 'g');
    f1.FaceColor = [0.90 0.90 0.90];
    f1.EdgeColor = 'none';
    hold on
    
    for i = 1:length(tempResult)
        legNames{i} = mat2str(keqvalsTested(i));
% % % %         if i == 3
        if((i == 1)||(i == 3)||(i == 18))
            plot(tempResult{i}.t,tempResult{i}.y(:,7),'-','color','k','linewidth',1.2)
        else
            plot(tempResult{i}.t,tempResult{i}.y(:,7),'-','color',greyTriplet)
%             plot(tempResult{i}.t,tempResult{i}.y(:,7),'-','color','k')
        end
%         plot(tempResult{i}.t,tempResult{i}.y(:,2),'.-','color',c(i,:))
        hold on
        xval = tempResult{i}.t(end-3); yval = tempResult{i}.y(end-3,7);% find the point at time 250, gettings its idxs. Get that as location.
        tempText = erase(sprintf('k_{eq} = %d',keqvalsTested(i)),".000000"); % prepare test keq
        if((i == 1)||(i == 3)||(i == 18)) % check when to stop
            if (i == 18)
                ht = text(xval,yval+0.005,tempText);% print there the k-value
            else
%                 ht = text(xval-10,yval+0.005,tempText);% print there the k-value
                ht = text(xval-10,yval-0.005,tempText);% print there the k-value
            end
            set(ht,'Rotation',rotationValue(i)); % determine rotation
        end
        
    end
    ylabel('NADPH concentration (mM)');xlabel('assay time (s)');

% subplot(2,2,2); % vmHXK
subplot(2,2,2); % vmHXK
    load('new_tempRes_figure_Keq_hxk.mat');
% % % %     load('new_tempRes_figure_Keq_hxk_mod2.mat');
    % highlight in grey the meaningful area
    x2 = [output_hxk_changing_pH.pHarray', fliplr(output_hxk_changing_pH.pHarray')];
    inBetween = [output_hxk_changing_pH.vm_uChange', fliplr(output_hxk_constant_pH.vm_uChange')];
    f1 = fill(x2, inBetween, 'g');
    f1.FaceColor = [0.90 0.90 0.90];
    f1.EdgeColor = 'none';
    hold on
% % % %     plot(output_hxk_changing_pH.pHarray, output_hxk_changing_pH.vm_uChange,'k.-','Linewidth',1.2,'MarkerSize',10)
    errorbar(output_hxk_changing_pH.pHarray, output_hxk_changing_pH.vm_uChange, 0.025*ones(size(output_hxk_changing_pH.vm_uChange)),'k.-','Linewidth',1.2,'MarkerSize',10)
    hold on
% % % %     plot(output_hxk_changing_pH.pHarray, output_hxk_constant_pH.vm_uChange,'.-','color',[0.5 0.5 0.5],'Linewidth',1.2,'MarkerSize',10)
    errorbar(output_hxk_changing_pH.pHarray, output_hxk_constant_pH.vm_uChange, 0.025*ones(size(output_hxk_constant_pH.vm_uChange)),'.-','color',[0.5 0.5 0.5],'Linewidth',1.2,'MarkerSize',10)
    ylabel('Enzyme capacity (umol mgP^{-1} min^{-1})'); xlabel('pH');
    
% subplot(2,2,3); % vmPYK
subplot(2,2,3); % vmPYK
    load('new_tempRes_figure_Keq_pyk.mat');
    legNames = cell(length(tempResult),1);
    c = cool(length(tempResult));
    c(8,:) = [0 0 0];
    
    % highlight in grey the meaningful area
    x2 = [tempResult{10}.t', fliplr(tempResult{11}.t')];
    inBetween = [tempResult{10}.y(:,2)', fliplr(tempResult{11}.y(:,2)')];
    f1 = fill(x2, inBetween, 'g');
    f1.FaceColor = [0.90 0.90 0.90];
    f1.EdgeColor = 'none';
    hold on
    
    % each plot and its style
    for i = 1:length(tempResult)
        legNames{i} = mat2str(keqvalsTested(i));
        
        if((i == 1)||(i == 6)||(i == 15))
% % % %         plot(tempResult{i}.t,tempResult{i}.y(:,2),'.-','color',c(i,:))
            plot(tempResult{i}.t,tempResult{i}.y(:,2),'-','color','k','linewidth',1.2)
        else
            plot(tempResult{i}.t,tempResult{i}.y(:,2),'-','color',greyTriplet)
        end
        
        hold on
        
        xval = tempResult{i}.t(end-3); yval = tempResult{i}.y(end-3,2);% find the point at time 250, gettings its idxs. Get that as location.
        tempText = erase(sprintf('k_{eq} = %d',keqvalsTested(i)),".000000"); % prepare test keq
        if((i == 1)||(i == 6)||(i == 15))
            if(i == 15)
                ht = text(xval-50,yval+0.005,'k_{eq} = 1e-9');% print there the k-value
            elseif(i == 1)
                ht = text(xval-50,yval+0.005,tempText);% print there the k-value
            else
                ht = text(xval-50,yval+0.005,tempText);% print there the k-value
            end
        end
        
    end
    annotation('textarrow',[0.225 0.275],[0.15 0.125],'String','grey area ','color',greyTriplet)
    ylim([0 0.08])
    ylabel('NADH concentration (mM)');xlabel('assay time (s)');

subplot(2,2,4);% vmPYK
    load('new_tempRes_figure_Keq_pyk.mat');
    % highlight in grey the meaningful area
    x2 = [output_pyk_changeKeq.pHarray', fliplr(output_pyk_constantKeq.pHarray')];
    inBetween = [output_pyk_changeKeq.vm_uChange', fliplr(output_pyk_constantKeq.vm_uChange')];
    f1 = fill(x2, inBetween, 'g');
    f1.FaceColor = [0.90 0.90 0.90];
    f1.EdgeColor = 'none';
    hold on
%     plot(output_pyk_changeKeq.pHarray, output_pyk_changeKeq.vm_uChange,'color',blueTriplet)
    errorbar(output_pyk_changeKeq.pHarray, output_pyk_changeKeq.vm_uChange, 0.125*ones(size(output_pyk_changeKeq.vm_uChange)),'k.-','Linewidth',1.2,'MarkerSize',10)
    hold on
%     plot(output_pyk_constantKeq.pHarray, output_pyk_constantKeq.vm_uChange,'color',orangeTriplet)
    errorbar(output_pyk_constantKeq.pHarray, output_pyk_constantKeq.vm_uChange, 0.125*ones(size(output_pyk_constantKeq.vm_uChange)),'.-','color',[0.5 0.5 0.5],'Linewidth',1.2,'MarkerSize',10)
    ylabel('Enzyme capacity (umol mgP^{-1} min^{-1})'); xlabel('pH');
    ylim([0 8])

% % % % [ax1,h1]=suplabel('Parameter Sensitivity Analysis                                           Parameter estimation');
% % % % [ax2,h2]=suplabel('PYK                                                                          HXK','y');
% [ax3,h2]=suplabel('super Y label (right)','yy');
% [ax4,h3]=suplabel('super Title'  ,'t');
% % % % set(h3,'FontSize',30)
orient portrait   
    
% suptitle('Final plot')
set(gcf,'color','w');


%% final plot
blueTriplet = [0, 0.4470, 0.7410];
orangeTriplet = [0.8500, 0.3250, 0.0980];
greyTriplet = [0.3 0.3 0.3];
rotationValue = 40 * ones(1,18);
    rotationValue(1) = 15;
    rotationValue(2) = 20;
    rotationValue(3) = 25;
colSteelBlue = [70/255 130/255 180/255]; % pH dependet
colLightBlue = [173/255 216/255 203/255]; % pH independent
    
    
figure(200)
% highlight in grey the meaningful area
% subplot(2,2,2); % vmHXK
load('new_tempRes_figure_Keq_hxk.mat');
eb1 = errorbar(output_hxk_changing_pH.pHarray, output_hxk_changing_pH.vm_uChange, 0.025*ones(size(output_hxk_changing_pH.vm_uChange)),'.-');
eb1.Color = colSteelBlue;
eb1.LineWidth = 1.2;
eb1.MarkerSize = 10;
eb1.CapSize = 2;
eb1.MarkerSize = 12;
hold on
eb2 = errorbar(output_hxk_changing_pH.pHarray, output_hxk_constant_pH.vm_uChange, 0.025*ones(size(output_hxk_constant_pH.vm_uChange)),'.-');
eb2.Color = colLightBlue;
eb2.LineWidth = 1.2;
eb2.MarkerSize = 10;
eb2.CapSize = 2;
eb2.MarkerSize = 12;
ylabel('Enzyme capacity (umol mgP^{-1} min^{-1})'); xlabel('pH');
hold off
set(gcf,'color','w');

