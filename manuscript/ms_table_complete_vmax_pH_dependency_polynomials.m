% % vmCompvsExp
% Generation of the manuscript figure that plots the experimental and
% computational estimates for the maximum reaction rates.

enzymeList = {'hxk_parEst';... %1 hxk
    'pgi_parEst';... %2 pgi
    'pfk_parEst';... %3 pfk
    'ald_parEst';... %4 ald
    'tpi_parEst';... %5 tpi
    'gapdh_parEst';... %6 gapdh_fwd
    'gapdhr_parEst';... %7 gapdh_rev
    'pgm_parEst';... %8 pgm
    'eno_parEst';... %9 eno_kmfixed
    'pyk_parEst';... %10 pyk
    'pdc_parEst'}; %11 pdc
numEnz = length(enzymeList);
enzymeName = {'hxk';... %1 hxk
    'pgi';... %2 pgi
    'pfk';... %3 pfk
    'ald';... %4 ald
    'tpi';... %5 tpi
    'gapdh_{fwd}';... %6 gapdh_fwd
    'gapdh_{rev}';... %7 gapdh_rev
    'pgm';... %8 pgm
    'eno';... %9 eno_kmfixed
    'pyk';... %10 pyk
    'pdc'}; %11 pdc
% % values that need to be multiplied by (1/2) (NOW IN LOOP)
% names2change = {'pgi';'pfk';'ald';'tpi'};
% val2multiply = (1/2);
% for i = 1:length(names2change)
%     tempName = ['output_',names2change{i}];
%     eval(['dataset = ',tempName,';']); % select the data
%     dataset.vm_uChange = dataset.vm_uChange * val2multiply;
%     dataset.vm_up_uChange = dataset.vm_up_uChange * val2multiply;
%     dataset.vm_down_uChange = dataset.vm_down_uChange * val2multiply;
%     dataset.Vmax_experimental = dataset.Vmax_experimental * val2multiply;
%     dataset.stDev_experimental = dataset.stDev_experimental * val2multiply;
%     % save back
%     eval([tempName,' = dataset;']);
% end

xlims_vals = ones(numEnz,2);
xticks_vals = ones(numEnz,5);
for i = 1:numEnz
    xlims_vals(i,:) = [6 8];
    xticks_vals(i,:) = [6 6.5 7 7.5 8];
end

% % 
% col1 = [189,215,231]/255;
% col2 = [107,174,214]/255;
% col3 = [49,130,189]/255;
col4 = [8,81,156]/255;

figure(400)
for i =1:numEnz
    % load and select the data
    load(enzymeList{i});
    subplot(4,3,i) % select subplot
    tempName = ['output_',extractBefore(enzymeList{i},"_parEst")];
    eval(['dataset = ',tempName,';']); % select the data
    
    % experimental or naïve parameter estimates
    errorbar(dataset.pHarray, dataset.Vmax_experimental, dataset.stDev_experimental,...
        'o-','MarkerSize',3,...%'LineWidth',1.2,...
        'MarkerFaceColor',col4,...
        'Color',col4,'CapSize',0)% plot errorbar experimentel estimation
    hold on % hold on
    
    % computational parameter estimates
    errorbar(dataset.pHarray, dataset.vm_uChange, dataset.vm_up_uChange-dataset.vm_uChange,...
        'ro:','MarkerSize',3,'MarkerFaceColor','r','CapSize',0)% plot errorbar computational estimation
    
    % readjust Y-axis
    ax = gca;
    ax.YLim(1) = 0;
    ax.YLim(2) = max([dataset.Vmax_experimental;dataset.vm_uChange])*1.1;
    
    % labels
    if((i == 1)||(i == 4)||(i == 7)||(i == 10))
        ylabel({'Vmax_{corrected}';'[\mumol mg_{P}^{-1} min^{-1}]'});
    end
    if((i == 10)||(i == 11))
        xlabel('pH value []');
    end
    title(enzymeName{i});
    
% %     if(i == 10)
% %         xlabel('pH','FontSize',16)
% %     end
% %     % ylabel
% % %     if((i == 1)||(i == 4)||(i == 7)||(i == 10))
% %     if(i == 1)
% %         ylabel('Vmax [\mumol mg_{P}^{-1} min^{-1}]','FontSize',16)
% %     end
    
end
set(gcf,'color','w');


% %%
% figure(200)%(1)
% for i =1:numEnz
%     load(enzymeList{i});
%     subplot(4,3,i) % select subplot
%     tempName = ['output_',extractBefore(enzymeList{i},"_parEst")];
%     eval(['dataset = ',tempName,';']); % select the data
%     % correcting units for the 1/2 factor
% % % % %     val2multiply = (1/2);           % NEW CORRECTION METHOD (2020 - 10 - 29)
%     val2multiply = (1);           % NEW CORRECTION METHOD (2020 - 10 - 29)
%     if((i==2)||(i==3)||(i==4)||(i==5))
%         tempName2 = ['output_',enzymeName{i}];
%         eval(['dataset = ',tempName2,';']); % select the data
%         dataset.vm_uChange = dataset.vm_uChange * val2multiply;
%         dataset.vm_up_uChange = dataset.vm_up_uChange * val2multiply;
%         dataset.vm_down_uChange = dataset.vm_down_uChange * val2multiply;
%         dataset.Vmax_experimental = dataset.Vmax_experimental * val2multiply;
%         dataset.stDev_experimental = dataset.stDev_experimental * val2multiply;
%         % save back
%         eval([tempName2,' = dataset;']);
%     end
%     % plotting
%     if i == i
%         plot(dataset.pHarray, dataset.vm_uChange,...
%             'bo:','MarkerSize',5,'MarkerFaceColor','b')% plot errorbar computational estimation
%     else
%         errorbar(dataset.pHarray, dataset.vm_uChange, dataset.vm_up_uChange-dataset.vm_uChange,...
%             'bo:','MarkerSize',5,'MarkerFaceColor','b','CapSize',0)% plot errorbar computational estimation
%     end
%     hold on % hold on
%     % tempName = extractAfter(enzymeList{i},"output_");
%     errorbar(dataset.pHarray, dataset.Vmax_experimental, dataset.stDev_experimental,...
%         'rs','MarkerSize',5,'LineWidth',1.2,'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',0)% plot errorbar experimentel estimation
%     ax = gca;
%     ax.XLim = xlims_vals(i,:); % xlimits
%     ax.XTick = xticks_vals(i,:); % xticklabels
%     ax.YLim(1) = 0; % ylim(1) = 0
%     tempYLim = [max(dataset.Vmax_experimental+dataset.stDev_experimental) max(dataset.vm_up_uChange)];
%     if i ~= 1
%         ax.YLim(2) = max(tempYLim)*1.1; % ylim(1) = 0
%     end
%     %     tempName2 = [strrep(tempName,'_','_{'),'}'];
%     tempName2 = extractBefore(enzymeList{i},"_parEst");
%     text(ax.XLim(1)+0.1,ax.YLim(2)*0.9,tempName2)% tag
% %     if((i ~= 10)||(i ~= 11))
%     if(i <= 9)
%         set(gca,'xticklabel',{[]});
%     end
%     if i == 11
% % % % %         lgd = legend('computational','experimental');
% % % % %         lgd.Position = lgd.Position + [0.3 0 0 0];
%     end
%     % xlabel
% %     if((i == 10)||(i == 11))
%     if(i == 10)
%         xlabel('pH','FontSize',16)
%     end
%     % ylabel
% %     if((i == 1)||(i == 4)||(i == 7)||(i == 10))
%     if(i == 1)
%         ylabel('Vmax [\mumol mg_{P}^{-1} min^{-1}]','FontSize',16)
%     end
%     
% end
% set(gcf,'color','w');


%%
% recall and structure the data (converting for the  '* (1/2)' here)
fullpHarray = output_eno.pHarray;
for i = 1 % write here the empty arrays and close it
    emptyAarray = zeros(size(fullpHarray));
    % ald
    ald_sim_vmax = emptyAarray;
    ald_sim_std = emptyAarray;
    ald_exp_vmax = emptyAarray;
    ald_exp_std = emptyAarray;
    % eno
    eno_sim_vmax = emptyAarray;
    eno_sim_std = emptyAarray;
    eno_exp_vmax = emptyAarray;
    eno_exp_std = emptyAarray;
    % gapdh
    gapdh_sim_vmax = emptyAarray;
    gapdh_sim_std = emptyAarray;
    gapdh_exp_vmax = emptyAarray;
    gapdh_exp_std = emptyAarray;
    % gapdhr
    gapdhr_sim_vmax = emptyAarray;
    gapdhr_sim_std = emptyAarray;
    gapdhr_exp_vmax = emptyAarray;
    gapdhr_exp_std = emptyAarray;
    % hxk
    hxk_sim_vmax = emptyAarray;
    hxk_sim_std = emptyAarray;
    hxk_exp_vmax = emptyAarray;
    hxk_exp_std = emptyAarray;
    % pdc
    pdc_sim_vmax = emptyAarray;
    pdc_sim_std = emptyAarray;
    pdc_exp_vmax = emptyAarray;
    pdc_exp_std = emptyAarray;
    % pfk
    pfk_sim_vmax = emptyAarray;
    pfk_sim_std = emptyAarray;
    pfk_exp_vmax = emptyAarray;
    pfk_exp_std = emptyAarray;
    % pgi
    pgi_sim_vmax = emptyAarray;
    pgi_sim_std = emptyAarray;
    pgi_exp_vmax = emptyAarray;
    pgi_exp_std = emptyAarray;
    % pgm
    pgm_sim_vmax = emptyAarray;
    pgm_sim_std = emptyAarray;
    pgm_exp_vmax = emptyAarray;
    pgm_exp_std = emptyAarray;
    % pyk
    pyk_sim_vmax = emptyAarray;
    pyk_sim_std = emptyAarray;
    pyk_exp_vmax = emptyAarray;
    pyk_exp_std = emptyAarray;
    % tpi
    tpi_sim_vmax = emptyAarray;
    tpi_sim_std = emptyAarray;
    tpi_exp_vmax = emptyAarray;
    tpi_exp_std = emptyAarray;
end
for i = 1:length(fullpHarray)
    pHtemp = fullpHarray(i);
    % ald
    idx = find(output_ald.pHarray == fullpHarray(i));
    if idx >= 1
        ald_sim_vmax(i) = output_ald.vm_uChange(idx) * (1); % not corrected before
        ald_sim_std(i) = output_eno.vm_up_uChange(idx) * (1)-output_eno.vm_uChange(idx) * (1); % not corrected before
        ald_exp_vmax(i) = output_ald.Vmax_experimental(idx) * (1); % not corrected before
        ald_exp_std(i) = output_ald.stDev_experimental(idx) * (1); % not corrected before
    else
        ald_sim_vmax(i) = 0;
        ald_sim_std(i) = 0;
        ald_exp_vmax(i) = 0;
        ald_exp_std(i) = 0;
    end
    % eno
    idx = find(output_eno.pHarray == fullpHarray(i));
    if idx >= 1
        eno_sim_vmax(i) = output_eno.vm_uChange(idx);
        eno_sim_std(i) = output_eno.vm_up_uChange(idx)-output_eno.vm_uChange(idx);
        eno_exp_vmax(i) = output_eno.Vmax_experimental(idx);
        eno_exp_std(i) = output_eno.stDev_experimental(idx);
    else
        eno_sim_vmax(i) = 0;
        eno_sim_std(i) = 0;
        eno_exp_vmax(i) = 0;
        eno_exp_std(i) = 0;
    end
    % gapdh
    idx = find(output_gapdh.pHarray == fullpHarray(i));
    if idx >= 1
        gapdh_sim_vmax(i) = output_gapdh.vm_uChange(idx);
        gapdh_sim_std(i) = output_eno.vm_up_uChange(idx)-output_eno.vm_uChange(idx);
        gapdh_exp_vmax(i) = output_gapdh.Vmax_experimental(idx);
        gapdh_exp_std(i) = output_gapdh.stDev_experimental(idx);
    else
        gapdh_sim_vmax(i) = 0;
        gapdh_sim_std(i) = 0;
        gapdh_exp_vmax(i) = 0;
        gapdh_exp_std(i) = 0;
    end
    % gapdhr
    idx = find(output_gapdhr.pHarray == fullpHarray(i));
    if idx >= 1
        gapdhr_sim_vmax(i) = output_gapdhr.vm_uChange(idx);
        gapdhr_sim_std(i) = output_gapdhr.vm_up_uChange(idx)-output_gapdhr.vm_uChange(idx);
        gapdhr_exp_vmax(i) = output_gapdhr.Vmax_experimental(idx);
        gapdhr_exp_std(i) = output_gapdhr.stDev_experimental(idx);
    else
        gapdhr_sim_vmax(i) = 0;
        gapdhr_sim_std(i) = 0;
        gapdhr_exp_vmax(i) = 0;
        gapdhr_exp_std(i) = 0;
    end
    % hxk
    idx = find(output_hxk.pHarray == fullpHarray(i));
    if idx >= 1
        hxk_sim_vmax(i) = output_hxk.vm_uChange(idx);
        hxk_sim_std(i) = output_hxk.vm_up_uChange(idx)-output_hxk.vm_uChange(idx);
        hxk_exp_vmax(i) = output_hxk.Vmax_experimental(idx);
        hxk_exp_std(i) = output_hxk.stDev_experimental(idx);
    else
        hxk_sim_vmax(i) = 0;
        hxk_sim_std(i) = 0;
        hxk_exp_vmax(i) = 0;
        hxk_exp_std(i) = 0;
    end
    % pdc
    idx = find(output_pdc.pHarray == fullpHarray(i));
    if idx >= 1
        pdc_sim_vmax(i) = output_pdc.vm_uChange(idx);
        pdc_sim_std(i) = output_pdc.vm_up_uChange(idx)-output_pdc.vm_uChange(idx);
        pdc_exp_vmax(i) = output_pdc.Vmax_experimental(idx);
        pdc_exp_std(i) = output_pdc.stDev_experimental(idx);
    else
        pdc_sim_vmax(i) = 0;
        pdc_sim_std(i) = 0;
        pdc_exp_vmax(i) = 0;
        pdc_exp_std(i) = 0;
    end
    % pfk
    idx = find(output_pfk.pHarray == fullpHarray(i));
    if idx >= 1
        pfk_sim_vmax(i) = output_pfk.vm_uChange(idx) * (1); % not corrected before
        pfk_sim_std(i) = output_pfk.vm_up_uChange(idx) * (1)-output_pfk.vm_uChange(idx) * (1); % not corrected before
        pfk_exp_vmax(i) = output_pfk.Vmax_experimental(idx) * (1); % not corrected before
        pfk_exp_std(i) = output_pfk.stDev_experimental(idx) * (1); % not corrected before
    else
        pfk_sim_vmax(i) = 0;
        pfk_sim_std(i) = 0;
        pfk_exp_vmax(i) = 0;
        pfk_exp_std(i) = 0;
    end
    % pgi
    idx = find(output_pgi.pHarray == fullpHarray(i));
    if idx >= 1
        pgi_sim_vmax(i) = output_pgi.vm_uChange(idx) * (1); % not corrected before
        pgi_sim_std(i) = output_pgi.vm_up_uChange(idx) * (1)-output_pgi.vm_uChange(idx) * (1); % not corrected before
        pgi_exp_vmax(i) = output_pgi.Vmax_experimental(idx) * (1); % not corrected before
        pgi_exp_std(i) = output_pgi.stDev_experimental(idx) * (1); % not corrected before
    else
        pgi_sim_vmax(i) = 0;
        pgi_sim_std(i) = 0;
        pgi_exp_vmax(i) = 0;
        pgi_exp_std(i) = 0;
    end
    % pgm
    idx = find(output_pgm.pHarray == fullpHarray(i));
    if idx >= 1
        pgm_sim_vmax(i) = output_pgm.vm_uChange(idx);
        pgm_sim_std(i) = output_pgm.vm_up_uChange(idx)-output_pgm.vm_uChange(idx);
        pgm_exp_vmax(i) = output_pgm.Vmax_experimental(idx);
        pgm_exp_std(i) = output_pgm.stDev_experimental(idx);
    else
        pgm_sim_vmax(i) = 0;
        pgm_sim_std(i) = 0;
        pgm_exp_vmax(i) = 0;
        pgm_exp_std(i) = 0;
    end
    % pyk
    idx = find(output_pyk.pHarray == fullpHarray(i));
    if idx >= 1
        pyk_sim_vmax(i) = output_pyk.vm_uChange(idx);
        pyk_sim_std(i) = output_pyk.vm_up_uChange(idx)-output_pyk.vm_uChange(idx);
        pyk_exp_vmax(i) = output_pyk.Vmax_experimental(idx);
        pyk_exp_std(i) = output_pyk.stDev_experimental(idx);
    else
        pyk_sim_vmax(i) = 0;
        pyk_sim_std(i) = 0;
        pyk_exp_vmax(i) = 0;
        pyk_exp_std(i) = 0;
    end
    % tpi
    idx = find(output_tpi.pHarray == fullpHarray(i));
    if idx >= 1
        tpi_sim_vmax(i) = output_tpi.vm_uChange(idx) * (1); % not corrected before
        tpi_sim_std(i) = output_tpi.vm_up_uChange(idx) * (1)-output_tpi.vm_uChange(idx) * (1); % not corrected before
        tpi_exp_vmax(i) = output_tpi.Vmax_experimental(idx) * (1); % not corrected before
        tpi_exp_std(i) = output_tpi.stDev_experimental(idx) * (1); % not corrected before
    else
        tpi_sim_vmax(i) = 0;
        tpi_sim_std(i) = 0;
        tpi_exp_vmax(i) = 0;
        tpi_exp_std(i) = 0;
    end
end % fill in the arrays

% create the table and excel file
T = table(fullpHarray,ald_sim_vmax,ald_sim_std,ald_exp_vmax,ald_exp_std,eno_sim_vmax,eno_sim_std,eno_exp_vmax,eno_exp_std,gapdh_sim_vmax,gapdh_sim_std,gapdh_exp_vmax,gapdh_exp_std,gapdhr_sim_vmax,gapdhr_sim_std,gapdhr_exp_vmax,gapdhr_exp_std,hxk_sim_vmax,hxk_sim_std,hxk_exp_vmax,hxk_exp_std,pdc_sim_vmax,pdc_sim_std,pdc_exp_vmax,pdc_exp_std,pfk_sim_vmax,pfk_sim_std,pfk_exp_vmax,pfk_exp_std,pgi_sim_vmax,pgi_sim_std,pgi_exp_vmax,pgi_exp_std,pgm_sim_vmax,pgm_sim_std,pgm_exp_vmax,pgm_exp_std,pyk_sim_vmax,pyk_sim_std,pyk_exp_vmax,pyk_exp_std,tpi_sim_vmax,tpi_sim_std,tpi_exp_vmax,tpi_exp_std);
T.Properties.VariableNames = {'fullpHarray','ald_sim_vmax','ald_sim_std','ald_exp_vmax','ald_exp_std','eno_sim_vmax','eno_sim_std','eno_exp_vmax','eno_exp_std','gapdh_sim_vmax','gapdh_sim_std','gapdh_exp_vmax','gapdh_exp_std','gapdhr_sim_vmax','gapdhr_sim_std','gapdhr_exp_vmax','gapdhr_exp_std','hxk_sim_vmax','hxk_sim_std','hxk_exp_vmax','hxk_exp_std','pdc_sim_vmax','pdc_sim_std','pdc_exp_vmax','pdc_exp_std','pfk_sim_vmax','pfk_sim_std','pfk_exp_vmax','pfk_exp_std','pgi_sim_vmax','pgi_sim_std','pgi_exp_vmax','pgi_exp_std','pgm_sim_vmax','pgm_sim_std','pgm_exp_vmax','pgm_exp_std','pyk_sim_vmax','pyk_sim_std','pyk_exp_vmax','pyk_exp_std','tpi_sim_vmax','tpi_sim_std','tpi_exp_vmax','tpi_exp_std'};

% %% % In case these figures are to be saved
% filename = '20201217_parameterValues.xlsx';
% writetable(T,filename,'Sheet',1,'Range','A1');

%%
% Easier visualization
% hxk
hxk_diff_sim_exp = hxk_sim_vmax - hxk_exp_vmax;
hxk_percent_exp = hxk_diff_sim_exp./hxk_exp_vmax;
hxk_percent_exp_abs_100 = abs(hxk_percent_exp) * 100;
hxk_percent_exp_abs_100(isnan(hxk_percent_exp_abs_100))=0;
% pgi
pgi_diff_sim_exp = pgi_sim_vmax - pgi_exp_vmax;
pgi_percent_exp = pgi_diff_sim_exp./pgi_exp_vmax;
pgi_percent_exp_abs_100 = abs(pgi_percent_exp) * 100;
pgi_percent_exp_abs_100(isnan(pgi_percent_exp_abs_100))=0;
% pfk
pfk_diff_sim_exp = pfk_sim_vmax - pfk_exp_vmax;
pfk_percent_exp = pfk_diff_sim_exp./pfk_exp_vmax;
pfk_percent_exp_abs_100 = abs(pfk_percent_exp) * 100;
pfk_percent_exp_abs_100(isnan(pfk_percent_exp_abs_100))=0;
% ald
ald_diff_sim_exp = ald_sim_vmax - ald_exp_vmax;
ald_percent_exp = ald_diff_sim_exp./ald_exp_vmax;
ald_percent_exp_abs_100 = abs(ald_percent_exp) * 100;
ald_percent_exp_abs_100(isnan(ald_percent_exp_abs_100))=0;
% tpi
tpi_diff_sim_exp = tpi_sim_vmax - tpi_exp_vmax;
tpi_percent_exp = tpi_diff_sim_exp./tpi_exp_vmax;
tpi_percent_exp_abs_100 = abs(tpi_percent_exp) * 100;
tpi_percent_exp_abs_100(isnan(tpi_percent_exp_abs_100))=0;
% gapdh_fwd
gapdh_diff_sim_exp = gapdh_sim_vmax - gapdh_exp_vmax;
gapdh_percent_exp = gapdh_diff_sim_exp./gapdh_exp_vmax;
gapdh_percent_exp_abs_100 = abs(gapdh_percent_exp) * 100;
gapdh_percent_exp_abs_100(isnan(gapdh_percent_exp_abs_100))=0;
% gapdh_rev
gapdhr_diff_sim_exp = gapdhr_sim_vmax - gapdhr_exp_vmax;
gapdhr_percent_exp = gapdhr_diff_sim_exp./gapdhr_exp_vmax;
gapdhr_percent_exp_abs_100 = abs(gapdhr_percent_exp) * 100;
gapdhr_percent_exp_abs_100(isnan(gapdhr_percent_exp_abs_100))=0;
% pgm
pgm_diff_sim_exp = pgm_sim_vmax - pgm_exp_vmax;
pgm_percent_exp = pgm_diff_sim_exp./pgm_exp_vmax;
pgm_percent_exp_abs_100 = abs(pgm_percent_exp) * 100;
pgm_percent_exp_abs_100(isnan(pgm_percent_exp_abs_100))=0;
% eno
eno_diff_sim_exp = eno_sim_vmax - eno_exp_vmax;
eno_percent_exp = eno_diff_sim_exp./eno_exp_vmax;
eno_percent_exp_abs_100 = abs(eno_percent_exp) * 100;
eno_percent_exp_abs_100(isnan(eno_percent_exp_abs_100))=0;
% pyk
pyk_diff_sim_exp = pyk_sim_vmax - pyk_exp_vmax;
pyk_percent_exp = pyk_diff_sim_exp./pyk_exp_vmax;
pyk_percent_exp_abs_100 = abs(pyk_percent_exp) * 100;
pyk_percent_exp_abs_100(isnan(pyk_percent_exp_abs_100))=0;
% pdc
pdc_diff_sim_exp = pdc_sim_vmax - pdc_exp_vmax;
pdc_percent_exp = pdc_diff_sim_exp./pdc_exp_vmax;
pdc_percent_exp_abs_100 = abs(pdc_percent_exp) * 100;
pdc_percent_exp_abs_100(isnan(pdc_percent_exp_abs_100))=0;
% T = table(fullpHarray,ald_sim_vmax,ald_sim_std,ald_exp_vmax,ald_exp_std,eno_sim_vmax,eno_sim_std,eno_exp_vmax,eno_exp_std,gapdh_sim_vmax,gapdh_sim_std,gapdh_exp_vmax,gapdh_exp_std,gapdhr_sim_vmax,gapdhr_sim_std,gapdhr_exp_vmax,gapdhr_exp_std,hxk_sim_vmax,hxk_sim_std,hxk_exp_vmax,hxk_exp_std,pdc_sim_vmax,pdc_sim_std,pdc_exp_vmax,pdc_exp_std,pfk_sim_vmax,pfk_sim_std,pfk_exp_vmax,pfk_exp_std,pgi_sim_vmax,pgi_sim_std,pgi_exp_vmax,pgi_exp_std,pgm_sim_vmax,pgm_sim_std,pgm_exp_vmax,pgm_exp_std,pyk_sim_vmax,pyk_sim_std,pyk_exp_vmax,pyk_exp_std,tpi_sim_vmax,tpi_sim_std,tpi_exp_vmax,tpi_exp_std);
T2 = table(fullpHarray,...
    hxk_sim_vmax,hxk_sim_std,hxk_exp_vmax,hxk_exp_std,hxk_percent_exp_abs_100,...
    pgi_sim_vmax,pgi_sim_std,pgi_exp_vmax,pgi_exp_std,pgi_percent_exp_abs_100,...
    pfk_sim_vmax,pfk_sim_std,pfk_exp_vmax,pfk_exp_std,pfk_percent_exp_abs_100,...
    ald_sim_vmax,ald_sim_std,ald_exp_vmax,ald_exp_std,ald_percent_exp_abs_100,...
    tpi_sim_vmax,tpi_sim_std,tpi_exp_vmax,tpi_exp_std,tpi_percent_exp_abs_100,...
    gapdh_sim_vmax,gapdh_sim_std,gapdh_exp_vmax,gapdh_exp_std,gapdh_percent_exp_abs_100,...
    gapdhr_sim_vmax,gapdhr_sim_std,gapdhr_exp_vmax,gapdhr_exp_std,gapdhr_percent_exp_abs_100,...
    pgm_sim_vmax,pgm_sim_std,pgm_exp_vmax,pgm_exp_std,pgm_percent_exp_abs_100,...
    eno_sim_vmax,eno_sim_std,eno_exp_vmax,eno_exp_std,eno_percent_exp_abs_100,...
    pyk_sim_vmax,pyk_sim_std,pyk_exp_vmax,pyk_exp_std,pyk_percent_exp_abs_100,...
    pdc_sim_vmax,pdc_sim_std,pdc_exp_vmax,pdc_exp_std,pdc_percent_exp_abs_100);
T2.Properties.VariableNames = {'fullpHarray',...
    'hxk_sim_vmax','hxk_sim_std','hxk_exp_vmax','hxk_exp_std','hxk_percent_diff',...
    'pgi_sim_vmax','pgi_sim_std','pgi_exp_vmax','pgi_exp_std','pgi_percent_diff',...
    'pfk_sim_vmax','pfk_sim_std','pfk_exp_vmax','pfk_exp_std','pfk_percent_diff',...
    'ald_sim_vmax','ald_sim_std','ald_exp_vmax','ald_exp_std','ald_percent_diff',...
    'tpi_sim_vmax','tpi_sim_std','tpi_exp_vmax','tpi_exp_std','tpi_percent_diff',...
    'gapdh_sim_vmax','gapdh_sim_std','gapdh_exp_vmax','gapdh_exp_std','gapdh_percent_diff',...
    'gapdhr_sim_vmax','gapdhr_sim_std','gapdhr_exp_vmax','gapdhr_exp_std','gapdhr_percent_diff',...
    'pgm_sim_vmax','pgm_sim_std','pgm_exp_vmax','pgm_exp_std','pgm_percent_diff',...
    'eno_sim_vmax','eno_sim_std','eno_exp_vmax','eno_exp_std','eno_percent_diff',...
    'pyk_sim_vmax','pyk_sim_std','pyk_exp_vmax','pyk_exp_std','pyk_percent_diff',...
    'pdc_sim_vmax','pdc_sim_std','pdc_exp_vmax','pdc_exp_std','pdc_percent_diff',...
    };

%%
filename = '20210224_parameterValues_differences.xlsx';
% writetable(T2,filename,'Sheet',1,'Range','A1');

%% calculate maximum diferences between experimental and computational
% hxk
hxk_diff_sim_exp = hxk_sim_vmax - hxk_exp_vmax;
hxk_percent_exp = hxk_diff_sim_exp./hxk_exp_vmax;
T1 = table(hxk_sim_vmax,hxk_exp_vmax,hxk_diff_sim_exp,hxk_percent_exp);
T1.Properties.VariableNames = {'hxk_sim_vmax','hxk_exp_vmax','hxk_diff_sim_exp','hxk_percent_exp'};
sprintf('Max. deviation hxk: %d percent',max(abs(hxk_percent_exp))*100)
% pgi
pgi_diff_sim_exp = pgi_sim_vmax - pgi_exp_vmax;
pgi_percent_exp = pgi_diff_sim_exp./pgi_exp_vmax;
T1 = table(pgi_sim_vmax,pgi_exp_vmax,pgi_diff_sim_exp,pgi_percent_exp);
T1.Properties.VariableNames = {'pgi_sim_vmax','pgi_exp_vmax','pgi_diff_sim_exp','pgi_percent_exp'};
sprintf('Max. deviation pgi: %d percent',max(abs(pgi_percent_exp))*100)
% pfk
pfk_diff_sim_exp = pfk_sim_vmax - pfk_exp_vmax;
pfk_percent_exp = pfk_diff_sim_exp./pfk_exp_vmax;
T1 = table(pfk_sim_vmax,pfk_exp_vmax,pfk_diff_sim_exp,pfk_percent_exp);
T1.Properties.VariableNames = {'pfk_sim_vmax','pfk_exp_vmax','pfk_diff_sim_exp','pfk_percent_exp'};
sprintf('Max. deviation pfk: %d percent',max(abs(pfk_percent_exp))*100)
% ald
ald_diff_sim_exp = ald_sim_vmax - ald_exp_vmax;
ald_percent_exp = ald_diff_sim_exp./ald_exp_vmax;
T1 = table(ald_sim_vmax,ald_exp_vmax,ald_diff_sim_exp,ald_percent_exp);
T1.Properties.VariableNames = {'ald_sim_vmax','ald_exp_vmax','ald_diff_sim_exp','ald_percent_exp'};
sprintf('Max. deviation ALD: %d percent',max(abs(ald_percent_exp))*100)
% tpi
tpi_diff_sim_exp = tpi_sim_vmax - tpi_exp_vmax;
tpi_percent_exp = tpi_diff_sim_exp./tpi_exp_vmax;
T1 = table(tpi_sim_vmax,tpi_exp_vmax,tpi_diff_sim_exp,tpi_percent_exp);
T1.Properties.VariableNames = {'tpi_sim_vmax','tpi_exp_vmax','tpi_diff_sim_exp','tpi_percent_exp'};
sprintf('Max. deviation tpi: %d percent',max(abs(tpi_percent_exp))*100)
% gapdh_fwd
gapdh_diff_sim_exp = gapdh_sim_vmax - gapdh_exp_vmax;
gapdh_percent_exp = gapdh_diff_sim_exp./gapdh_exp_vmax;
T1 = table(gapdh_sim_vmax,gapdh_exp_vmax,gapdh_diff_sim_exp,gapdh_percent_exp);
T1.Properties.VariableNames = {'gapdh_sim_vmax','gapdh_exp_vmax','gapdh_diff_sim_exp','gapdh_percent_exp'};
sprintf('Max. deviation gapdh: %d percent',max(abs(gapdh_percent_exp))*100)
% gapdh_rev
gapdhr_diff_sim_exp = gapdhr_sim_vmax - gapdhr_exp_vmax;
gapdhr_percent_exp = gapdhr_diff_sim_exp./gapdhr_exp_vmax;
T1 = table(gapdhr_sim_vmax,gapdhr_exp_vmax,gapdhr_diff_sim_exp,gapdhr_percent_exp);
T1.Properties.VariableNames = {'gapdhr_sim_vmax','gapdhr_exp_vmax','gapdhr_diff_sim_exp','gapdhr_percent_exp'};
sprintf('Max. deviation gapdhr: %d percent',max(abs(gapdhr_percent_exp))*100)
% pgm
pgm_diff_sim_exp = pgm_sim_vmax - pgm_exp_vmax;
pgm_percent_exp = pgm_diff_sim_exp./pgm_exp_vmax;
T1 = table(pgm_sim_vmax,pgm_exp_vmax,pgm_diff_sim_exp,pgm_percent_exp);
T1.Properties.VariableNames = {'pgm_sim_vmax','pgm_exp_vmax','pgm_diff_sim_exp','pgm_percent_exp'};
sprintf('Max. deviation pgm: %d percent',max(abs(pgm_percent_exp))*100)
% eno
eno_diff_sim_exp = eno_sim_vmax - eno_exp_vmax;
eno_percent_exp = eno_diff_sim_exp./eno_exp_vmax;
T1 = table(eno_sim_vmax,eno_exp_vmax,eno_diff_sim_exp,eno_percent_exp);
T1.Properties.VariableNames = {'eno_sim_vmax','eno_exp_vmax','eno_diff_sim_exp','eno_percent_exp'};
sprintf('Max. deviation eno: %d percent',max(abs(eno_percent_exp))*100)
% pyk
pyk_diff_sim_exp = pyk_sim_vmax - pyk_exp_vmax;
pyk_percent_exp = pyk_diff_sim_exp./pyk_exp_vmax;
T1 = table(pyk_sim_vmax,pyk_exp_vmax,pyk_diff_sim_exp,pyk_percent_exp);
T1.Properties.VariableNames = {'pyk_sim_vmax','pyk_exp_vmax','pyk_diff_sim_exp','pyk_percent_exp'};
sprintf('Max. deviation pyk: %d percent',max(abs(pyk_percent_exp))*100)
% pdc
pdc_diff_sim_exp = pdc_sim_vmax - pdc_exp_vmax;
pdc_percent_exp = pdc_diff_sim_exp./pdc_exp_vmax;
T1 = table(pdc_sim_vmax,pdc_exp_vmax,pdc_diff_sim_exp,pdc_percent_exp);
T1.Properties.VariableNames = {'pdc_sim_vmax','pdc_exp_vmax','pdc_diff_sim_exp','pdc_percent_exp'};
sprintf('Max. deviation pdc: %d percent',max(abs(pdc_percent_exp))*100)


%% memoryDump
% % old style figure(200)

% figure(200)%(1)
% for i =1:numEnz
%     load(enzymeList{i});
%     subplot(4,3,i) % select subplot
%     tempName = ['output_',extractBefore(enzymeList{i},"_parEst")];
%     eval(['dataset = ',tempName,';']); % select the data
%     % correcting units for the 1/2 factor
% % % % %     val2multiply = (1/2);           % NEW CORRECTION METHOD (2020 - 10 - 29)
%     val2multiply = (1);           % NEW CORRECTION METHOD (2020 - 10 - 29)
%     if((i==2)||(i==3)||(i==4)||(i==5))
%         tempName2 = ['output_',enzymeName{i}];
%         eval(['dataset = ',tempName2,';']); % select the data
%         dataset.vm_uChange = dataset.vm_uChange * val2multiply;
%         dataset.vm_up_uChange = dataset.vm_up_uChange * val2multiply;
%         dataset.vm_down_uChange = dataset.vm_down_uChange * val2multiply;
%         dataset.Vmax_experimental = dataset.Vmax_experimental * val2multiply;
%         dataset.stDev_experimental = dataset.stDev_experimental * val2multiply;
%         % save back
%         eval([tempName2,' = dataset;']);
%     end
%     % plotting
%     if i == i
%         plot(dataset.pHarray, dataset.vm_uChange,...
%             'bo:','MarkerSize',5,'MarkerFaceColor','b')% plot errorbar computational estimation
%     else
%         errorbar(dataset.pHarray, dataset.vm_uChange, dataset.vm_up_uChange-dataset.vm_uChange,...
%             'bo:','MarkerSize',5,'MarkerFaceColor','b','CapSize',0)% plot errorbar computational estimation
%     end
%     hold on % hold on
%     % tempName = extractAfter(enzymeList{i},"output_");
%     errorbar(dataset.pHarray, dataset.Vmax_experimental, dataset.stDev_experimental,...
%         'rs','MarkerSize',5,'LineWidth',1.2,'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',0)% plot errorbar experimentel estimation
%     ax = gca;
%     ax.XLim = xlims_vals(i,:); % xlimits
%     ax.XTick = xticks_vals(i,:); % xticklabels
%     ax.YLim(1) = 0; % ylim(1) = 0
%     tempYLim = [max(dataset.Vmax_experimental+dataset.stDev_experimental) max(dataset.vm_up_uChange)];
%     if i ~= 1
%         ax.YLim(2) = max(tempYLim)*1.1; % ylim(1) = 0
%     end
%     %     tempName2 = [strrep(tempName,'_','_{'),'}'];
%     tempName2 = extractBefore(enzymeList{i},"_parEst");
%     text(ax.XLim(1)+0.1,ax.YLim(2)*0.9,tempName2)% tag
% %     if((i ~= 10)||(i ~= 11))
%     if(i <= 9)
%         set(gca,'xticklabel',{[]});
%     end
%     if i == 11
% % % % %         lgd = legend('computational','experimental');
% % % % %         lgd.Position = lgd.Position + [0.3 0 0 0];
%     end
%     % xlabel
% %     if((i == 10)||(i == 11))
%     if(i == 10)
%         xlabel('pH','FontSize',16)
%     end
%     % ylabel
% %     if((i == 1)||(i == 4)||(i == 7)||(i == 10))
%     if(i == 1)
%         ylabel('Vmax [\mumol mg_{P}^{-1} min^{-1}]','FontSize',16)
%     end
%     
% end
% set(gcf,'color','w');


%% Final figure
colSteelBlue = [70/255 130/255 180/255]; % pH independet
colLightBlue = [173/255 216/255 203/255]; 
colLightSkyBlue = [135/255	206/255	250/255]; % pH dependent

figure(500)
for i =1:numEnz
    % load and select the data
    load(enzymeList{i});
    subplot(4,3,i) % select subplot
    tempName = ['output_',extractBefore(enzymeList{i},"_parEst")];
    eval(['dataset = ',tempName,';']); % select the data
    
    % experimental or naïve parameter estimates
    eb1 = errorbar(dataset.pHarray, dataset.Vmax_experimental, dataset.stDev_experimental,...
        'o-','MarkerSize',3,...%'LineWidth',1.2,...
        'MarkerFaceColor',col4,...
        'CapSize',0);% plot errorbar experimentel estimation
    hold on % hold on
    eb1.Color = colSteelBlue;
    eb1.MarkerFaceColor = colSteelBlue;
    
    % computational parameter estimates
%     eb2 = errorbar(dataset.pHarray, dataset.vm_uChange, dataset.vm_up_uChange-dataset.vm_uChange,...
%         'ro:','MarkerSize',3,'MarkerFaceColor','r','CapSize',0);% plot errorbar computational estimation
    eb2 = plot(dataset.pHarray, dataset.vm_uChange,':.');
    eb2.MarkerSize = 10;
    eb2.Color = colLightSkyBlue;
    
    % readjust Y-axis
    ax = gca;
    ax.YLim(1) = 0;
    ax.YLim(2) = max([dataset.Vmax_experimental;dataset.vm_uChange])*1.1;
    
% %     % labels
% %     if((i == 1)||(i == 4)||(i == 7)||(i == 10))
% %         ylabel({'Vmax_{corrected}';'[\mumol mg_{P}^{-1} min^{-1}]'});
% %     end
% %     if((i == 10)||(i == 11))
% %         xlabel('pH value []');
% %     end
% %     title(enzymeName{i});
        % labels
    if i == 1
        ylabel_h = ylabel('Enzyme capacity (\mumol mg_{P}^{-1} min^{-1})','FontSize',14);
        ylabel_h.Position(2) = -1.175;
    end
    % xlabel
    if i == 11
        xlabel_h = xlabel('pH','FontSize',14);
        xlabel_h.Position(2) = -1;
    end
%     title(enzymeName{i});
%     text(ax.XLim(1)*1.1,ax.YLim(2)*0.9,dilutionsConsidered{i})
    xloc = (ax.XLim(2) - ax.XLim(1))*0.035 + ax.XLim(1);
    if((i == 6)||(i == 7))
        yloc = (ax.YLim(2) - ax.YLim(1))*0.85 + ax.YLim(1);
    else
        yloc = (ax.YLim(2) - ax.YLim(1))*0.87 + ax.YLim(1);
    end
%     text(ax.XLim(1)*1.1,ax.YLim(2)*0.9,enzymeName{i})
    tbox = text(xloc,yloc,enzymeName{i},'FontSize',10.5,'BackgroundColor',[0.9 0.9 0.9]);
    
    [~,idx] = max(dataset.Vmax_experimental);
    pH_max = dataset.pHarray(idx);
    yloc2 = (ax.YLim(2) - ax.YLim(1))*0.05 + ax.YLim(1);
%     tbox = text(pH_max-0.04,yloc2,'*','FontSize',14);
    
    % vertical line pointing at ideal pH
    line([pH_max pH_max],ax.YLim,'Color','black','LineStyle',':')
    
    hold off
end
set(gcf,'color','w');


%% 2021/07/13 get the data + polynomial fits

% recallung original data in sincle matrix
original_data = cell(numEnz,1);
original_data{1} = T2.hxk_sim_vmax;
original_data{2} = T2.pgi_sim_vmax;
original_data{3} = T2.pfk_sim_vmax;
original_data{4} = T2.ald_sim_vmax;
original_data{5} = T2.tpi_sim_vmax;
original_data{6} = T2.gapdh_sim_vmax;
original_data{7} = T2.gapdhr_sim_vmax;
original_data{8} = T2.pgm_sim_vmax;
original_data{9} = T2.eno_sim_vmax;
original_data{10} = T2.pyk_sim_vmax;
original_data{11} = T2.pdc_sim_vmax;

% data to consider
dtps2consider = cell(numEnz,1);
for i = 1:numEnz
    idxs_zeros = find(original_data{i} == 0);
    idxs_ones = 1:12; idxs_ones(idxs_zeros) = [];
    dtps2consider{i} = idxs_ones;
end

% creating polynomials
n_orders = 5;
polynomial_fits = cell(numEnz,1);
polynomial_fits{1} = polyfitn(T2.fullpHarray(dtps2consider{1}), T2.hxk_sim_vmax(dtps2consider{1}), n_orders); % hxk
polynomial_fits{2} = polyfitn(T2.fullpHarray(dtps2consider{2}), T2.pgi_sim_vmax(dtps2consider{2}), n_orders); % pgi
polynomial_fits{3} = polyfitn(T2.fullpHarray(dtps2consider{3}), T2.pfk_sim_vmax(dtps2consider{3}), n_orders); % pfk
polynomial_fits{4} = polyfitn(T2.fullpHarray(dtps2consider{4}), T2.ald_sim_vmax(dtps2consider{4}), n_orders); % ald
polynomial_fits{5} = polyfitn(T2.fullpHarray(dtps2consider{5}), T2.tpi_sim_vmax(dtps2consider{5}), n_orders); % tpi
polynomial_fits{6} = polyfitn(T2.fullpHarray(dtps2consider{6}), T2.gapdh_sim_vmax(dtps2consider{6}), n_orders); % gapdh
polynomial_fits{7} = polyfitn(T2.fullpHarray(dtps2consider{7}), T2.gapdhr_sim_vmax(dtps2consider{7}), n_orders); % gapdhr
polynomial_fits{8} = polyfitn(T2.fullpHarray(dtps2consider{8}), T2.pgm_sim_vmax(dtps2consider{8}), n_orders); % pgm
polynomial_fits{9} = polyfitn(T2.fullpHarray(dtps2consider{9}), T2.eno_sim_vmax(dtps2consider{9}), n_orders); % eno
polynomial_fits{10} = polyfitn(T2.fullpHarray(dtps2consider{10}), T2.pyk_sim_vmax(dtps2consider{10}), n_orders); % pyk
polynomial_fits{11} = polyfitn(T2.fullpHarray(dtps2consider{11}), T2.pdc_sim_vmax(dtps2consider{11}), n_orders); % pdc
% polynomial_fits{1} = polyfitn(T2.fullpHarray,T2.hxk_sim_vmax,n_orders); % hxk
% polynomial_fits{2} = polyfitn(T2.fullpHarray,T2.pgi_sim_vmax,n_orders); % pgi
% polynomial_fits{3} = polyfitn(T2.fullpHarray,T2.pfk_sim_vmax,n_orders); % pfk
% polynomial_fits{4} = polyfitn(T2.fullpHarray,T2.ald_sim_vmax,n_orders); % ald
% polynomial_fits{5} = polyfitn(T2.fullpHarray,T2.tpi_sim_vmax,n_orders); % tpi
% polynomial_fits{6} = polyfitn(T2.fullpHarray,T2.gapdh_sim_vmax,n_orders); % gapdh
% polynomial_fits{7} = polyfitn(T2.fullpHarray,T2.gapdhr_sim_vmax,n_orders); % gapdhr
% polynomial_fits{8} = polyfitn(T2.fullpHarray,T2.pgm_sim_vmax,n_orders); % pgm
% polynomial_fits{9} = polyfitn(T2.fullpHarray,T2.eno_sim_vmax,n_orders); % eno
% polynomial_fits{10} = polyfitn(T2.fullpHarray,T2.pyk_sim_vmax,n_orders); % pyk
% polynomial_fits{11} = polyfitn(T2.fullpHarray,T2.pdc_sim_vmax,n_orders); % pdc

%% plotting
polynomial_texts = cell(numEnz,1);
it_vals = 6.19:0.01:7.90;
figure(601)
for i = 1:numEnz
    subplot(4,3,i)
    % plot polytif
    plot(it_vals, polyval(polynomial_fits{i}.Coefficients,it_vals), 'k-', 'LineWidth', 2)
    hold on
    % plot experimental data
%     plot(T2.fullpHarray, original_data{i}, 'o')
    idxs_zeros = find(original_data{i} == 0);
    idxs_ones = 1:12; idxs_ones(idxs_zeros) = [];
    plot(T2.fullpHarray(idxs_ones), original_data{i}(idxs_ones), 'ko', 'MarkerFaceColor', 'red', 'MarkerSize', 5)
    % title = enzymeName + polyfit equation + R2
    temp1 = char(vpa(poly2sym(polynomial_fits{i}.Coefficients),5));
    temp2 = num2str(polynomial_fits{i}.R2);
%     title({enzymeName{i}; ['y = ', temp1]; ['R2 = ', temp2]})
    title([enzymeName{i}, ', R2 = ', temp2])
    polynomial_texts{i} = temp1;
end
suptitle('Polynomial fit for each parameter estimates. NOT normalized')
set(gcf,'color','w');


%% saving to send
% save figure
if setup.saveOutput == 1
    savefig(601,'polynomialFits_Vmax_vs_pH')
end

%% save table data + polynomial
% create last table
T3 = table(enzymeName,polynomial_texts);
T3.Properties.VariableNames = {'enzyme_names','polynomial_fits'};
%%
% 
filename_T2 = 'Vmax_vs_pH.xlsx';
if setup.saveOutput == 1
    writetable(T2,filename_T2,'Sheet',1,'Range','A1');
end
%
filename_T3 = 'polynomialFits.xlsx';
if setup.saveOutput == 1
    writetable(T3,filename_T3,'Sheet',1,'Range','A1');
end

%% Polynomial 1
% creating polynomials
n_orders = 1;
polynomial_fits = cell(numEnz,1);
polynomial_fits{1} = polyfitn(T2.fullpHarray(dtps2consider{1}), T2.hxk_sim_vmax(dtps2consider{1}), n_orders); % hxk
polynomial_fits{2} = polyfitn(T2.fullpHarray(dtps2consider{2}), T2.pgi_sim_vmax(dtps2consider{2}), n_orders); % pgi
polynomial_fits{3} = polyfitn(T2.fullpHarray(dtps2consider{3}), T2.pfk_sim_vmax(dtps2consider{3}), n_orders); % pfk
polynomial_fits{4} = polyfitn(T2.fullpHarray(dtps2consider{4}), T2.ald_sim_vmax(dtps2consider{4}), n_orders); % ald
polynomial_fits{5} = polyfitn(T2.fullpHarray(dtps2consider{5}), T2.tpi_sim_vmax(dtps2consider{5}), n_orders); % tpi
polynomial_fits{6} = polyfitn(T2.fullpHarray(dtps2consider{6}), T2.gapdh_sim_vmax(dtps2consider{6}), n_orders); % gapdh
polynomial_fits{7} = polyfitn(T2.fullpHarray(dtps2consider{7}), T2.gapdhr_sim_vmax(dtps2consider{7}), n_orders); % gapdhr
polynomial_fits{8} = polyfitn(T2.fullpHarray(dtps2consider{8}), T2.pgm_sim_vmax(dtps2consider{8}), n_orders); % pgm
polynomial_fits{9} = polyfitn(T2.fullpHarray(dtps2consider{9}), T2.eno_sim_vmax(dtps2consider{9}), n_orders); % eno
polynomial_fits{10} = polyfitn(T2.fullpHarray(dtps2consider{10}), T2.pyk_sim_vmax(dtps2consider{10}), n_orders); % pyk
polynomial_fits{11} = polyfitn(T2.fullpHarray(dtps2consider{11}), T2.pdc_sim_vmax(dtps2consider{11}), n_orders); % pdc

% Plotting
polynomial_texts = cell(numEnz,1);
it_vals = 6.19:0.01:7.90;
figure(701)
for i = 1:numEnz
    subplot(4,3,i)
    % plot polytif
    plot(it_vals, polyval(polynomial_fits{i}.Coefficients,it_vals), 'k-', 'LineWidth', 2)
    hold on
    % plot experimental data
%     plot(T2.fullpHarray, original_data{i}, 'o')
    idxs_zeros = find(original_data{i} == 0);
    idxs_ones = 1:12; idxs_ones(idxs_zeros) = [];
    plot(T2.fullpHarray(idxs_ones), original_data{i}(idxs_ones), 'ko', 'MarkerFaceColor', 'red', 'MarkerSize', 5)
    % title = enzymeName + polyfit equation + R2
    temp1 = char(vpa(poly2sym(polynomial_fits{i}.Coefficients),5));
    temp2 = num2str(polynomial_fits{i}.R2);
%     title({enzymeName{i}; ['y = ', temp1]; ['R2 = ', temp2]})
    if((i == 4)||(i == 5)||(i == 7)||(i == 11))
        title(['\color{blue}', enzymeName{i}, ', R2 = ', temp2])
    else
        title([enzymeName{i}, ', R2 = ', temp2])
    end
    polynomial_texts{i} = temp1;
end
suptitle('Polynomial (order 1) fit for each parameter estimates. NOT normalized')
set(gcf,'color','w');

% table
T_01 = table(enzymeName,polynomial_texts);
T_01.Properties.VariableNames = {'enzyme_names','polynomial_fits'};

% saving output
filename_T_01 = 'polynomialFits_order1.xlsx';
if setup.saveOutput == 1
    writetable(T_01,filename_T_01,'Sheet',1,'Range','A1');
    savefig(701,'polynomialFits_Vmax_vs_pH_order1')
end

%% Polynomial 2
% creating polynomials
n_orders = 2;
polynomial_fits = cell(numEnz,1);
polynomial_fits{1} = polyfitn(T2.fullpHarray(dtps2consider{1}), T2.hxk_sim_vmax(dtps2consider{1}), n_orders); % hxk
polynomial_fits{2} = polyfitn(T2.fullpHarray(dtps2consider{2}), T2.pgi_sim_vmax(dtps2consider{2}), n_orders); % pgi
polynomial_fits{3} = polyfitn(T2.fullpHarray(dtps2consider{3}), T2.pfk_sim_vmax(dtps2consider{3}), n_orders); % pfk
polynomial_fits{4} = polyfitn(T2.fullpHarray(dtps2consider{4}), T2.ald_sim_vmax(dtps2consider{4}), n_orders); % ald
polynomial_fits{5} = polyfitn(T2.fullpHarray(dtps2consider{5}), T2.tpi_sim_vmax(dtps2consider{5}), n_orders); % tpi
polynomial_fits{6} = polyfitn(T2.fullpHarray(dtps2consider{6}), T2.gapdh_sim_vmax(dtps2consider{6}), n_orders); % gapdh
polynomial_fits{7} = polyfitn(T2.fullpHarray(dtps2consider{7}), T2.gapdhr_sim_vmax(dtps2consider{7}), n_orders); % gapdhr
polynomial_fits{8} = polyfitn(T2.fullpHarray(dtps2consider{8}), T2.pgm_sim_vmax(dtps2consider{8}), n_orders); % pgm
polynomial_fits{9} = polyfitn(T2.fullpHarray(dtps2consider{9}), T2.eno_sim_vmax(dtps2consider{9}), n_orders); % eno
polynomial_fits{10} = polyfitn(T2.fullpHarray(dtps2consider{10}), T2.pyk_sim_vmax(dtps2consider{10}), n_orders); % pyk
polynomial_fits{11} = polyfitn(T2.fullpHarray(dtps2consider{11}), T2.pdc_sim_vmax(dtps2consider{11}), n_orders); % pdc

% Plotting
polynomial_texts = cell(numEnz,1);
it_vals = 6.19:0.01:7.90;
figure(702)
for i = 1:numEnz
    subplot(4,3,i)
    % plot polytif
    plot(it_vals, polyval(polynomial_fits{i}.Coefficients,it_vals), 'k-', 'LineWidth', 2)
    hold on
    % plot experimental data
%     plot(T2.fullpHarray, original_data{i}, 'o')
    idxs_zeros = find(original_data{i} == 0);
    idxs_ones = 1:12; idxs_ones(idxs_zeros) = [];
    plot(T2.fullpHarray(idxs_ones), original_data{i}(idxs_ones), 'ko', 'MarkerFaceColor', 'red', 'MarkerSize', 5)
    % title = enzymeName + polyfit equation + R2
    temp1 = char(vpa(poly2sym(polynomial_fits{i}.Coefficients),5));
    temp2 = num2str(polynomial_fits{i}.R2);
%     title({enzymeName{i}; ['y = ', temp1]; ['R2 = ', temp2]})
    if((i == 4)||(i == 5)||(i == 7)||(i == 11))
        title(['\color{blue}', enzymeName{i}, ', R2 = ', temp2])
    else
        title([enzymeName{i}, ', R2 = ', temp2])
    end
    polynomial_texts{i} = temp1;
end
suptitle('Polynomial (order 2) fit for each parameter estimates. NOT normalized')
set(gcf,'color','w');

% table
T_02 = table(enzymeName,polynomial_texts);
T_02.Properties.VariableNames = {'enzyme_names','polynomial_fits'};

%% saving output
filename_T_02 = 'polynomialFits_order2.xlsx';
if setup.saveOutput == 1
    writetable(T_02,filename_T_02,'Sheet',1,'Range','A1');
    savefig(702,'polynomialFits_Vmax_vs_pH_order2')
end

%% Polynomial 3
% creating polynomials
n_orders = 3;
polynomial_fits = cell(numEnz,1);
polynomial_fits{1} = polyfitn(T2.fullpHarray(dtps2consider{1}), T2.hxk_sim_vmax(dtps2consider{1}), n_orders); % hxk
polynomial_fits{2} = polyfitn(T2.fullpHarray(dtps2consider{2}), T2.pgi_sim_vmax(dtps2consider{2}), n_orders); % pgi
polynomial_fits{3} = polyfitn(T2.fullpHarray(dtps2consider{3}), T2.pfk_sim_vmax(dtps2consider{3}), n_orders); % pfk
polynomial_fits{4} = polyfitn(T2.fullpHarray(dtps2consider{4}), T2.ald_sim_vmax(dtps2consider{4}), n_orders); % ald
polynomial_fits{5} = polyfitn(T2.fullpHarray(dtps2consider{5}), T2.tpi_sim_vmax(dtps2consider{5}), n_orders); % tpi
polynomial_fits{6} = polyfitn(T2.fullpHarray(dtps2consider{6}), T2.gapdh_sim_vmax(dtps2consider{6}), n_orders); % gapdh
polynomial_fits{7} = polyfitn(T2.fullpHarray(dtps2consider{7}), T2.gapdhr_sim_vmax(dtps2consider{7}), n_orders); % gapdhr
polynomial_fits{8} = polyfitn(T2.fullpHarray(dtps2consider{8}), T2.pgm_sim_vmax(dtps2consider{8}), n_orders); % pgm
polynomial_fits{9} = polyfitn(T2.fullpHarray(dtps2consider{9}), T2.eno_sim_vmax(dtps2consider{9}), n_orders); % eno
polynomial_fits{10} = polyfitn(T2.fullpHarray(dtps2consider{10}), T2.pyk_sim_vmax(dtps2consider{10}), n_orders); % pyk
polynomial_fits{11} = polyfitn(T2.fullpHarray(dtps2consider{11}), T2.pdc_sim_vmax(dtps2consider{11}), n_orders); % pdc

% Plotting
polynomial_texts = cell(numEnz,1);
it_vals = 6.19:0.01:7.90;
figure(703)
for i = 1:numEnz
    subplot(4,3,i)
    % plot polytif
    plot(it_vals, polyval(polynomial_fits{i}.Coefficients,it_vals), 'k-', 'LineWidth', 2)
    hold on
    % plot experimental data
%     plot(T2.fullpHarray, original_data{i}, 'o')
    idxs_zeros = find(original_data{i} == 0);
    idxs_ones = 1:12; idxs_ones(idxs_zeros) = [];
    plot(T2.fullpHarray(idxs_ones), original_data{i}(idxs_ones), 'ko', 'MarkerFaceColor', 'red', 'MarkerSize', 5)
    % title = enzymeName + polyfit equation + R2
    temp1 = char(vpa(poly2sym(polynomial_fits{i}.Coefficients),5));
    temp2 = num2str(polynomial_fits{i}.R2);
%     title({enzymeName{i}; ['y = ', temp1]; ['R2 = ', temp2]})
    if((i == 4)||(i == 5)||(i == 7)||(i == 11))
        title(['\color{blue}', enzymeName{i}, ', R2 = ', temp2])
    else
        title([enzymeName{i}, ', R2 = ', temp2])
    end
    polynomial_texts{i} = temp1;
end
suptitle('Polynomial (order 3) fit for each parameter estimates. NOT normalized')
set(gcf,'color','w');

% table
T_03 = table(enzymeName,polynomial_texts);
T_03.Properties.VariableNames = {'enzyme_names','polynomial_fits'};

%% saving output
filename_T_03 = 'polynomialFits_order3.xlsx';
if setup.saveOutput == 1
    writetable(T_03,filename_T_03,'Sheet',1,'Range','A1');
    savefig(703,'polynomialFits_Vmax_vs_pH_order3')
end

%% Polynomial 4
% creating polynomials
n_orders = 4;
polynomial_fits = cell(numEnz,1);
polynomial_fits{1} = polyfitn(T2.fullpHarray(dtps2consider{1}), T2.hxk_sim_vmax(dtps2consider{1}), n_orders); % hxk
polynomial_fits{2} = polyfitn(T2.fullpHarray(dtps2consider{2}), T2.pgi_sim_vmax(dtps2consider{2}), n_orders); % pgi
polynomial_fits{3} = polyfitn(T2.fullpHarray(dtps2consider{3}), T2.pfk_sim_vmax(dtps2consider{3}), n_orders); % pfk
polynomial_fits{4} = polyfitn(T2.fullpHarray(dtps2consider{4}), T2.ald_sim_vmax(dtps2consider{4}), n_orders); % ald
polynomial_fits{5} = polyfitn(T2.fullpHarray(dtps2consider{5}), T2.tpi_sim_vmax(dtps2consider{5}), n_orders); % tpi
polynomial_fits{6} = polyfitn(T2.fullpHarray(dtps2consider{6}), T2.gapdh_sim_vmax(dtps2consider{6}), n_orders); % gapdh
polynomial_fits{7} = polyfitn(T2.fullpHarray(dtps2consider{7}), T2.gapdhr_sim_vmax(dtps2consider{7}), n_orders); % gapdhr
polynomial_fits{8} = polyfitn(T2.fullpHarray(dtps2consider{8}), T2.pgm_sim_vmax(dtps2consider{8}), n_orders); % pgm
polynomial_fits{9} = polyfitn(T2.fullpHarray(dtps2consider{9}), T2.eno_sim_vmax(dtps2consider{9}), n_orders); % eno
polynomial_fits{10} = polyfitn(T2.fullpHarray(dtps2consider{10}), T2.pyk_sim_vmax(dtps2consider{10}), n_orders); % pyk
polynomial_fits{11} = polyfitn(T2.fullpHarray(dtps2consider{11}), T2.pdc_sim_vmax(dtps2consider{11}), n_orders); % pdc

% Plotting
polynomial_texts = cell(numEnz,1);
it_vals = 6.19:0.01:7.90;
figure(704)
for i = 1:numEnz
    subplot(4,3,i)
    % plot polytif
    plot(it_vals, polyval(polynomial_fits{i}.Coefficients,it_vals), 'k-', 'LineWidth', 2)
    hold on
    % plot experimental data
%     plot(T2.fullpHarray, original_data{i}, 'o')
    idxs_zeros = find(original_data{i} == 0);
    idxs_ones = 1:12; idxs_ones(idxs_zeros) = [];
    plot(T2.fullpHarray(idxs_ones), original_data{i}(idxs_ones), 'ko', 'MarkerFaceColor', 'red', 'MarkerSize', 5)
    % title = enzymeName + polyfit equation + R2
    temp1 = char(vpa(poly2sym(polynomial_fits{i}.Coefficients),5));
    temp2 = num2str(polynomial_fits{i}.R2);
%     title({enzymeName{i}; ['y = ', temp1]; ['R2 = ', temp2]})
    if((i == 4)||(i == 5)||(i == 7)||(i == 11))
        title(['\color{blue}', enzymeName{i}, ', R2 = ', temp2])
    else
        title([enzymeName{i}, ', R2 = ', temp2])
    end
    polynomial_texts{i} = temp1;
end
suptitle('Polynomial (order 4) fit for each parameter estimates. NOT normalized')
set(gcf,'color','w');

% table
T_04 = table(enzymeName,polynomial_texts);
T_04.Properties.VariableNames = {'enzyme_names','polynomial_fits'};

%% saving output
filename_T_04 = 'polynomialFits_order4.xlsx';
if setup.saveOutput == 1
    writetable(T_04,filename_T_04,'Sheet',1,'Range','A1');
    savefig(704,'polynomialFits_Vmax_vs_pH_order4')
end

