% % vmExp
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
enzymeName2 = {'hxk';... %1 hxk
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
dilutionsConsidered = {'DF 1 2 4 8';... %1 hxk
    'DF 1 2';... %2 pgi
    'DF 8 16 32';... %3 pfk
    'DF 1 2 4 8';... %4 ald
    'DF 8 16';... %5 tpi
    'DF 1 2 4 8 16';... %6 gapdh_fwd
    'DF 1 2 4 8';... %7 gapdh_rev
    'DF 2 4 8 16';... %8 pgm
    'DF 1 2 4 8';... %9 eno_kmfixed
    'DF 2 4 8 16 32';... %10 pyk
    'DF 1 2 4'}; %11 pdc

xlims_vals = ones(numEnz,2);
xticks_vals = ones(numEnz,5);
for i = 1:numEnz
    xlims_vals(i,:) = [6 8];
    xticks_vals(i,:) = [6 6.5 7 7.5 8];
end

% added colors
col1 = [189,215,231]/255; col2 = [107,174,214]/255; col3 = [49,130,189]/255; col4 = [8,81,156]/255;


%%
figure(300)
for i =1:numEnz
    load(enzymeList{i});
    subplot(4,3,i) % select subplot
    tempName = ['output_',extractBefore(enzymeList{i},"_parEst")];
    eval(['dataset = ',tempName,';']); % select the data
    % correcting units for the 1/2 factor
    val2multiply = (1);           % NEW CORRECTION METHOD (2020 - 10 - 29)
    if((i==2)||(i==3)||(i==4)||(i==5))
        tempName2 = ['output_',enzymeName{i}];
        eval(['dataset = ',tempName2,';']); % select the data
        dataset.vm_uChange = dataset.vm_uChange * val2multiply;
        dataset.vm_up_uChange = dataset.vm_up_uChange * val2multiply;
        dataset.vm_down_uChange = dataset.vm_down_uChange * val2multiply;
        dataset.Vmax_experimental = dataset.Vmax_experimental * val2multiply;
        dataset.stDev_experimental = dataset.stDev_experimental * val2multiply;
        % save back
        eval([tempName2,' = dataset;']);
    end
    % plotting
    errorbar(dataset.pHarray, dataset.Vmax_experimental, dataset.stDev_experimental,...
        'o-','MarkerSize',3,...%'LineWidth',1.2,...
        'MarkerFaceColor',col4,...
        'Color',col4,'CapSize',0)% plot errorbar experimentel estimation
    hold on % hold on
    % readjust Y-axis
    ax = gca;
    ax.YLim(1) = 0;
    ax.YLim(2) = max(dataset.Vmax_experimental)*1.5;
    
    % labels
    if((i == 1)||(i == 4)||(i == 7)||(i == 10))
        ylabel({'Vmax_{corrected}';'[\mumol mg_{P}^{-1} min^{-1}]'});
    end
    if((i == 10)||(i == 11))
        xlabel('pH value []');
    end
    title(enzymeName{i});
    text(ax.XLim(1)*1.1,ax.YLim(2)*0.9,dilutionsConsidered{i})
end
set(gcf,'color','w');


%%
figure(301)
for i =1:numEnz
    load(enzymeList{i});
    subplot(4,3,i) % select subplot
    tempName = ['output_',extractBefore(enzymeList{i},"_parEst")];
    eval(['dataset = ',tempName,';']); % select the data
    % correcting units for the 1/2 factor
    val2multiply = (1);           % NEW CORRECTION METHOD (2020 - 10 - 29)
    if((i==2)||(i==3)||(i==4)||(i==5))
        tempName2 = ['output_',enzymeName{i}];
        eval(['dataset = ',tempName2,';']); % select the data
        dataset.vm_uChange = dataset.vm_uChange * val2multiply;
        dataset.vm_up_uChange = dataset.vm_up_uChange * val2multiply;
        dataset.vm_down_uChange = dataset.vm_down_uChange * val2multiply;
        dataset.Vmax_experimental = dataset.Vmax_experimental * val2multiply;
        dataset.stDev_experimental = dataset.stDev_experimental * val2multiply;
        % save back
        eval([tempName2,' = dataset;']);
    end
    % plotting
    fill([dataset.pHarray' fliplr(dataset.pHarray')],...
        [dataset.Vmax_experimental'+dataset.stDev_experimental' fliplr(dataset.Vmax_experimental'-dataset.stDev_experimental')],...
        col2,...
        'linestyle','none');
    hold on
    plot(dataset.pHarray, dataset.Vmax_experimental,...
        'o-','MarkerSize',3,...
        'MarkerFaceColor',col4,...
        'Color',col4)
    
    
%     errorbar(dataset.pHarray, dataset.Vmax_experimental, dataset.stDev_experimental,...
%         'o-','MarkerSize',3,...%'LineWidth',1.2,...
%         'MarkerFaceColor',col4,...
%         'Color',col4,'CapSize',0)% plot errorbar experimentel estimation
%     hold on % hold on
    % readjust Y-axis
    ax = gca;
    ax.YLim(1) = 0;
    ax.YLim(2) = max(dataset.Vmax_experimental)*1.5;
    
    % labels
    if((i == 1)||(i == 4)||(i == 7)||(i == 10))
        ylabel({'Vmax_{corrected}';'[\mumol mg_{P}^{-1} min^{-1}]'});
    end
    if((i == 10)||(i == 11))
        xlabel('pH value []');
    end
    title(enzymeName{i});
    text(ax.XLim(1)*1.1,ax.YLim(2)*0.9,dilutionsConsidered{i})
end
set(gcf,'color','w');


%%
col1 = [189,215,231]/255; col2 = [107,174,214]/255; col3 = [49,130,189]/255; col4 = [8,81,156]/255;
colArray = [col1; col2; col3; col4];
figure(302)
for i = 1:length(enzymeName)
    % load
    loadName = [enzymeName2{i},'_initial_variables.mat'];
    load(loadName);
    [ij,ik] = size(DF);
    
    % plot
    subplot(4,3,i)
    for k = 1:ik
        plot(pH(:,k), Vmax_mw_opt_corr(:,k), 'o-','color',colArray(k,:),'markersize',4,'markerfacecolor',colArray(k,:))
        hold on
    end
    
    % labels
    if((i == 1)||(i == 4)||(i == 7)||(i == 10))
        ylabel({'Vmax_{corrected}';'[\mumol mg_{P}^{-1} min^{-1}]'});
    end
    if((i == 10)||(i == 11))
        xlabel('pH value []');
    end
    title(enzymeName{i});
end
set(302,'color','white')


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
filename = '20201102_parameterValues.xlsx';
% writetable(T,filename,'Sheet',1,'Range','A1');


%% Final figure (main)
colSteelBlue = [70/255 130/255 180/255]; % pH independent
colLightBlue = [173/255 216/255 203/255];
colLightSkyBlue = [135/255	206/255	250/255];

figure(400)
for i =1:numEnz
    load(enzymeList{i});
    subplot(4,3,i) % select subplot
    tempName = ['output_',extractBefore(enzymeList{i},"_parEst")];
    eval(['dataset = ',tempName,';']); % select the data
    % correcting units for the 1/2 factor
    val2multiply = (1);           % NEW CORRECTION METHOD (2020 - 10 - 29)
    if((i==2)||(i==3)||(i==4)||(i==5))
        tempName2 = ['output_',enzymeName{i}];
        eval(['dataset = ',tempName2,';']); % select the data
        dataset.vm_uChange = dataset.vm_uChange * val2multiply;
        dataset.vm_up_uChange = dataset.vm_up_uChange * val2multiply;
        dataset.vm_down_uChange = dataset.vm_down_uChange * val2multiply;
        dataset.Vmax_experimental = dataset.Vmax_experimental * val2multiply;
        dataset.stDev_experimental = dataset.stDev_experimental * val2multiply;
        % save back
        eval([tempName2,' = dataset;']);
    end
    % plotting
    errorbar(dataset.pHarray, dataset.Vmax_experimental, dataset.stDev_experimental,...
        'o-','MarkerSize',3,...%'LineWidth',1.2,...
        'MarkerFaceColor',colSteelBlue,...
        'Color',colSteelBlue,'CapSize',1)% plot errorbar experimentel estimation
    hold on % hold on
    % readjust Y-axis
    ax = gca;
    ax.YLim(1) = 0;
    ax.YLim(2) = max(dataset.Vmax_experimental)*1.5;
    
    % labels
    if i == 1
        ylabel_h = ylabel('Enzyme capacity (\mumol mg_{P}^{-1} min^{-1})','FontSize',14);
        ylabel_h.Position(2) = -1.5;
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


%% Final figure (appendix)

col1 = [189,215,231]/255; col2 = [107,174,214]/255; col3 = [49,130,189]/255; col4 = [8,81,156]/255;
colArray = [col1; col2; col3; col4];
figure(302)
for i = 1:length(enzymeName)
    % load
    loadName = [enzymeName2{i},'_initial_variables.mat'];
    load(loadName);
    [ij,ik] = size(DF);
    
    % plot
    subplot(4,3,i)
    for k = 1:ik
        plot(pH(:,k), Vmax_mw_opt_corr(:,k), 'o-','color',colArray(k,:),'markersize',4,'markerfacecolor',colArray(k,:))
        hold on
    end
    ax = gca;
    
    % labels
    if i == 1
        ylabel_h = ylabel('Enzyme capacity (\mumol mg_{P}^{-1} min^{-1})','FontSize',14);
        ylabel_h.Position(2) = -1.5;
    end
    % xlabel
    if i == 11
        xlabel_h = xlabel('pH','FontSize',14);
        xlabel_h.Position(2) = -1;
    end
    xloc = (ax.XLim(2) - ax.XLim(1))*0.035 + ax.XLim(1);
    if((i == 6)||(i == 7))
        yloc = (ax.YLim(2) - ax.YLim(1))*0.85 + ax.YLim(1);
    else
        yloc = (ax.YLim(2) - ax.YLim(1))*0.87 + ax.YLim(1);
    end
%     text(ax.XLim(1)*1.1,ax.YLim(2)*0.9,enzymeName{i})
    tbox = text(xloc,yloc,enzymeName{i},'FontSize',10.5,'BackgroundColor',[0.9 0.9 0.9]);
%     title(enzymeName{i});

    % legenda
    if i == length(enzymeName)
        leg_h = legend('DF x8','DF x4','DF x2','DF x1');
        leg_h.Position(1) = leg_h.Position(1) + 0.25;
        leg_h.Position(2) = leg_h.Position(2) - 0.015;
    end

    hold off
end
set(302,'color','white')


