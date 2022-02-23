% % vmExp
% Generation of the manuscript figure that plots the experimental and
% computational estimates for the maximum reaction rates.
% clear, close all
% 
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


%% Final figure (appendix)
col1 = [189,215,231]/255; col2 = [107,174,214]/255; col3 = [49,130,189]/255; col4 = [8,81,156]/255;
colArray = [col1; col2; col3; col4];
fh1 = figure(302);
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

%% Change 

% titles text, position and background
% 
fh1.Children(2).Children(1).String = 'PDC';
fh1.Children(3).Children(1).String = 'PYK';
fh1.Children(4).Children(1).String = 'ENO';
fh1.Children(5).Children(1).String = 'PGM';
fh1.Children(6).Children(1).String = 'GAPDHR';
fh1.Children(7).Children(1).String = 'GAPDH';
fh1.Children(8).Children(1).String = 'TPI';
fh1.Children(9).Children(1).String = 'ALD';
fh1.Children(10).Children(1).String = 'PFK';
fh1.Children(11).Children(1).String = 'PGI';
fh1.Children(12).Children(1).String = 'HXK';
% specify y-lims (needed for later)
fh1.Children(2).YLim = [0 3.5]; %'PDC';
fh1.Children(3).YLim = [0 8]; %'PYK';
fh1.Children(4).YLim = [0 2]; %'ENO';
fh1.Children(5).YLim = [0 30]; %'PGM';
fh1.Children(6).YLim = [0 5]; %'GAPDHR';
fh1.Children(7).YLim = [0 3]; %'GAPDH';
fh1.Children(8).YLim = [0 60]; %'TPI';
fh1.Children(9).YLim = [0 3]; %'ALD';
fh1.Children(10).YLim = [0 0.6]; %'PFK';
fh1.Children(11).YLim = [0 3]; %'PGI';
fh1.Children(12).YLim = [0 1]; %'HXK';
% 
for i = 2:12
    fh1.Children(i).FontSize  = 11;
    temp = fh1.Children(i).YLim(1) + 0.9*(fh1.Children(i).YLim(2) - fh1.Children(i).YLim(1));
    fh1.Children(i).Children(1).Position = [6.1 temp 0];
    fh1.Children(i).Children(1).BackgroundColor = 'none';
end

% legenda
fh1.Children(1).FontSize = 12;
fh1.Children(1).Position(1) = 0.75;
fh1.Children(1).Position(2) = 0.14;

%%
fh1.Position = [100 100 800 900]; % = [0.025    0.075    0.3    0.8225];
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%     % Saving: save them back to .png in location
%     savefig(1,'1appendix_ald_vm_dil');   
%     saveLoc = 'D:\OneDrive - TU Eindhoven\Documents\ch3_pHkinetics\results\manuscriptFigures\';
%     saveName_train_pdf = [saveLoc,'1appendix_ald_vm_dil.pdf'];
%     saveName_train_png = [saveLoc,'1appendix_ald_vm_dil.png'];
%     saveas(1,saveName_train_pdf);
%     saveas(1,saveName_train_png);


%%
% save
if setup.saveOutput == 1
    savefig(302,'1appendix_VmExp_differentDFs'); 
    % specs printing (method 3)
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters 1appendix_VmExp_differentDFs
    print -dpng -painters 1appendix_VmExp_differentDFs
end

