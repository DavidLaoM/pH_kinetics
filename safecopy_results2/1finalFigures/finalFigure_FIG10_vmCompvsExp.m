% % fig10_vmCompvsExp
% In this file, the vm estimated via the computational method are compared
% to Laura's value

% step 1. load
% step 2. plot
% (if needed) step 3. analyze

%% step 1. load
clear
load('output_hxk.mat');
load('output_pgi.mat');
load('output_ald.mat');
% load('output_gapdhr.mat');
load('output_eno_kmfixed.mat');
load('output_pyk.mat');
%     % to re run and correct
%     output_pyk.xres_selected(13:14) = [output_pyk.xres_selected(15), output_pyk.xres_selected(15)];
%     output_pyk.vm(7:8) = [output_pyk.vm(9); output_pyk.vm(9)];
%     output_pyk.vm_up(7:8) = [output_pyk.vm_up(9); output_pyk.vm_up(9)];
%     output_pyk.vm_down(7:8) = [output_pyk.vm_down(9); output_pyk.vm_down(9)];
%     output_pyk.vm_uChange(7:8) = [output_pyk.vm_uChange(9); output_pyk.vm_uChange(9)];
%     output_pyk.vm_up_uChange(7:8) = [output_pyk.vm_up_uChange(9); output_pyk.vm_up_uChange(9)];
%     output_pyk.vm_down_uChange(7:8) = [output_pyk.vm_down_uChange(9); output_pyk.vm_down_uChange(9)];
load('output_pfk.mat');
load('output_tpi.mat');
load('output_pgm.mat');
load('output_pdc.mat');
load('experimentalRates.mat');

%% step 2.1. plot not-normalized
plotORerrorbar = 2;
figure(31)
%   #1     HXK
%   #2     PGI
%   #3     PFK
%   #4     ALD
%   #5     TPI
%   #6     GAPDH
%   #7     GAPDHr
%   #8     PGM
%   #9     ENO
%   #10    PYK
%   #11    PDC
if plotORerrorbar == 1 % make plots
    subplot(4,3,1) % vmHXK
    plot(output_hxk.pHarray,output_hxk.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_hxk.pHarray,output_hxk.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_hxk.pHarray,output_hxk.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.hxk(:,1), ratesStruct.hxk(:,2), ratesStruct.hxk(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{HXK}, NADPH')

    subplot(4,3,2) % vmPGI
    plot(output_pgi.pHarray,output_pgi.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pgi.pHarray,output_pgi.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pgi.pHarray,output_pgi.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pgi(:,1), ratesStruct.pgi(:,2), ratesStruct.pgi(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PGI}, NADH')

    subplot(4,3,3) % vmPFK
    plot(output_pfk.pHarray,output_pfk.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pfk.pHarray,output_pfk.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pfk.pHarray,output_pfk.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pfk(2:end,1), ratesStruct.pfk(2:end,2), ratesStruct.pfk(2:end,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PFK}, NADH')

    subplot(4,3,4) % vmALD
    plot(output_ald.pHarray,output_ald.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_ald.pHarray,output_ald.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_ald.pHarray,output_ald.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.ald(:,1), ratesStruct.ald(:,2), ratesStruct.ald(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{ALD}, NADH')

    subplot(4,3,5) %   #5     TPI
    plot(output_tpi.pHarray,output_tpi.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_tpi.pHarray,output_tpi.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_tpi.pHarray,output_tpi.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.tpi([1,3:end],1), ratesStruct.tpi([1,3:end],2), ratesStruct.tpi([1,3:end],4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{TPI}, NADH')

    subplot(4,3,6) %   #6     GAPDH
    % plot(output_gapdh.pHarray,output_gapdh.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    % plot(output_gapdh.pHarray,output_gapdh.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    % plot(output_gapdh.pHarray,output_gapdh.vm_uChange,'.-','color','black')
    % hold on
    errorbar(ratesStruct.gapdh(:,1), ratesStruct.gapdh(:,2), ratesStruct.gapdh(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{GAPDH}, NADH')

    subplot(4,3,7) %   #7     GAPDHr
    % plot(output_gapdhr.pHarray,output_gapdhr.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    % plot(output_gapdhr.pHarray,output_gapdhr.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    % plot(output_gapdhr.pHarray,output_gapdhr.vm_uChange,'.-','color','black')
    % hold on
    errorbar(ratesStruct.gapdhr(:,1), ratesStruct.gapdhr(:,2), ratesStruct.gapdhr(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{GAPDHr}, NADH')

    subplot(4,3,8) %   #8     PGM
    plot(output_pgm.pHarray,output_pgm.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pgm.pHarray,output_pgm.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pgm.pHarray,output_pgm.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pgm(:,1), ratesStruct.pgm(:,2), ratesStruct.pgm(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PGM}, NADH')

    subplot(4,3,9) % vmENO
    plot(output_eno_kmfixed.pHarray,output_eno_kmfixed.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_eno_kmfixed.pHarray,output_eno_kmfixed.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_eno_kmfixed.pHarray,output_eno_kmfixed.vm_uChange,'.-','color','black')
    hold on
    % errorbar(ratesStruct.eno(:,1), ratesStruct.eno(:,2)*5.5, ratesStruct.eno(:,4)*5.5,...
    errorbar(ratesStruct.eno(:,1), ratesStruct.eno(:,2), ratesStruct.eno(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{ENO}, PEP')

    subplot(4,3,10) % vmPYK
    plot(output_pyk.pHarray,output_pyk.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pyk.pHarray,output_pyk.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pyk.pHarray,output_pyk.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pyk(:,1), ratesStruct.pyk(:,2), ratesStruct.pyk(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PYK}, NADH')

    subplot(4,3,11) %   #11    PDC
    % plot(output_pdc.pHarray,output_pdc.vm_up_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    % plot(output_pdc.pHarray,output_pdc.vm_down_uChange,'.-','color',[0.5 0.5 0.5]), hold on, 
    plot(output_pdc.pHarray,output_pdc.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pdc(:,1), ratesStruct.pdc(:,2), ratesStruct.pdc(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PDC}, NADH')

    suptitleText = {'v_{m} [umol_{metaboliteAssay} mg_{P}^{-1} min^{-1}]';'Not normalized'};
    suptitle(suptitleText);
    set(gcf,'color','w');
elseif plotORerrorbar == 2 % make errorbars
    subplot(4,3,1) % vmHXK
    errorbar(output_hxk.pHarray,output_hxk.vm_uChange,output_hxk.vm_up_uChange-output_hxk.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.hxk(:,1), ratesStruct.hxk(:,2), ratesStruct.hxk(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{HXK}, NADPH')

    subplot(4,3,2) % vmPGI
    errorbar(output_pgi.pHarray,output_pgi.vm_uChange,output_pgi.vm_up_uChange-output_pgi.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pgi(:,1), ratesStruct.pgi(:,2), ratesStruct.pgi(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PGI}, NADH')

    subplot(4,3,3) % vmPFK
    errorbar(output_pfk.pHarray,output_pfk.vm_uChange,output_pfk.vm_up_uChange-output_pfk.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pfk(2:end,1), ratesStruct.pfk(2:end,2), ratesStruct.pfk(2:end,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PFK}, NADH \color{red}{*solvePoint}')

    subplot(4,3,4) % vmALD
    errorbar(output_ald.pHarray,output_ald.vm_uChange,output_ald.vm_up_uChange-output_ald.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.ald(:,1), ratesStruct.ald(:,2), ratesStruct.ald(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{ALD}, NADH')

    subplot(4,3,5) %   #5     TPI
    errorbar(output_tpi.pHarray,output_tpi.vm_uChange,output_tpi.vm_up_uChange-output_tpi.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.tpi([1,3:end],1), ratesStruct.tpi([1,3:end],2), ratesStruct.tpi([1,3:end],4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{TPI}, NADH')

    subplot(4,3,6) %   #6     GAPDH
    % errorbar(output_gapdh.pHarray,output_gapdh.vm_uChange,output_gapdh.vm_up_uChange-output_gapdh.vm_uChange,'.-','color','black')
    % hold on
    errorbar(ratesStruct.gapdh(:,1), ratesStruct.gapdh(:,2), ratesStruct.gapdh(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{GAPDH}, NADH \color{red}{*simplify?}')

    subplot(4,3,7) %   #7     GAPDHr
    % errorbar(output_gapdhr.pHarray,output_gapdhr.vm_uChange,output_gapdhr.vm_up_uChange-output_gapdhr.vm_uChange,'.-','color','black')
    % hold on
    errorbar(ratesStruct.gapdhr(:,1), ratesStruct.gapdhr(:,2), ratesStruct.gapdhr(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{GAPDHr}, NADH \color{red}{*simplify?}')

    subplot(4,3,8) %   #8     PGM 
    errorbar(output_pgm.pHarray,output_pgm.vm_uChange,output_pgm.vm_up_uChange-output_pgm.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pgm(:,1), ratesStruct.pgm(:,2), ratesStruct.pgm(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PGM}, NADH')

    subplot(4,3,9) % vmENO
    errorbar(output_eno_kmfixed.pHarray,output_eno_kmfixed.vm_uChange,output_eno_kmfixed.vm_up_uChange-output_eno_kmfixed.vm_uChange,'.-','color','black')
    hold on
    % errorbar(ratesStruct.eno(:,1), ratesStruct.eno(:,2)*5.5, ratesStruct.eno(:,4)*5.5,...
    errorbar(ratesStruct.eno(:,1), ratesStruct.eno(:,2), ratesStruct.eno(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{ENO}, PEP \color{red}{*expData?}')

    subplot(4,3,10) % vmPYK
    errorbar(output_pyk.pHarray,output_pyk.vm_uChange,output_pyk.vm_up_uChange-output_pyk.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pyk(:,1), ratesStruct.pyk(:,2), ratesStruct.pyk(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PYK}, NADH \color{red}{*solvePoint}')

    subplot(4,3,11) %   #11    PDC
%     errorbar(output_pdc.pHarray,output_pdc.vm_uChange,output_pdc.vm_up_uChange-output_pdc.vm_uChange,'.-','color','black')
    plot(output_pdc.pHarray,output_pdc.vm_uChange,'.-','color','black')
    hold on
    errorbar(ratesStruct.pdc(:,1), ratesStruct.pdc(:,2), ratesStruct.pdc(:,4),...
        's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
    title('v_{PDC}, NADH \color{red}{*reg,0.1}')

    suptitleText = {'v_{m} [umol_{metaboliteAssay} mg_{P}^{-1} min^{-1}]';'Not normalized'};
    suptitle(suptitleText);
    set(gcf,'color','w');    
end


%% step 2.2. plot normalized @pH6.81
figure(32)
%   #1     HXK
%   #2     PGI
%   #3     PFK
%   #4     ALD
%   #5     TPI
%   #6     GAPDH
%   #7     GAPDHr
%   #8     PGM
%   #9     ENO
%   #10    PYK
%   #11    PDC

subplot(4,3,1) % vmHXK
plot(output_hxk.pHarray,output_hxk.vm_up_uChange./output_hxk.vm_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_hxk.pHarray,output_hxk.vm_down_uChange./output_hxk.vm_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_hxk.pHarray,output_hxk.vm_uChange./output_hxk.vm_uChange(6),'.-','color','black')
hold on
plot(ratesStruct.hxk(:,1), ratesStruct.hxk(:,2)./ratesStruct.hxk(6,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{HXK}, NADPH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,2) % vmPGI
plot(output_pgi.pHarray,output_pgi.vm_up_uChange./output_pgi.vm_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pgi.pHarray,output_pgi.vm_down_uChange./output_pgi.vm_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pgi.pHarray,output_pgi.vm_uChange./output_pgi.vm_uChange(6),'.-','color','black')
hold on
plot(ratesStruct.pgi(:,1), ratesStruct.pgi(:,2)./ratesStruct.pgi(6,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{PGI}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,3) % vmPFK
plot(output_pfk.pHarray,output_pfk.vm_up_uChange./output_pfk.vm_up_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pfk.pHarray,output_pfk.vm_down_uChange./output_pfk.vm_down_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pfk.pHarray,output_pfk.vm_uChange./output_pfk.vm_uChange(6),'.-','color','black')
hold on
plot(ratesStruct.pfk(2:end,1), ratesStruct.pfk(2:end,2)./ratesStruct.pfk(7,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{PFK}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,4) % vmALD
plot(output_ald.pHarray,output_ald.vm_up_uChange./output_ald.vm_uChange(2),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_ald.pHarray,output_ald.vm_down_uChange./output_ald.vm_uChange(2),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_ald.pHarray,output_ald.vm_uChange./output_ald.vm_uChange(2),'.-','color','black')
hold on
plot(ratesStruct.ald(:,1), ratesStruct.ald(:,2)./ratesStruct.ald(6,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{ALD}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,5) %   #5     TPI
plot(output_tpi.pHarray,output_tpi.vm_up_uChange./output_tpi.vm_up_uChange(4),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_tpi.pHarray,output_tpi.vm_down_uChange./output_tpi.vm_down_uChange(4),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_tpi.pHarray,output_tpi.vm_uChange./output_tpi.vm_uChange(4),'.-','color','black')
hold on
plot(ratesStruct.tpi([1,3:end],1), ratesStruct.tpi([1,3:end],2)./ratesStruct.tpi(5,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{TPI}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,6) %   #6     GAPDH
% plot(output_gapdh.pHarray,output_gapdh.vm_up_uChange./output_gapdh.vm_up_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
% plot(output_gapdh.pHarray,output_gapdh.vm_down_uChange./output_gapdh.vm_down_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
% plot(output_gapdh.pHarray,output_gapdh.vm_uChange./output_gapdh.vm_uChange(6),'.-','color','black')
% hold on
plot(ratesStruct.gapdh(:,1), ratesStruct.gapdh(:,2)./ratesStruct.gapdh(6,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{GAPDH}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 10],'color','blue','LineStyle','--') %vertical line
ylim([0 10])

subplot(4,3,7) %   #7     GAPDHr
% plot(output_gapdhr.pHarray,output_gapdhr.vm_up_uChange./output_gapdhr.vm_up_uChange(5),'.-','color',[0.5 0.5 0.5]), hold on, 
% plot(output_gapdhr.pHarray,output_gapdhr.vm_down_uChange./output_gapdhr.vm_down_uChange(5),'.-','color',[0.5 0.5 0.5]), hold on, 
% plot(output_gapdhr.pHarray,output_gapdhr.vm_uChange./output_gapdhr.vm_uChange(5),'.-','color','black')
% hold on
plot(ratesStruct.gapdhr(:,1), ratesStruct.gapdhr(:,2)./ratesStruct.gapdhr(5,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{GAPDHr}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,8) %   #8     PGM
plot(output_pgm.pHarray,output_pgm.vm_up_uChange./output_pgm.vm_up_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pgm.pHarray,output_pgm.vm_down_uChange./output_pgm.vm_down_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pgm.pHarray,output_pgm.vm_uChange./output_pgm.vm_uChange(6),'.-','color','black')
hold on
plot(ratesStruct.pgm(:,1), ratesStruct.pgm(:,2)./ratesStruct.pgm(6,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{PGM}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,9) % vmENO
plot(output_eno_kmfixed.pHarray,output_eno_kmfixed.vm_up_uChange./output_eno_kmfixed.vm_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_eno_kmfixed.pHarray,output_eno_kmfixed.vm_down_uChange./output_eno_kmfixed.vm_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_eno_kmfixed.pHarray,output_eno_kmfixed.vm_uChange./output_eno_kmfixed.vm_uChange(6),'.-','color','black')
hold on
plot(ratesStruct.eno(:,1), ratesStruct.eno(:,2)./ratesStruct.eno(6,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{ENO}, PEP')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,10) % vmPYK
plot(output_pyk.pHarray,output_pyk.vm_up_uChange./output_pyk.vm_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pyk.pHarray,output_pyk.vm_down_uChange./output_pyk.vm_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pyk.pHarray,output_pyk.vm_uChange./output_pyk.vm_uChange(6),'.-','color','black')
hold on
plot(ratesStruct.pyk(:,1), ratesStruct.pyk(:,2)./ratesStruct.pyk(6,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{PYK}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

subplot(4,3,11) %   #11    PDC
% plot(output_pdc.pHarray,output_pdc.vm_up_uChange./output_pdc.vm_up_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
% plot(output_pdc.pHarray,output_pdc.vm_down_uChange./output_pdc.vm_down_uChange(6),'.-','color',[0.5 0.5 0.5]), hold on, 
plot(output_pdc.pHarray,output_pdc.vm_uChange./output_pdc.vm_uChange(6),'.-','color','black')
hold on
plot(ratesStruct.pdc(:,1), ratesStruct.pdc(:,2)./ratesStruct.pdc(6,2),...
    's','MarkerSize',5,'LineWidth',1,'MarkerEdgeColor','blue','MarkerFaceColor','cyan')
title('v_{PDC}, NADH')
line([6 8], [1 1],'color','blue','LineStyle','--') %horizontal line
line([6.81 6.81], [0 2],'color','blue','LineStyle','--') %vertical line
ylim([0 2])

suptitleText = {'v_{m} [umol_{metaboliteAssay} mg_{P}^{-1} min^{-1}]';'Normalized \color{blue}@pH6.81'};
suptitle(suptitleText);
set(gcf,'color','w');


% %% montage pics from Laura
% fileNames = cell(6,1);
% fileNames{1} = 'HXK_3.png';
% fileNames{2} = 'PGI_3.png';
% fileNames{3} = 'PFK_3.png';
% fileNames{4} = 'ALD_3.png';
% fileNames{5} = 'TPI_3.png';
% fileNames{6} = 'GAPDH_3.png';
% fileNames{7} = 'GAPDHr_3.png';
% fileNames{8} = 'PGM_3.png';
% fileNames{9} = 'ENO_3.png';
% fileNames{10} = 'PYK_3.png';
% fileNames{11} = 'PDC_3.png';
% 
% figure(33)
% %   #1     HXK
% %   #2     PGI
% %   #3     PFK
% %   #4     ALD
% %   #5     TPI
% %   #6     GAPDH
% %   #7     GAPDHr
% %   #8     PGM
% %   #9     ENO
% %   #10    PYK
% %   #11    PDC
% montage(fileNames, 'Size', [3 2])
% set(gcf,'color','w')


%% saving plots
% once finished


%% write ouput vmaxs not normalized
% hxk
output.hxk.pH = output_hxk.pHarray;
output.hxk.vm = output_hxk.vm_uChange;
output.hxk.vm_norm = output_hxk.vm_up_uChange./output_hxk.vm_uChange(6);
% pgi
output.pgi.pH = output_pgi.pHarray;
output.pgi.vm = output_pgi.vm_uChange;
output.pgi.vm_norm = output_pgi.vm_up_uChange./output_pgi.vm_uChange(6);
% pfk
output.pfk.pH = output_pfk.pHarray;
output.pfk.vm = output_pfk.vm_uChange;
output.pfk.vm_norm = output_pfk.vm_up_uChange./output_pfk.vm_up_uChange(6);
% ald
output.ald.pH = output_ald.pHarray;
output.ald.vm = output_ald.vm_uChange;
output.ald.vm_norm = output_ald.vm_up_uChange./output_ald.vm_uChange(2);
% tpi
output.tpi.pH = output_tpi.pHarray;
output.tpi.vm = output_tpi.vm_uChange;
output.tpi.vm_norm = output_tpi.vm_up_uChange./output_tpi.vm_up_uChange(4);
% pgm
output.pgm.pH = output_pgm.pHarray;
output.pgm.vm = output_pgm.vm_uChange;
output.pgm.vm_norm = output_pgm.vm_up_uChange./output_pgm.vm_up_uChange(6);
% eno
output.eno.pH = output_eno_kmfixed.pHarray;
output.eno.vm = output_eno_kmfixed.vm_uChange;
output.eno.vm_norm = output_eno_kmfixed.vm_up_uChange./output_eno_kmfixed.vm_uChange(6);
% pyk
output.pyk.pH = output_pyk.pHarray;
output.pyk.vm = output_pyk.vm_uChange;
output.pyk.vm_norm = output_pyk.vm_up_uChange./output_pyk.vm_uChange(6);
% pdc
output.pdc.pH = output_pdc.pHarray;
output.pdc.vm = output_pdc.vm_uChange;
output.pdc.vm_norm = output_pdc.vm_uChange./output_pdc.vm_uChange(6);
% gapdh_fwd
output.gapdh_fwd.pH = ratesStruct.gapdh(:,1);
output.gapdh_fwd.vm = ratesStruct.gapdh(:,2);
output.gapdh_fwd.vm_norm = ratesStruct.gapdh(:,2)./ratesStruct.gapdh(6,2);
% gapdh_rev
output.gapdh_rev.pH = ratesStruct.gapdhr(:,1);
output.gapdh_rev.vm = ratesStruct.gapdhr(:,2);
output.gapdh_rev.vm_norm = ratesStruct.gapdhr(:,2)./ratesStruct.gapdhr(5,2);

save('20200928_enzAssay_vm_vs_pH.mat','output')


%% write into a table
writetable(struct2table(output.hxk), 'hxk_vm_vs_pH.xlsx')
writetable(struct2table(output.pgi), 'pgi_vm_vs_pH.xlsx')
writetable(struct2table(output.pfk), 'pfk_vm_vs_pH.xlsx')
writetable(struct2table(output.ald), 'ald_vm_vs_pH.xlsx')
writetable(struct2table(output.tpi), 'tpi_vm_vs_pH.xlsx')
writetable(struct2table(output.pgm), 'pgm_vm_vs_pH.xlsx')
writetable(struct2table(output.eno), 'eno_vm_vs_pH.xlsx')
writetable(struct2table(output.pyk), 'pyk_vm_vs_pH.xlsx')
writetable(struct2table(output.pdc), 'pdc_vm_vs_pH.xlsx')
writetable(struct2table(output.gapdh_fwd), 'gapdh_fwd_vm_vs_pH.xlsx')
writetable(struct2table(output.gapdh_rev), 'gapdh_rev_vm_vs_pH.xlsx')






