% % finalFigure_fig11_keq_effect
% 1 load data
% 2 prepare empty plots
% 3 fill info


%% Prepare plots: pEst estimates
figure(201)


load('tempRes_figure_Keq_hxk.mat');
subplot(4,3,1) % vmHXK
plot(output_hxk_changing_pH.pHarray, output_hxk_changing_pH.vm_uChange,'color','blue')
hold on
plot(output_hxk_changing_pH.pHarray, output_hxk_constant_pH.vm_uChange,'color','red')
title('v_{HXK}, NADPH')


load('output_pgi.mat');
subplot(4,3,2) % vmPGI
plot(output_pgi.pHarray, output_pgi.vm_uChange,'color','black')
title('v_{PGI}, NADH')


load('output_pfk.mat');
load('output_pfk_keq_constant.mat');
subplot(4,3,3) % vmPFK
plot(output_pfk.pHarray, output_pfk.vm_uChange,'color','blue')
hold on
plot(output_pfk_keq_constant.pHarray, output_pfk_keq_constant.vm_uChange,'color','red')
title('v_{PFK}, NADH')


load('output_ald.mat');
load('output_ald_keq_constant.mat');
subplot(4,3,4) % vmALD
plot(output_ald.pHarray, output_ald.vm_uChange,'color','blue')
hold on
plot(output_ald_keq_constant.pHarray, output_ald_keq_constant.vm_uChange,'color','red')
title('v_{ALD}, NADH')


load('output_tpi.mat');
subplot(4,3,5) %   #5     TPI
plot(output_tpi.pHarray, output_tpi.vm_uChange,'color','black')
title('v_{TPI}, NADH')


subplot(4,3,6) %   #6     GAPDH
title('v_{GAPDH}, NADH')


subplot(4,3,7) %   #7     GAPDHr
title('v_{GAPDHr}, NADH')


load('output_pgm.mat');
subplot(4,3,8) %   #8     PGM
plot(output_pgm.pHarray, output_pgm.vm_uChange,'color','black')
title('v_{PGM}, NADH')


load('output_eno_kmfixed.mat');
subplot(4,3,9) % vmENO
plot(output_eno_kmfixed.pHarray, output_eno_kmfixed.vm_uChange,'color','black')
title('v_{ENO}, PEP')


load('tempRes_figure_Keq_pyk.mat');
subplot(4,3,10) % vmPYK
plot(output_pyk_changeKeq.pHarray, output_pyk_changeKeq.vm_uChange,'color','blue')
hold on
plot(output_pyk_constantKeq.pHarray, output_pyk_constantKeq.vm_uChange,'color','red')
title('v_{PYK}, NADH')


load('output_pdc.mat');
load('output_pdc_keq_constant.mat');
subplot(4,3,11) %   #11    PDC
plot(output_pdc.pHarray, output_pdc.vm_uChange,'color','blue')
hold on
plot(output_pdc_keq_constant.pHarray, output_pdc_keq_constant.vm_uChange,'color','red')
title('v_{PDC}, NADH')


suptitleText = {'Parameter estimates, constant and changing Keq',...
    'umol_{metabolite} mg_{P}^{-1} min^{-1} vs pH',...
    '\color{blue}{k_{eq.Changing}} \color{red}{k_{eq.Constant}} \color{black}{k_{eq.unaffected}}'};
suptitle(suptitleText);
set(gcf,'color','w');


%% Prepare plots: pEst PSA
figure(202)

load('tempRes_figure_Keq_hxk.mat');
subplot(4,3,1) % vmHXK
legNames = cell(length(tempResult),1);
c = cool(length(tempResult));
c(3,:) = [0 0 0];
for i = 1:length(tempResult)
    legNames{i} = mat2str(keqvalsTested(i));
    plot(tempResult{i}.t,tempResult{i}.y(:,2),'.-','color',c(i,:))
    hold on
end
% colormap('cool')
% % % % lgd = legend(legNames,'location','EastOutside','Orientation','vertical');
% ylabel('NADPH @pH7.90')
% % % % xlabel('time [s]')
% % % % title(lgd,'K_{Eq} values tested')
title('v_{HXK}, NADPH')


subplot(4,3,2) % vmPGI
title('v_{PGI}, NADH')


subplot(4,3,3) % vmPFK
title('v_{PFK}, NADH')


subplot(4,3,4) % vmALD
title('v_{ALD}, NADH')


subplot(4,3,5) %   #5     TPI
title('v_{TPI}, NADH')


subplot(4,3,6) %   #6     GAPDH
title('v_{GAPDH}, NADH')


subplot(4,3,7) %   #7     GAPDHr
title('v_{GAPDHr}, NADH')


subplot(4,3,8) %   #8     PGM
title('v_{PGM}, NADH')


subplot(4,3,9) % vmENO
title('v_{ENO}, PEP')


load('tempRes_figure_Keq_pyk.mat');
subplot(4,3,10) % vmPYK
legNames = cell(length(tempResult),1);
c = cool(length(tempResult));
c(8,:) = [0 0 0];
for i = 1:length(tempResult)
    legNames{i} = mat2str(keqvalsTested(i));
    plot(tempResult{i}.t,tempResult{i}.y(:,2),'.-','color',c(i,:))
    hold on
end
% colormap('cool')
% % % % lgd = legend(legNames,'location','EastOutside','Orientation','vertical');
% ylabel('NADH @pH7.90')
% % % % xlabel('time [s]')
% % % % title(lgd,'K_{Eq} values tested')
title('v_{PYK}, NADH')


subplot(4,3,11) %   #11    PDC
title('v_{PDC}, NADH')


suptitleText = {'PSA Keq, metabolite concentration [mM] @pH7.90 vs time [s]'};
suptitle(suptitleText);
set(gcf,'color','w');

