%% 1) load and rename matrixes
clear
% 
load('ref_ald_parEst.mat'); ref_ald = output_ald; clear output_ald
load('ref_pdc_parEst.mat'); ref_pdc = output_pdc; clear output_pdc
load('ref_pfk_parEst.mat'); ref_pfk = output_pfk; clear output_pfk
load('ref_pyk_parEst.mat'); ref_pyk = output_pyk; clear output_pyk
% 
load('old_keqconst_ald.mat'); keqconst_ald = output_ald_keq_constant; clear output_ald_keq_constant
load('old_keqconst_pdc.mat'); keqconst_pdc = output_pdc_keq_constant; clear output_pdc_keq_constant
load('old_keqconst_pfk.mat'); keqconst_pfk = output_pfk_keq_constant; clear output_pfk_keq_constant
load('old_keqconst_pyk.mat','output_pyk_constantKeq'); keqconst_pyk = output_pyk_constantKeq; clear output_pyk_constantKeq

%% 2) plot it
fig_h101 = figure(101);
% ALD
sp1 = subplot(2,2,1);
plot(ref_ald.pHarray, ref_ald.vm_uChange, 'b.-')
hold on
plot(keqconst_ald.pHarray, keqconst_ald.vm_uChange, 'r.-')
legend('Keq pH-dependent', 'Keq pH-independent')
title('ALD')
hold off
% PDC
sp2 = subplot(2,2,2);
plot(ref_pdc.pHarray, ref_pdc.vm_uChange, 'b.-')
hold on
plot(keqconst_pdc.pHarray, keqconst_pdc.vm_uChange, 'r.-')
legend('Keq pH-dependent', 'Keq pH-independent')
title('PDC')
hold off
% PFK
sp3 = subplot(2,2,3);
plot(ref_pfk.pHarray, ref_pfk.vm_uChange, 'b.-')
hold on
plot(keqconst_pfk.pHarray, keqconst_pfk.vm_uChange, 'r.-')
legend('Keq pH-dependent', 'Keq pH-independent')
title('PFK')
hold off
% PYK
sp4 = subplot(2,2,4);
plot(ref_pyk.pHarray, ref_pyk.vm_uChange, 'b.-')
hold on
plot(keqconst_pyk.pHarray, keqconst_pyk.vm_uChange, 'r.-')
legend('Keq pH-dependent', 'Keq pH-independent')
title('PYK')
hold off


% set position
set(101, 'Position', [100 100 750 750])


%% needed stop not to get truncated... #matlabUselessSecrets
% % save
% % specs printing (method 3)
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters fits

% %%
% ref_ald.vm
% ref_ald.vm_uChange