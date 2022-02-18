% % temp.m
% 1. Run parameter estimation in prallel
set_paths_pHstudy
% % % % delete(gcp('nocreate'))
% % % % parpool(4)
% % % % parfor i = 1:4
% % % %     if i == 1
% % % % dev_pEst_ald;
% % % % dev_pEst_hxk;
% % % % dev_pEst_tpi;
% % % %     elseif i == 2
% % % % dev_pEst_pgm;
% % % % dev_pEst_pfk;
% % % % dev_pEst_eno;
% % % %     elseif i == 3
dev_pEst_pdc;
dev_pEst_pyk;
dev_pEst_gapdh_fwd;
% % % %     elseif i == 4
% % % % dev_pEst_pgi;
% % % % dev_pEst_gapdh_rev;
% % % %     end
% % % % end
% % % % 
