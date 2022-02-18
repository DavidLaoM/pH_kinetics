% % temp.m
% 1. Run parameter estimation in prallel
set_paths_pHstudy
% % % % delete(gcp('nocreate'))
% % % % parpool(4)
% % % % parfor i = 1:4
% % % %     if i == 1
dev_pEst_ald;
dev_pEst_hxk;
dev_pEst_tpi;
% % % %     elseif i == 2
dev_pEst_eno;
dev_pEst_pfk;
% % % %     elseif i == 3
dev_pEst_gapdh_fwd;
dev_pEst_pdc;
dev_pEst_pyk;
% % % %     elseif i == 4
dev_pEst_gapdh_rev;
dev_pEst_pgi;
dev_pEst_pgm;
% % % %     end
% % % % end

