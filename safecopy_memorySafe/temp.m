% % temp.m
% 1. Run parameter estimation in prallel
set_paths_pHstudy
delete(gcp('nocreate'))
parpool(4)
parfor i = 1:4
    if i == 1
        x = dev_pEst_ald;
        x = dev_pEst_hxk;
        x = dev_pEst_tpi;
    elseif i == 2
        x = dev_pEst_eno;
        x = dev_pEst_pfk;
    elseif i == 3
        x = dev_pEst_gapdh_fwd;
        x = dev_pEst_pdc;
        x = dev_pEst_pyk;
    elseif i == 4
        x = dev_pEst_gapdh_rev;
        x = dev_pEst_pgi;
        x = dev_pEst_pgm;
    end
end

