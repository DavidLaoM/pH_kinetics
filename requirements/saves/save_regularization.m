% % SAVE_REGULARIZATION
if(setup.saveOutput == 1)
    saveName = ['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_regularizationResults.mat'];
    save(saveName,'array_xres','lambdalist','eParameters','eData','array_eParams','array_eData','xres_selected','array_resnorm','array_residual','array_Jacobian','array_raw_error')
end
if((setup.plotOutput == 1)&&(setup.saveOutput == 1))
    savefig(101,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_trainData_fit_metabolites_reg.fig']);
    savefig(102,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_testData_fit_fluxes_reg.fig']);
    savefig(103,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_regularization.fig']);
end