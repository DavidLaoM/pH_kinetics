% % SAVE_INITIAL_RECALL
if(setup.saveOutput == 1)
    tempFolderName = ['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName];
    if(~exist(tempFolderName,'dir'))
        mkdir('manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName)
    end
    save(['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_initial_variables.mat'],'Vmax_mw_opt_corr','idxs2consider','DF','pH');
end
if((setup.plotOutput == 1)&&(setup.saveOutput == 1))
    savefig(1,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_concentrations_basezero.fig']);
    savefig(2,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_df.fig']);
    savefig(3,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_mw_vmax_vs_movingWindow.fig']);
    savefig(4,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_mw_R2_vs_movingWindow.fig']);
    savefig(5,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_mw_iniPoint_vs_movingWindow.fig']);
    savefig(6,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_experimental_vmax_vs_pH.fig']);
end