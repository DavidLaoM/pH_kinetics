% % SAVE_OUTPUT
if(setup.saveOutput == 1)
    saveName = ['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_parEst.mat'];
    outputName = ['output_',setup.enzymeName];
    save(saveName,outputName);
    if(setup.plotOutput == 1)
        savefig(105,['manuscript/supplementary_enzyme_by_enzyme/',setup.enzymeName,'/',setup.enzymeName, '_parEstimates.fig']);
    end
end