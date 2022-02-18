% control initial variables
% control parameter estimation
% control regularization

% loadName
loadName_ini_var = [setup.enzymeName, '_initial_variables.mat'];
%     loadName_par_est = [setup.enzymeName, '_parEst.mat'];
loadName_reg = [setup.enzymeName, '_regularizationResults.mat'];
loadName_con_ini_var = ['control_', setup.enzymeName, '_initial_variables.mat'];
%     loadName_con_par_est = ['control_', setup.enzymeName, '_parEst.mat'];
loadName_con_reg = ['control_', setup.enzymeName, '_regularizationResults.mat'];
% load files
load(loadName_ini_var)
    ini_var.DF = DF;
    ini_var.idxs2consider = idxs2consider;
    ini_var.pH = pH;
    ini_var.Vmax_mw_opt_corr = Vmax_mw_opt_corr;
load(loadName_reg)
    reg.error_sum = sum(abs(array_eData{1}));
load(loadName_con_ini_var)
    ini_var_con.DF = DF;
    ini_var_con.idxs2consider = idxs2consider;
    ini_var_con.pH = pH;
    ini_var_con.Vmax_mw_opt_corr = Vmax_mw_opt_corr;
load(loadName_con_reg)
    reg_con.error_sum = sum(abs(array_eData{1}));
if((ini_var.DF(end,end) == ini_var_con.DF(end,end)) &&(ini_var.idxs2consider(end,end) == ini_var_con.idxs2consider(end,end))&&(ini_var.pH(end,end) == ini_var_con.pH(end,end))&&(ini_var.Vmax_mw_opt_corr(end,end) == ini_var_con.Vmax_mw_opt_corr(end,end))&&(reg.error_sum == reg_con.error_sum))
    disp(['Does ',setup.enzymeName,' estimation match control safecopy? YES.'])
else
    disp(['Does ',setup.enzymeName,' estimation match control safecopy? NO.'])
end
