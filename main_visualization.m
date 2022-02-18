% % MAIL_VISUALIZATION.m
% The plots obtained in the paper are generated
close all
dbstop if error
set_paths_pHstudy;
setup.saveOutput = 0;
    % 0: do not save plots
    % 1: asve plots in 'data/processed_data' folder ('setup.plotOutput' must be equal to 1).
setup.fast_option = 0;
    % 0: the entire code is run    
    % 1: some materials that were saved are loaded, to obtain the plots quickly

% which figure to run? (select one at a time)
% core article
setup.sel.m_figure_enolase = 1;
setup.sel.m_figure_fits = 0;
setup.sel.m_figure_hexokinase = 0;
% supplementary materials
setup.sel.ms_figure_ald_r2_vs_mw = 0;
setup.sel.ms_figure_ald_vm_vs_df = 0;
setup.sel.ms_figure_gapdh_direct_vs_curve_fit = 0;
setup.sel.ms_figure_gapdh_literature = 0;
setup.sel.ms_figure_keq_pH_dependency = 0;
setup.sel.ms_figure_style_train_test = 0;
setup.sel.ms_figure_vm_exp_vs_df = 0;
setup.sel.ms_table_complete_vmax_pH_dependency_polynomials = 0;
    
% https://nl.mathworks.com/matlabcentral/answers/231799-exported-pdf-figures-have-inaccurate-dimensions
% core article
if setup.sel.m_figure_enolase == 1
    tic
    m_figure_enolase, %end
    toc
elseif setup.sel.m_figure_fits == 1
    tic
    m_figure_fits, %end
    toc
elseif setup.sel.m_figure_hexokinase == 1
    tic
    m_figure_hexokinase, %end
    toc
% supplementary materials
elseif setup.sel.ms_figure_ald_r2_vs_mw == 1
    tic
    ms_figure_ald_r2_vs_mw, %end
    toc
elseif setup.sel.ms_figure_ald_vm_vs_df == 1
    tic
    ms_figure_ald_vm_vs_df, %end
    toc
elseif setup.sel.ms_figure_gapdh_direct_vs_curve_fit == 1
    tic
    ms_figure_gapdh_direct_vs_curve_fit, %end
    toc
elseif setup.sel.ms_figure_gapdh_literature == 1
    tic
    ms_figure_gapdh_literature, %end
    toc
elseif setup.sel.ms_figure_keq_pH_dependency == 1
    tic
    ms_figure_keq_pH_dependency, %end
    toc
elseif setup.sel.ms_figure_style_train_test == 1
    tic
    ms_figure_style_train_test, %end
    toc
elseif setup.sel.ms_figure_vm_exp_vs_df == 1
    tic
    ms_figure_vm_exp_vs_df, %end
    toc
elseif setup.sel.ms_table_complete_vmax_pH_dependency_polynomials == 1
    tic
    ms_table_complete_vmax_pH_dependency_polynomials, %end
    toc
end
