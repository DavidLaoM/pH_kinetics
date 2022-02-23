# pH_kinetics
 Descrition:
 pH dependency of glycolytic enzymes. Progression curve analysis.
 
 Citation:
 TBA
 
 Keywords:
 Saccharomyces cerevisiae, glycolysis, enzymes, pH dependency, progression curve analysis
 
 Model overview:
 The models consist of ordinary differential equations (ODE) representing the mass balances involved in each reaction. 
 Reactions covered in yeast glycolysis:
 - HXK, E.C.: 2.7.1.1
 - PGI, E.C.: 5.3.1.9
 - PFK, E.C.: 2.7.1.11
 - ALD, E.C.: 4.1.2.13
 - TPI, E.C.: 5.3.1.1
 - GAPDH, E.C.: 1.2.1.12
 - PGM, E.C.: 5.4.2.11
 - ENO, E.C.: 4.2.1.11
 - PYK, E.C.: 2.7.1.40
 - PDC, E.C.: 4.1.1.1
 pH dependency was checked for the pH range [6.19:7.91]
 
 Required software:
 - A functional matlab installation (MATLAB 9.3 or higher).

 Reproducibility:
 - The entire computational work, including article figures and supplementary figures, can be reproduced using matlab. For this we recommend runnign the 'main_(...).m' files in the root folder.
 - Models were also generated in python and SBML format. These can be found in the folder /models, /python and /sbml, respectively.

 Usage:
 Three main scripts have been created to drive the user along the pipeline.
 - main_import_data: imports and organises the raw data
 - main_estimate_parameters: estimates the paraemters in the different progression curve analysis experiments
 - main_manuscript_visualizations: produces the plots found in the article

 Contributors:
 TBA
 
 Folder structure:
 - data
	- processed data
	- raw data
 - manuscript: including supplementary information.
	- supplementary
	- supplementary_enzyme_by_enzyme
 - models
	- models
 - parameter_estimation
	- parameter_estimation
 - requirements
	- additional_assays
	- control
	- cost_functions
	- experimental_data_import
	- file_exchange
	- saves
	- simulations
 - set_paths_pHstudy.m
 - main_import_data.m
 - main_estimate_parameters.m
 - main_manuscript_visualizations.m
