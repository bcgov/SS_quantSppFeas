Current modeling/analyses workflows are as follows:

Species suitability & abundance modeling- with Bayesian analysis 
1. climr_getdata_plots.R- run climR downscale query for all BEC plot locations and merge with BEC plot abundance data  (outpt= clim_dat.plots.Rdata)
2. explore_veg_data.R- Compile/filter BEC plot data for tree species (output= tree_data_cleaned.Rdata)
3. feas_tables.R- a) pull in spp suitabiliy tables and merge with BEC data (output= feas_abundance_data.Rdata) b) joins w/ plot level climate data (output=feas_abund_clim_data.Rdata)
4. climPCAs.R- runs PCA on plot level climate parameters and merge back with feas/abundance data- sourced in feas.tables.R (output=feas_abund_clim_data.Rdata)
5. explore_tree_models_bayesian.R- Models BEC plot data using expert informed suitability ratings as Bayesian priors

Species Abundance modeling- with random forest 
1. climr_getdata_plots.R- run climR downscale query for all BEC plot locations and merge with BEC plot data (sourced in script 2) 
2. explore_veg_data.R- Compile/filter BEC plot data for tree species abundance modeling- Save output 
3. explore_tree_models.R- train random forest and ordinal forest models with climate and tree species data from steps 1 & 2- save trained models (in development)
4. climr_getdata_projections.R- run climR downscale query for 800m PRISM DECM for BC. Aggregate to 2km resolution. Save output 
5. project_tree_models.R- predict from trained models (step 3) onto full BC climate surface by tree species and edatope. Save output (in development)
6. plot_model_projections.R- spatially map model predictions across BC for each tree species x edatope combination. 