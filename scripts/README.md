Current modeling/analyses workflows are as follows:

Species suitability & abundance modeling- with Bayesian analysis 
1. climr_getdata_plots.R- run climR downscale query for all BEC plot locations and merge with BEC plot abundance data  (outpt= clim_dat.plots.Rdata)
2. explore_veg_data.R- Compile/filter BEC plot data for tree species (output= tree_data_cleaned.Rdata)
3. Bayesian_analysis/feas_tables.R- a) pull in spp suitabiliy tables and merge with BEC data (output= feas_abundance_data.Rdata) b) joins w/ plot level climate data (output=feas_abund_clim_data.Rdata)
4. Bayesian_analysis/climPCAs.R- runs PCA on plot level climate parameters and merge back with feas/abundance data- sourced in feas.tables.R (output=feas_abund_clim_data.Rdata)
5. Bayesian_analysis/explore_tree_models_bayesian.R- Models BEC plot data using expert informed suitability ratings as Bayesian priors
