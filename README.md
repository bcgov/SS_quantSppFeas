# SS_quantSppFeas
Quantitative analyses on current species/site feasibility

---
## Readme 
---
This repository holds data, scripts, model outputs and figures of preliminary analyses of species abundance and climatic suitability predictions for BC tree species at the BEC site-series level.

## Features
TBD

## Usage
Current scripting workflow for species feasibility modeling is as follows 
1. climr_getdata_plots.R- run climR downscale query for all BEC plot locations and merge with BEC plot data (sourced in script 2) 
2. explore_veg_data.R- Compile/filter BEC plot data for tree species abundance modeling- Save output 
3. explore_tree_models.R- train random forest and ordinal forest models with climate and tree species data from steps 1 & 2- save trained models (in development)
4. climr_getdata_projections.R- run climR downscale query for 800m PRISM DECM for BC. Aggregate to 2km resolution. Save output 
5. project_tree_models.R- predict from trained models (step 3) onto full BC climate surface by tree species and edatope. Save output (in development)
6. plot_model_projections.R- spatially map model predictions across BC for each tree species x edatope combination. 

## Requirements
TBD

## Installation
TBD

## Project Status
Initializing 

## Goals/Roadmap
Inform and refine current BC species feasibility tables including CFRG and ecological feasibility tables for downstream use in Climate Change Informed Species Selection (CCISS) tools

Develop site-specific (micro-climatically relevant) species distribution/suitability models 

## Getting Help or Reporting an Issue
send an email to Courtney.Collins[at]gov.bc.ca

## How to Contribute
Government employees, public and members of the private sector are encouraged to contribute to the repository by forking and submitting a pull request.

(If you are new to GitHub, you might start with a basic tutorial and check out a more detailed guide to pull requests.)

Pull requests will be evaluated by the repository guardians on a schedule and if deemed beneficial will be committed to the master.

All contributors retain the original copyright to their stuff, but by contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the license under which this project is distributed.


"Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms."

## License

  Copyright 2019 Province of British Columbia

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at 

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

Detailed guidance around licenses is available 
[here](LICENSE.md)
