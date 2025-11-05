# mostly-true-short
Code to accompany "Most claimed statistical findings in cross-sectional return predictability are likely true" (Journal of Finance: Insights and Perspectives). 

This code is a shorter and cleaner version of the `mostly-true` repo, which had results for the longer previous drafts of the paper. For simulation evidence that the FDR estimates work, see the `mostly-true` repo.

# Exhibits in the paper
- Table 1: Easy Bounds on the FDR
  - `ez-combined.tex`
  - Created from `ez-table.xlsx` w/ excel2latex plugin (manually)
- Figure 1: Conservative Extrapolation from Published t-stats
  - Created by `02-easy-bounds.r`
  - Output name: `hlz-intuition.pdf`
- Figure 2: A Visual Bound on the FDR
  - Created by `03-visual-bounds.r`
  - `dm-viz-storey-err.pdf` 
  - `dm-viz-ez-err.pdf`
- Table 2: FDR Estimates Controlling for Value Weighting and Factor Adjustments
  - Created by `04-factor-adjustments-yz.r`
  - Output name: `yz-fdr.tex`
- Figure 3: Two Interpretations of Harvey, Liu, and Zhu's (2016) SMM Estimates 
  - Created by: `05-hlz-simpler.r`
  - Output names: `hlz-simp-standard.pdf` and `hlz-simp-post-truth.pdf`


## R Scripts

### Main Workflow
- **main.r**: Coordinates the full analysis workflow by sourcing project scripts in sequence
- **config-and-functions.r**: Provides project configuration utilities and shared helper functions

### Data Preparation
- **01-prep-data.r**: Downloads and processes predictor return data from Google Drive and Dropbox; creates cleaned datasets for Chen et al., CLZ, and Yan-Zheng predictors; also runs the CLZ bootstraps
  - Creates: `data/emp_data.Rdata` (required by most other scripts)
  - Creates: `data/bootact.Rdata` (bootstrap samples for CLZ returns)

### Main Analysis
- **02-easy-bounds.r**: Visualizes easy- and visual-bound intuition figures used in the paper
  - Requires: `data/emp_data.Rdata`
- **03-visual-bounds.r**: Builds diagnostic figures and tables contrasting empirical distribution with null models
  - Requires: `data/emp_data.Rdata`, optionally `data/bootact.Rdata` (regenerates if missing)
- **04-factor-adjustments-yz.r**: Estimates factor-adjusted performance for Yan-Zheng predictors and exports LaTeX tables
  - Requires: `data/emp_data.Rdata`
- **04alt-factor-adjustments-clz.r**: Computes factor-adjusted statistics for CLZ predictors and produces LaTeX output
  - Requires: `data/emp_data.Rdata`
- **05-hlz-simpler.r**: Provides a simplified HLZ simulation for two key comparison figures
  - Requires: None (fully independent, generates own simulations)

## R Environment

- R version 4.5.0 (2025-04-11) — "How About a Twenty-Six"
- Package versions used:
  - here 1.0.1
  - tidyverse 2.0.0
  - data.table 1.17.8
  - ggplot2 4.0.0
  - extrafont 0.19
  - latex2exp 0.9.6
  - kableExtra 1.4.0
  - googledrive 2.1.2
  - haven 2.5.5  
  - foreach 1.5.2
  - doParallel 1.0.17
