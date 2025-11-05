# mostly-true-short
Code to accompany "Most claimed statistical findings in cross-sectional return predictability are likely true" (Journal of Finance: Insights and Perspectives). 

This code is a shorter and cleaner version of the `mostly-true` repo, which had results for the longer previous drafts of the paper. For simulation evidence that the FDR estimates work, see the `mostly-true` repo.

# Exhibits in the paper
- Table 1: Easy Bounds on the FDR
  - `ez-combined.tex`
  - Created with `ez-table.xlsx`
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
  - Created by: `hlz-simpler.r`
  - Output names: `hlz-simp-standard.pdf` and `hlz-simp-post-truth.pdf`


## R Scripts

### Main Workflow
- **main.r**: Coordinates the full analysis workflow by sourcing project scripts in sequence
- **config-and-functions.r**: Provides project configuration utilities and shared helper functions

### Data Preparation
- **01-prep-data.r**: Downloads and processes predictor return data from Google Drive and Dropbox; creates cleaned datasets for Chen et al., CLZ, and Yan-Zheng predictors
  - Creates: `data/emp_data.Rdata` (required by most other scripts)

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

## Dependency Structure

**To run everything**: Simply run `Rscript main.r` — it executes all scripts in the correct order.

**For individual scripts**, most can be run independently after data preparation:

1. **Core dependency**: Run `01-prep-data.r` first (downloads data, creates `emp_data.Rdata`)
2. **Fully independent**: Script 05 generates its own simulations and has no data dependencies

### Option 1: Run everything (recommended)
```r
Rscript main.r
```
This runs all scripts in the proper dependency order:
1. Data preparation (01)
2. Main analyses (02-05)

### Option 2: Run scripts individually
```r
# 1. Prepare data (required for most analyses)
source("01-prep-data.r")

# 2. Run main analyses (can run in any order)
source("02-easy-bounds.r")
source("03-visual-bounds.r")  # auto-generates bootstrap if needed
source("04-factor-adjustments-yz.r")
source("04alt-factor-adjustments-clz.r")
source("05-hlz-simpler.r")  # independent
```  


