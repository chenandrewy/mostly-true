# mostly-true
Code to accompany "Most claimed statistical findings in cross-sectional return predictability are likely true"

Code performs non- and semi-parametric FDR estimates a la Benjamini-Hochberg-Storey on the anomaly zoo, as well as simulation verifications.

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
- **05-factor-adjustments-clz.r**: Computes factor-adjusted statistics for CLZ predictors and produces LaTeX output
  - Requires: `data/emp_data.Rdata`
- **06-literature-reconciliation-hlz.r**: Simulates and visualizes Harvey-Liu-Zhu style factor selection to reconcile literature findings
  - Requires: None (fully independent, generates own simulations)
- **07-hlz-simpler.r**: Provides a simplified HLZ simulation for two key comparison figures
  - Requires: None (fully independent, generates own simulations)

### Simulations
- **91-run-bootstraps.r**: Runs bootstrap simulations for CLZ returns to support downstream inference
  - Requires: `data/emp_data.Rdata`
  - Creates: `data/bootnull.Rdata`, `data/bootact.Rdata`
- **92-bootstrap-validation.r**: Validates bootstrap null distribution against a standard normal benchmark
  - Requires: `data/emp_data.Rdata`, `data/bootnull.Rdata`
- **93-simulations-theory-free.r**: Runs theory-free simulations comparing empirical and simulated discovery rates
  - Requires: `data/emp_data.Rdata`, `data/bootnull.Rdata`
- **94-simulations-publication-bias.r**: Simulates publication-bias scenarios to benchmark CZ extrapolation
  - Requires: `data/emp_data.Rdata`
- **95-simulations-standard-errors.r**: Simulates standard-error estimation accuracy for large-scale inference histograms
  - Requires: None (fully independent, generates own simulations)

## Dependency Structure

**To run everything**: Simply run `Rscript main.r` — it executes all scripts in the correct order.

**For individual scripts**, most can be run independently after data preparation:

1. **Core dependency**: Run `01-prep-data.r` first (downloads data, creates `emp_data.Rdata`)
2. **Bootstrap dependency**: Run `91-run-bootstraps.r` before 92 and 93 (creates bootstrap samples)
3. **Fully independent**: Scripts 06, 07, and 95 generate their own simulations and have no data dependencies

### Option 1: Run everything (recommended)
```r
Rscript main.r
```
This runs all scripts in the proper dependency order:
1. Data preparation (01)
2. Main analyses (02-07)
3. Simulations (91-95)

### Option 2: Run scripts individually
```r
# 1. Prepare data (required for most analyses)
source("01-prep-data.r")

# 2. Run main analyses (can run in any order)
source("02-easy-bounds.r")
source("03-visual-bounds.r")  # auto-generates bootstrap if needed
source("04-factor-adjustments-yz.r")
source("05-factor-adjustments-clz.r")
source("06-literature-reconciliation-hlz.r")  # independent
source("07-hlz-simpler.r")  # independent

# 3. Run bootstrap simulations (only needed for 92, 93)
source("91-run-bootstraps.r")
source("92-bootstrap-validation.r")
source("93-simulations-theory-free.r")

# 4. Run other simulations (independent or depend on 01 only)
source("94-simulations-publication-bias.r")
source("95-simulations-standard-errors.r")  # independent
```  


