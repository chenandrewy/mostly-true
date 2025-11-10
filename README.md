# mostly-true
Code to accompany "Most claimed statistical findings in cross-sectional return predictability are likely true" (Journal of Finance: Insights and Perspectives). Working paper version: https://arxiv.org/abs/2206.15365.

This code omits the bootstrap simulation "proof" that standard FDR estimates work in cross-sectional asset pricing. For that evidence, see the [long-version branch](https://github.com/chenandrewy/mostly-true/tree/long-version) of the code and the [previous draft](https://arxiv.org/abs/2206.15365v6) of the paper.

The code downloads all the data you need, including the Yan-Zheng (2017, RFS) data-mined strategies. Special thanks to [Sterling Yan](https://business.lehigh.edu/directory/xuemin-sterling-yan) and [Lingling Zheng](https://linglingzheng.com/research) for sharing their data.

# Repo Contents

## Excel Table
- `ez-table.xlsx`: used to create Table 1: Easy Bounds on the FDR
  - The bounds are easy, after all. No need for real code.
  - The tex input is created manually using excel2latex plugin.

## R Scripts

### Main Workflow
- **main.r**: Runs all scripts in order. 
- **config-and-functions.r**: Provides project configuration utilities and shared helper functions

### Data Preparation
- **01-prep-data.r**: Downloads and processes predictor return data from Google Drive and Dropbox; creates cleaned datasets for Chen et al., CLZ, and Yan-Zheng predictors; also runs the CLZ bootstraps
  - Creates: `data/emp_data.Rdata` (required by most other scripts)
  - Creates: `data/bootact.Rdata` (bootstrap samples for CLZ returns)

### Make Exhibits
- **02-easy-bounds.r**: Makes Figure 1: Conservative Extrapolation from Published t-stats
  - Requires: `data/emp_data.Rdata`
  - Output name: `hlz-intuition.pdf`
- **03-visual-bounds.r**: Makes Figure 2: A Visual Bound on the FDR
  - Requires: `data/emp_data.Rdata` and `data/bootact.Rdata` 
  - Output names: `dm-viz-storey-err.pdf` and `dm-viz-ez-err.pdf`
- **04-factor-adjustments-yz.r**: Makes Table 2: FDR Estimates Controlling for Value Weighting and Factor Adjustments
  - Requires: `data/emp_data.Rdata`
  - Output name: `yz-fdr.tex`
- **04alt-factor-adjustments-clz.r**: Makes an alternative version of Table 2, but is not used in paper.
  - Requires: `data/emp_data.Rdata`
  - Output name: `vw-ffn-visual-raw.tex`
- **05-hlz-simpler.r**: Makes Figure 3: Two Interpretations of Harvey, Liu, and Zhu's (2016) SMM Estimates 
  - Requires: None (fully independent, generates own simulations)
  - Output names: `hlz-simp-standard.pdf` and `hlz-simp-post-truth.pdf`

## Environment

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
- Run on MacBook M4 Pro
  - 24 GB RAM
- Run time is very roughly 5 minutes
