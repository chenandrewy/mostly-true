# ABOUTME: Runs bootstrap simulations for CLZ returns to support downstream inference.
# Inputs:
#   - functions.r (bootstrap helpers)
#   - data/emp_data.Rdata (predictor returns from 01-prep-data.r)
# Outputs:
#   - data/bootnull.Rdata
#   - data/bootact.Rdata
# How to run:
#   Rscript 91-run-bootstraps.r
#   Rscript 91-run-bootstraps.r --vanilla

# setup ---------------------------
rm(list = ls())

library(here)
here::i_am("91-run-bootstraps.r")

source(here("functions.r"))

data_dir <- here("data")
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

bootnull_path <- file.path(data_dir, "bootnull.Rdata")
bootact_path <- file.path(data_dir, "bootact.Rdata")

load(file.path(data_dir, "emp_data.Rdata"))

set.boot = list(
  nboot = 1000,
  min_obs_pct = 50,
  ncore = 8
)

ret = clz_ret

# bootstrap residuals ------------------

tic = Sys.time()
bootnull = bootstrap_flex(ret, nboot = set.boot$nboot, min_obs_pct = set.boot$min_obs_pct, demean = TRUE, ncore = set.boot$ncore, output_cor = TRUE)
toc = Sys.time()
print(paste0('bootstrap done in min: ', difftime(toc, tic, units='mins')))

# save
save(bootnull, file = bootnull_path)

# bootstrap returns -----------------
tic = Sys.time()
bootact = bootstrap_flex(ret, nboot = set.boot$nboot, min_obs_pct = set.boot$min_obs_pct, demean = FALSE, ncore = set.boot$ncore)
toc = Sys.time()
print(paste0('bootstrap done in min: ', difftime(toc, tic, units='mins')))

# save
save(bootact, file = bootact_path)
