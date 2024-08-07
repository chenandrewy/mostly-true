# runs all bootstraps and saves to disk
# created 2024 07 to clean things up

# setup ---------------------------
rm(list = ls())
source('0-functions.r')
load('../data/emp_data.Rdata')

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
save(bootnull, file = '../data/bootnull.Rdata')

# bootstrap returns -----------------
tic = Sys.time()
bootact = bootstrap_flex(ret, nboot = set.boot$nboot, min_obs_pct = set.boot$min_obs_pct, demean = FALSE, ncore = set.boot$ncore)
toc = Sys.time()
print(paste0('bootstrap done in min: ', difftime(toc, tic, units='mins')))

# save
save(bootact, file = '../data/bootact.Rdata')
