# ABOUTME: Simulates standard-error estimation accuracy for large-scale inference histograms.
# Inputs:
#   - config-and-functions.r (utility functions including bootstrap_flex)
# Outputs:
#   - results/sim-se.pdf
#   - results/sim-se-count.pdf
# How to run:
#   Rscript 95-simulations-standard-errors.r
#   Rscript 95-simulations-standard-errors.r --vanilla

# Consistent with Efron's intuition, if T is small compared to N, then the nonpar bootstrap SE is too large

# N = 10000 # number of stocks / genes
# T =  1*12 # number of months / patients
# set.cor = list(type = "indep")

# On the other hand, it does seem that when T and N are comparable to the anomaly returns, the bias is very small. For example, if N = 1000 and T = 120 months then the bias is mostly gone for tabs < 3

# you really need a huge rho, like 0.9999 to get a difference in RMS correlations. But then the sims can be badly behaved.
# For rho = 0.999, the bootstrap actually works pretty well. The efron SE is too small, at least the way I calculated it (which could be wrong)

# For realistic N and T realistic, the bootstrap actually works pretty well.

# Right now, the Efron spline SE is incomplete. It doesn't find the central disperson using a quadratic approx and instead just uses the actual central dispersion.

# Setup -------------------------------------------

rm(list = ls())

library(here)
here::i_am("95-simulations-standard-errors.r")

source(here("config-and-functions.r"))

paths <- project_paths()
results_dir <- paths$results

start_time <- Sys.time()

set.seed(123) # for reproducibility
nboot <- 1000
nsim <- 1000
edge <- c(-Inf, seq(-20, 20, 0.5), Inf)
nround <- 3
sig_mu <- 0.50 # sd of mu (see Chen-Zim RAPS Table 3)
sig_noise <- 7

N <- 29300 # number of stocks / genes
T <-  600 # number of months / patients
# set.cor = list(type = "indep")
set.cor <- list(type = "fac", n_fac = 20, w_fac = 0.5) # super fast
# set.cor = list(type = "ar1", rho = 0.99) # can be slow

# check cores: detectCores()
nsamp_for_cor <- 1000
ncore <- 8
targetvar <- "tabs"

# Function for simulating data ------------------

generate_panel_data <- function() {
  # simulate noise component based on set.cor
  if (set.cor$type == "ar1") {
    panel <- foreach(i = 1:N, .combine = rbind) %do% {
      if (i == 1) {
        noise <- rnorm(T, 0, sig_noise)
      } else {
        noise <- rnorm(T, set.cor$rho * lastnoise, sqrt((1 - set.cor$rho^2) * sig_noise^2))
      }
      lastnoise <- noise
      data.table(stock = i, month = 1:T, noise = noise)
    }

    # simulate noise with n_fac factors
  } else if (set.cor$type == "fac") {
    # generate eigenvectors for each time period(T x n_fac)
    eigvec <- matrix(rnorm(T * set.cor$n_fac), nrow = T, ncol = set.cor$n_fac)
    eigvec <- apply(eigvec, 2, function(x) x / sqrt(sum(x^2)))

    # generate normalized eigenvalues
    eigval <- abs(rnorm(set.cor$n_fac))
    eigval <- eigval / sum(eigval) * T %>% sort(decreasing = TRUE)

    # generate factor realizations (T x n_fac)
    fac <- eigvec %*% diag(sqrt(eigval))

    # stock loadings on factors (N x n_fac)
    beta <- matrix(rnorm(N * set.cor$n_fac), nrow = N, ncol = set.cor$n_fac)

    # generate matrix of noise realizations (T x N)
    # noisemat <- fac %*% t(beta) # old    
    ep = rnorm(N*T) %>% matrix(T, N)
    noisemat = sqrt(set.cor$w_fac) * fac %*% t(beta) + sqrt(1-set.cor$w_fac) * ep       

    # rescale variance to match sig_noise
    tempvar <- mean(colMeans(noisemat^2))
    noisemat <- noisemat * sqrt(sig_noise^2 / tempvar)

    # # debug check
    # sqrt(colMeans(noisemat^2)) %>% summary()
    # sqrt(colMeans((fac %*% t(beta))^2)) %>% summary()

    # convert to long format data table
    panel <- data.table(noisemat) %>%
      mutate(month = 1:T) %>%
      melt(id.vars = "month", variable.name = "stock", value.name = "noise") %>%
      .[, stock := as.integer(stock)] %>%
      setorder(stock, month)
  } else if (set.cor$type == "indep") {
    panel <- expand.grid(stock = 1:N, month = 1:T) %>%
      as.data.table() %>%
      mutate(noise = rnorm(N * T, 0, sig_noise))
  } # end if set.cor$type

  cross <- data.table(stock = 1:N, mu = rnorm(N, 0, sig_mu))
  panel <- panel[cross, on = .(stock), mu := i.mu] %>%
    .[, r := mu + noise]

  sumstat <- panel[, .(tstat = mean(r) / sd(r) * sqrt(T), se = sd(r) / sqrt(T)), by = stock] %>%
    mutate(tabs = abs(tstat))

  return(list(panel = panel, sumstat = sumstat))
}

# Simulate empirical data -----------------------

emp <- generate_panel_data()

# check volatility (should match sig_noise)
emp$panel[ , .(vol = sd(r)), by = stock] %>% summary()


# Find actual standard error by simulation ---------------

# simulate true model nsim times

# Set up parallel backend
registerDoParallel(cores = ncore)

sim <- foreach(i = 1:nsim, .combine = rbind, .packages = c("data.table", "dplyr", "foreach")) %dopar% {
  temp <- generate_panel_data()
  temp$sumstat$sim <- i
  temp$sumstat
}

# Stop parallel backend
stopImplicitCluster()

# find SE of histogram
sim_hist <- histogram_by_group(sim, edge, varname = targetvar, group = "sim")
se_sim <- sim_hist %>%
  group_by(mids) %>%
  summarize(
    se_sim = sd(count / ntotal),
    fhist = mean(count / ntotal),
    fhist05 = quantile(count / ntotal, 0.05),
    fhist95 = quantile(count / ntotal, 0.95), .groups = "drop"
  )

# Bootstrap SE -------------------------------------------
bootdat <- bootstrap_flex(ret = emp$panel, nboot = nboot, coli = "stock", colt = "month", colr = "r", demean = FALSE, ncore = ncore)

# find boot SE
bootdat <- bootdat %>% mutate(tstat = mean / vol * sqrt(nmonth), tabs = abs(tstat))
boothist <- histogram_by_group(bootdat, edge, varname = targetvar, group = "booti")
se_boot <- boothist %>%
  group_by(mids) %>%
  summarize(
    se_boot = sd(count / ntotal),
    fhist = mean(count / ntotal),
    fhist05 = quantile(count / ntotal, 0.05),
    fhist95 = quantile(count / ntotal, 0.95), .groups = "drop"
  )

# Find Efron SE --------------------------

# Efron 2009 "Correlated z-values..." Equation 1.4
bin_width <- mean(diff(edge[is.finite(edge)]))

# find rms correlation
tempi <- sample(1:N, nsamp_for_cor)
tempwide <- dcast(emp$panel[stock %in% tempi], month ~ stock, value.var = "r")
cmat <- cor(tempwide[, -1])
alpha <- sqrt(mean(cmat[lower.tri(cmat)]^2))

# fit histogram
temp <- hist(emp$sumstat %>% pull(targetvar), breaks = edge, plot = FALSE)
emp$hist <- data.table(mids = temp$mids, fhist = temp$counts / N) %>%
  mutate(mids = round(mids, nround)) %>%
  filter(!is.infinite(fhist + mids), fhist > 0)

# second derivative (Efron book (7.46) cov_1)
sig_act <- sqrt((sig_mu / sig_noise * sqrt(T))^2 + 1)

f2actual <- ((emp$hist$mids - 0)^2 - sig_act^2) / sig_act^4 * dnorm(emp$hist$mids, 0, sig_act)

se_efron <- emp$hist %>%
  mutate(se_bern = sqrt(fhist * (1 - fhist) / N),
    Vpen = (bin_width * sig_act^2 * alpha / sqrt(2) * f2actual)^2,
    se_closed = sqrt(se_bern^2 + Vpen)
  )

# DEBUG ====
if (FALSE) {
se_efron[mids > -1 & mids < 1] 
sd(emp$sumstat$tstat)
sig_act

# the difference is in Vpen, let's break out the components
bin_width
sig_act
alpha
tibble(
  mids = emp$hist$mids,
  f2 = ((mids - 0)^2 - sig_act^2) / sig_act^4 * dnorm(mids, 0, sig_act),
  Vpen6 = (bin_width * sig_act^2 * alpha / sqrt(2) * f2)^2*1e6,
  Vpen_p1 = (bin_width * sig_act^2 * alpha / sqrt(2))^2*1e6
) %>%
  filter(mids > -1 & mids < 1)

# omg the alpha is everything. "Small" errors in alpha multiply everything. The alpha is estimated with noise, which seems to bias it upward when all correlations are truly zero

}

# Plots -------------------------------------------

## Plot Freq ====
groupdat <- tibble(
  levels = c("se_sim", "se_boot", "se_bern", "se_closed"),
  labels = c("Simulation", "Bootstrap", "Bernoulli", "Efron Closed Form")
)

# Omits Efron for now
plotme <- data.table(mids = emp$hist$mids, freq = emp$hist$fhist) %>%
  mutate(mids = round(mids, nround)) %>%
  filter(freq > 0) %>%
  merge(se_sim, by = "mids", all.x = TRUE) %>%
  mutate(se_bern = sqrt(freq * (1 - freq) / N)) %>%
  merge(se_boot %>% select(mids, se_boot), by = "mids", all.x = TRUE) %>%
  merge(se_efron %>% select(mids, se_closed), by = "mids", all.x = TRUE) %>%
  filter(is.finite(mids)) %>%
  pivot_longer(cols = c(se_sim, se_bern, se_boot, se_closed), names_to = "method", values_to = "se") %>%
  mutate(method = factor(method, levels = groupdat$levels, labels = groupdat$labels)) %>%
  filter(!is.na(se), method != "Efron Closed Form")

plt <- plotme %>%
  mutate(se = se * 100) %>%
  ggplot(aes(x = mids, y = se, color = method, linetype = method)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("Simulation" = MATBLUE, "Bernoulli" = MATRED, "Bootstrap" = MATYELLOW)) +
  theme(legend.position = c(0.8, 0.8), text = element_text(size = 22), legend.title = element_blank()) +
  xlab("Absolute t-statistic") +
  ylab("SE of Hist Freq (%)") +
  coord_cartesian(xlim = c(0,8))

ggsave(file.path(results_dir, "sim-se.pdf"), device = cairo_pdf, scale = 1.4, height = 2.5, width = 5)


## Plot Counts ====

plt <- plotme %>%
  mutate(se = se * N) %>%
  ggplot(aes(x = mids, y = se, color = method, linetype = method)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("Simulation" = MATBLUE, "Bernoulli" = MATRED, "Bootstrap" = MATYELLOW)) +
  theme(legend.position = c(0.8, 0.8), text = element_text(size = 22), legend.title = element_blank()) +
  xlab("Absolute t-statistic") +
  ylab("SE of Count") +
  coord_cartesian(xlim = c(0,8))

ggsave(file.path(results_dir, "sim-se-count.pdf"), device = cairo_pdf, scale = 1.4, height = 2.5, width = 5)



# Time taken -------------------------------------------
end_time <- Sys.time()
time_diff <- difftime(end_time, start_time, units = "mins")
print(time_diff)
