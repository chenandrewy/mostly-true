# ABOUTME: Builds diagnostic figures contrasting empirical distribution with null models.
# ABOUTME: Creates dm-viz-storey-err.pdf and dm-viz-ez-err.pdf with bootstrap error bars.
# Inputs:
#   - config-and-functions.r (project utilities)
#   - data/emp_data.Rdata (summaries for CLZ signals)
#   - data/bootact.Rdata (bootstrap samples saved by 91-run-bootstraps.r)
# Outputs:
#   - results/dm-viz-storey-err.pdf
#   - results/dm-viz-ez-err.pdf
# How to run:
#   Rscript 03-visual-bounds.r
#   Rscript 03-visual-bounds.r --vanilla

# Setup -------------------------------------

rm(list = ls())

library(here)
here::i_am("03-visual-bounds.r")

source(here("config-and-functions.r"))

paths <- project_paths()
data_dir <- paths$data
results_dir <- paths$results

load(file.path(data_dir, "emp_data.Rdata"))
load(file.path(data_dir, "bootact.Rdata"))

# ensure data.table format
bootact <- as.data.table(bootact)

# plot settings
color_emp <- "gray50"
color_null <- MATRED
bar_alpha <- 0.65
ylimnum <- c(0, 12000)
discovery_y <- 11000
intuition_y <- 3000
yticks <- seq(0, 12000, 2000)

# add t-stats to bootact
bootact[, tstat := mean / vol * sqrt(nmonth)] %>%
  .[, tabs := abs(tstat)]

# Calculations --------------------------

# distributions
F_null <- function(tabs) { (2 * (pnorm(tabs) - 0.5)) }
F_dm <- ecdf(clz_sum$tabs)

# FDR ez (using rounding)
h_disc <- 2
Pr_disc <- mean(clz_sum$tabs > h_disc)
FDRmaxez <- 0.05 / Pr_disc
n_dm <- length(clz_sum$tabs)

# storey
tmax_small <- 0.5
Pr_small <- mean(clz_sum$tabs < tmax_small)
pFmax <- F_dm(tmax_small) / F_null(tmax_small)

# bootstrap SE
bootFDR0 <- bootact[, .(
  Pr_disc = mean(tabs > h_disc),
  Pr_small = mean(tabs < tmax_small),
  FDRmaxez = 0.05 / mean(tabs > h_disc),
  pFmax = mean(tabs < tmax_small) / F_null(tmax_small)
), by = booti] %>%
  mutate(FDRmaxviz = FDRmaxez * pFmax) %>%
  melt(id.vars = "booti", variable.name = "stat") %>%
  group_by(stat) %>%
  summarize(
    mean = mean(value), sd = sd(value),
    p01 = quantile(value, 0.01), p99 = quantile(value, 0.99),
    p05 = quantile(value, 0.05), p95 = quantile(value, 0.95), .groups = "drop"
  )

bootFDR <- tibble(
  stat = c("Pr_disc", "Pr_small", "FDRmaxez", "pFmax", "FDRmaxviz"),
  point = c(Pr_disc, Pr_small, FDRmaxez, pFmax, FDRmaxez * pFmax)
) %>%
  left_join(bootFDR0, by = "stat") %>%
  rename(boot_mean = mean, boot_sd = sd)

# make data for plotting
edge <- seq(0, 10, 0.5)
temp1 <- make_dist_dat(F_dm, edge, n_dm, F_null, edge, n_dm * pFmax, x_match = NULL) %>%
  mutate(group = factor(group, levels = c(1, 2), labels = c("emp", "null")))
temp2 <- make_dist_dat(F_dm, edge, n_dm, F_null, edge, n_dm * 1, x_match = NULL) %>%
  mutate(group = factor(group, levels = c(1, 2), labels = c("emp", "null_ez")))
plotme <- temp1 %>% rbind(temp2 %>% filter(group == "null_ez"))

## add bootstrap histogram ----------------------
boothist0 <- bootact %>%
  # histogram counts within each bootstrap
  mutate(bin = cut(tabs, breaks = edge, include.lowest = TRUE)) %>%
  group_by(booti, bin) %>%
  summarize(count = n(), .groups = "drop") %>%
  # normalize by total signals
  left_join(bootact %>% group_by(booti) %>% summarize(ntotal = n()), by = "booti") %>%
  # select definition of dF
  mutate(dF = count) %>%
  # find bin midpoints
  mutate(
    bin = str_remove_all(bin, "[\\(\\)\\[\\]]"),
    left = str_split(bin, ","),
    left = sapply(left, function(x) as.numeric(x[1])),
    right = str_split(bin, ","),
    right = sapply(right, function(x) as.numeric(x[2])),
    mids = (left + right) / 2
  )

boothist <- boothist0 %>%
  group_by(mids) %>%
  summarize(
    boot05 = quantile(dF, 0.05),
    boot95 = quantile(dF, 0.95),
    bootsd = sd(dF),
    .groups = "drop"
  ) %>%
  mutate(group = "emp") %>%
  left_join(plotme %>% select(mids, group, dF),
    by = c("group", "mids")
  )

## define error bars ====
boothist <- boothist %>%
  mutate(dF_lo = dF - 1 * bootsd, dF_hi = dF + 1 * bootsd)
# boothist = boothist  %>%
#   mutate(dF_lo = boot05, dF_hi = boot95)

plotme_err <- plotme %>%
  left_join(boothist %>% select(group, mids, starts_with("dF_")),
    by = c("group", "mids")
  )

# plot with error bars -----------------------------------------

# plot storey null
plt <- ggplot(plotme_err[group != "null_ez"], aes(x = mids, y = dF)) +
  coord_cartesian(xlim = c(0, 8), ylim = ylimnum) +
  scale_y_continuous(breaks = yticks) +
  xlab(TeX("Absolute t-statistic ($|t_i|$)")) +
  ylab("Number of Signals") +
  # bars
  geom_bar(
    stat = "identity", position = "identity", alpha = bar_alpha,
    aes(fill = group)
  ) +
  scale_fill_manual(
    values = c(color_emp, color_null)
    # , labels=c('Data', paste0('False = ', round(pFmax*100,0), '% of Data'))
    , labels = c("Data: CLZ EW", paste0("Null Component Bound")),
    name = ""
  ) +
  # error bars
  geom_errorbar(aes(ymin = dF_lo, ymax = dF_hi), width = 0.15, color = "grey30", size = 0.3) +
  # discovery line
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(
    geom = "text", x = 2.1, y = discovery_y, hjust = 0,
    label = "Discoveries ->", color = MATRED
  ) +
  theme(legend.position = c(80, 80) / 100) +
  # write out intuition
  geom_segment(aes(xend = 22 / 10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +
  # annotate with SE
  annotate(
    "text", x = 33 / 10, y = intuition_y, hjust = 0,
    label = TeX(paste0(
      "FDR $\\leq$ "
      , round(100*bootFDR$point[bootFDR$stat == "FDRmaxviz"], 1)
      , " $\\pm$ "
      , round(100*bootFDR$boot_sd[bootFDR$stat == "FDRmaxviz"], 1)
      , "\\%"
    ), output = "character"),
    parse = TRUE
  )

ggsave(file.path(results_dir, "dm-viz-storey-err.pdf"), scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# plot ez null
plt <- ggplot(plotme_err[group != "null"], aes(x = mids, y = dF)) +
  coord_cartesian(xlim = c(0, 8), ylim = ylimnum) +
  scale_y_continuous(breaks = yticks) +
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX("Absolute t-statistic ($|t_i|$)")) +
  ylab("Number of Signals") +
  # bars
  geom_bar(
    stat = "identity", position = "identity", alpha = bar_alpha,
    aes(fill = group)
  ) +
  scale_fill_manual(
    values = c(color_emp, color_null)
    # , labels=c('Data', paste0('False = ', round(pFmax*100,0), '% of Data'))
    , labels = c("Data: CLZ EW", paste0("Null Component Bound")),
    name = ""
  ) +
  # error bars
  geom_errorbar(aes(ymin = dF_lo, ymax = dF_hi), width = 0.15, color = "grey30", size = 0.3) +
  # discovery line
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(
    geom = "text", x = 2.1, y = discovery_y, hjust = 0,
    label = "Discoveries ->", color = MATRED
  ) +
  theme(legend.position = c(80, 80) / 100) +
  # write out intuition
  geom_segment(aes(xend = 22 / 10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +
  # annotate with SE
  annotate(
    "text", x = 33 / 10, y = intuition_y, hjust = 0,
    label = TeX(paste0(
      "FDR $\\leq$ "
      , round(100*bootFDR$point[bootFDR$stat == "FDRmaxez"], 1)
      , " $\\pm$ "
      , round(100*bootFDR$boot_sd[bootFDR$stat == "FDRmaxez"], 1)
      , "\\%"
    ), output = "character"),
    parse = TRUE
  )

ggsave(file.path(results_dir, "dm-viz-ez-err.pdf"), scale = 1, height = 2.5, width = 5, device = cairo_pdf)
