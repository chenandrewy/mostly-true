# 2022 05 10: simulation for dm data

# it's quite fast, the bootstrap is pre-computed

# Setup -----------------------------------------------------------------------
rm(list = ls())

# Set working directory to unbreakable-bh folder
if (basename(getwd()) != "unbreakable-bh") {
  # Try to find unbreakable-bh directory
  if (dir.exists("unbreakable-bh")) {
    setwd("unbreakable-bh")
  } else if (dir.exists("../unbreakable-bh")) {
    setwd("../unbreakable-bh")  
  } else {
    stop("Please run this script from the unbreakable-bh directory or its parent directory.")
  }
}

source("0-functions.r")
load("../data/emp_data.Rdata")
load('../data/bootnull.Rdata')

## User Settings ====

# select data
ret <- clz_ret # yz_ret or clz_ret
retsum <- clz_sum

# parameters
pF_list <- c(0.01, seq(0.05, 0.95, 0.05), 0.99)
mutrue_list <- c(0.25, 0.5, 0.75)

parlist <- expand_grid(pF = pF_list, mutrue = mutrue_list)

# seed
set.seed(1120)

# fdr estimation
h_disc <- 2 # cutoff for a discovery
fdrhat_numer <- 0.05 # plug in this for Pr(F|disc)

# bootstrap residuals
# these residuals get reused for different expected return parameter values

bootcor <- bootnull[!is.na(corsamp), .(corsamp, booti)]
bootresid <- bootnull %>% transmute(booti, ebar = mean, vol, nmonth)

# Simulate various parameters -------------------------------------
# sim cross ====
fdrlist <- tibble()

library(doParallel)
library(foreach)
library(dplyr)
library(data.table)

tic = Sys.time()
ncore <- 3 # set the number of cores
cl <- makeCluster(ncore)
registerDoParallel(cl)

famstat = foreach(pari = 1:dim(parlist)[1], .combine = rbind, .packages = c("dplyr", "data.table")) %dopar% {
  print(paste0("Simulating pari = ", pari, " out of ", dim(parlist)[1]))

  # load par
  par <- parlist[pari, ]

  # simulate mu, rbar, tstat
  ntot <- nrow(bootresid)
  bootresid[, verity := runif(ntot) > par$pF] %>%
    .[, mu := verity * par$mutrue + (1 - verity) * 0] %>% 
    .[, tabs := abs((mu+ebar)/vol*sqrt(nmonth))]

  # = find fdr (true and estimated) = #
  h_storey <- 0.5

  famstat <- bootresid[, .(
    fdp = sum(!verity & tabs > h_disc) / sum(tabs > h_disc),
    dr = sum(tabs > h_disc) / .N,
    pFmax = mean(tabs < h_storey) / (2 * (pnorm(h_storey) - 0.5))
  ), by = .(booti)] %>% 
  mutate(fdrmax = fdrhat_numer/dr, fdrmax2=fdrmax*pFmax) %>% 
  # add parameter values
  cbind(par)

  return(famstat)
} # for pari

stopCluster(cl)
toc = Sys.time()
print(paste0("Elapsed minutes: ", difftime(toc, tic, units = "mins")))


# FIG: COR DIST --------------------------------------------

# find empirical correlations
nsignalemp <- 1000 # otherwise it takes forever
sampid <- sample(unique(ret$signalname), nsignalemp)
rsampemp <- ret[signalname %in% sampid, ] %>%
  dcast(date ~ signalname, value.var = "ret") %>%
  select(-date) %>%
  as.matrix()
corsampemp <- cor(rsampemp, use = "pairwise.complete.obs")
corsampemp <- corsampemp[lower.tri(corsampemp)] %>% as.vector()

cdat <- rbind(
  data.frame(c = bootcor$corsamp, group = "sim", color = NICEBLUE),
  data.frame(c = corsampemp, group = "emp", color = "gray")
)

edge <- seq(-1, 1, 0.05)
plotme <- cdat %>%
  group_by(group) %>%
  summarise(
    cmid = hist(c, edge, plot = F)$mids,
    density = hist(c, edge, plot = F)$density
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("sim", "emp"),
      labels = c("Simulated", "Empirical")
    )
  )

plt <- ggplot(
  plotme, aes(x = cmid, y = density, group = group)
) +
  geom_line(
    aes(linetype = group, color = group),
    size = 2
  ) +
  theme_economist_white(gray_bg = F) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.1, "cm"),
    legend.position = c(80, 80) / 100,
    legend.key.width = unit(1, "cm"),
    legend.spacing.y = unit(0.000001, "cm"),
    legend.background = element_rect(colour = "black", fill = "white")
  ) +
  labs(
    x = "Pairwise Correlation",
    y = "Density"
  ) +
  scale_color_manual(
    values = c(NICEBLUE, "gray")
  ) +
  scale_linetype_manual(values = c("solid", "31"))
ggsave(filename = "../results/cor-theory-free.pdf", width = 5, height = 4)


# Plot for Paper --------------------------------------------------------

plotme <- famstat %>%
  group_by(pF, mutrue) %>%
  mutate(pF = 100 * pF) %>%
  summarize(
    fdr_act = 100 * mean(fdp), fdr_max = 100 * mean(fdrmax),
    fdr_max2 = 100 * mean(fdrmax2)
  ) %>%
  pivot_longer(cols = starts_with("fdr")) %>%
  mutate(name = factor(
    name,
    levels = c("fdr_act", "fdr_max2", "fdr_max"),
    labels = c("Actual", "Visual Bound", "Easy Bound")
  ))

# reordering legend (but preserving plotting order)
# https://stackoverflow.com/questions/50420205/how-to-reorder-overlaying-order-of-geoms-but-keep-legend-order-intact-in-ggplot
plotme$id2 <- factor(plotme$name, levels = c("Easy Bound", "Visual Bound", "Actual"))

# loop over mutrue values
mutrue_list <- unique(plotme$mutrue)
for (mutruei in mutrue_list) {
  plt <- plotme %>%
    filter(mutrue == mutruei) %>%
    ggplot(aes(x = pF, y = value, group = name)) +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_hline(yintercept = 100, color = "gray50") +
    # plot FDR and bounds
    geom_line(aes(linetype = name, color = name), size = 1.2) +
    scale_color_manual(
      values = c(
        "Easy Bound" = "gray60",
        "Visual Bound" = MATBLUE, "Actual" = MATRED
      ),
      breaks = levels(plotme$id2),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c(
        "Easy Bound" = "dotdash",
        "Visual Bound" = "31", "Actual" = "solid"
      ),
      breaks = levels(plotme$id2),
      name = NULL
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype"),
      legend.key.size = unit(0.1, "cm"),
      legend.position = c(20, 70) / 100,
      legend.key.width = unit(1, "cm"),
      legend.spacing.y = unit(0.001, "cm"),
      legend.background = element_rect(colour = "white", fill = "white"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = TeX("Proportion Null Overall $Pr(null_i)$ (%)"),
      y = TeX("$FDR_{|t|>2}$ (%)")
    ) +
    coord_cartesian(ylim = c(0, 100))

  ggsave(
    paste0("../results/sim-dm-visual-", mutruei, ".pdf"), plt,
    width = 5, height = 2.5,
    scale = 1, device = cairo_pdf
  )
} # for mutruei


# Plot for Short Slides --------------------------------------------------------
# Omit Visual Bound for 25 min or less

plotme <- famstat %>%
  group_by(pF, mutrue) %>%
  mutate(pF = 100 * pF) %>%
  summarize(
    fdr_act = 100 * mean(fdp), fdr_max = 100 * mean(fdrmax),
    fdr_max2 = 100 * mean(fdrmax2)
  ) %>%
  pivot_longer(cols = starts_with("fdr")) %>%
  mutate(name = factor(
    name,
    levels = c("fdr_act", "fdr_max2", "fdr_max"),
    labels = c("Actual", "Visual Bound", "Easy Bound")
  )) %>% 
  filter(name != "Visual Bound")

# reordering legend (but preserving plotting order)
# https://stackoverflow.com/questions/50420205/how-to-reorder-overlaying-order-of-geoms-but-keep-legend-order-intact-in-ggplot
plotme$id2 <- factor(plotme$name, levels = c("Easy Bound", "Visual Bound", "Actual"))

# loop over mutrue values
mutrue_list <- unique(plotme$mutrue)
for (mutruei in mutrue_list) {
  plt <- plotme %>%
    filter(mutrue == mutruei) %>%
    ggplot(aes(x = pF, y = value, group = name)) +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_hline(yintercept = 100, color = "gray50") +
    # plot FDR and bounds
    geom_line(aes(linetype = name, color = name), size = 1.2) +
    scale_color_manual(
      values = c(
        "Easy Bound" = "gray60",
        "Visual Bound" = MATBLUE, "Actual" = MATRED
      ),
      breaks = levels(plotme$id2),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c(
        "Easy Bound" = "dotdash",
        "Visual Bound" = "31", "Actual" = "solid"
      ),
      breaks = levels(plotme$id2),
      name = NULL
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype"),
      legend.key.size = unit(0.1, "cm"),
      legend.position = c(20, 70) / 100,
      legend.key.width = unit(1, "cm"),
      legend.spacing.y = unit(0.001, "cm"),
      legend.background = element_rect(colour = "white", fill = "white"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = TeX("Proportion Null Overall $Pr(null_i)$ (%)"),
      y = TeX("$FDR_{|t|>2}$ (%)")
    ) +
    coord_cartesian(ylim = c(0, 100))

  ggsave(
    paste0("../results/sim-dm-ez-only-", mutruei, ".pdf"), plt,
    width = 5, height = 2.5,
    scale = 1, device = cairo_pdf
  )
} # for mutruei
