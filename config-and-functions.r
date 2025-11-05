# ABOUTME: Provides project configuration utilities alongside shared helper functions.
# Inputs:
#   - Sourced by scripts (e.g., 01-prep-data.r) after calling here::i_am(<script>).
#   - Requires here, data.table, tidyverse, ggplot2, ggthemes, gridExtra, latex2exp, foreach, doParallel, extrafont.
# Outputs:
#   - Ensures project directories exist and exposes reusable helper utilities.
# How to run:
#   Rscript -e 'source(here::here("config-and-functions.r"))'
#   """
#   library(here)
#   source(here::here("config-and-functions.r"))
#   paths <- project_paths()
#   """
# Folder layout:
#   """
#   project-root/
#   ├─ main.r
#   ├─ config-and-functions.r
#   ├─ data/              # created if missing; RData/CSV outputs from prep scripts
#   ├─ results/           # created if missing; PDFs/TEX/plots from analysis scripts
#   ├─ temp/              # scratch space for intermediates (optional; created on demand)
#   ├─ 0X-*.r             # numbered analysis scripts
#   └─ …
#   """

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

project_paths <- function() {
  results_root <- ensure_dir(here::here("results"))
  list(
    data = ensure_dir(here::here("data")),
    temp = ensure_dir(here::here("temp")),
    results = results_root
  )
}

library(here)

paths <- project_paths()
DATA_DIR <- paths$data
RESULTS_DIR <- paths$results

library(data.table)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(foreach)
library(doParallel)
library(extrafont)

# Global settings -----------------------------

# AESTHETICS ====

library(latex2exp)
library(extrafont)

MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)
MATPURPLE = rgb(0.4940, 0.1840, 0.5560)
MATGREEN = rgb(0.4660, 0.6740, 0.1880)

NICEBLUE = "#619CFF"
NICEGREEN = "#00BA38"
NICERED = "#F8766D"

theme_set(
  theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype")
      # text = element_text(family = "ArialMT")
    )
)

# histogram data prep function -----------------------------------------------

# creates data for comparing cdf F1 to cdf F2 in a plot
# automatically adjusts for different x-binning
make_dist_dat = function(F1, edge1, N1, F2, edge2, N2
  , x_match = c(-Inf,Inf), showplot = F){
  
  # adjust for different x-binning
  if (!is.null(x_match)){
    rescale_fac = diff(F1(x_match))/diff(F2(x_match)) * diff(edge1)[1] /diff(edge2)[1]
  } else {
    rescale_fac = 1
  }
  
  # make histogram counts, with normalization adjustments
  dat = tibble(
    edge = edge1, F = N1*F1(edge1), group = 1
  ) %>% 
    rbind(
      tibble(
        edge = edge2, F = N2*F2(edge2)*rescale_fac, group = 2
      )
    ) %>% 
    # take first differences, find midpoints
    group_by(group) %>% 
    mutate(
      F = F
      , dF = F - lag(F)
      , mids = 0.5*(edge + lag(edge))
    ) %>% 
    filter(!is.na(dF)) %>% 
    setDT()
  
  if (showplot) {
    dat %>% 
      ggplot(aes(x=edge, y=dF)) +
      geom_line(aes(color = group))
  }
  
  return(dat)
  
} # make_dist_dat

# define bootstrap function -------------------------------
bootstrap_flex <- function(ret, nboot, coli = "signalname", colt = "date", colr = "ret", min_obs_pct = 50, demean = TRUE, ncore = 1, output_cor = FALSE) {
    # min_obs_pct = 80 starts in 1986, 50 starts in 1972
    
    # convert to wide matrix
    ret = ret %>% select(all_of(c(coli, colt, colr)))
    colnames(ret) = c('signalname','date','ret')
    rmat <- dcast(ret, signalname ~ date, value.var = "ret") %>%
        as.matrix()
    row.names(rmat) <- rmat[, "signalname"]
    rmat <- rmat[, -which(colnames(rmat) == "signalname")]

    # define sample period ------------------------------------
    obs_by_yearm <- colSums(!is.na(rmat), na.rm = TRUE) / nrow(rmat) * 100
    col_ok <- which(obs_by_yearm > min_obs_pct)
    col_ok <- min(col_ok):max(col_ok)
    rmat <- rmat[, col_ok]

    # de-mean -----------------------------------------------
    if (demean){rmat <- rmat - rowMeans(rmat, na.rm=TRUE)}

    # inner function for a single boot
    bootstrap_once <- function() {
        # sample dates (make this flexible)
        tempdate <- sample(colnames(rmat), ncol(rmat), replace = TRUE) %>% sort()

        # make bootstrapped panel
        tempmat <- rmat[, tempdate]

        # summary stats
        tempmean <- rowMeans(tempmat, na.rm=T)
        tempsd <- sqrt(rowMeans(tempmat^2, na.rm=T)-tempmean^2)
        tempnmonth <- rowSums(!is.na(tempmat))
        temptstat <- tempmean / tempsd * sqrt(tempnmonth)

        # add correlations (in a slightly janky way)
        tempcor <- rep(NA, length(temptstat))
        if (output_cor) {
          # if we sample 200 signals, we end up with 200*199/2 = 19900 correlations
          tempid = sample(1:nrow(tempmat), 200, replace = FALSE)
          tempcmat = cor(tempmat[tempid,] %>% t(), use='pairwise.complete.obs')
          tempcorvec <- tempcmat[lower.tri(tempcmat)]
          tempcor[1:length(tempcorvec)] <- tempcorvec
        }
        return(data.table(mean = tempmean, vol = tempsd, nmonth = tempnmonth,  corsamp = tempcor))
    }

    # bootstrap nboot times (in parallel or not)
    if (ncore > 1) {
        print(paste0('bootstrapping with ncore = ', ncore))

        # set up cluster
        cl <- makeCluster(ncore)
        registerDoParallel(cl)
        on.exit(stopCluster(cl))

        bootdat <- foreach(booti = 1:nboot, .combine = rbind, .packages = c("data.table", "dplyr")) %dopar% {
            print(paste0("bootstrapping ", booti, " of ", nboot))
            bootstrap_once() %>% mutate(booti = booti)
        }
    } else {
        bootdat <- foreach(booti = 1:nboot, .combine = rbind, .packages = c("data.table", "dplyr")) %do% {
            print(paste0("bootstrapping ", booti, " of ", nboot))
            bootstrap_once() %>% mutate(booti = booti)
        }
    }

    return(bootdat)
} # end bootstrap_flex