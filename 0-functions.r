# 2021 08 Andrew
# frequently used functions for bh with pub bias

# PATHS AND LIBRARIES ====

dir.create('../data/', showWarnings = F)
dir.create('../results/', showWarnings = F)
dir.create('../results/sim-theory-free/', showWarnings = F)
dir.create('../results/sim-extrap-pub/', showWarnings = F)

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(latex2exp)
library(foreach)
library(doParallel)

# Global settings -----------------------------


# STATS ====

# truncated gamma moments
mom_trunc_gamma = function(shape,scale,tmin,ord=1){
  
  tempd = function(t){dgamma(t,shape,scale = scale)}
  
  Fh = integrate(tempd, 0, tmin)$value
  
  intme = function(t){
    t^ord*tempd(t)/(1-Fh)
  }
  temp = integrate(intme,tmin,Inf)
  return = temp$value
}

# fit truncated gamma
est_trunc_gamma = function(tabs, tgood, shape, ord=1){
  # tabs = rgamma(n = 1e3, shape = 0.5, scale = 1) 
  # shape = 0.5
  # ord = 1
  # tgood = 2
  
  minme = function(scale){
    (
      mom_trunc_gamma(shape,scale,tgood,ord) - mean(tabs[tabs>tgood]^ord)
    )^2
  }
  
  temp = optimize(minme, c(0.1,6)/shape)
  est = tibble(shape = shape, scale = temp$minimum, obj = temp$objective)
  
  return = est
  
} # est_trunc_gamma

# TOOLS ====

detach_all = function(){
  invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
}

histcomp = function(
  dat1,dat2
  ,label1 = '1',label2 = '2'
  ,edge=seq(-3,15,0.5),tgood=-Inf
){
  t = dat1
  t = t[t>min(edge) & t<max(edge)]
  hemp = hist(t,edge)
  
  t = dat2
  t = t[t>min(edge) & t<max(edge)]
  hsim = hist(t,edge)
  
  plotme = rbind(
    data.frame(
      t = hemp$mids
      , f = hemp$density / sum(dat1>tgood) * length(dat1)
      , group = label1
    )
    , data.frame(
      t = hsim$mids
      , f = hsim$density / sum(dat2>tgood) * length(dat2)
      , group = label2
    )
  )
  
  p1 = ggplot(plotme, aes(x=t, y=f, fill=group)) +
    geom_bar(stat='identity', position='identity',alpha=0.6, show.legend = F) 
  
  p2 = ggplot(
    plotme %>% filter(t>tgood)
    , aes(x=t, y=f, fill=group)
  ) +
    geom_bar(stat='identity', position='identity',alpha=0.6, show.legend = F) 
  
  
  grid.arrange(p1, p2, nrow=1)  
} # end function

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

## set theme ====
theme_set(
  theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype")
    )
)

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

# function for turning bootstrap results into a histogram
histogram_by_group = function(bootact, edge, varname = 'tstat', group = 'booti') {

  bootact = bootact %>% select(all_of(c(group, varname)))
  colnames(bootact) = c('booti', 'xvar')

  boothist = bootact %>%
    # histogram counts within each bootstrap
    mutate(bin = cut(xvar, breaks = edge, include.lowest = TRUE)) %>%
    group_by(booti, bin) %>%
    summarize(count = n(), .groups = 'drop')  %>% 
    # normalize by total signals
    left_join(bootact %>% group_by(booti) %>% summarize(ntotal=n()), by = "booti") %>% 
    # find bin midpoints 
    mutate(bin = str_remove_all(bin, "[\\(\\)\\[\\]]")
        , left = str_split(bin, ",")
        , left = sapply(left, function(x) as.numeric(x[1]))
        , right = str_split(bin, ",")
        , right = sapply(right, function(x) as.numeric(x[2]))
        , mids = (left + right) / 2
    ) %>% 
    setDT()
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
