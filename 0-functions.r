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