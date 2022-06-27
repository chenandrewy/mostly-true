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

NICEBLUE = "#619CFF"
NICEGREEN = "#00BA38"
NICERED = "#F8766D"

