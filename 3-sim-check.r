# 2021 08 Andrew
# Simulation verification

# ENVIRONMENT ====
rm(list=ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(latex2exp)
source('0-functions.r')

# load data
load('../data/emp_data.Rdata')

# balanced data
retwide = pivot_wider(
  emp_ret
  , c(signalname,date,ret), names_from = signalname, values_from = ret
) %>% 
  select(-date) %>% 
  filter(complete.cases(.))
retmat = as.matrix(retwide) 

# generate residuals
rbar = colMeans(retmat) %>% as.matrix() %>% t()
l = 1+integer(nrow(retmat)) %>% as.matrix 
ematemp = retmat - l %*% rbar
Nemp = dim(ematemp)[2]
Temp = dim(ematemp)[1]


## FUNCTIONS USED THROUGHOUT SECTIONS
# tdat is a dataframe with columns (t, other stuff)
tdat_to_tdatpub = function(tdat, tbad=1.96, tgood=2.6, smarg=0.5, sbar=1){
  tdat$u = runif(dim(tdat)[1]) 
  tdat = tdat %>% 
    mutate(
      pub = F
      , pub = case_when(
        t > tbad & t <= tgood & u < smarg ~ T
        , t > tgood & u < sbar ~ T
      )
    ) %>% 
    filter(pub) %>% 
    select(-pub)
} # end function 


estatsim = function(N, T_){
  # emat is a global
  
  Nperblock = dim(emat)[2]
  nblock = floor(N/Nperblock)+1
  
  # first draw a ton of residuals
  imonth  = sample(1:Temat, nblock*T_ , replace = T) # draw list of months
  esim = emat[imonth, ] %>% as.matrix
  
  # average within blocks 
  tic = Sys.time()
  ebar = numeric(nblock*Nemat)
  evol = ebar
  for (blocki in 1:nblock){
    tstart = (blocki-1)*T_ + 1
    tend = tstart + T_ - 1
    
    Nstart = (blocki-1)*Nemat + 1
    Nend = Nstart + Nemat - 1
    ebar[Nstart:Nend] = colMeans(esim[tstart:tend, ])
    evol[Nstart:Nend] = sqrt(colMeans(esim[tstart:tend, ]^2))
  }
  
  # clean and output
  estat = data.frame(
    bar = ebar[1:N]
    , vol = evol[1:N]
  )
} # end ebarsim


# function for generating observables and nulls
estat_to_tdat = function(N,pnull,Emualt,estat){
  
  # adjust for true
  nnull = sum(runif(N) < pnull)
  null = logical(N)
  null[1:nnull] = T
  mu = matrix(0, N)
  mu[!null] = Emualt
  
  # observables
  t = (mu + estat$bar)/estat$vol*sqrt(T_)
  
  # pack and output
  data.frame(t = abs(t), null = null, mu = mu, traw = t)
  
} # end e_to_t


average_many_sims = function(
  N, T_, pnull, Emualt, tgood, smarg
  , nsim 
  , tgoodhat, pnullhat, shapehat
){
  simmany = data.frame()
  for (simi in 1:nsim){
    print(paste0('simulation number ', simi))
    
    estat = estatsim(N, T_)
    tdat = estat_to_tdat(N, pnull, Emualt, estat)
    tdatpub = tdat_to_tdatpub(tdat, tgood = tgood, smarg = smarg )
    
    # actual fdr
    fdr_actual = estimate_fdr(
      tdat$t, tbarlist = tbarlist, null = tdat$null
    ) %>% 
      transmute(tbar, dr_actual = dr, fdr_actual)
    
    # fdrhat using exp
    est_bias = estimate_exponential(tdatpub$t, tgoodhat)
    fdr_exp = estimate_fdr(
      tdatpub$t, tbarlist = tbarlist, C = est_bias$C, nulldf = nulldf
    ) %>% 
      transmute(tbar, fdrhat_exp = fdrhat)
    
    # fdrhat using mix
    est_bias = estimate_mixture(tdatpub$t, tgood = tgoodhat, pnull = pnullhat, shape = shapehat, 1)
    fdr_mix = estimate_fdr(
      tdatpub$t, tbarlist = tbarlist, C = est_bias$C, nulldf = nulldf
    ) %>% 
      transmute(tbar, fdrhat_mix = fdrhat)
    
    est_all = fdr_actual %>% 
      left_join(fdr_exp, by = 'tbar') %>% 
      left_join(fdr_mix, by = 'tbar')
    
    est_all$simi = simi
    
    simmany = rbind(simmany, est_all)  
    
  } # end for simi
  
  # average across simulations
  manysum = simmany %>% 
    group_by(tbar) %>% 
    summarise_at(
      vars(2:(dim(simmany)[2])-2), mean, na.rm=F
    ) %>% 
    mutate(
      Ndisc = dr_actual*N
    )
  
} # end average_many_sims


# SIMULATE ILLUSTRATIONS ====



## BASELINE SETTINGS ====

# residuals
emat = ematemp 
Nemat = dim(emat)[2]
Temat = dim(emat)[1]

# modelS for simulation
N = 1e4
T_ = 200

# fdr parameters
pnull_small  = 0.5
Emualt_small = 0.5

pnull_large = 0.99
Emualt_large = 0.25

# bias parameters
tgood_real = 2.6
smarg_real = 0.5

tgood_cray = 5.0
smarg_cray = 0.25

# number of sims, fdr estimates, other
nsim = 100
nsim_large = 100
nulldf = 100
tbarlist = seq(0,6,0.2)
tgoodhat = 2.6
pnullhat = 0.9
shapehat = 2

# plot settings
custom_plot = function(manysum){
  manysum_plot = manysum %>% 
    select(tbar, starts_with('fdr')) %>% 
    mutate_at(
      .vars = vars(c(-tbar))
      , .funs = ~round(.*100,1)
    ) %>% 
    pivot_longer(-tbar, names_to = 'type', values_to = 'fdr') %>% 
    mutate(
      type = factor(
        type
        , levels = c("fdr_actual","fdrhat_mix","fdrhat_exp")
        , labels = c('Actual','Est: Conservative','Est: Exponential')
      )
    )
  
  ggplot(
    manysum_plot
    , aes(x=tbar, y=fdr, group = type)
  ) +
    geom_line(aes(linetype=type, color = type),size = 2.5) +
    coord_cartesian(ylim = c(0,100)) +
    scale_linetype_manual(values=c("dashed", "solid", "dotted")) +
    scale_color_manual(values=c("#619CFF","#00BA38", "#F8766D")) +
    theme_economist_white(gray_bg = FALSE) + 
    theme(
      axis.title = element_text(size = 40)
      , axis.text = element_text(size = 30)      
      , legend.title = element_blank()
      , legend.text = element_text(size = 30)
      , legend.background = element_rect(colour = 'black', fill = 'white', linetype = 'solid')
    )  +
    labs(
      x = TeX('\\bar{t}')
      , y = TeX('FDR for $|t_i|>\\bar{t}$ (\\%)')
    )  +
    theme(
      legend.position = c(70,80)/100
      , legend.key.width = unit(3,'cm')
    )
  
} # end custom_plot

## SIMULATE ====



# realistic bias, small fdr
# show off why fdrhat_exp is relevant
manysum_real_small = average_many_sims(
  N, T_
  , pnull = pnull_small, Emualt = Emualt_small,
  tgood = tgood_real, smarg = smarg_real
  , nsim = 10
  , tgoodhat, pnullhat, shapehat
)

# realistic bias, large fdr
# show off how fdr_hat mix works even when most discoveries are false
manysum_real_large = average_many_sims(
  N, T_
  , pnull = pnull_large, Emualt = Emualt_large
  , tgood = tgood_real, smarg = smarg_real
  , nsim = nsim_large
  , tgoodhat, pnullhat, shapehat
)


# cray bias, small fdr
# show off why fdrhat_exp is relevant
manysum_cray_small = average_many_sims(
  N, T_
  , pnull = pnull_small, Emualt = Emualt_small,
  tgood = tgood_cray, smarg = smarg_cray
  , nsim = 10
  , tgoodhat, pnullhat, shapehat
)

# cray bias, large fdr
# show off how fdr_hat mix works even when most discoveries are false
manysum_cray_large = average_many_sims(
  N, T_
  , pnull = pnull_large, Emualt = Emualt_large
  , tgood = tgood_cray, smarg = smarg_cray
  , nsim = nsim_large
  , tgoodhat, pnullhat, shapehat
)


## PLOT ====

# realistic bias
custom_plot(manysum_real_small)
ggsave(filename = '../results/sim_ex_smallfdr.pdf', width = 10, height = 8)

custom_plot(manysum_real_large) +
  theme(
    legend.position = 'none'
  )
ggsave(filename = '../results/sim_ex_largefdr.pdf', width = 10, height = 8)

# cray bias
custom_plot(manysum_cray_small)
ggsave(filename = '../results/sim_ex_mis_smallfdr.pdf', width = 10, height = 8)

custom_plot(manysum_cray_large) +
  theme(
    legend.position = 'none'
  )

ggsave(filename = '../results/sim_ex_mis_largefdr.pdf', width = 10, height = 8)


# A CLOSER LOOK ====

# requires sim from previous section

##  FDRHATS ARE DIFFERENT ====
temp = manysum_real_small %>% 
  transmute(
    tbar, fdrhat_exp, group = 'real'
  ) %>% 
  rbind(
    manysum_cray_small %>% 
      transmute(
        tbar, fdrhat_exp, group = 'cray'
      ) 
  )

# visual table check
temp %>% 
  pivot_wider(
    names_prefix = 'fdrhat_', names_from = group, values_from = fdrhat_exp
  )

ggplot(temp, aes(x=tbar,y=fdrhat_exp)) +
  geom_line(aes(color = group))


ggsave(filename = '../results/detail_fdrhat.pdf', width = 10, height = 8)

## XX DETAILED SIM ====

# notes:
# offsetting effects: cray => C is smaller => fdrhat smaller
# but cray => t-obs pushed to right => fdrhat is larger

# choose spec
nsim = 1e5
pnull_det = pnull_small
Emualt_det = Emualt_small

# simulate all
estat = estatsim(N, T_)
tdat = estat_to_tdat(N, pnull=pnull_det,  Emualt=Emualt_det, estat)
t_all = tdat$t

# realistic bias
tdatpub = tdat_to_tdatpub(tdat, tgood = tgood_real, smarg = smarg_real )
t_real = tdatpub$t
fit_real = estimate_mixture(t_real, tgood = tgoodhat, pnull = pnullhat, shape = shapehat, 1)
tfit_real = simmix(nsim,pnullhat,shape=shapehat,fit_real$scalehat,1)$t

fdr_real = estimate_fdr(t_real, tbarlist, fit_real$C, nulldf=nulldf)

# crazy bias
tdatpub = tdat_to_tdatpub(tdat, tgood = tgood_cray, smarg = smarg_cray )
t_cray = tdatpub$t
fit_cray = estimate_mixture(t_cray, tgood = tgoodhat, pnull = pnullhat, shape = shapehat, 1)
fdr_cray = estimate_fdr(t_cray, tbarlist, fit_real$C, nulldf=nulldf)
tfit_cray = simmix(nsim,pnullhat,shape=shapehat,fit_cray$scalehat,1)$t

# plot ====

# settings
edge = seq(0,8,0.2)
density_max = 1.5

temp = data.frame(
  t = t_real, group = 'obs'
) %>% 
  rbind(
    data.frame(
      t = tfit_real, group = 'fit'
    )
  ) %>%
  rbind(
    data.frame(
      t = t_all, group = 'all'
    ) 
  ) %>% 
  filter(t>edge[1],t<edge[length(edge)]) %>% 
  group_by(group) %>% 
  summarise(
    density = hist(t,edge)$density / sum(t>tgoodhat)*n()
    , t = hist(t,edge)$mids
  ) 

p_real = ggplot(
  temp %>% filter(group %in% c('obs','all'))
  , aes(x=t, y=density, fill  = group)
) +
  geom_bar(stat = 'identity', alpha = 0.6, position = 'dodge') +
  geom_line(
    data=temp %>% filter(group=='fit') 
  ) +
  coord_cartesian(ylim =c(0,density_max))
  
temp = data.frame(
  t = t_cray, group = 'obs'
) %>% 
  rbind(
    data.frame(
      t = tfit_cray, group = 'fit'
    )
  ) %>% 
  rbind(
    data.frame(
      t = t_all, group = 'all'
    ) 
  ) %>%   
  filter(t>edge[1],t<edge[length(edge)]) %>% 
  group_by(group) %>% 
  summarise(
    density = hist(t,edge)$density / sum(t>tgoodhat)*n()
    , t = hist(t,edge)$mids
  ) 

p_cray = ggplot(
  temp %>% filter(group %in% c('obs','all'))
  , aes(x=t, y=density, fill  = group)
) +
  geom_bar(stat = 'identity', alpha = 0.6, position = 'dodge') +
  geom_line(
    data=temp %>% filter(group=='fit') 
  ) +
  coord_cartesian(ylim =c(0,density_max))


grid.arrange(p_real, p_cray)


# save to disk
ggsave(p_real, filename = '../results/simdet_fitreal.pdf', width = 10, height = 4)

ggsave(p_cray, filename = '../results/simdet_fitcray.pdf', width = 10, height = 4)
