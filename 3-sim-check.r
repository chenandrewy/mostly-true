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
debugSource('0-functions.r')

# load data
load('../data/emp_data.Rdata')

# balanced data
retwide = pivot_wider(
  emp_retbal
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


# SIMULATION VERIFICATION ====



## SETTINGS ====

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
Emualt_large = 0.12

# bias parameters
tgood_real = 2.6
smarg_real = 0.5

tgood_cray = 5.0
smarg_cray = 0.25

# number of sims, 
nsim = 200
nsim_large = 200
nulldf = 100
tbarlist = seq(0,12,0.2)

# tuning parameters
tgoodhat = 2.6
pnullhat = 0
shapehat = 0.5



## SIMULATE ====

# realistic bias, small fdr
# show off why fdrhat_exp is relevant
manysum_real_small = average_many_sims(
  ematemp, N, T_
  , pnull = pnull_small, Emualt = Emualt_small,
  tgood = tgood_real, smarg = smarg_real
  , nsim = nsim
  , tgoodhat, pnullhat, shapehat, nulldf, tbarlist
)

# realistic bias, large fdr
manysum_real_large = average_many_sims(
  ematemp, N, T_
  , pnull = pnull_large, Emualt = Emualt_large
  , tgood = tgood_real, smarg = smarg_real
  , nsim = nsim_large
  , tgoodhat, pnullhat, shapehat, nulldf, tbarlist
)


# cray bias, small fdr
# show off why fdrhat_exp is relevant
manysum_cray_small = average_many_sims(
  ematemp, N, T_
  , pnull = pnull_small, Emualt = Emualt_small,
  tgood = tgood_cray, smarg = smarg_cray
  , nsim = nsim
  , tgoodhat, pnullhat, shapehat, nulldf, tbarlist
)

# cray bias, large fdr
# show off how fdr_hat mix works even when most discoveries are false
manysum_cray_large = average_many_sims(
  ematemp, N, T_
  , pnull = pnull_large, Emualt = Emualt_large
  , tgood = tgood_cray, smarg = smarg_cray
  , nsim = nsim_large
  , tgoodhat, pnullhat, shapehat, nulldf, tbarlist
)


## PLOT ====
pdf_w = 10/2.5
pdf_h = 8/2.5

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
        , levels = c("fdrhat_exp", "fdrhat_mix", "fdr_actual")
        , labels = c('Exponential','Conservative', 'Actual FDR')
      )
    )
  
      
  ggplot(
    manysum_plot
    , aes(x=tbar, y=fdr, group = type)
  ) +
    geom_line(aes(linetype=type, color = type),size = 2.5) +
    coord_cartesian(ylim = c(0,100)) +
    scale_linetype_manual(values = c('solid','31','11')) +
    scale_color_manual(values=c(niceblue, nicegreen, 'gray')) +
    theme_economist_white(gray_bg = FALSE) + 
    theme(
      axis.title = element_text(size = 12)
      , axis.text = element_text(size = 10)      
      , legend.background = element_rect(colour = 'black', fill = 'white', linetype = 'solid')      
      , legend.title = element_blank()
      , legend.text = element_text(size = 10)
      , legend.spacing.y = unit(0.000001, 'cm')      
      , legend.position = c(75,80)/100
      , legend.key.width = unit(1.25,'cm')      
    )  +
    labs(
      x = TeX('\\bar{t}')
      , y = TeX('FDR Upper Bound for $|t_i|>\\bar{t}$ (\\%)')
    )  
} # end custom_plot




# realistic bias
p1 = custom_plot(manysum_real_small) +
  coord_cartesian(
    xlim = c(0,6)
  )
ggsave(filename = '../results/sim_ex_smallfdr.pdf', width = pdf_w, height = pdf_h)

p2 = custom_plot(manysum_real_large) +
  theme(
    legend.position = 'none'
  ) +
  coord_cartesian(
    xlim = c(0,10)
  )
ggsave(filename = '../results/sim_ex_largefdr.pdf', width = pdf_w, height = pdf_h)

# cray bias
p3 = custom_plot(manysum_cray_small) +
  coord_cartesian(
    xlim = c(0,6)
  )
ggsave(filename = '../results/sim_ex_mis_smallfdr.pdf', width = pdf_w, height = pdf_h)

p4 = custom_plot(manysum_cray_large) +
  theme(
    legend.position = 'none'
  )  +
  coord_cartesian(
    xlim = c(0,10)
  )
ggsave(filename = '../results/sim_ex_mis_largefdr.pdf', width = pdf_w, height = pdf_h)


grid.arrange(p1,p2,p3,p4)


## CHECKS ====

emat = ematemp

Nperblock = dim(emat)[2]
nblock = floor(N/Nperblock)+1
Nemat = dim(emat)[2]


# first draw a ton of residuals
imonth  = sample(1:dim(emat)[1], nblock*T_ , replace = T) # draw list of months
esim = emat[imonth, ] %>% as.matrix
esim = matrix(esim, nrow = T_ )

c = cor(esim[, 1:200])

hist(c,40)

hist(cor(emat),40)











