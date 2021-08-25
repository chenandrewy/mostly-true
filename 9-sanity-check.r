# Andrew 2021 08
# for sanity checks
# had sanity checks in a scratch script but having a clean
# sanity check script helps with sanity

# ENVIRONMENT ====
rm(list=ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(latex2exp)
library(truncdist)
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
        , labels = c('Actual','Mix','Exp')
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
      axis.title = element_text(size = 12)
      , axis.text = element_text(size = 10)      
      , legend.title = element_blank()
      , legend.text = element_text(size = 10)
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


# AVERAGE MANY SIMS CHECK  ====
# the exp estimation fails(somewhat) for selected pnull in (0.9, 0.95)
# with Emualt_large = 0.50... but then works if pnull is higher

debugSource('0-functions.r')
manysum_real_large = average_many_sims(
  ematemp
  , N = 1e4
  , T_ = 200
  , pnull = 0.5
  , Emualt = 0.5
  , tgood = 2.6
  , smarg = 0.5
  , nsim = 10
  , tgoodhat = 2.6
  , pnullhat = 0
  , shapehat = 1/8
  , nulldf = 100
  , tbarlist = seq(0,6,0.5)
)


custom_plot(manysum_real_large)


# SINGLE SIM CHECK ====
# it seems that the exponential estimator can get thrown off
# if you have both (a) a ton of nulls and (b)
# non nulls are really really distinct,
# it seems like in this case, the "approximation" that both
# nulls and non nulls are drawn from the same distrbituion
# kinda falls apart?
# but empirically, we don't really see data like this

N = 1e5
T_ = 200
emat = ematemp
pnull = 0.5
Emualt = 0.12
tgood = 2.6
smarg = 0.5

tbarlist = seq(0,6,0.5)

estat = estatsim(emat, N, T_)
tdat = estat_to_tdat(N, T_, pnull, Emualt, estat)
tdatpub = tdat_to_tdatpub(tdat, tgood = tgood, smarg = smarg )

# actual fdr
fdr_actual = estimate_fdr(
  tdat$t, tbarlist = tbarlist, null = tdat$null
) %>% 
  transmute(tbar, dr_actual = dr, fdr_actual)

fdr_actual

# hand estimate exponential
scalehat = mean(tdatpub$t[tdatpub$t>2.6]) - 2.6
 
# estiamte gamma 
gshapehat = 1/8
minme = function(scale){
  (mean(tdatpub$t[tdatpub$t>2.6])
   - extrunc('gamma', a=2.6, b = Inf, shape = gshapehat, rate = 1/scale))^2
}
temp = optimise(minme, c(0.1, 10 ) )
gscalehat = temp$minimum


fdr =  fdr_actual %>% 
  mutate(
    dr_hat = exp(-tbar/scalehat)
    , fdr_hat = 2*pt(-tbar, df=100)/dr_hat
    , fdr_actual_hat = fdr_actual/fdr_hat
    , fdr_ghat = 2*pt(-tbar, df=100)/(1-pgamma(tbar,gshapehat,1/gscalehat))
    , fdr_actual_ghat = fdr_actual/fdr_ghat
  ) %>% 
  mutate_at(
    .vars = vars(starts_with('fdr'))
    , .funs = ~pmin(. , 1.0)
  ) %>% 
  mutate_all(
    .funs = ~round(.,3)
  )


# plot check
# plot fdr
p_fdr = ggplot(
  fdr %>%
    select(tbar, fdr_actual, fdr_hat, fdr_ghat) %>%
    pivot_longer(-tbar, names_to = 'group', values_to = 'fdr')
  , aes(x=tbar, y=fdr, group = group) 
  ) +
  geom_line(
    aes(color = group)
    , size = 2
  )

# check fit 
edge = seq(0,8,0.1)
t_fit = rexp(1e5, 1/scalehat)
t_gfit = rgamma(1e5, gshapehat, 1/gscalehat)
tgoodhat = 2.6
temp = data.frame(
  t = tdatpub$t, group = 'obs'
) %>% 
  rbind(
    data.frame(t = t_fit, group = 'fit')
  ) %>%
  rbind(
    data.frame(t = tdat$t, group = 'all') 
  ) %>% 
  rbind(
    data.frame(t=t_gfit, group = 'gfit')
  ) %>% 
  filter(t>edge[1],t<edge[length(edge)]) %>% 
  group_by(group) %>% 
  summarise(
    density = hist(t,edge)$density / sum(t>tgoodhat)*n()
    , t = hist(t,edge)$mids
  ) 

p_fit = ggplot(
  temp, aes(x=t,y=density,group=group)
) +
  geom_line(aes(color=group), size = 1.6) +
  coord_cartesian(
    ylim = c(0,30)
  )


fdr

grid.arrange(p_fdr,p_fit)
