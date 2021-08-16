# 2021 08 Andrew
# Simulation verification

# ENVIRONMENT ====
rm(list=ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(gridExtra)
source('0-functions.r')

# full data stats
ret = fread('../data/clean_ret.csv')
signalsum = ret %>% 
  filter(insamp) %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , nmonth = n()
    , t = rbar/vol*sqrt(nmonth)
  )

# balanced data
retwide = fread('../data/balanced_ret.csv')
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



# SIM1: RESIDUAL BOOTSTRAP ====
nsim = 100
pnull = 0.9
Emualt = 0.5
T_ = 200
tbarlist = seq(0,5,0.25)
emat = ematemp

## functions 
esimbootc = function(T_){
  imonth  = sample(1:Temp,T_, replace = T) # draw list of months
  esim = emat[imonth, ]
}

e_to_tdat = function(pnull,Emualt,emat){
  T_ = dim(emat)[1]
  N = dim(emat)[2]
  
  # adjust for true
  nnull = sum(runif(N) < pnull)
  null = logical(N)
  null[1:nnull] = T
  mumat = matrix(0, T_, N)
  mumat[ , !null ] = Emualt
  
  # observables
  rsim = emat + mumat
  t = colMeans(rsim)/apply(rsim,2,sd)*sqrt(T_)
  
  # pack and output
  data.frame(t = abs(t), null = null, mu = mumat[1, ], traw = t)
  
} # end e_to_t

average_many_sims = function(){
  # simulate many times
  simmany = data.frame()
  for (simi in 1:nsim){
    e = esimbootc(T_)
    tdat = e_to_tdat(pnull, Emualt, e)
    est = estimate_fdr(t = tdat$t, tbarlist = tbarlist, nulldf = 60, null = tdat$null)  
    est$simi = simi
    
    simmany = rbind(simmany, est)  
    
  } # end for simi
  
  # average across simulations
  manysum = simmany %>% 
    group_by(tbar) %>% 
    summarize_at(
      vars(2:(dim(simmany)[2])-2), mean, na.rm=F
    )
} # end average_many_sims

# output
manysum = average_many_sims()

manysum %>% as.data.frame()


ggplot(
  manysum %>% 
    select(tbar, starts_with('fdr')) %>% 
    pivot_longer(-tbar, names_to = 'type', values_to = 'fdr')
  , aes(x=tbar, y=fdr, group = type)
) +
  geom_line(aes(linetype=type, color = type)) 

# SIM 2: DIRECTLY SIM T DISTS BASED ON BLOCK EXTRAP OF EMP COR ====
library(mvtnorm)

# read data
retwide = fread('../data/balanced_ret.csv')
retmat = as.matrix(retwide) 
colnames(retmat) = NULL
coremp = cor(retmat) 
Nperblock = dim(coremp)[1]

# sim blocks of mv t

## settings
# model for simulation
N = 1e5
pnull = 0.9
Emualt = 0.5
vol = 5
df = 200
T_ = 200

# number of sims, fdr estimates, other
nsim = 10
nulldf = df
tbarlist = seq(2,6,0.2)

## FUNCTIONS
# function for simulating residuals' sample mean
ebarsim = function(){
  nblock = floor(N/Nperblock)+1
  # ebar0 = rmvt(nblock, coremp, df)
  # ebar0 = rmvnorm(nblock, sigma = coremp)
  ebar0 = rnorm(N)
  ebar = as.vector(ebar0)
  ebar = ebar[1:N]
} # end function

# function for generating observables and nulls
ebar_to_tdat = function(pnull,Emualt,ebar){

  # adjust for true
  nnull = sum(runif(N) < pnull)
  null = logical(N)
  null[1:nnull] = T
  mu = matrix(0, N)
  mu[!null] = Emualt
  
  # observables
  t = mu/vol*sqrt(T_) + ebar
  
  # pack and output
  data.frame(t = abs(t), null = null, mu = mu, traw = t)
  
} # end e_to_t

tic = Sys.time()

# test sim
ebar = ebarsim()
tdat = ebar_to_tdat(pnull, Emualt, ebar)
tdat = tdat_to_tdatpub(tdat)
# est_bias = estimate_exponential(tdat$t, 2.6)
est_bias = estimate_mixture(tdat$t, tgood = 2.6, pnull = 0.9, shape = 2, 1)
est = estimate_fdr(tdat$t, nulldf = 50, null = tdat$null, C = est_bias$C)

Sys.time() - tic


# simulate many times ====
simmany = data.frame()
for (simi in 1:nsim){
  print(paste0('simulation number ', simi))
  
  ebar = ebarsim()
  tdat = ebar_to_tdat(pnull, Emualt, ebar)
  tdatpub = tdat_to_tdatpub(tdat, tgood = 2.6, smarg = 0.5 )
  
  # actual fdr
  fdr_actual = estimate_fdr(
    tdat$t, tbarlist = tbarlist, null = tdat$null
  ) %>% 
  transmute(tbar, dr_actual = dr, fdr_actual)
  
  # fdrhat using exp
  est_bias = estimate_exponential(tdat$t, 2.6)
  fdr_exp = estimate_fdr(
    tdatpub$t, tbarlist = tbarlist, C = est_bias$C
  ) %>% 
  transmute(tbar, fdrhat_exp = fdrhat)
    
  # fdrhat using mix
  est_bias = estimate_mixture(tdatpub$t, tgood = 2.6, pnull = 0.9, shape = 2, 1)
  fdr_mix = estimate_fdr(
    tdatpub$t, tbarlist = tbarlist, C = est_bias$C
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
  summarize_at(
    vars(2:(dim(simmany)[2])-2), mean, na.rm=F
) %>% 
  mutate(
    Ndisc = dr_actual*N
  )


# output 
p1 = ggplot(
  manysum %>% 
    filter(Ndisc > 5) %>% 
    select(tbar, starts_with('fdr')) %>% 
    pivot_longer(-tbar, names_to = 'type', values_to = 'fdr')
  , aes(x=tbar, y=fdr, group = type)
) +
  geom_line(aes(linetype=type, color = type)) +
  coord_cartesian(ylim=c(0,0.5))

p2 = ggplot(
  manysum %>% select(tbar, Ndisc)
  , aes(x=tbar, y=Ndisc)
) +
  geom_line() +
  scale_y_continuous(trans='log')

grid.arrange(p1,p2,nrow=1)

# SIM3: BLOCKS OF RESIDUAL BOOTSTRAP ====


emat = ematemp 
Nemat = dim(emat)[2]
Temat = dim(emat)[1]

## settings
# model for simulation
N = 1e4
pnull = 0.5
Emualt = 0.5
T_ = 200
Nperblock = Nemat
nblock = floor(N/Nperblock)+1

# number of sims, fdr estimates, other
nsim = 10
nulldf = 100
tbarlist = seq(0,6,0.2)


## FUNCTIONS 

estatsim = function(){
  
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
estat_to_tdat = function(pnull,Emualt,estat){
  
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


average_many_sims = function(){
  simmany = data.frame()
  for (simi in 1:nsim){
    print(paste0('simulation number ', simi))
    
    estat = estatsim()
    tdat = estat_to_tdat(pnull, Emualt, estat)
    tdatpub = tdat_to_tdatpub(tdat, tgood = 2.6, smarg = 0.5 )
    
    # actual fdr
    fdr_actual = estimate_fdr(
      tdat$t, tbarlist = tbarlist, null = tdat$null
    ) %>% 
      transmute(tbar, dr_actual = dr, fdr_actual)
    
    # fdrhat using exp
    est_bias = estimate_exponential(tdat$t, 2.6)
    fdr_exp = estimate_fdr(
      tdatpub$t, tbarlist = tbarlist, C = est_bias$C
    ) %>% 
      transmute(tbar, fdrhat_exp = fdrhat)
    
    # fdrhat using mix
    est_bias = estimate_mixture(tdatpub$t, tgood = 2.6, pnull = 0.95, shape = 2, 1)
    fdr_mix = estimate_fdr(
      tdatpub$t, tbarlist = tbarlist, C = est_bias$C
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
    summarize_at(
      vars(2:(dim(simmany)[2])-2), mean, na.rm=F
    ) %>% 
    mutate(
      Ndisc = dr_actual*N
    )
  
} # end average_many_sims


## DO STUFF 
manysum = average_many_sims()


ggplot(
  manysum %>% 
    select(tbar, starts_with('fdr')) %>% 
    pivot_longer(-tbar, names_to = 'type', values_to = 'fdr')
  , aes(x=tbar, y=fdr, group = type)
) +
  geom_line(aes(linetype=type, color = type)) +
  coord_cartesian(ylim = c(0,1))



# SIM 4: SUPER FAST SUPER SIMPLE ====

# for fast simple testing
N = 1e7
pnull = 0.95
Et_alt = 0.25/2.7*sqrt(184)
Nalt = floor(N*(1-pnull))+1

# simulate
e = rnorm(N)
null = runif(N) < pnull
Et = numeric(N)
Et[!null] = Et_alt
t = Et + e
t = abs(t)


est = estimate_fdr(t=t, tbarlist = tbarlist, null = null)


ggplot(
  est %>% 
    select(tbar, starts_with('fdr')) %>% 
    pivot_longer(-tbar, names_to = 'type', values_to = 'fdr')
  , aes(x=tbar, y=fdr, group = type)
) +
  geom_line(aes(linetype=type, color = type)) +
  coord_cartesian(ylim = c(0,1))
