# 2021 08 Andrew
# Simulation verification


# ENVIRONMENT ====
rm(list=ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(gridExtra)
source('0-functions.r')

# read data
retwide = fread('../data/balanced_ret.csv')
retmat = as.matrix(retwide)

# generate residuals
rbar = colMeans(retmat) %>% as.matrix() %>% t()
l = 1+integer(nrow(retmat)) %>% as.matrix 
emat = retmat - l %*% rbar

# SIMULATE ====

# It seems like no matter what kind of simulation you set up,
# the FDR is always tiny for t>3.0.  
# This seems like a numerical issue, since theoreitcally
# FDR(3.0) can be as high as pnull, and will be pnull
# if pnull = 1 or the alternative dist -> N(0,1)

# thing is, when I crank pnull -> 1 or alt dist -> N(0,1)
# I get no observatiosn for t > 3.0.  
# Intuitively, pnorm(-3) is about 1 in 1000, so 
# 

# I feel like something is up with my sim

set.seed(1)
source('0-functions.r')
N = 1e5;
pnull = 0.9
shape = 1
scale = 1
loc   = 0.5
Nemp = dim(retmat)[2]

simulate = function(N,pnull,shape,scale,loc){
  
  # simulate null
  Nalt = sum(runif(N) > pnull)
  null = integer(N)+1
  null[1:Nalt] = 0
  null = as.logical(null)
  
  # simulate mu
  mu = numeric(N)
  # mu[!null] = rgamma(Nalt, shape, 1/scale)
  mu[!null] = rnorm(Nalt, loc, scale)
  mu[1:round(Nalt/2)] = -mu[1:round(Nalt/2)]
  
  # simulate returns
  # i = sample(1:Nemp,N,replace=T)
  # tau = sample(1:nrow(emat), 500, replace = T)
  # esim1 = emat[tau, ] # check me
  # esim  = esim1[   ,i]
  
  # debug
  esim = rnorm(N*500,0,5)
  esim = matrix(as.matrix(esim), nrow = 500)

  musim = as.matrix(1+integer(nrow(esim))) %*% t(as.matrix(mu))
  rsim = musim + esim
  
  # find t-stats
  tsim = colMeans(rsim)/apply(rsim,2,sd)*sqrt(nrow(rsim))
  tsim = abs(tsim)
  
  # output sorted dataset
  isort = order(-tsim)
  simsort = data.frame(
    t = tsim[isort]
    , mu = mu[isort]
    , null = null[isort]
    , p = 2*pnorm(-tsim[isort])
  )
  
  simsort = simsort %>% mutate(fdrdat = cummean(null))

} # end function
select_for_pub = function(simsort, tbad, tgood, smarg, sbar){
  simsort$u = runif(dim(simsort)[1]) 
  simsort = simsort %>% 
    mutate(
      pub = F
      , pub = case_when(
         t > tbad & t <= tgood & u < smarg ~ T
         , t > tgood & u < sbar ~ T
      )
    )
}




# ====




# plot simulation distributions
p1 = ggplot(
  sim  %>% filter(abs(t)<10)
  ,aes(x=t,fill=null)) +
  geom_histogram(position = 'identity', alpha = 0.6)
p2 = ggplot(
  sim  %>% filter(pub, abs(t)<10)
  ,aes(x=t,fill=null)) +
  geom_histogram(position = 'identity', alpha = 0.6)
grid.arrange(p1,p2)




## estimate 
source('0-functions.r')
nsim = 1e4

temp = sim %>% 
  arrange(p) %>% 
  mutate(
    fdrexp = 2*pnorm(-t)/(row_number()/N)*fitexp$C
    , fdrmix = 2*pnorm(-t)/(row_number()/N)*fitmix$C
  )




# ====

debugSource('0-functions.r')


bh = run_bh_plus(sim$raw$t, C=1, qlist = c(0.05, 0.10, 0.25, 0.50))
by = run_bh_plus(sim$raw$t, C=6, qlist = c(0.05, 0.10, 0.25, 0.50))


ggplot(data=bh$fdr, aes(x=tsort, y=fdrhat)) + geom_line()

plotme = bh$fdr %>% transmute(
  t = tsort, fdr = fdrhat, name = 'bh'
) %>% rbind(
  sim$sort %>% transmute(t, fdr, name = 'oracle')
) %>% rbind(
  by$fdr %>% transmute(
    t = tsort, fdr = fdrhat, name = 'by'
  ) 
) %>% 
filter( t<3)

ggplot(data=plotme, aes(x=t, y=fdr, group=name)) +
  geom_line(aes(linetype=name, color = name))

# actual fdrs ====

thurdle = seq(0,6,0.5)

# add actuals
hurdle_actual = bh$hurdle

t = sim$raw$t
null = sim$raw$null

fdr_actual = numeric(length(thurdle))
for (ti in 1:length(thurdle)){
  idisc =t>thurdle[ti]
  fdr_actual[ti] = mean(as.numeric(null[idisc]))
}

hurdle_actual$fdr_actual = fdr_actual

hurdle_actual

# ====

tbarlist = c(2.0, 2.6, 3.0)
C = 1

Ns = length(tsim)
fdis = 1:length(tbarlist)
pval_dis = fdis
for (tbari in seq(1,length(tbarlist))){
  tbar = tbarlist[tbari]
  
  fdis[tbari] = sum(tsim>tbar)/Ns
  pval_dis[tbari] = 2*pnorm(-tbar)
}


fdrhat = pval_dis/fdis*C
fdrhat
fdis
pval_dis

# ====



# BKY ==== 

q = 0.05
q2 = q/(1+q)

# step 1: use linear step up using q2
Ns = length(tsim)
psort = sort(2*pnorm(-tsim))
fdrhat = psort/(1:Ns/Ns)
istar2 = sum(fdrhat < q2)

# step 2: estimate n null
Nnull = Ns - istar2
pnullhat = Nnull/Ns

# step 3: use linear step up with qstar
qstar = q2/pnullhat
istarstar = sum(fdrhat < qstar)

psort[istarstar]


pnullhat


