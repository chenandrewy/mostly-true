# 2021 08 Andrew
# Shows validity of BH under correlations


# ==== ENVIRONMENT ====
rm(list=ls())
library(data.table)
library(tidyverse)
library(ggplot2)
source('0-functions.r')

# read data
ret0 = fread('../data/clean_ret.csv')

# make a balanced matrix
retwide0 = pivot_wider(
  ret0, c(signalname,date,ret), names_from = signalname, values_from = ret
)
retmat = retwide0 %>% select(-date) %>% as.matrix()
colnames(retmat) = NULL
numt = colSums(!is.na(retmat))
igood = numt > 500
retmat = retmat[ , igood]
retmat = retmat[complete.cases(retmat) , ]

# generate residuals
rbar = colMeans(retmat) %>% as.matrix() %>% t()
l = 1+integer(nrow(retmat)) %>% as.matrix 
emat = retmat - l %*% rbar

# ==== SIMULATE ====
source('0-functions.r')
N = 1000;
pnull = 0.5
shape = 1
scale = 0.5
Nemp = dim(retmat)[2]

simulate = function(N,pnull,shape,scale){
  
  # simulate null
  Nalt = sum(runif(N) > pnull)
  null = integer(N)+1
  null[1:Nalt] = 0
  null = as.logical(null)
  
  # simulate mu
  mu = numeric(N)
  mu[!null] = rgamma(Nalt, shape, 1/scale)
  # mu[!null] = scale
  mu[1:round(Nalt/2)] = -mu[1:round(Nalt/2)]
  
  # simulate returns
  i = sample(1:Nemp,N,replace=T)
  tau = sample(1:nrow(emat), 500, replace = T)
  esim = emat[tau,i] # check me
  musim = as.matrix(1+integer(nrow(esim))) %*% t(as.matrix(mu))
  rsim = musim + esim
  
  # find t-stats
  tsim = colMeans(rsim)/apply(rsim,2,sd)*sqrt(nrow(rsim))
  tsim = abs(tsim)
  
  # output selected
  simtab = data.table(
    t = tsim
    , mu = mu
    , null = null
  )
  
  # sorted dataset
  isort = order(-tsim)
  simsort = data.frame(
    t = tsim[isort]
    , mu = mu[isort]
    , null = null[isort]
    , p = 2*pnorm(-tsim[isort])
  )
  
  simsort = simsort %>% mutate(fdr = cummean(null))
  
  
  sim = list(
    rmat = rsim
    , raw = simtab
    , sort = simsort
  )
  
} # end function

sim = simulate(N,pnull,shape,scale)

# ==== ====

# simple check better ==== 

dat = sim$sort %>%  as.data.table()

tlist = seq(0,3,0.1)

fdr = numeric(length(tlist))
fdrhat = fdr
for (ti in 1:length(tlist)){
  i = which(dat$t>tlist[ti])
  fdr[ti] = mean(as.numeric(dat$null[i]))
}
fdr

# ez bh calc ====


temp = sim$sort %>% 
  mutate(
    fdrhat = 2*pnorm(-t)/(row_number()/N)
    , fdr/fdrhat
  )


plotme = temp %>% 
  select(c(t, fdr,fdrhat)) %>% 
  filter( t < 10) %>% 
  pivot_longer(
    c(fdr,fdrhat)
  , names_to = 'type', values_to = 'fdr'
)
ggplot(data=plotme, aes(x=t, y=fdr, group = type)) + 
  geom_line(aes(linetype=type)) + 
  ylim(0,1.0)


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

# ====

thurdle = bh$hurdle$thurdle

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


