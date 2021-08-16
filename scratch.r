# GETTING BOOTSTRAP SIMULATION TO WORK ====
# want to get FDR > 20 for t > 4

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


set.seed(1)
source('0-functions.r')
N = 1e5;
nmonth = 20
pnull = 0.99
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
  mu[!null] = loc
  # mu[!null] = rgamma(Nalt, shape, 1/scale)
  # mu[!null] = rnorm(Nalt, loc, scale)
  # mu[1:round(Nalt/2)] = -mu[1:round(Nalt/2)]
  
  # simulate returns
  i = sample(1:Nemp,N,replace=T_)
  tau = sample(1:nrow(emat), nmonth, replace = T_)
  esim1 = emat[tau, ] # check me
  esim  = esim1[   ,i]
  
  # # debug
  # esim = rnorm(N*500,0,5)
  # esim = matrix(as.matrix(esim), nrow = 500)
  
  musim = as.matrix(1+integer(nrow(esim))) %*% t(as.matrix(mu))
  rsim = musim + esim
  
  # find t-stats
  tsim = colMeans(rsim)/apply(rsim,2,sd)*sqrt(nrow(rsim))
  
  simsort = data.frame(
    traw = tsim
    , t = abs(tsim)
    , null = null
    , p = 2*pnorm(-abs(tsim))
  )
  
  
} # end function

sim = simulate(N,pnull,shape,scale,loc)

## fdr calculations

# find actual fdrs
tbarlist = seq(0,6,0.25)
fdrlist = numeric(length(tbarlist))*NA
for (ti in 1:length(tbarlist)){
  temp = sim %>% filter(t>tbarlist[ti]) 
  fdrlist[ti] = mean(temp$null)
}

# estimated FDRs
Fhatall = ecdf(sim$t)

tbardat = data.frame(
  tbar = tbarlist
  , fdr = fdrlist
) %>% 
  mutate(
    fdrhat = if_else(
      Fhatall(tbar) < 1
      , 2*pnorm(-tbar)/(1-Fhatall(tbar))
      , NaN
    )
  ) 

plotme = tbardat %>% 
  pivot_longer(
    c(starts_with('fdr'))
    , names_to = 'type', values_to = 'fdr'
  )
ggplot(data=plotme, aes(x=tbar, y=fdr, group = type, color = type)) + 
  geom_line(aes(linetype=type)) 


# checking for normality of errors ====
N = 1e4;
nmonth = 50
pnull = 0.99
shape = 1
scale = 1
loc   = 0
Nemp = dim(retmat)[2]

sim = simulate(N,pnull,shape,scale,loc)


hist(sim$traw,40)
sd(sim$traw)

qqnorm(sim$traw, pch =1, frame = F)
qqline(sim$traw, lwd = 3)

# what the plot should look like ====
esim = rnorm(N*nmonth,0,5)
esim = matrix(as.matrix(esim), nrow = nmonth)
rsim = esim

# find t-stats
tsim = colMeans(rsim)/apply(rsim,2,sd)*sqrt(nrow(rsim))

qqnorm(tsim, pch =1, frame = F)
qqline(tsim, lwd = 3)


# trying to get bootstrap to be normal try 1 ====
# this one: random strat, then random months
# basically a full bootstrap
# it looks normal

Nemp = dim(emat)[2]

t = numeric(N)*NA_real_
for (booti in 1:N){
  iemp  = sample(1:Nemp,1) # pick a random empirical strategy
  r = sample(emat[,iemp],nmonth,T_) # pick some random months
  t[booti] = mean(r)/sd(r)*sqrt(nmonth)
}

qqnorm(t, pch =1, frame = F)
qqline(t, lwd = 3)


# trying to get cluster bootstrap to be normal ====
# this does really not look normal 
# for nmonth > 100, tails are too thin
# for nmonth = 10, tails are too fat
# for nmonth = 20, you get variable results

# perhaps what's going on is that, as N->infty, 
# Fhat(t_i)->F(t), and F(t) is not quite normal
# For small nmonth, F(t) is t-distributd, intuitively, and has fat tails
# But for large nmonth, something weird is happening that 
# seems to be related to the fact that the Temp is on ly 184
# so I have to redraw the some obs.  That seems to limit the
# fatness of the tails?
# could be realted to how I construct the residuals


library(e1071)    

N = 1e5
nmonth = 20
Nemp = dim(emat)[2]
Temp = dim(emat)[1]

# draw random set of months
tboot = sample(1:Temp, nmonth, replace = T_)

t = numeric(N)*NA_real_
for (booti in 1:N){
  iemp  = sample(1:Nemp,1) # pick a random empirical strategy
  r = emat[tboot,iemp]  # use previously drawn random monthhs (fixed)
  t[booti] = mean(r)/sd(r)*sqrt(nmonth)
}

# qqnorm(t, pch =1, frame = F)
# qqline(t, lwd = 3)

hist(t,40)
kurtosis(t)
kurtosis(rnorm(N))


# alternativ residual calc ====
# trying to subtract grand mean instead of col means
# this looks better, Kurtosis \approx 0.27
# not sure this is legit though, since 
# each t_i has a different mean.

# perhaps one alternative is to use the correlation matrix
# from this estimation only, to ensure
# that each t_i is marginally N(0,1)
# and that correlations match.
# the bootstrap only ensure correlations match,

# under what conditions is a bootstrap marginally N(0,1)?
# fundamentally the components are a dirac distribution
# under what conditions does the average draw from a dirac distribution
# normalized, follow N(0,1)?  

# I think Gosset's theorem has one answer

retwide = fread('../data/balanced_ret.csv')
retmat = as.matrix(retwide)
emat = retmat - mean(as.vector(retmat))

N = 1e5
nmonth = 100
Nemp = dim(emat)[2]
Temp = dim(emat)[1]

# draw random set of months
tboot = sample(1:Temp, nmonth, replace = T_)

t = numeric(N)*NA_real_
for (booti in 1:N){
  iemp  = sample(1:Nemp,1) # pick a random empirical strategy
  r = emat[tboot,iemp]  # use previously drawn random monthhs (fixed)
  t[booti] = mean(r)/sd(r)*sqrt(nmonth)
}


hist(t,40)
kurtosis(t)
kurtosis(rnorm(N))


# qqnorm(t, pch =1, frame = F)
# qqline(t, lwd = 3)

# mvnormal simulateion ====
# let's try just using empirical correlation matrix
# and then 

library(MASS)
cormat = cor(emat)
muvec = numeric(dim(cormat)[1])

r = mvrnorm(n=10, muvec, cormat)
r = c(r)

kurtosis(r)

hist(r)

# independent block cluster bootstrap attempt ====
# I realized the problem is that I was having residuals that are 
# perfectly correlated because I had only one draw of months
# The simplest way to fix this is to have independent draws of months

# shoot now I'm getting tails that are consistently fatter than normal
# could this be because there is still some residual correlation 
# across blocks that are far apart?  Since the 
# bootstrap is just reusing the data over and over, it seems there
# must be far horizon correlatino effects, especially when 
# a lot of months are drawn

nblock = 1000
nmonth = 30

tblock = matrix(NA , nrow = nblock, ncol = Nemp)
for (blocki in 1:nblock){
  imonth  = sample(1:Temp,nmonth, replace = T_) # draw list of months
  eboot = emat[imonth, ]  # find returns
  tblock[blocki, ] = colMeans(eboot)/apply(eboot,2,sd)*sqrt(nmonth) 
}

t = as.vector(tblock)

hist(abs(t),40)

sum(abs(t)>4)/length(t) / (2*pnorm(-4))

kurtosis(t)

# qqnorm(t)
# qqline(t)

# simple block bootstrap fdr ====
# no selection bias, just get something simple without distributional
# assumptions

rm(list=ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(MASS)
source('0-functions.r')

# read data
retwide = fread('../data/balanced_ret.csv')
retmat = as.matrix(retwide)

# generate residuals
rbar = colMeans(retmat) %>% as.matrix() %>% t()
l = 1+integer(nrow(retmat)) %>% as.matrix
emat = retmat - l %*% rbar


Nemp = dim(emat)[2]
Temp = dim(emat)[1]

Tsim = 100
pnull = 0.9
mualt = 0.5

## bootstrap residuals
# for pnull > 0.8 this becomes badly behved because there are so few
# discoveries for t>3.  note this is not an issue with empirical data
# because we have a lot of discoveries for t>3
# Nsim = Nemp
# imonth  = sample(1:Temp,Tsim, replace = T_) # draw list of months
# esim = emat[imonth, ]

## bootstrap block mvnormal
Nsim = 1000
rho = cor(emat)
esim = mvrnorm(n=Nsim*Tsim , numeric(Nemp), rho)
esim = matrix(esim, Tsim, Nsim, byrow = T_)


nnull = sum(runif(Nsim) < pnull)
null = logical(Nsim)
null[1:nnull] = T_
musim = matrix(0, Tsim, Nsim)
musim[ , !null ] = mualt


rsim = esim + musim
tsim = colMeans(rsim)/apply(rsim,2,sd)*sqrt(Tsim)
tsim = abs(tsim)



# estimated FDRs
Fhatall = ecdf(tsim)

tbardat = data.frame(
  tbar = tbarlist
  , fdractual = fdrlist
) %>% 
  mutate(
    fdrhat = if_else(
      Fhatall(tbar) < 1
      , 2*pnorm(-tbar)/(1-Fhatall(tbar))
      , NaN
    )
  ) 

ggplot(
  tbardat %>% pivot_longer(-tbar, names_to = 'type', values_to = 'fdr')
  , aes(x=tbar, y=fdr, group = type)
) +
  geom_line(aes(linetype=type))



# ====
N = 1e4
T_ = 100
p = 0.5

x = runif(N*T_) < p
x = matrix(x, nrow = T_)

e = x - p
ebar = colMeans(e) 
evol = sqrt(colMeans(e^2))


t = ebar/evol*sqrt(T_)

hist(t)


qqnorm(t)
qqline(t)