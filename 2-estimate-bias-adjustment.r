# 2021 08 estimate bias adjustment 

# ==== ENVIRONMENT ====

rm(list = ls())
library(tidyverse)
library(data.table)
library(googledrive)
library(readxl)
library(RColorBrewer)
library(lubridate)
library(boot)
library(e1071)
library(truncdist)
library(gridExtra)

ret = fread('../data/clean_ret.csv')

# univariate
signalsum = ret %>% 
  filter(insamp) %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , nmonth = n()
    , t = rbar/vol*sqrt(nmonth)
  )


# ==== SYMMETRIC MIXTURE ESTIMATION ====
pnulllist = seq(0.05,0.95,0.05)
pnulllist = 0.90
shape = 2
scale = 1/shape
n = 10000
tgood = 2.6

# t_emp = sample( signalsum$t,length(signalsum$t),replace=T)
t_emp = signalsum$t


qlist = seq(0.1,0.9,0.1)
# qlist = seq(0.25,0.75,0.25)
momfun = function(t,i){
  quantile(t[i],qlist)
}
momboot = boot(t_emp, momfun, R=100)


randt = function(pnull,scale,tmin){
  altp  = rgamma(round(n*(1-pnull)/2),shape,1/scale)
  altn = -rgamma(round(n*(1-pnull)/2),shape,1/scale)
  tnull = rnorm(round(n*pnull))

  t = abs(c(tnull,altp,altn))
  tlist = list(
    t = t[t>tmin]
    ,talt  = c(altp,altn)
    ,tnull = tnull  
  )
}

# check t assumption
tsim = randt(pnull, scale, -Inf)
plotme = data.frame(t = tsim$talt, group = 'alt') %>%
  rbind(
    data.frame(t=tsim$tnull, group = 'null')
  )
ggplot(plotme, aes(x=t, fill=group)) +
  geom_histogram(position = 'identity', alpha = 0.6, bins=40)


w = 1/diag(var(momboot$t))
obj = function(par){
  # par is (pnull,scale)
  x = randt(par[1],par[2],tgood)
  mean( w*(quantile(x$t, qlist) 
           - quantile(t_emp[t_emp>tgood], qlist))^2 )  
}


# optimize in 1D many times
scalelist = numeric(length(pnulllist))
objlist = numeric(length(pnulllist))
for (pi in 1:length(pnulllist)){
  obj1 = function(scale){obj(c(pnulllist[pi],scale))}
  temp = optimize(obj1, c(0.1/shape, 20.0/shape))    
  scalelist[pi] = temp$minimum
  objlist[pi] = temp$objective
}

istar = which(objlist == min(objlist))
pnullhat = pnulllist[istar]
scalehat = scalelist[istar]
objhat = objlist[istar]

# plot
edge = seq(-3,15,0.5)

t = t_emp
t = t[t>min(edge) & t<max(edge)]
hemp = hist(t,edge)

tsim = randt(pnullhat,scalehat,-Inf)
t = tsim$t
t = t[t>min(edge) & t<max(edge)]
hsim = hist(t,edge)


plotme = rbind(
  data.frame(
    t = hemp$mids
    , f = hemp$density / sum(t_emp>tgood) * length(t_emp)
    , group = 'emp'
  )
  , data.frame(
    t = hsim$mids
    , f = hsim$density / sum(tsim$t>tgood) * length(tsim$t)
    , group = 'fit'
  )
)

p1 = ggplot(plotme, aes(x=t, y=f, fill=group)) +
  geom_bar(stat='identity', position='identity',alpha=0.6, show.legend = F) 

p2 = ggplot(
  plotme %>% filter(t>1.0)
  , aes(x=t, y=f, fill=group)
) +
  geom_bar(stat='identity', position='identity',alpha=0.6, show.legend = F) 

# plot t assumption
plotme = data.frame(t = tsim$talt, group = 'alt') %>%
rbind(
  data.frame(t=tsim$tnull, group = 'null')
)
p0 = ggplot(plotme, aes(x=t, fill=group)) +
  geom_histogram( aes(y=..density..), position = 'identity', bins=40, alpha = 0.6, show.legend = F ) 

grid.arrange(p1, p2, p0, nrow=1)

C = length(tsim$t)/sum(tsim$t>2.6)
print(C)

print(pnullhat)
print(scalehat)
print(scalehat*shape)
print(objhat)


nalt = sum(tsim$talt > 2.6)
nnull = sum(tsim$tnull > 2.6)

nnull/(nnull+nalt)

# ==== FDR BOUNDS ====

tbarlist = c(2.0, 2.6, 3.0)
C = 17

Ns = length(t_emp)
fdis = 1:length(tbarlist)
pval_dis = fdis
for (tbari in seq(1,length(tbarlist))){
  tbar = tbarlist[tbari]
  
  fdis[tbari] = sum(t_emp>tbar)/Ns
  pval_dis[tbari] = 2*pnorm(-tbar)
}


fdrhat = pval_dis/fdis*C
fdrhat
fdis
pval_dis

