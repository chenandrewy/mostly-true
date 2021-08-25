# 2021 08 andrew testing bh style analysis in r


# ENVIRONMENT ====

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

### USER ENTRY
# root of April 2021 release on Gdrive
pathRelease = 'https://drive.google.com/drive/folders/1I6nMmo8k_zGCcp9tUvmMedKTAkb9734R'
url_prefix = 'https://drive.google.com/uc?export=download&id='


# login to gdrive
# this prompts a login
pathRelease %>% drive_ls()



# DOWNLOAD DATA =====

# summary
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='Portfolios') %>% drive_ls() %>% 
  filter(name=='Full Sets OP') %>% drive_ls() %>% 
  filter(name=='PredictorSummary.xlsx')

# (make sure file is closed)
drive_download(target_dribble, path = '../data/deleteme.xlsx', overwrite = T)

cz = read_excel('../data/deleteme.xlsx',sheet='short')


# GAMMA ESTIMATION ====
scale = 2
n = 5000
tgood = 2.6

t_emp = sample(cz$tstat,length(cz$tstat),replace=T)
t_emp = cz$tstat

shape = 0.25


randt = function(scale,tmin){
  rtrunc(n, 'gamma', a=tmin, b=Inf, shape, scale)
}

obj = function(scale){
  x = randt(scale,tgood)
  error = (
    mean(x)-mean(t_emp[t_emp>tgood])
  )^2 
  + (
    sd(x)-sd(t_emp[t_emp>tgood])
  )^2
}


# optimize
temp = optimize(obj, c(0.1,4))
scalehat = temp$minimum
objhat = temp$objective



# plot
edge = seq(-3,15,0.5)

t = t_emp
t = t[t>min(edge) & t<max(edge)]
hemp = hist(t,edge)

tsim = randt(scalehat,0)
t = tsim
t = t[t>min(edge) & t<max(edge)]
hsim = hist(t,edge)


plot(hsim$density)

plotme = rbind(
  data.frame(
    t = hemp$mids
    , f = hemp$density / sum(t_emp>tgood) * length(t_emp)
    , group = 'emp'
  )
  , data.frame(
    t = hsim$mids
    , f = hsim$density / sum(tsim>tgood) * length(tsim)
    , group = 'fit'
  )
)

p1 = ggplot(plotme, aes(x=t, y=f, fill=group)) +
  geom_bar(stat='identity', position='identity',alpha=0.6) 

p2 = ggplot(
    plotme %>% filter(t>1.0)
    , aes(x=t, y=f, fill=group)
    ) +
  geom_bar(stat='identity', position='identity',alpha=0.6) +
  xlim(0,15)


grid.arrange(p1, p2, nrow=1)

C = length(tsim)/sum(tsim>2)
print(C)
print(scalehat)
print(objhat)


# MIXTURE ESTIMATION ====
shape = 3
n = 1e6
tgood = 2.6
pnull = 0.95

t_emp = sample(cz$tstat,length(cz$tstat),replace=T)

t_emp = cz$tstat


qlist = c(0.25, 0.5, 0.75)

momfun = function(t,i){
  quantile(t[i],qlist)
}
momboot = boot(t_emp, momfun, R=100)


randt = function(pnull,scale,tmin){
  temp1 = rtrunc(round(n*(1-pnull)), 'gamma', a=tmin, b=Inf, shape, scale)
  temp2 = rtrunc(round(n*pnull), 'norm', a=tmin, b=Inf, 0, 1)
  t = c(temp1,temp2) 
}


w = 1/diag(var(momboot$t))
obj = function(par){
  # par is (pnull,scale)
  x = randt(par[1],par[2],tgood)
  mean( w*(quantile(x, qlist) 
  - quantile(t_emp[t_emp>tgood], qlist))^2 )  
}
obj1 = function(par1){obj(c(pnull,par1))}

# optimize

temp = optimize(obj1, c(0.1,10))


scalehat = temp$minimum
objhat = temp$objective



# plot
edge = seq(-3,15,0.5)

t = t_emp
t = t[t>min(edge) & t<max(edge)]
hemp = hist(t,edge)

tsim = randt(pnull,scalehat,-Inf)
t = tsim
t = t[t>min(edge) & t<max(edge)]
hsim = hist(t,edge)


plot(hsim$density)

plotme = rbind(
  data.frame(
    t = hemp$mids
    , f = hemp$density / sum(t_emp>tgood) * length(t_emp)
    , group = 'emp'
  )
  , data.frame(
    t = hsim$mids
    , f = hsim$density / sum(tsim>tgood) * length(tsim)
    , group = 'fit'
  )
)

p1 = ggplot(plotme, aes(x=t, y=f, fill=group)) +
  geom_bar(stat='identity', position='identity',alpha=0.6) 

p2 = ggplot(
  plotme %>% filter(t>1.0)
  , aes(x=t, y=f, fill=group)
) +
  geom_bar(stat='identity', position='identity',alpha=0.6) +
  xlim(0,15)


grid.arrange(p1, p2, nrow=1)

C = length(tsim)/sum(tsim>2)
print(C)
print(scalehat)
print(objhat)



# SYMMETRIC MIXTURE ESTIMATION ====
pnull = 0.5
shape = 2
scale = 1/shape
n = 5000
tgood = 2.6


t_emp = sample(cz$tstat,length(cz$tstat),replace=T)

t_emp = cz$tstat


qlist = c(0.5)

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



# optimize in 2D
# fit = optim(c(pnull,scale), obj)
# pnullhat = fit$par[1]
# scalehat = fit$par[2]

# optimize in 1D many times
pnulllist = seq(0.1,0.9,0.07)
scalelist = numeric(length(pnulllist))
objlist = numeric(length(pnulllist))
for (pi in 1:length(pnulllist)){
  obj1 = function(scale){obj(c(pnulllist[pi],scale))}
  temp = optimize(obj1, c(0.1/shape, 10.0/shape))    
  scalelist[pi] = temp$minimum
  objlist[pi] = temp$objective
}

istar = which(objlist == min(objlist))
pnullhat = pnulllist[istar]
scalehat = scalelist[istar]

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
  geom_histogram( position = 'identity', bins=40, alpha = 0.6, show.legend = F ) 

grid.arrange(p1, p2, p0, nrow=1)

C = length(tsim$t)/sum(tsim$t>2.6)
print(C)

print(pnullhat)
print(scalehat)
print(scalehat*shape)

# FDR BOUNDS ====

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


# ====
x = fread('d:/Google Drive/Work/Chen and Velikov/JFQA Revision/Inputs/hf-spreads-all/data-from-wrds/hf_monthly.csv')


dup = x %>% group_by(permno,yearm) %>% filter( n() > 1)
dup = dup %>% arrange(yearm,permno)

View(dup)

# ====

retsum = ret0 %>% filter(insamp) %>% group_by(signalname) %>% summarize( t = mean(ret)/sd(ret)*sqrt(n()))

t20 = quantile(retsum$t, 0.20)

bh = run_bh_plus(retsum$t, 3)

i20 = sum(bh$fdr$tsort > t20)

bh$fdr$fdrhat[i20]

# extrunc stuff ====

shape = 0.5
scale = 1/shape
extrunc('gamma', a=2.6, b=Inf, shape, 1/scale)


t = rgamma(1e5, shape, 1/scale)

hist(t,30)

C = 1/(sum(t>2.6)/length(t))
print(C)


# playing with pareto ====
library(EnvStats)
library(ggplot2)
library(dplyr)
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

t = signalsum$t
t = t[t>2]



hist(rpareto(1e3,1,100)-1,40)
hist(t,40)


mean(t)

loc = 1
shapehat = length(t)/sum(log(t/loc))


plotme = data.frame(
  t = rpareto(1e3,loc,shapehat), group = 'fit'
) %>% rbind(
  data.frame(
    t = t, group = 'emp'
  )
)

ggplot(plotme, aes(x=t,fill=group)) +
  geom_histogram(aes(y=..density..), bins=40, position = 'identity', alpha = 0.6) +
  xlim(0,10)
  
  
# plot
edge = seq(-3,15,0.5)

ttemp = t_emp
ttemp = ttemp[ttemp>min(edge) & ttemp<max(edge)]
hemp = hist(ttemp,edge)

tsim = rpareto(1e3,loc,shapehat)
ttemp = tsim
ttemp = ttemp[ttemp>min(edge) & ttemp<max(edge)]
hsim = hist(ttemp,edge)


plotme = rbind(
  data.frame(
    t = hemp$mids
    , f = hemp$density / sum(t_emp>tgood) * length(t_emp)
    , group = 'emp'
  )
  , data.frame(
    t = hsim$mids
    , f = hsim$density / sum(tsim>tgood) * length(tsim)
    , group = 'fit'
  )
)

ggplot(plotme, aes(x=t, y=f, fill=group)) +
  geom_bar(stat='identity', position='identity',alpha=0.6, show.legend = F) 


# Closed form estimation of bias ====
library(data.table)

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

t_emp = signalsum$t
t_emp = t_emp[t_emp>1.96]





# ====

# add actuals
thurdle = seq(0,6,0.5)
fdr_actual = numeric(length(thurdle)) + NA_real_
t = sim$t
null = sim$null

fdr_actual = numeric(length(thurdle))
for (ti in 1:length(thurdle)){
  idisc =t>thurdle[ti]
  fdr_actual[ti] = mean(as.numeric(null[idisc]))
}


fdr = data.frame(
  t = t, fdr_actual
)



# check by hand robustness of extrapolations! ====
# 08 23

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

emat = ematemp
Temat = Temp
Nemat = Nemp


N = 1e4
T_ = 200
pnull = 0.9
Emualt = 0.5
tgood = 2.6
smarg = 0.5

# simulate  
estat = estatsim(N, T_)
tdat = estat_to_tdat(N, pnull, Emualt, estat)
tdatpub = tdat_to_tdatpub(tdat, tgood = tgood, smarg = smarg )

tbarlist = seq(0,6,0.5)
fdr_actual = estimate_fdr(
  tdat$t, tbarlist = tbarlist, null = tdat$null
) %>% 
  transmute(tbar, dr_actual = dr, fdr_actual)

# hand estimate
shapehat = mean(tdatpub$t[tdatpub$t>2.6]) - tgood
tsim = rexp(1e5, 1/shapehat)

fdr_actual %>% 
  mutate(
    # dr_hat = round(exp(-tbar/shapehat) , 4)
    dr_hat = round( 1-pexp(tbar, 1/shapehat) , 4)
    , fdr_hat = round( 2*pt(-tbar,100)/dr_hat, 4)
    , actual_hat = fdr_actual/fdr_hat
  ) 


sum(tsim > 1.0) / length(tsim)

hist(tdatpub$t)





# bootstrapping with correlated blocks ====


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

nblock = 10
T_ = 100

N = 2000
pnull = 0.5
Emualt = 0.5
tgood = 2.6
smarg = 0.5

emat = ematemp
rho = 0.5
i = sample(1:dim(emat)[2], N, replace = T)
t = sample(1:dim(emat)[1], T_, replace = T)


noise = matrix(rnorm(N*T_,0,5), N, T_)
emat2 = rho*emat[t,i]+(1-rho)*noise
ebar = colMeans(emat2)
evol = sqrt(colMeans(emat2^2))

temp = cor(emat2)

# clean and output
estat = data.frame(
  bar = ebar
  , vol = evol
)


tdat = estat_to_tdat(N, T_, pnull, Emualt, estat)
tdatpub = tdat_to_tdatpub(tdat, tgood = tgood, smarg = smarg )

tbarlist = seq(0,6,0.5)
fdr_actual = estimate_fdr(
  tdat$t, tbarlist = tbarlist, null = tdat$null
) %>% 
  transmute(tbar, dr_actual = dr, fdr_actual)


# hand estimate
shapehat = mean(tdatpub$t[tdatpub$t>2.6]) - tgood
tsim = rexp(1e5, 1/shapehat)

fdr_actual %>% 
  mutate(
    # dr_hat = round(exp(-tbar/shapehat) , 4)
    dr_hat = round( 1-pexp(tbar, 1/shapehat) , 4)
    , fdr_hat = round( 2*pt(-tbar,100)/dr_hat, 4)
    , actual_hat = fdr_actual/fdr_hat
  ) 
