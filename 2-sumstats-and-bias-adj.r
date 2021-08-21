# 2021 08 estimate bias adjustment 

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
source('0-functions.r')

load('../data/emp_data.Rdata')


# SUMMARY STATS ====

## univariate
tab_univar_quantiles = emp_sum %>% select(-signalname) %>% 
  apply(2, quantile, probs=qlist) %>% 
  t()


## correlations

# first do pairwise complete
retwideall = pivot_wider(
  emp_ret
  , c(signalname,date,ret), names_from = signalname, values_from = ret
) %>% 
  select(-date)
corpairwise = cor(retwideall, use = 'pairwise.complete.obs')

# then do balanced matrix with minimum missing
retwidebal = pivot_wider(
  emp_ret
  , c(signalname,date,ret), names_from = signalname, values_from = ret
) %>% 
  select(-date)
corbalanced = cor(retwidebal, use = 'complete.obs')

tab_corr_quantiles = rbind(
  quantile(
    corpairwise[lower.tri(corpairwise)], probs = qlist
  ) 
  , quantile(
    corbalanced[lower.tri(corbalanced)], probs = qlist
  ) 
) 
row.names(tab_corr_quantiles) = c('pairwise','subset complete')


## table prep
qlist = seq(0.1,0.9,0.1)

# 
finalMatrix <- rbind(tab_univar_quantiles, tab_corr_quantiles)

# Need to make these numbers shorter. Otherwise, they won't fit in one page.
finalMatrix <- round(finalMatrix, 1)
finalMatrix[3,] = as.integer(round(finalMatrix[3,],0))
finalMatrix

# Take transpose
finalMatrix <- t(finalMatrix)
# Rearrange columns
order <- c("t", "rbar", "vol", "nmonth", "pairwise", "subset complete")
finalMatrix <- finalMatrix[,order]
#Take transpose again
finalMatrix <- t(finalMatrix)

# % character in column names will cause issues with latex. Remove them.
colnames(finalMatrix) <- gsub("\\%", "", colnames(finalMatrix))

# Specify rownames
rownames(finalMatrix) <- c(Hmisc::latexTranslate("| t |"),
                           "Mean Return",
                           "Volatility",
                           "Num of Obs",
                           "Pairwise",
                           "Subset Complete")


## table 

# Produces latex table code
capture.output(
  Hmisc::latex(finalMatrix,
               file = "",
               title = '',
               table.env=F,
               cgroup = "Percentile",
               n.cgroup = 9,
               rgroup = c("Section 1", "Section 2"),
               n.rgroup = c(length(row.names(tab_univar_quantiles)),
                            length(row.names(tab_corr_quantiles))),
               cgroupTexCmd = "normalfont",
               rgroupTexCmd = "normalfont",
               col.just = c('c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'), 
               colhead = colnames(finalMatrix),
               math.row.names = F,
               first.hline.double = FALSE,
               insert.bottom = "",
               rdec = c(2, 2, 2, 0, 2, 2)),
  file = "../results/mattab-sum.tex")


# ESTIMATE AND SIMULATE ====
t_emp = emp_sum$t

# FIXED PARAMETERS
tgood = 2.6

pnull = 0.9
shape = 2
sigma = 1

# simulation
nsim = 1e6

fitexp = estimate_exponential(t_emp,tgood)
fitmix = estimate_mixture(t_emp,tgood,pnull,shape,sigma)

# output to console
fitexp %>% t()
fitmix %>% t()

# simulate to plot dist
# eventually: just plot the distribution function
datmix = simmix(nsim,fitmix$pnull,shape=fitmix$shape,fitmix$scalehat,sigma=1)
datexp = simmix(nsim,pnull=0,shape=1,fitexp$scalehat,sigma=1)

# PLOT TO ILLUSTRATE MODEL ====

## Plot mixture ====
edge = seq(0,8,0.1)
Gscale = 5
legtitle = 'Drawn from'

summix = datmix %>% 
  filter(t>min(edge),t<max(edge)) %>% 
  group_by(null) %>% 
  summarise(
    f = hist(t,edge)$count/nsim
    , t = hist(t,edge)$mids
  ) 

plotme = summix %>%
  mutate(
    group = factor(
      null
      , levels = c(F,T)
      , labels = c(
        'Gamma (right axis, rescaled)'
        ,'|Z| (left axis)'
        )
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    f = if_else(null,f,f*Gscale)
  )

p_mix = ggplot(plotme, aes(x=t,y=f,group=group)) +
  geom_line(aes(linetype=group, color=group), size = 2.5) +
  scale_y_continuous(
    'density (drawn from folded normal)' 
    , sec.axis = sec_axis(~.*Gscale, name= 'density (drawn from Gamma)')
  ) +
  labs(
    x = '|t|'
    , linetype = legtitle, color = legtitle
  )

p_mix

## Plot exponential ====

edge = seq(0,8,0.1)
Gscale = 5

plotme = datexp %>% 
  filter(t>min(edge),t<max(edge)) %>% 
  summarise(
    f = hist(t,edge)$count/nsim
    , t = hist(t,edge)$mids
  ) 

p_exp = ggplot(plotme, aes(x=t,y=f)) +
  geom_line(size = 2.5) +
  labs(
    x = '|t|', y = 'density'
    , linetype = legtitle, color = legtitle
  )

grid.arrange(p_mix, p_exp)

ggsave(p_mix, filename = '../results/illus_mix.pdf', width = 10, height = 8)
ggsave(p_exp, filename = '../results/illus_exp.pdf', width = 10, height = 8)

# PLOT TO SHOW EMPIRICAL FIT ====
source('0-functions.r')


## create data frame with all groups
t_exp = datemp$t
t_mix = datmix$t


datall = data.frame(t = t_emp, group = 'emp') %>% 
  rbind(
    data.frame(t = t_exp, group = 'exp')
  ) %>% 
  rbind(
    data.frame(t = t_mix, group = 'mix')
  )


edge = seq(0,10,0.25)
hall = datall %>% 
  filter(t>min(edge), t<max(edge)) %>% 
  group_by(group) %>% 
  summarise(
    tmid = hist(t,edge)$mid
    , density = hist(t,edge)$density
  ) %>% 
  left_join(
    datall %>% group_by(group) %>% summarise(Pr_good = sum(t>tgood)/n())
  ) %>% 
  mutate(
    density_good = density/Pr_good
  )

## plot and save
p_fit = ggplot(
  hall %>% filter(group != 'emp')
  , aes(x=tmid, y=density_good, color = group)
) +
  geom_line(
    aes(linetype = group, color = group) 
  ) +
  geom_bar(
    data = hall %>% filter(group=='emp')
    , stat = 'identity', alpha = 0.6
  ) 


ggsave(p_fit, filename = '../results/fitboth.pdf', width = 10, height = 8)

# ==== MIX GAMMA SMM ON DECILES (FOR APPENDIX)  ====
pnulllist = seq(0.05,0.95,0.05)
pnulllist = 0.5
shape = 2
scale = 1/shape
n = 1e5
tgood = 2.6

# t_emp = sample( signalsum$t,length(signalsum$t),replace=T)
t_emp = signalsum$t


qlist = seq(0.1,0.9,0.1)
# qlist = seq(0.5)
momfun = function(t,i){
  quantile(t[i],qlist)
}
momboot = boot(t_emp, momfun, R=100)

randt = function(pnull,scale,tmin){
  talt  = rgamma(round(n*(1-pnull)),shape,1/scale)
  tnull = rnorm(round(n*pnull))
  
  t = abs(c(tnull,talt))
  tlist = list(
    t = t[t>tmin]
    ,talt  = talt
    ,tnull = tnull  
  )
}


w = 1/diag(var(momboot$t))
# obj = function(par){
#   # par is (pnull,scale)
#   x = randt(par[1],par[2],tgood)
#   mean( w*(quantile(x$t, qlist) 
#            - quantile(t_emp[t_emp>tgood], qlist))^2 )  
# }

# testing mean
obj = function(par){
  # par is (pnull,scale)
  x = randt(par[1],par[2],tgood)
  ( mean(x$t) - mean(t_emp[t_emp>tgood]) )^2  
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
  plotme %>% filter(t>2)
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

grid.arrange(p1, p2, nrow=1)

C = length(tsim$t)/sum(tsim$t>2.6)

print(C)
print(pnullhat)
print(scalehat)
print(scalehat*shape)
print(objhat)


