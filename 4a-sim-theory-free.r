# 2022 05 10: simulation for yz data only

# about 10 min for nsim = 1000

# Setup -----------------------------------------------------------------------
rm(list=ls())
source('0-functions.r')

load('../data/emp_data.Rdata')

## User Settings ====

# select data
dm_ret = clz_ret # yz_ret or clz_ret
dm_sum = clz_sum

# data cleaning 
min_nmonth = 200
min_nsignal = 100

# dimensions
ndate = 500
nsim = 1000

# parameters
pF_list     = c(0.01, seq(0.05, 0.95, 0.05), 0.99)
mutrue_list = c(0.25, 0.5, 0.75)

parlist = expand_grid(pF = pF_list, mutrue = mutrue_list)

# seed
set.seed(1120)

# fdr estimation
h_disc = 2 # cutoff for a discovery
fdrhat_numer = 0.05 # plug in this for Pr(F|disc)
# PREP RESIDUALS  ====
dmmonthsum = dm_ret %>%
  group_by(date) %>%
  summarize(nsignal = sum(!is.na(ret)))

# dm residuals
resid = dm_ret %>% 
    left_join(dm_sum, by = 'signalname') %>%
    left_join(dmmonthsum, by = 'date') %>%
    filter(
      nmonth >= min_nmonth, nsignal >= min_nsignal, !is.na(ret)
    ) %>% 
    group_by(signalname) %>% 
    mutate(
      e = as.vector(scale(ret, center = T, scale = F))
    )

emat = resid %>% 
  pivot_wider(
    id_cols = c(date)
    , names_from = signalname
    , values_from = e
    ) %>% 
  select(-date) %>% 
  as.matrix() 

N = dim(emat)[2]

# FIG: COR DIST ====

# simulate residuals once
set.seed(339)
dateselect = sample(1:dim(emat)[1], ndate, replace = T)
eboot = emat[dateselect, ]

# select subset to plot (otherwise it takes forever)
nplot = 1000

# find correlation matricies
isim = sample(1:N, nplot, replace = F)
csim = cor(eboot[ , isim], use = 'pairwise.complete.obs')
csim2 = csim[lower.tri(csim)]

temp = dm_ret %>% 
    left_join(dm_sum, by = 'signalname') %>%
    left_join(dmmonthsum, by = 'date') %>%
    filter(
      nmonth >= min_nmonth, nsignal >= min_nsignal, !is.na(ret)
    ) %>% 
    select(date, signalname, ret) %>%
    pivot_wider(
      names_from = signalname, values_from = ret
    ) %>% 
    select(-date)
iemp = sample(1:ncol(temp), nplot, replace = F)    
cemp = cor(temp[ , iemp], use = 'pairwise.complete.obs')
cemp2 = cemp[lower.tri(cemp)]

cdat = data.frame(
  c = csim2, group = 'sim', color = NICEBLUE
) %>% rbind(
  data.frame(
    c = cemp2 , group = 'emp', color = 'gray'
  ) 
)

edge = seq(-1,1,0.05)
plotme = cdat %>% 
  group_by(group) %>% 
  summarise(
    cmid = hist(c,edge,plot=F)$mids
    , density = hist(c,edge,plot=F)$density
  ) %>% 
  mutate(
    group = factor(
      group
      , levels = c('sim','emp')
      , labels = c('Simulated','Empirical')
    )
  )

plt = ggplot(
  plotme, aes(x=cmid, y=density, group = group)
) +
  geom_line(
    aes(linetype = group, color = group), size = 2
  ) +
  theme_economist_white(gray_bg = F) +
  theme(
    axis.title = element_text(size = 12)
    , axis.text = element_text(size = 10)      
    , legend.title = element_blank()
    , legend.text = element_text(size = 10)
    , legend.key.size = unit(0.1, 'cm')
    , legend.position = c(80,80)/100
    , legend.key.width = unit(1,'cm')    
    , legend.spacing.y = unit(0.000001, 'cm')
    , legend.background = element_rect(colour = 'black', fill = 'white')    
  ) +
  labs(
    x = 'Pairwise Correlation'
    , y = 'Density'
  ) +
  scale_color_manual(
    values=c(NICEBLUE, 'gray')
  ) +
  scale_linetype_manual(values = c('solid','31'))
ggsave(filename = '../results/cor-theory-free.pdf', width = 5, height = 4)

# Simulate nsim times --------------------------------------------------------

simdat = tibble()

for (simi in 1:nsim){
  tic = Sys.time()
  
  print(
    paste0('sim-theory-free: simi = ', simi, ' nsim = ', nsim)
  )
  
  ## sim noise mat style ====
  dateselect = sample(1:dim(emat)[1], ndate, replace = T)
  eboot = emat[dateselect, ]
  
  cross0 = tibble(
    signalid = 1:N
    , ebar = apply(eboot, 2, mean, na.rm=T)
    , vol  = apply(eboot, 2, sd, na.rm=T)
    , ndate  = apply(eboot, 2, function(x) sum(!is.na(x)))
  )

  # sim cross ====
  fdrlist = tibble()
  
  for (pari in 1:dim(parlist)[1]){
    
    # load par
    par = parlist[pari,]
    
    # simulate mu, rbar, tstat
    cross = cross0 %>% 
      mutate(
      signalid = 1:N
      , verity = runif(N) > par$pF
      , mu = verity*par$mutrue + (1-verity)*0
    ) %>% 
      mutate(
        rbar = mu + ebar, tstat = rbar/vol*sqrt(ndate), tabs = abs(tstat)
      )

    # = find fdr (true and estimated) = #
    h_storey = 0.5
    famstat = cross %>% 
      summarize(
        fdp = sum(!verity & tabs > h_disc)/sum(tabs > h_disc)
        , dr = sum(tabs > h_disc)/N
        , fdrmax = fdrhat_numer/dr
        , pFmax = mean(tabs<h_storey)/(2*(pnorm(h_storey)-0.5))
        , fdrmax2 = fdrmax*pFmax
      )

    # save 
    temp = cbind(par, famstat) %>% mutate(pari = pari, simi = simi) %>% 
      select(pari, simi, everything())
    
    simdat = rbind(simdat, temp)
    
  } # for pari

  toc = Sys.time()
  print(toc - tic)
  
} # for simi
# calculate fdr
temp = simdat %>% 
  group_by(pari) %>% 
  summarize(
    fdr = mean(fdp)
  )
simdat = simdat %>% 
  left_join(
    temp, by = 'pari'
  )

# make dt
setDT(simdat)

# Convenience Save --------------------------------------------------------
save.image('../data/deleteme-sim-theory-free.RData')

# Convenience Load --------------------------------------------------------
load('../data/deleteme-sim-theory-free.RData')

# Exhibits --------------------------------------------------------

plotme = simdat %>% group_by(pF,mutrue) %>% 
  mutate(pF = 100*pF) %>%
  summarize(fdr_act = 100*mean(fdp), fdr_max = 100*mean(fdrmax)
            , fdr_max2 = 100*mean(fdrmax2)) %>% 
  pivot_longer(cols = starts_with('fdr')) %>% 
  mutate(name=factor(
    name
    , levels = c('fdr_act', 'fdr_max2', 'fdr_max')
    , labels = c('Actual', 'Visual Bound', 'Easy Bound')
  ))

# reordering legend (but preserving plotting order)
# https://stackoverflow.com/questions/50420205/how-to-reorder-overlaying-order-of-geoms-but-keep-legend-order-intact-in-ggplot
plotme$id2 = factor(plotme$name, levels = c('Easy Bound', 'Visual Bound', 'Actual'))

# loop over mutrue values
mutrue_list = unique(plotme$mutrue)
for (mutruei in mutrue_list){  
  
  plt = plotme %>% 
    filter(mutrue == mutruei) %>% 
    ggplot(aes(x=pF, y=value, group=name)) +
    geom_hline(yintercept=0, color='gray50') + 
    geom_hline(yintercept=100, color='gray50') +    
    # plot FDR and bounds
    geom_line(aes(linetype = name, color=name), size = 1.2) +
    scale_color_manual(values=c('Easy Bound'='gray60'
                                , 'Visual Bound'=MATBLUE, 'Actual'=MATRED)
                       , breaks=levels(plotme$id2)
                       , name=NULL) +
    scale_linetype_manual(values = c('Easy Bound'='dotdash'
                                     , 'Visual Bound'='31', 'Actual'='solid')
                          , breaks=levels(plotme$id2)
                          , name=NULL) +
    theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype")
      # , axis.title = element_text(size = 12)
      # , axis.text = element_text(size = 10)      
      # , legend.title = element_blank()
      # , legend.text = element_text(size = 8)
      , legend.key.size = unit(0.1, 'cm')
      , legend.position = c(20,80)/100
      , legend.key.width = unit(1,'cm')    
      , legend.spacing.y = unit(0.001, 'cm')
      , legend.background = element_rect(colour = 'white', fill = 'white')    
      , panel.grid.minor = element_blank()
    ) +
    labs(
      x = TeX('Proportion False Overall $Pr(F_d)$ (%)')
      , y = TeX('$FDR_{|t|>2}$ (%)')
    ) +
    coord_cartesian(ylim = c(0, 100)) 
  
  ggsave(
    paste0('../results/sim-dm-visual-',mutruei,'.pdf'), plt, width = 5, height = 2.5
    , scale = 1, device = cairo_pdf
  )
} # for mutruei

