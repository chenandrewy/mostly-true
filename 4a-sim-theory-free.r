# 2022 05 10: simulation for yz data only

# SETUP ====
rm(list=ls())
source('0-functions.r')

load('../data/emp_data.Rdata')

## User Settings ====

# data cleaning 
min_nmonth = 200
min_nsignal = 100

# dimensions
ndate = 500
nsim = 1000

# parameters
pF_list     = c(seq(0.5, 0.9, 0.1), 0.95, 0.99)
mutrue_list = c(0.25, 0.5, 0.75)

parlist = expand_grid(pF = pF_list, mutrue = mutrue_list)

# seed
set.seed(1120)

# fdr estimation
h_disc = 2 # cutoff for a discovery
fdrhat_numer = 0.05 # plug in this for Pr(F|disc)


# PREP RESIDUALS  ====
yzmonthsum = yz_ret %>%
  group_by(date) %>%
  summarize(nsignal = sum(!is.na(ret)))

# yz residuals
resid = yz_ret %>% 
    left_join(
      yz_sum, by = 'signalname'
    ) %>%
    left_join(
      yzmonthsum, by = 'date'
    ) %>%
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
isim = sample(1:N, nplot, replace = F)
iemp = sample(1:length(unique(yz_ret$signalname)), nplot, replace = F)

# find correlation matricies
csim = cor(eboot[ , isim], use = 'pairwise.complete.obs')
csim2 = csim[lower.tri(csim)]

temp = yz_ret %>% 
  pivot_wider(
    names_from = signalname, values_from = ret
  ) %>% 
  select(-date)
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
    cmid = hist(c,edge)$mids
    , density = hist(c,edge)$density
  ) %>% 
  mutate(
    group = factor(
      group
      , levels = c('sim','emp')
      , labels = c('Simulated','Yan-Zheng Data')
    )
  )


ggplot(
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


ggsave(
  filename = '../results/cor-theory-free.pdf', width = 5, height = 4
)


# SIMULATE MANY TIMES ====

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
    famstat = cross %>% 
      filter( tabs > h_disc ) %>% 
      summarize(
        fdp = mean(!verity), dr = n()/N, fdrmax = fdrhat_numer/dr
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

# TABLE: BOUND CHECK ====


## table for output ====
tabout = simdat %>% 
  group_by(pF, mutrue) %>%
  summarize(
    x1_fdr = mean(fdr)*100
    , x2_ave = mean(fdrmax/fdr)
    , x4_pct_ok =mean(fdrmax > fdr) *100
  ) %>% 
  pivot_longer(
    cols = starts_with('x'), names_to = 'stat', values_to = 'value'
  ) %>% 
  pivot_wider(
    names_from = pF, names_prefix = 'pF_', values_from = 'value'
  ) %>% 
  arrange(
    -mutrue, stat
  )


write_csv(tabout, '../results/tab-sim-theory-free.csv')

