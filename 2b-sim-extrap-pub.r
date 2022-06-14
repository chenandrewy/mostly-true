# 2022 05 10: simulation for cz data 

# SETUP ====
rm(list=ls())
source('0-functions.r')

load('../data/emp_data.Rdata')

## settings ====

# data cleaning 
min_nmonth = 200
min_nsignal = 100

# dimensions
N = 1e4
ndate = 200
nsim = 20
weight_emp = 0.9
vol_noise = mean(cz_sum$vol)

# parameters
pF_list     = c(seq(0.5, 0.9, 0.1), 0.95, 0.99)
mutrue_list = c(0.5, 0.25)
# mutrue_list = seq(0.1, 0.5, 0.1)
tgood_list = c(2.6, 5.0)

tbad = 1.96
smarg = 0.5

parlist = expand_grid(
  pF = pF_list
  , mutrue = mutrue_list
  , tgood = tgood_list
  , tbad
  , smarg
)

# seed
set.seed(1120)

# fdr estimation
# h_disc = 2.6 # cutoff for a discovery
# fdrhat_numer = 2*pt(-h_disc, 100) # plug in this for Pr(F|disc)

h_disc = 2.0 # cutoff for a discovery
fdrhat_numer = 0.05

tgoodhat = 2.6


# PREP RESIDUALS  ====
monthsum = cz_ret %>%
  group_by(date) %>%
  summarize(nsignal = sum(!is.na(ret)))

# residuals
resid = cz_ret %>% 
  left_join(
    cz_sum, by = 'signalname'
  ) %>%
  left_join(
    monthsum, by = 'date'
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
    id_cols = c(signalname,date,e)
    , names_from = signalname
    , values_from = e
  ) %>% 
  select(-date) %>% 
  as.matrix() 

# FIG: COR DIST ====

# simulate residuals once
set.seed(339)

signalselect = sample(1:dim(emat)[2], N, replace = T)
dateselect = sample(1:dim(emat)[1], ndate, replace = T)
eboot = weight_emp* emat[dateselect, signalselect] + 
      (1-weight_emp)*matrix(rnorm(N*ndate, 0, vol_noise), nrow = ndate)

# select subset to plot (for sim only)
nplot = 1000
isim = sample(1:N, nplot, replace = F)

# find correlation matricies
csim = cor(eboot[ , isim], use = 'pairwise.complete.obs')
csim2 = csim[lower.tri(csim)]

temp = cz_ret %>% 
  pivot_wider(
    names_from = signalname, values_from = ret
  ) %>% 
  select(-date)
cemp = cor(temp, use = 'pairwise.complete.obs')
cemp2 = cemp[lower.tri(cemp)]

cdat = data.frame(
  c = csim2, group = 'sim', color = NICEBLUE
) %>% rbind(
  data.frame(
    c = cemp2 , group = 'emp', color = 'gray'
  ) 
)

edge = seq(-1,1,0.1)
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
      , labels = c('Simulated','Chen-Zimmermann Data')
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
  filename = '../results/cor-extrap-pub.pdf', width = 5, height = 4
)


# SIMULATE MANY TIMES ====

for (simi in 1:nsim){
  tic = Sys.time()
  
  print(
    paste0('sim-extrap-pub: simi = ', simi, ' nsim = ', nsim)
  )
  
  ## sim noise ====
  signalselect = sample(1:dim(emat)[2], N, replace = T)
  dateselect = sample(1:dim(emat)[1], ndate, replace = T)
  eboot = weight_emp* emat[dateselect, signalselect] + 
    (1-weight_emp)*matrix(rnorm(N*ndate, 0, vol_noise), nrow = ndate)
  
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
      ) %>% 
      # sim publication
      mutate(
        u = runif(N)
        , pub = case_when(
          tabs <= par$tbad ~ F
          , tabs > par$tbad & tabs <= par$tgood & u < par$smarg ~ T
          , tabs > par$tgood ~ T
        )
      )
    
    pubcross = cross %>% filter(pub)
    
    # estimate lambda
    mean_tabs_good = pubcross %>% filter(tabs > tgoodhat) %>% pull(tabs) %>% mean()
    lambda = mean_tabs_good - tgoodhat
    
    # q = 0.75
    # q_tabs_good = pubcross %>% filter(tabs > tgoodhat) %>% pull(tabs) %>% quantile()
    # lambda = (med_tabs_good - tgoodhat)/log(1/(1-q))
    
    
    # = find fdr (true and estimated) = #
    famstat = cross %>% 
      filter( tabs > h_disc ) %>% 
      summarize(
        fdp = mean(!verity)
        , share_disc_actual = n()/N
        , share_disc = exp(-1/lambda*h_disc)
        , fdrhat = fdrhat_numer/share_disc
        , bounded = fdrhat >= fdp
        , npub = dim(pubcross)[1]
      )
    

    # save and feedback
    fdrlist = rbind(fdrlist, cbind(par, famstat) )
    
  } # for pari
  
  # save to disk ====
  
  fdrlist = fdrlist %>% 
    mutate(simi = simi) %>% 
    select(simi, everything())
  
  
  
  write_csv(
    fdrlist
    , paste0('../results/sim-extrap-pub/fdrlist-sim', sprintf('%04d', simi), '.csv')
  )
  
  
  toc = Sys.time()
  print(toc - tic)
  print(fdrlist)
  
} # for simi


# TABLE: BOUND CHECK ====

## compile sims ====
csvlist = dir('../results/sim-extrap-pub/', full.names = T)

simdat = tibble()
for (csvname in csvlist){
  temp = fread(csvname)
  simdat = rbind(simdat, temp)
} # for i 


pardat = simdat %>% 
  group_by(pF, mutrue, tgood) %>% 
  summarize(
    fdr = mean(fdp)
  ) 

simdat = simdat %>% 
  left_join(
    pardat, by = c('pF','mutrue', 'tgood')
  )


## table for output ====
tabout = simdat %>% 
  group_by(pF, mutrue, tgood) %>%
  mutate(
    denom = fdr
  ) %>% 
  summarize(
    x1_fdract = mean(fdp)*100
    , x2_fdrhat = mean(fdrhat)*100
    , x2b_ave = mean(fdrhat/denom)
    , x3_sd = -sd(fdrhat/denom)
    , x4_pctok = mean(fdrhat>denom)*100
    , x5_npub = mean(npub)
  ) %>% 
  mutate_all(round, 2) %>% 
  pivot_longer(
    cols = starts_with('x'), names_to = 'stat', values_to = 'value'
  ) %>% 
  pivot_wider(
    names_from = pF, names_prefix = 'pF_', values_from = 'value'
  ) %>% 
  arrange(
    -mutrue, tgood, stat
  )


write_csv(tabout, '../results/tab-sim-extrap-pub.csv')

tabout %>% 
  as.data.frame()

tabout %>% 
  filter(
    grepl('_fdract', stat)
  )


tabout %>% 
  filter(
    grepl('_ave', stat)
  )


tabout %>% 
  filter(
    grepl('pctok', stat)
  )

tabout %>% 
  filter(stat != 'x5_npub', stat != 'x2_fdrhat') %>% 
  as.data.frame() %>% 
    arrange(tgood, -mutrue, stat)


# TEST ====

set.seed(934)

tgoodhat = 2

par = data.frame(
  pF = 0.95, mutrue = 0.75, tgood = 2.6, tbad = 1.96, smarg = 0.5
)


# sim noise ===
signalselect = sample(1:dim(emat)[2], N, replace = T)
dateselect = sample(1:dim(emat)[1], ndate, replace = T)
eboot = weight_emp* emat[dateselect, signalselect] + 
  (1-weight_emp)*matrix(rnorm(N*ndate, 0, vol_noise), nrow = ndate)

cross0 = tibble(
  signalid = 1:N
  , ebar = apply(eboot, 2, mean, na.rm=T)
  , vol  = apply(eboot, 2, sd, na.rm=T)
  , ndate  = apply(eboot, 2, function(x) sum(!is.na(x)))
)

# simulate mu, rbar, tstat
cross = cross0 %>% 
  mutate(
    signalid = 1:N
    , verity = runif(N) > par$pF
    , mu = verity*par$mutrue + (1-verity)*0
  ) %>% 
  mutate(
    rbar = mu + ebar, tstat = rbar/vol*sqrt(ndate), tabs = abs(tstat)
  ) %>% 
  # sim publication
  mutate(
    u = runif(N)
    , pub = case_when(
      tabs <= par$tbad ~ F
      , tabs > par$tbad & tabs <= par$tgood & u < par$smarg ~ T
      , tabs > par$tgood ~ T
    )
  )

pubcross = cross %>% filter(pub)


quantile(pubcross %>% filter(tabs>tgoodhat) %>% pull(tabs))


mean_tabs_good = pubcross %>% filter(tabs > tgoodhat) %>% pull(tabs) %>% mean()

q = 0.75
q_tabs_good = pubcross %>% filter(tabs > tgoodhat) %>% pull(tabs) %>% quantile(q)

lambda = mean_tabs_good - tgoodhat
lambda_med = (q_tabs_good - tgoodhat)/log(1/(1-q))

cross %>% 
  filter( tabs > h_disc ) %>% 
  summarize(
    fdp = mean(!verity)
    , share_disc_actual = n()/N
    , share_disc = exp(-1/lambda*h_disc)
    , fdrhat = fdrhat_numer/share_disc
    , bounded = fdrhat >= fdp
    , npub = dim(pubcross)[1]
    , fdrhat_alt =fdrhat_numer / exp(-1/lambda_med*h_disc)
    , b_alt = fdrhat_alt >= fdp
  )




hist(pubcross$tabs)
hist(pubcross %>% filter(tabs>tgoodhat) %>% pull(tabs))
hist(cross %>%  pull(tabs), breaks = 50)


