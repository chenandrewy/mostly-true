# 2022 05 10: simulation for cz data 

# SETUP ====
rm(list=ls())
source('0-functions.r')

load('../data/emp_data.Rdata')

# bootstraps residual effects om crpss=sectopm
sim_noise = function(emat, N, ndate, weight_emp, vol_noise){
  
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
  
  return = cross0
  
} # end sim noise

# turns residuals into a cross section 
sim_other_stuff = function(cross0, N, par){
  
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
        tabs > par$tbad & tabs <= par$tgood & u < par$smarg ~ T
        , tabs > par$tgood ~ T
        , T ~ F # otherwise false
      )
    )
  
  return = cross
  
} # end sim_other_stuff

# creates data for comparing cdf F1 to cdf F2 in a plot
make_dist_dat = function(F1, edge1, F2, edge2, x_match = c(-Inf,Inf), N1 = 1, showplot = F){
  
  rescale_fac = diff(F1(x_match))/diff(F2(x_match)) * diff(edge1)[1] /diff(edge2)[1]
  
  dat = tibble(
    edge = edge1, F = F1(edge1), group = 1
  ) %>% 
    rbind(
      tibble(
        edge = edge2, F = F2(edge2)*rescale_fac, group = 2
      )
    ) %>% 
    group_by(group) %>% 
    mutate(
      F = N1*F
      , dF = F - lag(F)
      , mids = 0.5*(edge + lag(edge))
      , group = factor(group, levels = c(1,2))
    ) %>% 
    filter(!is.na(dF))
  
  if (showplot) {
    dat %>% 
      ggplot(aes(x=edge, y=dF)) +
      geom_line(aes(color = group))
  }
  
  return(dat)
  
} # make_dist_dat


## User entry ====

# data cleaning 
min_nmonth = 200
min_nsignal = 100

# dimensions
N = 1e4
ndate = 200
nsim = 1000
weight_emp = 0.9
vol_noise = mean(cz_sum$vol)

# parameters
pF_list     = c(seq(0.5, 0.9, 0.1), 0.95, 0.99)
mutrue_list = c(0.25, 0.375, 0.5, 0.75)
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
parlist = parlist %>% 
  mutate(pari = 1:dim(parlist)[1]) %>% 
  select(pari, everything())


# seed
set.seed(1120)

# fdr estimation
h_disc = 2.0 # cutoff for a discovery
fdrhat_numer = 0.05
tgoodhat = 2.6 # 


## Prep Residuals  ====
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

truthlist = tibble()
estlist = tibble()
for (simi in 1:nsim){
  tic = Sys.time()
  
  print(
    paste0('sim-extrap-pub: simi = ', simi, ' nsim = ', nsim)
  )
  
  # simulate noise (reuse for each par)
  cross0 = sim_noise(emat, N, ndate, weight_emp, vol_noise)
  
  # loop over different par in parlist
  for (pari in 1:dim(parlist)[1]){
    
    # load par
    par = parlist[pari,]
    
    # simulate mu and add to noise to make rbar, tstat 
    cross = sim_other_stuff(cross0, N, par)
    pubcross = cross %>% filter(pub)
    
    # summarize truth
    truth = par %>% 
      mutate(
        npub = dim(pubcross)[1]
        , h_disc = h_disc
        , fdp = mean(!cross$verity[cross$tabs > h_disc])
        , dr  = sum(cross$tabs > h_disc)/N
        , simi = simi
        , pari = pari
      ) %>% 
      select(pari, simi, everything())
    
    # extrapolate and estimate fdrmax
    temp.exp = est_trunc_gamma(pubcross$tabs, tgoodhat, 1) %>% as_tibble()
    temp.alt = est_trunc_gamma(pubcross$tabs, tgoodhat, 0.5) %>% as_tibble()
    fit = rbind(temp.exp, temp.alt) %>% 
      mutate(
        h_disc = h_disc
        , dr_hat = 1-pgamma(h_disc, shape = shape, scale = scale)
        , fdrmax = fdrhat_numer/dr_hat
      ) %>% 
      mutate(simi = simi, pari = pari) %>% 
      select(pari, simi, everything())
    
    # store
    truthlist = rbind(truthlist, truth)
    estlist = rbind(estlist, fit)

  } # for pari
  
  toc = Sys.time()
  print(toc - tic)

} # for simi


# TABLE: BOUND CHECK ====

temp = truthlist %>% 
  group_by(pari) %>% 
  summarize(
    fdr = mean(fdp)
  )

# merge all sim data
simdat = estlist %>% 
  left_join(
    truthlist, by = c('pari','simi','h_disc')
  ) %>% 
  left_join(
    temp, by = 'pari'
  )

  
## table for output ====
tabout = simdat %>% 
  # drop sims with not enough pubs
  filter(npub > 50) %>% 
  # force fdr max <= 1 for ease of interpretation
  mutate(fdrmax = pmin(fdrmax,1)) %>% 
  group_by(tgood, pF, mutrue, shape) %>%
  summarize(
    x1_fdr = min(fdr)*100
    , x2_fdrmax   = mean(fdrmax)*100
    , x3_ratio    = mean(fdrmax/fdr)
    , x3b_ratio_med    = median(fdrmax/fdr)    
    , x3c_ratio_sd = -sd(fdrmax/fdr)
    , x4_pctok    = mean(fdrmax>=fdr)*100
    , x5_npub     = mean(npub)
  ) %>% 
  mutate_all(round, 2) %>% 
  pivot_longer(
    cols = starts_with('x'), names_to = 'stat', values_to = 'value'
  ) %>% 
  pivot_wider(
    names_from = pF, names_prefix = 'pF_', values_from = 'value'
  ) %>% 
  arrange(
    tgood, -shape, -mutrue, stat, 
  )


tabsmall = tabout %>% 
  filter(
    stat %in% c('x1_fdr', 'x3_ratio', 'x4_pctok'), tgood == 2.6, mutrue <= 0.5
  ) %>% 
  arrange(
    tgood, -mutrue,  -shape, stat
  ) %>% 
  print(n = 50)

write_csv(tabsmall, '../results/tab-sim-extrap-bench.csv')  


tabsmall2 = tabout %>% 
  filter(
    stat %in% c('x1_fdr', 'x3_ratio', 'x4_pctok'), tgood == 5, mutrue <= 0.5
  ) %>% 
  arrange(
    tgood, -mutrue,  -shape, stat
  ) %>% 
  print(n = 50)

write_csv(tabsmall2, '../results/tab-sim-extrap-alt.csv')  

# MANUALLY INSPECT ====

# set.seed(934)

## settings ====
tgoodhat = 2.6

# par = data.frame(pF = 0.8, mutrue = 0.25, tgood = 2.6, tbad = 1.96, smarg = 0.5)
par = data.frame(pF = 0.9, mutrue = 0.4, tgood = 2.6, tbad = 1.96, smarg = 0.5)

shape = 0.5

## sim ===
# simulate noise (reuse for each par)
cross0 = sim_noise(emat, N, ndate, weight_emp, vol_noise)

# simulate mu and add to noise to make rbar, tstat 
# cross = sim_other_stuff(cross0, N, par)
# temptrue = rexp(N, rate = 1/par$mutrue)
temptrue = par$mutrue
cross = cross0 %>% 
  mutate(
    signalid = 1:N
    , verity = runif(N) > par$pF
    , mu = verity*temptrue + (1-verity)*0
  ) %>% 
  mutate(
    rbar = mu + ebar, tstat = rbar/vol*sqrt(ndate), tabs = abs(tstat)
  ) %>% 
  # sim publication
  mutate(
    u = runif(N)
    , pub = case_when(
      tabs > par$tbad & tabs <= par$tgood & u < par$smarg ~ T
      , tabs > par$tgood ~ T
      , T ~ F # otherwise false
    )
  )


pubcross = cross %>% filter(pub)

# estimate 
fit = est_trunc_gamma(pubcross$tabs, tgoodhat, shape = shape) %>% 
  mutate(
    dr_hat = 1-pgamma(h_disc, shape = shape, scale = scale)
    , fdrmax = fdrhat_numer/dr_hat
  ) 

stat.fdr = cross %>% 
  filter( tabs > h_disc ) %>% 
  summarize(
    h_disc = h_disc            
    , fdp = mean(!verity)
    , dr = n()/N
  ) %>% 
  cbind(
    fit
  )

stat.fdr

## plot ====


edge = seq(0,10,0.2)
tabs_match = 3

edge_fit = seq(0,10,0.1)
F_fit = function(tabs) pgamma(tabs, shape = fit$shape, scale = fit$scale)

match_fac = sum(cross$tabs > tabs_match) / (1-F_fit(tabs_match)) *
      diff(edge)[1] /diff(edge_fit)[1]


tempdat = tibble(
  edge = edge_fit, Fdat = F_fit(edge_fit)*match_fac
) %>% 
  mutate(
    dF = Fdat - lag(Fdat), mids = 0.5*(edge + lag(edge)), pub = T, verity = T
  ) %>% 
  filter(!is.na(dF))

xlimnum = c(0, quantile(cross$tabs, 0.99)*1.2)
ylimnum = c(0, max(tempdat$dF)*1)


p1 = ggplot(cross, aes(x=tabs, group = pub)) +
  geom_histogram(aes(fill = pub), breaks = edge, position = 'stack') +
  geom_line(
    data = tempdat, aes(x=mids, y=dF)
  ) + 
  geom_vline(xintercept = tgoodhat) +
  theme(legend.position = c(7,7)/10) +
  coord_cartesian(xlim = xlimnum)


p2 = ggplot(cross, aes(x=tabs, group = verity)) +
  geom_histogram(aes(fill = verity), breaks = edge, position = 'stack') +
  geom_line(
    data = tempdat, aes(x=mids, y=dF)
  ) +
  scale_fill_brewer() +
  geom_vline(xintercept = tgoodhat) +
  theme(legend.position = c(7,7)/10) +
  coord_cartesian(xlim = xlimnum)

p3 = ggplot(pubcross, aes(x=tabs, group = verity)) +
  geom_histogram(aes(fill = verity), breaks = edge, position = 'stack') +
  geom_line(
    data = tempdat, aes(x=mids, y=dF)
  ) +
  scale_fill_economist() +
  geom_vline(xintercept = tgoodhat) +
  theme(legend.position = c(7,7)/10) +
  coord_cartesian(xlim = xlimnum) +
  ylab('count, pub only')

if (T){
  ylimmax = 500
  p1 = p1 + coord_cartesian(xlim = xlimnum, ylim = c(0,ylimmax))
  p2 = p2 + coord_cartesian(xlim = xlimnum, ylim = c(0,ylimmax))
  p3 = p3 + coord_cartesian(xlim = xlimnum, ylim = c(0,ylimmax))
}

grid.arrange(p1,p2,p3, nrow = 1)

# show stats to console
stat.fdr



# extra check ====
ggplot(pubcross, aes(x=tabs, group = verity)) +
  geom_histogram(aes(fill = verity), breaks = edge, position = 'stack')