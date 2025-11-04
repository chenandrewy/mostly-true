# ABOUTME: Simulates publication-bias scenarios to benchmark CZ extrapolation.
# Inputs:
#   - functions.r (project utilities)
#   - data/emp_data.Rdata (predictor returns from 01-prep-data.r)
# Outputs:
#   - data/deleteme-sim-extrap-pub.RData
#   - results/cor-extrap-pub.pdf
#   - results/sim-extrap-gamma-<mutrue>.pdf
#   - results/sim-extrap-<mutrue>.pdf
#   - results/deleteme.pdf
# How to run:
#   Rscript 94-simulations-publication-bias.r
#   Rscript 94-simulations-publication-bias.r --vanilla

# takes about 2 minutes for nsim = 200
# this is just for the appendix now

# Setup -----------------------------------------------------------------------
rm(list=ls())

library(here)
here::i_am("94-simulations-publication-bias.r")

source(here("functions.r"))

data_dir <- here("data")
results_dir <- here("results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

sim_rdata_path <- file.path(data_dir, "deleteme-sim-extrap-pub.RData")
cor_extrap_pdf_path <- file.path(results_dir, "cor-extrap-pub.pdf")
deleteme_pdf_path <- file.path(results_dir, "deleteme.pdf")

load(file.path(data_dir, "emp_data.Rdata"))

## User entry ====

# data cleaning 
min_nmonth = 200
min_nsignal = 100

# dimensions
N = 1e4
nsim = 1000
ndate = 200
weight_emp = 0.65 # 0.65  matches the cor dist well
vol_noise = mean(cz_sum$vol)

# parameters
pF_list     = c(0.01, seq(0.05, 0.95, 0.05), 0.99)
mutrue_list = c(0.25, 0.5, 0.75)
tgood_list = c(2.6) # should not be a list for now

tbad = 1.96
smarg = 0.5

# fdr estimation
h_disc = 2.0 # cutoff for a discovery
fdrhat_numer = 0.05
tgoodhat = 2.0

parlist = expand_grid(
  pF = pF_list, mutrue = mutrue_list, tgood = tgood_list, tbad, smarg
)
parlist = parlist %>% 
  mutate(pari = 1:dim(parlist)[1]) %>% 
  select(pari, everything())

# seed
set.seed(1120)

## Functions ====

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

## Prep Residuals  ====
monthsum = cz_ret %>%
  group_by(date) %>%
  summarize(nsignal = sum(!is.na(ret)))

# residuals
resid = cz_ret %>% 
  left_join(cz_sum, by = 'signalname') %>%
  left_join(monthsum, by = 'date') %>%
  filter(nmonth >= min_nmonth, nsignal >= min_nsignal, !is.na(ret)) %>% 
  group_by(signalname) %>% 
  mutate(e = as.vector(scale(ret, center = T, scale = F)))

emat = resid %>% 
  pivot_wider(
    id_cols = c(date)
    , names_from = signalname
    , values_from = e
  ) %>% 
  select(-date) %>% 
  as.matrix() 

# FIG: COR DIST ----------------------------------------------------------------

# simulate residuals once
set.seed(339)

signalselect = sample(1:dim(emat)[2], N, replace = T)
dateselect = sample(1:dim(emat)[1], ndate, replace = T)
eboot = weight_emp* emat[dateselect, signalselect] + 
      (1-weight_emp)*matrix(rnorm(N*ndate, 0, vol_noise), nrow = ndate)

# select subset to plot (for sim only)
nplot = 2000
isim = sample(1:N, nplot, replace = F)

# find correlation matricies
csim = cor(eboot[ , isim], use = 'pairwise.complete.obs')
csim2 = csim[lower.tri(csim)]

temp = cz_ret %>% 
  pivot_wider(names_from = signalname, values_from = ret) %>% 
  select(-date)
# cemp = cor(temp, use = 'pairwise.complete.obs')
  cemp = cor(emat, use = 'pairwise.complete.obs')
cemp2 = cemp[lower.tri(cemp)]

cdat = data.frame(c = csim2, group = 'sim', color = NICEBLUE) %>% 
rbind(
  data.frame(c = cemp2 , group = 'emp', color = 'gray') 
)

edge = seq(-1,1,0.05)
plotme = cdat %>% 
  group_by(group) %>% 
  summarise(cmid = hist(c,edge,plot=F)$mids, density = hist(c,edge,plot=F)$density) %>% 
  mutate(
    group = factor(
      group
      , levels = c('sim','emp')
      , labels = c('Simulated','Chen-Zimmermann Data')
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
ggsave(filename = cor_extrap_pdf_path, 
  , width = 5, height = 4)

# Simulate many times ----------------------------------------------------------

for (simi in 1:nsim){

  if (simi==1){truthlist = tibble(); estlist = tibble()}

  tic = Sys.time()  
  print(paste0('sim-extrap-pub: simi = ', simi, ' nsim = ', nsim))
  
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
        , fdp = mean(!pubcross$verity[pubcross$tabs > h_disc]) # unclear if should be cross or pubcross
        , simi = simi
        , pari = pari
      ) %>% 
      select(pari, simi, everything())
    
    # extrapolate and estimate fdrmax
    h_storey = 0.5
    temp.exp = est_trunc_gamma(pubcross$tabs, tgoodhat, 1) %>% as_tibble()
    temp.alt = est_trunc_gamma(pubcross$tabs, tgoodhat, 0.5) %>% as_tibble()    
    fit = rbind(temp.exp, temp.alt) %>% 
      mutate(h_disc = h_disc
        , dr_hat = 1-pgamma(h_disc, shape=shape, rate=1/scale)
        , fdrmax = pmin(fdrhat_numer/dr_hat, 1)
        , pFmax = pmin(pgamma(h_storey, shape-shape, rate=1/scale)/(2*(pnorm(h_storey) - 0.5)), 1)
        , fdrmax2 = fdrmax*pFmax
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

## Clean up, make simdat ----------------------------------------------------------

# calculate fdr
temp = truthlist %>% 
  group_by(pari) %>% 
  summarize(fdr = mean(fdp))

# merge all sim data
simdat = estlist %>% 
  left_join(
    truthlist, by = c('pari','simi','h_disc')
  ) %>% 
  left_join(
    temp, by = 'pari'
  )

# Convenience Save --------------------------------------------------------
save.image(sim_rdata_path)

# Convenience Load --------------------------------------------------------
load(sim_rdata_path)

# Exhibits --------------------------------------------------------

## Figure of fdr actual vs bound w/ gamma fit ====

truthdat = simdat %>% 
  group_by(pF,mutrue,tgood,tbad,smarg) %>% summarize(fdr=100*mean(fdr)) %>% ungroup() %>% 
  mutate(name='act')
fitdat = simdat %>% 
      group_by(pF,mutrue,tgood,tbad,smarg,shape) %>% summarize(fdr=100*mean(fdrmax)) %>% 
      ungroup() %>% 
      mutate(name=if_else(shape==1, 'exp', 'gamma')) %>% 
      select(-shape) 
plotme = rbind(truthdat, fitdat) %>% 
  mutate(pF=100*pF) %>% 
  mutate(name=factor(name
    , levels = c('act', 'exp', 'gamma')
    , labels = c('Actual', 'Exp Easy Bound', 'Gamma Easy Bound')))

# reordering legend (but preserving plotting order)
# https://stackoverflow.com/questions/50420205/how-to-reorder-overlaying-order-of-geoms-but-keep-legend-order-intact-in-ggplot
plotme$id2 = factor(plotme$name, levels = c('Gamma Easy Bound', 'Exp Easy Bound', 'Actual'))

# loop over mutrue values
mutrue_list = unique(plotme$mutrue)
for (mutruei in mutrue_list){  

  plt = plotme %>% 
    filter(mutrue == mutruei) %>% 
    ggplot(aes(x=pF, y=fdr, group=name)) +
    geom_hline(yintercept=0, color='gray50') + 
    geom_hline(yintercept=100, color='gray50') +    
    # plot FDR and bounds
    geom_line(aes(linetype = name, color=name), size = 1.2) +
    scale_color_manual(values=c('Exp Easy Bound'='gray60', 'Gamma Easy Bound'=MATBLUE, 'Actual'=MATRED)
      , breaks=levels(plotme$id2)) +
    scale_linetype_manual(values = c('Exp Easy Bound'='dotdash', 'Gamma Easy Bound'='31', 'Actual'='solid')
      , breaks=levels(plotme$id2)) +
    theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype")
      , axis.title = element_text(size = 12)
      , axis.text = element_text(size = 10)      
      , legend.title = element_blank()
      , legend.text = element_text(size = 10)
      , legend.key.size = unit(0.1, 'cm')
      , legend.position = c(30,75)/100
      , legend.key.width = unit(1,'cm')    
      # , legend.spacing.y = unit(0.5, 'cm')
      , legend.background = element_rect(colour = 'white', fill = 'white')    
      , panel.grid.minor = element_blank()
    ) +
    labs(
      x = TeX('Proportion Null Overall $Pr(null_i)$ (%)')
      , y = TeX('$FDR_{|t|>2}$ (%)')
    ) +
    coord_cartesian(ylim = c(0, 100)) 
  
  ggsave(
    file.path(results_dir, sprintf('sim-extrap-gamma-%s.pdf', mutruei)), plt, width = 8, height = 4
    , scale = 0.6, device = cairo_pdf
  )
} # for mutruei

## Figure of fdr actual vs bound exp only fit ====

truthdat = simdat %>% 
  group_by(pF,mutrue,tgood,tbad,smarg) %>% summarize(fdr=100*mean(fdr)) %>% ungroup() %>% 
  mutate(name='act')
fitdat = simdat %>% 
      group_by(pF,mutrue,tgood,tbad,smarg,shape) %>% summarize(fdr=100*mean(fdrmax)) %>% 
      ungroup() %>% 
      mutate(name=if_else(shape==1, 'exp', 'gamma')) %>% 
      select(-shape) 
plotme = rbind(truthdat, fitdat) %>% 
  mutate(pF=100*pF) %>% 
  mutate(name=factor(name
    , levels = c('act', 'exp', 'gamma')
    , labels = c('Actual', 'Exp Easy Bound', 'Gamma Easy Bound')))

# reordering legend (but preserving plotting order)
# https://stackoverflow.com/questions/50420205/how-to-reorder-overlaying-order-of-geoms-but-keep-legend-order-intact-in-ggplot
plotme$id2 = factor(plotme$name, levels = c('Gamma Easy Bound', 'Exp Easy Bound', 'Actual'))

# loop over mutrue values
mutrue_list = unique(plotme$mutrue)
for (mutruei in mutrue_list){  

  plt = plotme %>% 
    filter(mutrue == mutruei) %>% 
    filter(name != 'Gamma Easy Bound') %>% 
    ggplot(aes(x=pF, y=fdr, group=name)) +
    geom_hline(yintercept=0, color='gray50') + 
    geom_hline(yintercept=100, color='gray50') +    
    # plot FDR and bounds
    geom_line(aes(linetype = name, color=name), size = 1.2) +
    scale_color_manual(values=c('Exp Easy Bound'='gray60', 'Gamma Easy Bound'=MATBLUE, 'Actual'=MATRED)
      , breaks=levels(plotme$id2)) +
    scale_linetype_manual(values = c('Exp Easy Bound'='dotdash', 'Gamma Easy Bound'='31', 'Actual'='solid')
      , breaks=levels(plotme$id2)) +
    theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype")
      , axis.title = element_text(size = 12)
      , axis.text = element_text(size = 10)      
      , legend.title = element_blank()
      , legend.text = element_text(size = 10)
      , legend.key.size = unit(0.1, 'cm')
      , legend.position = c(30,75)/100
      , legend.key.width = unit(1,'cm')    
      # , legend.spacing.y = unit(0.5, 'cm')
      , legend.background = element_rect(colour = 'white', fill = 'white')    
      , panel.grid.minor = element_blank()
    ) +
    labs(
      x = TeX('Proportion Null Overall $Pr(null_i)$ (%)')
      , y = TeX('$FDR_{|t|>2}$ (%)')
    ) +
    coord_cartesian(ylim = c(0, 100)) 
  
  ggsave(
    file.path(results_dir, sprintf('sim-extrap-%s.pdf', mutruei)), plt, width = 8, height = 4
    , scale = 0.6, device = cairo_pdf
  )
} # for mutruei

## Illustrate fitting problems ====
# Fitting problems are consistent with "do t-stat hurdles need to be raised"
# In principle, the mass below t = 1.96 can be arbitrarily large
# and this mass is invisible.

# Previous drafts of "mostly true" dealt with this using a Gamma
# distribution with shape 1/2, which is just very conservative
# and does the for reasonable parameter values
# But the shape 1/2 is both arbitrary and does not allow for
# closed form estimation (memorylessness)

# Moreover, the data mining bounds show that the extremely large
# missing mass is not reasonable. Data mining doesn't generate a crazy
# large missing mass, so the research process should not either

# For these reasons, I'm not sure it's worth including
# all of these sims

# find par that violates the bounds
parbadlist = simdat %>% group_by(mutrue,pF,tgood,tbad,smarg,shape) %>% 
  summarize(fdr=mean(fdr), fdrmax=mean(fdrmax)) %>% 
  filter(fdr > fdrmax)
print(parbadlist)
parbad = parbadlist[1, ]

# simulate 
cross0 = sim_noise(emat, N, ndate, weight_emp, vol_noise)
cross = sim_other_stuff(cross0, N, parbad)
pubcross = cross %>% filter(pub)

# simulate fit and bind
# exponential
tgoodhattemp = 2
Etabs = mean(pubcross$tabs[pubcross$tabs>tgoodhattemp]) - tgoodhattemp

# gamma
fitgamma = est_trunc_gamma(pubcross$tabs, tgoodhattemp, 1/2, 1)

plotme = tibble(tabs=rexp(N, rate=1/Etabs), name='exp') %>%
  rbind(
    tibble(tabs=rgamma(N, shape=fitgamma$shape, rate=1/fitgamma$scale)
      , name='gamma')
  ) %>%
  rbind(
    tibble(tabs=cross$tabs, name='act')
  )

# plot
p = ggplot(plotme, aes(x=tabs)) +
  # geom_histogram(breaks=seq(0,20,0.2), position='identity', alpha=0.7
  #   , aes(fill=name)) +
  geom_density(aes(color=name, linetype=name), size=1.2, adjust=3) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 8)) +
  scale_x_continuous(breaks = seq(0, 15, 2))
ggsave(deleteme_pdf_path)
