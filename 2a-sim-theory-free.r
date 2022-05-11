# 2022 05 10: simulation for yz data only

# SETUP ====
rm(list=ls())
source('0-functions.r')

load('../data/emp_data.Rdata')

# MAKE RESIDUALS  ====

# parameters for balanced panel
min_nmonth = 200
min_nsignal = 150

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
    id_cols = c(signalname,date,e)
    , names_from = signalname
    , values_from = e
    )

# SIMULATE MANY TIMES ====

## settings ====

# dimensions and randomization
N = 1e4
ndate = 500
nsim = 100
signalreplace = T

# parameters
pF_list     = c(0.5, 0.9)
mutrue_list = c(0.2, 0.5)

parlist = expand_grid(pF = pF_list, mutrue = mutrue_list)

# seed
set.seed(1120)


for (simi in 1:nsim){
  tic = Sys.time()
  
  print(
    paste0('sim-theory-free: simi = ', simi, ' nsim = ', nsim)
  )
  
  ## sim noise ====
  # noise is reused for different mu parameters, which is fine
  # since I'm simulating nsim times anyway
  # takes about 5 seconds
  
  datelist = resid$date %>% unique()
  signallist = resid$signalname %>% unique()
  
  datesim = sample(datelist, ndate, replace = T)
  signalsim = sample(signallist, N, replace = signalreplace)
      # signalsim = signallist[1:N] # no random signal
  signallink = tibble(signalname = signalsim) %>% 
    mutate(signalid = row_number())
  
  esim0 = expand.grid(
    signalname = signalsim, date = datesim
  ) %>% as_tibble() %>% 
    left_join(signallink, by = 'signalname')
  
  esim = esim0 %>% inner_join(
    resid %>% select(signalname, date, e) %>% filter(!is.na(e))
    , by = c('signalname','date')
  ) 
  
  # calculate noise's contribution to t-stat (early for speed)
  ebarsim = esim %>% group_by(signalid) %>% summarize(
    ebar = mean(e), vol = sd(e), ndate = n()
  )

  
  ## sim cross ====
  h_disc = 2
  fdrhat_numer = 0.05
  
  fdrlist = tibble()
  
  for (pari in 1:dim(parlist)[1]){
    
    # load par
    par = parlist[pari,]
    
    # simulate mu, rbar, tstat
    cross = tibble(
      signalid = 1:N
      , verity = runif(N) > par$pF
      , mu = verity*par$mutrue + (1-verity)*0
    )  %>% 
      left_join(
        ebarsim, by = 'signalid'
      ) %>% 
      mutate(
        rbar = mu + ebar, tstat = rbar/vol*sqrt(ndate), tabs = abs(tstat)
      )

    # = find fdr (true and estimated) = #
    famstat = cross %>% 
      filter( tabs > h_disc ) %>% 
      summarize(
        fdp = mean(!verity), share_disc = n()/N, fdrhat = fdrhat_numer/share_disc
        , bounded = fdrhat >= fdp
      )
    
    
    # save and feedback
    fdrlist = rbind(fdrlist, cbind(par, famstat) )
    
  } # for pari
  
  ## save to disk ====
  
  fdrlist = fdrlist %>% 
    mutate(simi = simi) %>% 
    select(simi, everything())
  
  
  
  write_csv(
    fdrlist
    , paste0('../results/sim-theory-free/fdrlist-sim', sprintf('%04d', simi), '.csv')
  )
  
  
  toc = Sys.time()
  print(toc - tic)
  print(fdrlist)
  
} # for simi

# COMPUTE TABLE (TESTING) ====

csvlist = dir('../results/sim-theory-free/', full.names = T)

fdrall = tibble()
for (csvname in csvlist){
  temp = fread(csvname)
  fdrall = rbind(fdrall, temp)
} # for i 

temp = fdrall %>% 
  mutate(
    err = fdrhat - fdp
  )

hist(temp$err)

quantile(temp$err, seq(0,1,0.2))

parstat = fdrall %>% 
  group_by(pF, mutrue) %>% 
  summarize(
    fdr = mean(fdp), fdrhat_mean = mean(fdrhat)
  )

temp = fdrall %>% 
  left_join(parstat, by = c('pF','mutrue')) %>% 
  mutate(
    err = fdrhat - fdr
  ) 

hist(temp$err)

temp %>% 
  group_by(pF, mutrue) %>% 
  summarize(
    p05 = quantile(err, 0.05)
    , p25 = quantile(err, 0.25)
    , p50 = quantile(err, 0.50)
    , p75 = quantile(err, 0.75)
    , p95 = quantile(err, 0.95)
  )
