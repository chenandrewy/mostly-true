# 2024 05 bootstrap check of the standard normal null
# nboot = 1000 takes about 7 minutes on 4 cores

# Setup -----------------------------------------------------------------------
rm(list=ls())
source('0-functions.r')
load('../data/emp_data.Rdata')

## User ====
nboot = 1000 
ncores = 4
min_obs_pct = 80 # minimum observations define sample
histedge = c(-Inf, seq(-6, 6, 0.25), Inf)

ret = clz_ret %>% copy() # could use clz_vw or yz

# define sample period --------------------------------------------------------
# tbc: unify with sim-theory-free

ret[ , yr := year(date)]

obs_sum = ret[ , .(nobs = .N), by = c('yr', 'signalname')] %>% 
    group_by(yr) %>% 
    summarize(n_complete = sum(nobs==12)) %>% 
    mutate(pct_complete = n_complete/length(unique(ret$signalname))*100)

obs_sum %>% print(n=Inf)

samp = obs_sum %>% 
    filter(pct_complete > min_obs_pct) %>% 
    summarize(start = min(yr), end = max(yr))

# create panel of residuals ----------------------------------------

# filter to sample period
ret2 = ret[yr %in% samp$start:samp$end] %>% copy() %>% 
    select(-yr)

# calculate residuals
ret2[  
    , rbar := mean(ret), by = signalname
][
    , resid := ret - rbar
]

# order for fast merging
setkey(ret2, signalname, date)

# bootstrap -------------------------------------------------------------------
nmonth = length(unique(ret2$date))

datelist = unique(ret2$date) %>% sort()
signallist = unique(ret2$signalname)

# create histogram info
histinfo = data.table(left = histedge[-length(histedge)], right = histedge[-1]) %>% 
    mutate(mid = (left + right)/2, bin = 1:(length(histedge)-1))

# create bootstrap samples
tic = Sys.time()
set.seed(123)
cl = makeCluster(ncores)
registerDoParallel(cl)
boothist = foreach(i = 1:nboot, .combine = rbind
    , .packages = c('data.table', 'dplyr')
    ) %dopar% {

    print(paste0('bootstrapping ', i, ' of ', nboot))

    # draw bootstrap sample
    tempdate = sample(datelist, nmonth, replace = TRUE) %>% # draw nmonth dates w/ replacement
        sort()
    tempkey = expand.grid(
        signalname = signallist, date = tempdate
    ) %>% as.data.table() # combine w/ signalnames
    setkey(tempkey, signalname, date)
    tempret = ret2[tempkey, on = c('signalname', 'date')] # merge
    tempret = tempret[!is.na(ret), ]

    # find t-stats
    tempsum = tempret[
        , .(tstat = mean(resid)/sd(resid)*sqrt(.N), nmonth = .N)
        , by = 'signalname'
    ] %>% mutate(name = 'resid') %>% 
    rbind(
        tempret[
            , .(tstat = mean(ret)/sd(ret)*sqrt(.N), nmonth = .N)
            , by = 'signalname'
        ] %>% mutate(name = 'ret')
    ) %>% 
    mutate(tabs = abs(tstat)) %>% 
    melt(id.vars = c('signalname', 'name', 'nmonth'), variable.name = 'stat')     
    
    # compress to histogram frequencies by (name,stat)
    temphist = tempsum %>% 
        group_by(name,stat) %>% 
        mutate(bin = cut(value, histedge, labels = FALSE)) %>% 
        count(name,stat,bin) %>% mutate(f = n/sum(n)) %>%
        ungroup() %>% 
        left_join(histinfo, by = 'bin') %>% 
        setDT()

    temphist$booti = i

    # save
    return(temphist)
} # end foreach i = 1:nboot
stopCluster(cl)
print(paste0('done bootstrap in ', difftime(Sys.time(), tic, units = 'min')
    , ' minutes'))

# summarize, add standard normal ----------------------------------------------

# summarize
histdat = boothist %>% 
    group_by(mid, left, right, name, stat) %>% 
    summarize(f_mean = mean(f), f_low=quantile(f, 0.05), f_high=quantile(f, 0.95)) %>% 
    ungroup() %>% 
    arrange(name, desc(stat), mid) 

# add normal
histdat = histdat %>% mutate(fnorm = if_else(
    stat=='tstat', pnorm(right) - pnorm(left)
    , (mid>0)*2*(pnorm(right) - pnorm(left)))) %>% 
    setDT()

# plot simpler -----------------------------------------------------
lablist = c('Bootstrapped Null', 'Normal(0,1)')

p = histdat[name=='resid' & stat=='tabs' & mid > 0] %>%
    select(mid, f_mean, fnorm) %>%
    pivot_longer(cols=c('f_mean', 'fnorm')) %>% 
    ggplot(aes(x = mid, y = value), size = 2) +
    geom_line(aes(color = name), size=0.5) +
    geom_point(aes(color = name, shape = name), size=3) +
    labs(x = 'Absolute t-statistic', y = 'Probability') +
    scale_x_continuous(breaks = seq(-8, 20, 1)) +
    coord_cartesian(xlim=c(0,4)) +
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
        , axis.line.x = element_line(color = 'black')
        , axis.line.y = element_line(color = 'black')
    ) +    
    scale_color_manual(values=c(MATBLUE, MATRED), labels = lablist) +
    scale_shape_manual(values=c(1, 4), labels = lablist)

ggsave('../results/boot-null-tabs.pdf', p, width = 5, height = 2.5, scale = 1
    , device = cairo_pdf)

