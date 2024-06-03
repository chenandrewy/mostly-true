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

ret = clz_ret %>% copy() # could use clz_vw or yz
# ret = clzvw_ret %>% copy()

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

datelist = unique(ret2$date)
signallist = unique(ret2$signalname)
histedge = c(-Inf, seq(-6, 6, 0.50), Inf)

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
    tempdate = sample(datelist, nmonth, replace = TRUE) # draw nmonth dates w/ replacement
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
    ]
    
    # compress to histogram frequencies
    temp = hist(tempsum$tstat, histedge, plot = FALSE)
    temp2 = hist(abs(tempsum$tstat), histedge, plot = FALSE)
    temphist = data.table(
        tmid = temp$mids, tleft=histedge[-length(histedge)], tright=histedge[-1]
        , f_t = temp$counts / sum(temp$counts), f_tabs = temp2$counts / sum(temp2$counts)
    )
    temphist$booti = i

    # save
    return(temphist)
} # end foreach i = 1:nboot
stopCluster(cl)
print(paste0('done bootstrap in ', difftime(Sys.time(), tic, units = 'min')
    , ' minutes'))

# summarize, add standard normal ----------------------------------------------

# reshape
boothistlong = melt(boothist, id.vars = c('tmid', 'tleft','tright','booti')
    , measure.vars = c('f_t', 'f_tabs')
    , variable.name = 'name', value.name = 'f') %>% 
    mutate(name = ifelse(name == 'f_t', 'boot_t', 'boot_tabs'))

# summarize
histdat = boothistlong %>% 
    group_by(tmid, tleft, tright, name) %>% 
    summarize(f_mean = mean(f), f_low=quantile(f, 0.05), f_high=quantile(f, 0.95)) %>% 
    ungroup() %>% 
    arrange(desc(name), tmid) 

# add normal
histdat = histdat %>% mutate(fnorm = if_else(
    name=='boot_t', pnorm(tright) - pnorm(tleft)
    , (tmid>0)*2*(pnorm(tright) - pnorm(tleft))))

# Plot tabs errorbars ---------------------------------------------------------

# plot bar with errorbars
p = histdat %>% filter(name == 'boot_tabs', tmid > 0) %>% 
    ggplot(aes(x = tmid, y = f_mean)) +
    # plot bootstrap
    geom_bar(stat='identity', position='identity', aes(fill = name), alpha=1) +
    geom_errorbar(aes(ymin = f_low, ymax = f_high), width = 0.1, color='black') +
    scale_fill_manual(values = MATRED, labels = 'Bootstrap', name = NULL) +
    # add normal
    geom_line(aes(x = tmid, y = fnorm, color = name)) +
    scale_color_manual(values = MATBLUE, labels = 'Normal(0,1)', name = NULL) +
    labs(x = 'Absolute t-statistic', y = 'Probability') +
    theme(legend.position = c(80,70)/100
        , legend.title = element_blank()
        , legend.margin = margin(t = -15, r = 20, b = 0, l = 5)) +
    scale_x_continuous(breaks = seq(-6, 6, 1)) +
    coord_cartesian(xlim=c(0, 1)*4) 

ggsave('../results/boot-null-tabs-bar.pdf', p, width = 5, height = 2.5, scale = 1
    , device = cairo_pdf)

# plot simpler -----------------------------------------------------
lablist = c('Bootstrap', 'Normal(0,1)')

p = histdat %>% filter(name == 'boot_tabs', tmid > 0) %>%  
    select(tmid, f_mean, fnorm) %>%
    pivot_longer(cols=c('f_mean', 'fnorm')) %>% 
    ggplot(aes(x = tmid, y = value), size = 2) +
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

