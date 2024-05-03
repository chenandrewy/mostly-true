# 2024 05 bootstrap check of the standard normal null

# Setup -----------------------------------------------------------------------
rm(list=ls())
source('0-functions.r')
load('../data/emp_data.Rdata')

ret = clz_ret %>% copy()

# define sample period --------------------------------------------------------
# tbc: unify with sim-theory-free

# using minimum observations
min_obs_pct = 80

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
# for nboot = 1,000 takes about 7 minutes on 4 cores

nboot = 1000
nmonth = length(unique(ret2$date))
ncores = 4

datelist = unique(ret2$date)
signallist = unique(ret2$signalname)
histedge = c(-Inf, seq(-6, 6, 0.25), Inf)

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
    temphist = data.table(tmid = temp$mids, f = temp$counts / sum(temp$counts))
    temphist$booti = i

    # save
    return(temphist)
} # end foreach i = 1:nboot
stopCluster(cl)
print(paste0('done bootstrap in ', difftime(Sys.time(), tic, units = 'min')
    , ' minutes'))

# plot tstat dist -----------------------------------------------------

# add standard normal
histdat = boothist %>% 
    group_by(tmid) %>% 
    summarize(fboot = mean(f)) 
histdat$tleft = histedge[-length(histedge)]    
histdat$tright = histedge[-1]
histdat = histdat %>% mutate(fnorm = pnorm(tright) - pnorm(tleft))

lablist = c('Bootstrap', 'Normal(0,1)')

p = histdat %>% 
    pivot_longer(cols = c('fboot', 'fnorm')) %>%
    ggplot(aes(x = tmid, y = value)) +
    geom_line(aes(color = name), size=0.5) +
    geom_point(aes(color = name, shape = name), size=3) +
    scale_x_continuous(breaks = seq(-8, 8, 1)) +
    coord_cartesian(xlim=c(-1, 1)*4) +
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
    labs(x = 't-statistic', y = 'Probability') +
    scale_color_manual(values=c(MATBLUE, MATRED), labels = lablist) +
    scale_shape_manual(values=c(1, 4), labels = lablist)

ggsave('../results/boot-null-t.pdf', p, width = 5, height = 4, scale = 1, 
    device = cairo_pdf)

# plot abs(tstat) dist -----------------------------------------------------

# calculate tabs dist
histdat_abs = histdat %>% 
    mutate(tabsmid = abs(tmid)) %>% 
    arrange(tabsmid) %>% 
    group_by(tabsmid) %>%
    summarize(fboot = sum(fboot), fnorm = sum(fnorm)) 

p = histdat_abs %>%
    pivot_longer(cols = c('fboot', 'fnorm')) %>%
    ggplot(aes(x = tabsmid, y = value), size = 2) +
    geom_line(aes(color = name), size=0.5) +
    geom_point(aes(color = name, shape = name), size=3) +
    scale_x_continuous(breaks = seq(-8, 20, 1)) +
    coord_cartesian(xlim=c(0,5)) +
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
    labs(x = 't-statistic', y = 'Probability') +
    scale_color_manual(values=c(MATBLUE, MATRED), labels = lablist) +
    scale_shape_manual(values=c(1, 4), labels = lablist)


ggsave('../results/boot-null-tabs.pdf', p, width = 10, height = 5, scale = 0.5)
