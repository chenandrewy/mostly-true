# Setup --------------------------------------------------------

rm(list = ls())
source('0-functions.r')
load('../data/emp_data.Rdata')
load('../data/bootact.Rdata')

## Make data for plotting ------------------
# fdr calculations
h_disc = 2
F_dm = ecdf(clz_sum$tabs)
Pr_disc = 1-F_dm(h_disc)
Pr_disc_F = 0.05
fdrmax = Pr_disc_F/Pr_disc
n_dm = length(clz_sum$tabs)

# plot settings
color_emp = 'gray50'
color_null = MATRED
bar_alpha = 0.65
ylimnum = c(0, 12000)
discovery_y = 11000
intuition_y = 3000
yticks = seq(0,12000,2000)

# storey
F_null = function(tabs){(2*(pnorm(tabs)-0.5))}
t_null_max = 0.5
pF = F_dm(t_null_max)/F_null(t_null_max)
pF # check...

# make data for plotting 
edge = seq(0,10,0.5)
temp1 = make_dist_dat(F_dm, edge, n_dm, F_null, edge, n_dm*pF, x_match = NULL) %>% 
  mutate(group=factor(group,levels=c(1,2),labels=c('emp','null')))
temp2 = make_dist_dat(F_dm, edge, n_dm, F_null, edge, n_dm*1, x_match = NULL) %>% 
  mutate(group=factor(group,levels=c(1,2),labels=c('emp','null_ez')))
plotme = temp1 %>% rbind(temp2 %>% filter(group=='null_ez'))

## add bootstrap data ----------------------
bootact[ , tstat:=mean/vol*sqrt(nmonth)]

boothist = bootact %>%
    # histogram counts within each bootstrap
    mutate(bin = cut(abs(tstat), breaks = edge, include.lowest = TRUE)) %>%
    group_by(booti, bin) %>%
    summarize(count = n(), .groups = 'drop')  %>% 
    # normalize by total signals
    left_join(bootact %>% group_by(booti) %>% summarize(ntotal=n()), by = "booti") %>% 
    # select definition of dF
    mutate(dF=count)  %>% 
    # find bin midpoints 
    mutate(bin = str_remove_all(bin, "[\\(\\)\\[\\]]")
        , left = str_split(bin, ",")
        , left = sapply(left, function(x) as.numeric(x[1]))
        , right = str_split(bin, ",")
        , right = sapply(right, function(x) as.numeric(x[2]))
        , mids = (left + right) / 2
    )

bootstat = boothist  %>% 
  group_by(mids)  %>% 
  summarize(boot05 = quantile(dF, 0.05),
            boot95 = quantile(dF, 0.95),
            bootsd = sd(dF),
            .groups = 'drop')  %>% 
  mutate(group='emp') %>% 
  left_join(plotme %>% select(mids,group, dF),
  by=c('group','mids')) 

# define error bars now
bootstat = bootstat  %>% 
  mutate(dF_lo = dF-2*bootsd,
  dF_hi = dF+2*bootsd)

plotme_err = plotme %>% 
  left_join(bootstat  %>% select(group, mids, starts_with('dF_')),
  by = c('group','mids'))

# plots for slides ---------------------------------------

# label explicitly for clarity in talks
lab_dat = 'Data: CLZ\'s 29,000 EW Returns'

discovery_y = 11000

plotme20min = plotme_err

plotme20min = plotme20min %>% 
  rbind(
    plotme20min[group=='null'] %>% 
    mutate(dF = dF/2, group = 'null_small')
  )

# plot axes only 
plt = ggplot(plotme20min[group!='null_ez'] %>% mutate(dF = 0)
  , aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +  
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') 

ggsave('../results/dm-viz-axes-only.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# plot data only (need to introduce CLZ for short slides)
# plotme20min[group!='null_ez'] %>% mutate(dF = if_else(group=='null',0,dF))
plt = ggplot(plotme20min[group=='emp'], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +  
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    , labels=c(lab_dat, paste0('Null Component Bound: N(0,1)'))    
    , name='') 

ggsave('../results/dm-viz-data-only.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# plot a null that's too small
plt = ggplot(plotme20min[group%in%c('emp','null_small')], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +  
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    , labels=c(lab_dat, paste0('Null Component Bound: N(0,1)'))    
    , name='')
ggsave('../results/dm-viz-null-small.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)    

# plot storey null
plt = ggplot(plotme20min[group%in%c('emp','null')], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +  
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    , labels=c(lab_dat, paste0('Null Component Bound: N(0,1)'))    
    , name='') 

ggsave('../results/dm-viz-storey-color-0.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

plt2 = plt +     
  # discovery line 
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=discovery_y, hjust = 0
    , label='Discoveries ->', color = MATRED) +
  theme(legend.position = c(70,70)/100) +
  # write out intuition
  geom_segment(aes(xend = 22/10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  )   +
  annotate(geom = 'text', x = 33/10, y = intuition_y, hjust = 0,
    , label = TeX(paste0(
    '$FDR_{|t|>2}\\leq$ red / total = 9%'
    ))    
  ) 

ggsave('../results/dm-viz-storey-color.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# plot ez null
plt = ggplot(plotme20min[group%in%c('emp','null_ez')], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +  
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)

    , labels=c(lab_dat, paste0('Null Component Bound: N(0,1)'))    
    , name='') 

ggsave('../results/dm-viz-ez-color-0.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

plt2 = plt +     
  # discovery line 
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=discovery_y, hjust = 0
    , label='Discoveries ->', color = MATRED) +
  theme(legend.position = c(70,70)/100) +
  # write out intuition
  geom_segment(aes(xend = 22/10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  )   +
  annotate(geom = 'text', x = 33/10, y = intuition_y, hjust = 0,
    , label = TeX(paste0(
    '$FDR_{|t|>2}\\leq$ red / total = 15%'
    ))    
  ) 

ggsave('../results/dm-viz-ez-color.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)


# plot with error bars -----------------------------------------

# plot storey null
plt = ggplot(plotme_err[group!='null_ez'], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +  
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Predictors') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    # , labels=c('Data', paste0('False = ', round(pF*100,0), '% of Data'))
    , labels=c('Data: CLZ EW', paste0('Null Component Bound'))    
    , name='') +     
  # error bars
  geom_errorbar(aes(ymin = dF_lo, ymax=dF_hi), width = 0.15, color = 'grey30', size=0.3) +
  # discovery line 
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=discovery_y, hjust = 0
    , label='Discoveries ->', color = MATRED) +
  theme(legend.position = c(80,80)/100) +
  # write out intuition
  geom_segment(aes(xend = 22/10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 33/10, y = intuition_y, hjust = 0,
    , label = TeX(paste0(
    'FDR $\\leq \\frac{5\\%}{', round(Pr_disc,2), '}$ ('
      , round(pF,2), ') = '
      , round(fdrmax*pF*100, 0), '%'
    ))    
  )   

ggsave('../results/dm-viz-storey-err.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)



# plot with formulas ====

# plot storey null
plt = ggplot(plotme[group!='null_ez'], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +  
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Predictors') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    # , labels=c('Data', paste0('False = ', round(pF*100,0), '% of Data'))
    , labels=c('Data: CLZ EW', paste0('Null Component Bound'))    
    , name='') +     
  # discovery line 
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=discovery_y, hjust = 0
    , label='Discoveries ->', color = MATRED) +
  theme(legend.position = c(80,80)/100) +
  # write out intuition
  geom_segment(aes(xend = 22/10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 33/10, y = intuition_y, hjust = 0,
    , label = TeX(paste0(
    'FDR $\\leq \\frac{5\\%}{', round(Pr_disc,2), '}$ ('
      , round(pF,2), ') = '
      , round(fdrmax*pF*100, 0), '%'
    ))    
  )   

ggsave('../results/dm-viz-storey.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# plot easy null
plt = ggplot(plotme[group!='null'], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Predictors') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    # , labels=c('Data', paste0('False = ', round(1*100,0), '% of Data'))
    , labels=c('Data: CLZ EW', paste0('Null Component Bound'))    
    , name='') +     
  # discovery line 
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=discovery_y, hjust = 0
    , label='Discoveries ->', color = MATRED) +
  theme(legend.position = c(80,80)/100) +
  # write out intuition
  geom_segment(aes(xend = 22/10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 33/10, y = intuition_y, hjust = 0
    , label = TeX(paste0(
    'FDR $\\leq \\frac{5\\%}{', round(Pr_disc,2), '}$ ('
      ,'1.00) = '
      , round(fdrmax*1*100, 0), '%'
    ))    
  ) 

ggsave('../results/dm-viz-ez.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

