# Setup -------------------------------------------------------------------

rm(list = ls())
source('0-functions.r')
load('../data/emp_data.Rdata')

# Data-Mining with Storey w/ EZ benchmark -------------------------------

## plot setup ====
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
intuition_y = 3600
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

## plot start ====
# plot easy null
plt = ggplot(plotme[group!='null'], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    # , labels=c('Data', paste0('False = ', round(1*100,0), '% of Data'))
    , labels=c('Data', paste0('False Upper Bound'))    
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

# plot storey null
plt = ggplot(plotme[group!='null_ez'], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8), ylim=ylimnum) +
  scale_y_continuous(breaks = yticks) +  
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    # , labels=c('Data', paste0('False = ', round(pF*100,0), '% of Data'))
    , labels=c('Data', paste0('False Upper Bound'))    
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

## plots for slides ---------------------------------------
# (zoomed in)
# uses plot data from previous block
discovery_y = 6000

# plot data
plt = ggplot(plotme[group=='emp'], aes(x=mids, y=dF))  +
  coord_cartesian(xlim = c(0,8), ylim=c(0,7000)) +
  theme(legend.position = c(80,80)/100) +  
  scale_y_continuous(breaks = yticks) +  
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    , labels=c('CLZ Data-Mined EW', paste0('Normal(0,1)'))    
    , name='') 
ggsave('../results/storey-anim-1.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# plot draw a null that's too small
tempdat = plotme[group%in%c('emp','null')] %>% 
  mutate(dF=if_else(group=='null', 0.5*dF, dF))
plt = ggplot(tempdat, aes(x=mids, y=dF))  +
  coord_cartesian(xlim = c(0,8), ylim=c(0,7000)) +
  theme(legend.position = c(80,80)/100) +  
  scale_y_continuous(breaks = yticks) +  
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    , labels=c('CLZ Data-Mined EW', paste0('Normal(0,1)'))    
    , name='') 
ggsave('../results/storey-anim-2a.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# draw the ez null (too big)
plt_ez = ggplot(plotme[group%in%c('emp','null_ez')], aes(x=mids, y=dF))  +
  coord_cartesian(xlim = c(0,8), ylim=c(0,7000)) +
  theme(legend.position = c(80,80)/100) +  
  scale_y_continuous(breaks = yticks) +  
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    , labels=c('CLZ Data-Mined EW', paste0('Normal(0,1)'))    
    , name='') 
ggsave('../results/storey-anim-2b.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# draw the good null
plt_good = ggplot(plotme[group%in%c('emp','null')], aes(x=mids, y=dF))  +
  coord_cartesian(xlim = c(0,8), ylim=c(0,7000)) +
  theme(legend.position = c(80,80)/100) +  
  scale_y_continuous(breaks = yticks) +  
  xlab(TeX('Absolute t-statistic')) +
  ylab('Number of Strategies') +  
  # bars 
  geom_bar(stat='identity', position='identity', alpha=bar_alpha
    , aes(fill=group)) +
  scale_fill_manual(values=c(color_emp, color_null)
    , labels=c('CLZ Data-Mined EW', paste0('Normal(0,1)'))    
    , name='') 
ggsave('../results/storey-anim-2c.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# draw discoveries line
plt_good2 = plt_good +
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=discovery_y, hjust = 0
    , label='Discoveries ->', color = MATRED) 
ggsave('../results/storey-anim-3a.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# add FDR text annotation
plt_good3 = plt_good2 +
  geom_segment(aes(xend = 22/10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 33/10, y = 2500, hjust = 0,
    , label = TeX(paste0(
    'FDR $\\leq$ red / gray = 9%'
    ))    
  )   

ggsave('../results/storey-anim-3b.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# plot ez null with annotation
plt_good3c = plt_ez +
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=discovery_y, hjust = 0
           , label='Discoveries ->', color = MATRED)    +
  geom_segment(aes(xend = 24/10, yend = 380, x = 3.2, y = 2300),
               arrow = arrow(length = unit(0.03, "npc")),
               colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 33/10, y = 2500, hjust = 0,
           , label = TeX(paste0(
             'FDR $\\leq$ red / gray = 15%'
           ))    
  )   
ggsave('../results/storey-anim-3c.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# use FDR annotation that relates visual to EZ bound
plt_good2 +
  # write out intuition
  geom_segment(aes(xend = 22/10, yend = 250, x = 3.2, y = 2300),
               arrow = arrow(length = unit(0.03, "npc")),
               colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 33/10, y = 2800, hjust = 0,
           , label = TeX(paste0(
             'FDR $\\leq \\frac{5\\%}{', round(Pr_disc,2), '}$ ('
             , round(pF,2), ') = '
             , round(fdrmax*pF*100, 0), '%'
           ))    
  )   

ggsave('../results/storey-anim-3b-alt.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)


