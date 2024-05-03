# super simple now (2022 04)
# outputs figures to ../results/
#   dm-intuition.pdf and hlz-intuition.pdf
# Setup -------------------------------------------------------------------

rm(list = ls())
source('0-functions.r')
load('../data/emp_data.Rdata')

# Data-Mining -------------------------------------------------------------

# fdr calculations
h_disc = 2
F_dm = ecdf(clz_sum$tabs)
Pr_disc = 1-F_dm(h_disc)
Pr_disc_F = 0.05
fdrmax = Pr_disc_F/Pr_disc
n_dm = length(clz_sum$tabs)

# make data for plotting 
edge = seq(0,10,0.5)
plotme = make_dist_dat(F_dm, edge, n_dm, F_dm, edge, n_dm)

# plot 
plt = ggplot(plotme, aes(x=mids, y=dF)) +
  geom_bar(
    data = plotme %>% filter(group ==1 )
    , stat = 'identity', position = 'identity', fill = 'gray'
    , aes(fill = group)
  ) +
  coord_cartesian(xlim = c(0,8)) +
  theme(
    legend.position = c(0.7, 0.7)
  ) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') +
  # discovery line 
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=5400, hjust = 0
    , label='Discoveries ->', color = MATRED) +
  # write out intuition 
  geom_segment(aes(xend = 27/10, yend = 500, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 3.7, y = 2600,
    label = TeX(paste0('Pr($|t_i|>', h_disc, '$) = ', round(Pr_disc,2)))    
  ) +
  annotate(geom = 'text', x = 6.7, y = 1800
    , label = TeX(paste0(
    'FDR $\\leq \\frac{5\\%}{', round(Pr_disc,2), '}$ = ',  round(fdrmax, 3)*100, '%'
    ))    
  )   

ggsave('../results/dm-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

## plots for slides ====

# data only
p1 = ggplot(plotme, aes(x=mids, y=dF)) +
  geom_bar(
    data = plotme %>% filter(group ==1 )
    , stat = 'identity', position = 'identity', fill = 'gray'
    , aes(fill = group)
  ) +
  coord_cartesian(xlim = c(0,8)) +
  theme(
    legend.position = c(0.7, 0.7)
  ) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') 

ggsave('../results/dm-anim-1.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

p2 = p1 + 
  # discovery line 
  geom_vline(xintercept = h_disc, color = MATRED) +
  # annotate(geom='text', x=2.1, y=5400, hjust = 0
  #   , label='Discoveries ->', color = MATRED) +
  # write out intuition 
  geom_segment(aes(xend = 27/10, yend = 1000, x = 3.3, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 3.7, y = 2600,
    label = paste0(round(Pr_disc*100,0), '%')
  )   
ggsave('../results/dm-anim-2.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

# Extrapolated HLZ ------------------------------------------------------

## define distributions
# empirical dist
F_cz = ecdf(cz_sum$tabs)
# hlz extrap
F_fit = function(tabs) pgamma(tabs, shape = 1, scale = 2)
# data mined
F_dm = ecdf(clz_sum$tabs)

## make plotting data
# define scaling
n_cz = length(cz_sum$tabs)
edge_cz = seq(0,10,0.5)
edge_fit = seq(0,10,0.1)
edge_dm = seq(0,10,0.5)

# scaled data
plotme1 = make_dist_dat(
  F_cz, edge_cz, n_cz, F_fit, edge_fit, n_cz, x_match = c(3.0,Inf), showplot = T
) 
plotme2 = make_dist_dat(
  F_cz, edge_cz, n_cz, F_dm, edge_dm, n_cz, x_match = c(3.0,Inf), showplot = F 
) %>% 
  filter(group==2) %>% 
  mutate(group=3)
plotme = plotme1 %>% rbind(plotme2) %>% 
  mutate(group = factor(group, levels = c(1,2,3)
  ))

# numbers for annotations
h_disc = 2
Pr_disc = 1-F_fit(h_disc)
fdrmax = 0.05/Pr_disc

groupdat = tibble(
  group = c('Exponential w/ Mean = 2','Data-Mined (CLZ)')
  , color = c(MATBLUE, MATRED)
  , linetype = c('solid', 'dashed')
)

# actual plot
plt = ggplot(data = plotme, aes(x=mids, y=dF)) +
  geom_bar(
    data = plotme %>% filter(group == 1)
    , stat = 'identity', position = 'identity'
    , aes(fill = group)
  ) +
  geom_line(
    data = plotme %>% filter(group %in% c(2,3))
    , aes(color = group, linetype = group)
  ) +
  scale_fill_manual(values = 'gray', labels = 'Published                ', name = NULL) +  
  scale_color_manual(values = groupdat$color, labels = groupdat$group, name = NULL) +
  scale_linetype_manual(values = groupdat$linetype, labels = groupdat$group, name = NULL)  +
  geom_vline(xintercept = h_disc) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') +
  scale_x_continuous(
    breaks = seq(0,14,2)
  ) +
  # discovery line
  geom_vline(xintercept = h_disc, color = 'black') +
  theme(
    legend.position = c(70,60)/100
    , legend.margin = margin(t = -15, r = 20, b = 0, l = 5),
  ) +
  coord_cartesian(
    xlim = c(0,8), ylim = c(0, 150)
  ) 

ggsave('../results/hlz-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

cz_sum %>% filter(tabs>3.0) %>% summarize(mean_t_gt_2 = mean(tabs)) %>% 
  mutate(scale = mean_t_gt_2 - 3)
clz_sum %>% filter(tabs>2.0) %>% summarize(mean_t_gt_2 = mean(tabs)) %>% 
  mutate(scale = mean_t_gt_2 - 2)

## plots for slides ====

plt = ggplot(data = plotme, aes(x=mids, y=dF)) +
  geom_bar(
    data = plotme %>% filter(group == 1)
    , stat = 'identity', position = 'identity'
    , aes(fill = group)
  ) +
  scale_fill_manual(
    values = 'gray', labels = 'Published', name = NULL
  ) +  
  geom_line(
    data = plotme %>% filter(group %in% groupdat$group), aes(color = group)
  ) +
  scale_color_manual(
    values = MATBLUE, labels = 'Extrapolated', name = NULL
  ) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') +
  scale_x_continuous(
    breaks = seq(0,14,2)
  ) +
  theme(
    legend.position = c(75,75)/100
    , legend.margin = margin(t = -15, r = 20, b = 0, l = 5),
  ) +
  coord_cartesian(
    xlim = c(0,8), ylim = c(0, 120)
  )

ggsave('../results/hlz-intuition-1.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

plt2 = plt + 
  geom_vline(xintercept = h_disc, color = MATRED) +
  # write out intuition  
  geom_segment(aes(xend = 2.75, yend = 13
                   , x = 3.8, y = 60),
               arrow = arrow(length = unit(0.03, "npc")),
               colour = "black", size = 0.1
  ) +  
  annotate(
    geom = 'text'
    , label = paste0(round(Pr_disc*100,1), '%')
    , x = 42/10, y = 70
  )

ggsave('../results/hlz-intuition-2.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

