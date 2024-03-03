# super simple now (2022 04)
# outputs figures to ../results/
#   dm-intuition.pdf and hlz-intuition.pdf

# Setup -------------------------------------------------------------------

rm(list = ls())
source('0-functions.r')
load('../data/emp_data.Rdata')

# creates data for comparing cdf F1 to cdf F2 in a plot
# automatically adjusts for different x-binning
make_dist_dat = function(F1, edge1, N1, F2, edge2, N2
  , x_match = c(-Inf,Inf), showplot = F){
  
  # adjust for different x-binning
  if (!is.null(x_match)){
    rescale_fac = diff(F1(x_match))/diff(F2(x_match)) * diff(edge1)[1] /diff(edge2)[1]
  } else {
    rescale_fac = 1
  }
  
  # make histogram counts, with normalization adjustments
  dat = tibble(
    edge = edge1, F = N1*F1(edge1), group = 1
  ) %>% 
    rbind(
      tibble(
        edge = edge2, F = N2*F2(edge2)*rescale_fac, group = 2
      )
    ) %>% 
    # take first differences, find midpoints
    group_by(group) %>% 
    mutate(
      F = F
      , dF = F - lag(F)
      , mids = 0.5*(edge + lag(edge))
    ) %>% 
    filter(!is.na(dF)) %>% 
    setDT()
  
  if (showplot) {
    dat %>% 
      ggplot(aes(x=edge, y=dF)) +
      geom_line(aes(color = group))
  }
  
  return(dat)
  
} # make_dist_dat

## set theme ====
theme_set(
  theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype")
    )
)

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

# Extrapolated HLZ -----------------------------------------------------------

## define distributions
# empirical dist
F_cz = ecdf(cz_sum$tabs)
# hlz extrap
F_fit = function(tabs) pgamma(tabs, shape = 1, scale = 2)

## make plotting data 
# define scaling
n_cz = length(cz_sum$tabs)
edge_cz = seq(0,10,0.5)
edge_fit = seq(0,10,0.1)

# scaled data
plotme = make_dist_dat(
  F_cz, edge_cz, n_cz, F_fit, edge_fit, n_cz, x_match = c(3.0,Inf), showplot = T
) %>% 
  mutate(
    group = factor(group, levels =c(1,2,3))
  )

# numbers for annotations
Pr_disc = 1-F_fit(h_disc)
fdrmax = 0.05/Pr_disc

groupdat = tibble(
  group = c(2, 3)
  , color = c(MATBLUE, MATRED)
  , labels = c('Exp','Gamma')
)

# plot
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
  geom_vline(xintercept = h_disc) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') +
  scale_x_continuous(
    breaks = seq(0,14,2)
  ) +
  # discovery line
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=112, hjust = 0
    , label='Discoveries ->', color = MATRED) +
  # write out intuition  
  geom_segment(aes(xend = 2.75, yend = 13
                  , x = 3.0, y = 60),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(
    geom = 'text', label = TeX(paste0('Pr($|t_i|>', h_disc, '$) = ', round(Pr_disc,2)))
    , x = 42/10, y = 70
  ) +
  annotate(geom = 'text', x = 68/10, y = 40
      , label = TeX(paste0(
      'FDR $\\leq \\frac{', round(0.05, 3)*100
      , '\\%}{'
      , round(Pr_disc,2), '}$ = ',  round(fdrmax, 3)*100, '%'
    ))    
  ) + 
  theme(
    legend.position = c(75,75)/100
    , legend.margin = margin(t = -15, r = 20, b = 0, l = 5),
  ) +
  coord_cartesian(
    xlim = c(0,8), ylim = c(0, 120)
  )

ggsave('../results/hlz-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

cz_sum %>% filter(tabs>3.0) %>% summarize(mean_t_gt_2 = mean(tabs)) %>% 
  mutate(scale = mean_t_gt_2 - 3)
clz_sum %>% filter(tabs>2.0) %>% summarize(mean_t_gt_2 = mean(tabs)) %>% 
  mutate(scale = mean_t_gt_2 - 2)
# Data-Mining with Storey ----------------------------------------

# fdr calculations
h_disc = 2
F_dm = ecdf(clz_sum$tabs)
Pr_disc = 1-F_dm(h_disc)
Pr_disc_F = 0.05
fdrmax = Pr_disc_F/Pr_disc
n_dm = length(clz_sum$tabs)

# storey
F_null = function(tabs){(2*(pnorm(tabs)-0.5))}
t_null_max = 0.5
pF = F_dm(t_null_max)/F_null(t_null_max)
pF # check...

# make data for plotting 
edge = seq(0,10,0.5)
plotme = make_dist_dat(F_dm, edge, n_dm, F_null, edge, n_dm*pF, x_match = NULL) %>% 
  mutate(group=factor(group,levels=c(1,2),labels=c('emp','null')))

# plot try 2
plt = ggplot(plotme[group=='emp'], aes(x=mids, y=dF)) +
  coord_cartesian(xlim = c(0,8)) +
  theme(legend.position = c(0.7, 0.7)) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') +
  # emp bars 
  geom_bar(stat='identity', position='identity', fill='gray') +
  # discovery line 
  geom_vline(xintercept = h_disc, color = MATRED) +
  annotate(geom='text', x=2.1, y=5400, hjust = 0
    , label='Discoveries ->', color = MATRED) +
  # null bars
  geom_bar(data=plotme[group=='null']
    , stat='identity', position='identity', alpha = 0.7
    , aes(fill=group)) +
  scale_fill_manual(values=c('null'=MATRED), labels=c('null'='False (56% of Total)')
    , name='') +
  theme(legend.position = c(80,80)/100) +
  # write out intuition
  geom_segment(aes(xend = 22/10, yend = 250, x = 3.2, y = 2300),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(geom = 'text', x = 33/10, y = 2600, hjust = 0,
    , label = TeX(paste0(
    'FDR $\\leq \\frac{5\\%}{', round(Pr_disc,2), '}('
      , round(pF,2),')$ = '
      , round(fdrmax*pF*100, 0), '%'
    ))    
  )   

ggsave('../results/dm-storey.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

