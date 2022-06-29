# super simple now (2022 04)


# SETUP ====
rm(list = ls())
source('0-functions.r')
library(tidyverse)
library(data.table)
load('../data/emp_data.Rdata')


h_disc = 2

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
    ) %>% 
    filter(!is.na(dF))
  
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
F_yz = ecdf(yz_sum$tabs)
Pr_disc = 1-F_yz(h_disc)
# Pr_disc_F = 2*(1-pnorm(h_disc))
Pr_disc_F = 0.05
fdrmax = Pr_disc_F/Pr_disc
n_yz = length(yz_sum$tabs)

# make data for plotting 
edge = seq(0,10,0.5)
plotme = make_dist_dat(F_yz, edge, F_yz, edge, N1 = n_yz)

# plot 
ggplot(plotme, aes(x=mids, y=dF)) +
  geom_bar(
    data = plotme %>% filter(group == 1)
    , stat = 'identity', position = 'identity', fill = 'gray'
    , aes(fill = group)
  ) +
  coord_cartesian(xlim = c(0,8)) +
  theme(
    legend.position = c(0.7, 0.7)
  ) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') +
  # discovery line =
  geom_vline(xintercept = h_disc, color = NICERED) +
  # write out intuition =
  geom_segment(
    aes(
      xend = 27/10, yend = 500
      , x = 3.2, y = 2300
    ),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(
    geom = 'text', label = TeX(paste0('Pr($|t_i|>', h_disc, '$) = ', round(Pr_disc,2)))
    , x = 3.7, y = 2600
  ) +
  annotate(
    geom = 'text', label = TeX(paste0(
    'FDR $\\leq \\frac{5\\%}{', round(Pr_disc,2), '}$ = ',  round(fdrmax, 3)*100, '%'
    ))
    , x = 6.7, y = 1800
  )   

ggsave('../results/yz-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)




# Extrapolated -----------------------------------------------------------

# empirical dist
F_cz = ecdf(cz_sum$tabs)

# hlz extrap
F_hlz = function(tabs) pgamma(tabs, shape = 1, scale = 2)

# gamma extrap
fit = est_trunc_gamma(cz_sum$tabs, tgood = 2.6, shape = 0.5)
F_gamma = function(tabs) pgamma(tabs, shape = fit$shape, scale = fit$scale)


## make plotting data ====
# define scaling
n_cz = length(cz_sum$tabs)
edge_cz = seq(0,10,0.5)
edge_fit = seq(0,10,0.1)

# scaled data
plotme1 = make_dist_dat(
  F_cz, edge_cz, F_hlz, edge_fit, x_match = c(3.0,Inf), N1 = n_cz, showplot = T
)
plotme2 = make_dist_dat(
  F_cz, edge_cz, F_gamma, edge_fit, x_match = c(3.0,Inf), N1 = n_cz, showplot = T
) 

# merge and clean
plotmeboth = plotme1 %>% 
  rbind(
    plotme2 %>% 
      filter(group == 2) %>% 
      mutate(group = 3)
  ) %>% 
  mutate(
    group = factor(group, levels =c(1,2,3))
  )

# numbers for annotations
Pr_disc = 1-F_hlz(h_disc)
fdrmax = 0.05/Pr_disc


# plot
groupdat = tibble(
  group = c(2, 3)
  , color = c(MATBLUE, MATRED)
  , labels = c('Exp','Gamma')
)

ggplot(data = plotmeboth, aes(x=mids, y=dF)) +
  geom_bar(
    data = plotmeboth %>% filter(group == 1)
    , stat = 'identity', position = 'identity'
    , aes(fill = group)
  ) +
  scale_fill_manual(
    values = 'gray', labels = 'Published', name = NULL
  ) +  
  geom_line(
    data = plotmeboth %>% filter(group %in% groupdat$group), aes(color = group)
  ) +
  scale_color_manual(
    values = groupdat$color, labels = groupdat$labels, name = NULL
  ) +
  geom_vline(xintercept = h_disc) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') +
  scale_x_continuous(
    breaks = seq(0,14,2)
  ) +
  # discovery line
  geom_vline(xintercept = h_disc, color = NICERED) +
  # write out intuition  
  geom_segment(
    aes(
      xend = 2.75, yend = 13
      , x = 3.0, y = 60
    ),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(
    geom = 'text', label = TeX(paste0('Pr($|t_i|>', h_disc, '$) = ', round(Pr_disc,2)))
    , x = 42/10, y = 70
  ) +
  annotate(
    geom = 'text', label = TeX(paste0(
      'FDR $\\leq \\frac{', round(0.05, 3)*100
      , '\\%}{'
      , round(Pr_disc,2), '}$ = ',  round(fdrmax, 3)*100, '%'
    ))
    , x = 68/10, y = 40
  ) + 
  theme(
    legend.position = c(75,75)/100
    , legend.margin = margin(t = -15, r = 20, b = 0, l = 5),
  ) +
  coord_cartesian(
    xlim = c(0,8), ylim = c(0, 120)
  )


  

ggsave('../results/hlz-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)





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
  F_cz, edge_cz, F_fit, edge_fit, x_match = c(3.0,Inf), N1 = n_cz, showplot = T
) %>% 
  mutate(
    group = factor(group, levels =c(1,2,3))
  )

# numbers for annotations
Pr_disc = 1-F_fit(h_disc)
fdrmax = 0.05/Pr_disc

# plot
ggplot(data = plotme, aes(x=mids, y=dF)) +
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
  geom_vline(xintercept = h_disc, color = NICERED) +
  # write out intuition  
  geom_segment(
    aes(
      xend = 2.75, yend = 13
      , x = 3.0, y = 60
    ),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(
    geom = 'text', label = TeX(paste0('Pr($|t_i|>', h_disc, '$) = ', round(Pr_disc,2)))
    , x = 42/10, y = 70
  ) +
  annotate(
    geom = 'text', label = TeX(paste0(
      'FDR $\\leq \\frac{', round(0.05, 3)*100
      , '\\%}{'
      , round(Pr_disc,2), '}$ = ',  round(fdrmax, 3)*100, '%'
    ))
    , x = 68/10, y = 40
  ) + 
  theme(
    legend.position = c(75,75)/100
    , legend.margin = margin(t = -15, r = 20, b = 0, l = 5),
  ) +
  coord_cartesian(
    xlim = c(0,8), ylim = c(0, 120)
  )




ggsave('../results/hlz-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)







# Extrapolated Gamma -----------------------------------------------------------

## define distributions

# empirical dist
F_cz = ecdf(cz_sum$tabs)

# fit gamma
fit = est_trunc_gamma(cz_sum$tabs, tgood = 2.6, shape = 0.5)
F_fit = function(tabs) pgamma(tabs, shape = fit$shape, scale = fit$scale)

## make plotting data 
# define scaling
n_cz = length(cz_sum$tabs)
edge_cz = seq(0,10,0.5)
edge_fit = seq(0,10,0.1)

# scaled data
plotme = make_dist_dat(
  F_cz, edge_cz, F_fit, edge_fit, x_match = c(3.0,Inf), N1 = n_cz, showplot = T
) %>% 
  mutate(
    group = factor(group, levels =c(1,2,3))
  )

# numbers for annotations
Pr_disc = 1-F_fit(h_disc)
fdrmax = 0.05/Pr_disc

# plot
ggplot(data = plotme, aes(x=mids, y=dF)) +
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
  geom_vline(xintercept = h_disc, color = NICERED) +
  # write out intuition  
  geom_segment(
    aes(
      xend = 2.75, yend = 13
      , x = 3.0, y = 60
    ),
    arrow = arrow(length = unit(0.03, "npc")),
    colour = "black", size = 0.1
  ) +  
  annotate(
    geom = 'text', label = TeX(paste0('Pr($|t_i|>', h_disc, '$) = ', round(Pr_disc,2)))
    , x = 42/10, y = 70
  ) +
  annotate(
    geom = 'text', label = TeX(paste0(
      'FDR $\\leq \\frac{', round(0.05, 3)*100
      , '\\%}{'
      , round(Pr_disc,2), '}$ = ',  round(fdrmax, 3)*100, '%'
    ))
    , x = 68/10, y = 40
  ) + 
  theme(
    legend.position = c(75,75)/100
    , legend.margin = margin(t = -15, r = 20, b = 0, l = 5),
  ) +
  coord_cartesian(
    xlim = c(0,8), ylim = c(0, 200)
  )




ggsave('../results/gamma-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)



