# super simple now (2022 04)


# SETUP ====
rm(list = ls())
source('0-functions.r')
library(tidyverse)
library(data.table)
load('../data/emp_data.Rdata')


tabs_cut = 2

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
      , group = factor(group, levels = c(1,2))
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




# NON-PARAMETRIC ====

# fdr calculations
F_yz = ecdf(yz_sum$tabs)
Pr_disc = 1-F_yz(tabs_cut)
# Pr_disc_F = 2*(1-pnorm(tabs_cut))
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
  geom_vline(xintercept = tabs_cut, color = NICERED) +
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
    geom = 'text', label = TeX(paste0('Pr($|t_i|>', tabs_cut, '$) = ', round(Pr_disc,2)))
    , x = 3.7, y = 2600
  ) +
  annotate(
    geom = 'text', label = TeX(paste0(
    'FDR $\\leq \\frac{5\\%}{', round(Pr_disc,2), '}$ = ',  round(fdrmax, 3)*100, '%'
    ))
    , x = 6.7, y = 1800
  )   

ggsave('../results/yz-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)



# Extrap ====

# fdr calculations
mean_all = 2.0 # 2.0 is just easier
Pr_disc = exp(-1/mean_all*tabs_cut)
# Pr_disc_F = 2*(1-pnorm(tabs_cut))
Pr_disc_F = 0.05
fdrmax = Pr_disc_F/Pr_disc

# make data for plotting
F_cz = ecdf(cz_sum$tabs)
n_cz = length(cz_sum$tabs)
edge_cz = seq(0,10,0.5)
F_fit = function(tabs) pexp(tabs, rate = 1/mean_all)

edge_fit = seq(0,10,0.1)
plotme = make_dist_dat(
  F_cz, edge_cz, F_fit, edge_fit, x_match = c(3.0,Inf), N1 = n_cz, showplot = T
)


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
    data = plotme %>% filter(group == 2), aes(color = group)
  ) +
  scale_color_manual(
    values = MATBLUE, labels = 'Extrapolated', name = NULL
  ) +
  geom_vline(xintercept = tabs_cut) +
  xlab(TeX('Absolute t-statistic ($|t_i|$)')) +
  ylab('Number of Strategies') +
  scale_x_continuous(
    breaks = seq(0,14,2)
  ) +
  # discovery line
  geom_vline(xintercept = tabs_cut, color = NICERED) +
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
    geom = 'text', label = TeX(paste0('Pr($|t_i|>', tabs_cut, '$) = ', round(Pr_disc,2)))
    , x = 42/10, y = 70
  ) +
  annotate(
    geom = 'text', label = TeX(paste0(
      'FDR $\\leq \\frac{', round(Pr_disc_F, 3)*100
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
    xlim = c(0,8)
  )


  

ggsave('../results/hlz-intuition.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)

