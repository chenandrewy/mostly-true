# New in 2024
# Outputs ../results/hlz-pub-scatter.pdf
# doesn't require any data, can be run separately, in principle.

# Setup -------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(data.table)
library(ggplot2)
library(extrafont)
library(latex2exp)

dir.create('../results/')

MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)
MATPURPLE = rgb(0.4940, 0.1840, 0.5560)


chen_theme =   theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1)
    
    # Font sizes
    , axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 26),
    axis.text.y = element_text(size = 26),
    legend.text = element_text(size = 24),
    
    # Tweaking legend
    legend.position = c(.80, .15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 20, b = 5, l = 5), 
    legend.key.size = unit(1.5, "cm")
    , legend.title = element_blank()    
  ) 

# Simulate hlz ----------------------------------------------------------------

nfactorall = 1378 # Table 5, rho = 0.2, M
nsim = 1e4
n = nfactorall*nsim

set.seed(121)

# hlz baseline
v  = runif(n) > 0.444 # Table 5, rho = 0.2, p_0
se = 1500/sqrt(12)/sqrt(240) # page 29, N = 240 and page 30, annual volatilty of 15% 
mu = rexp(n, 1/55.5); mu[!v] = 0 # Table 5, rho = 0.2, \lambda
theta = mu / se
theta_scatter = theta; theta_scatter[!v] = rnorm(sum(!v), 0, 0.1)
mu_scatter = mu; mu_scatter[!v] = rnorm(sum(!v), 0, 0.1)*se
pubnoise = runif(n)

# Fabian Winkler's fast method for common correlations
rho = 0.2
sigc = sqrt(rho)
sige = sqrt(1-sigc^2)

# simulate common component (in blocks)
c = matrix(rnorm(nsim, 0, sigc), nrow = nsim, ncol = nfactorall) %>% t() %>% 
  as.vector
# add idiosyncratic noise
e = rnorm(nfactorall*nsim, 0, sige)
Z = c + e

# sanity check
# c2 = rnorm(n, 0, sigc)
# z1 = c2 + rnorm(n, 0, sige)
# z2 = c2 + rnorm(n, 0, sige)
# cor(z2,z1)
# var(z1)
# var(z2)

# assemble into data table
dat = data.table(
  Z, theta, v, pubnoise, theta_scatter
) %>% 
  mutate(
    t = theta + Z
    , tabs = abs(t)
    , pub = case_when(
      tabs <= 1.96 ~ FALSE
      , (1.96 < tabs) & (tabs <= 2.57) ~ pubnoise < 0.5
      , 2.57 < tabs  ~ TRUE
    )
    , mu = mu
    , mu_scatter = mu_scatter
  ) %>% 
  mutate(
    tselect = tabs
  )


# Significance tests using direct quotes ------------------------------------

# we only really need to compute the distribution of t-stats and nulls
# the significance tests we can pull the info from the text

h_bonf = 3.78 # page 24, bottom
h_by13 = 3.39 # page 26, top
h_by13_5 = 2.78 # page 26, top

hurdle_dat = tibble(
    name  = c('Bonferroni', 'BY13', 'BY13_5')
    , h = c(h_bonf, h_by13, h_by13_5)
)

npub = 316 # page 8 and Fig 3 (page 25) (I guess we don't use this one)
npubsig = 296 # conclusion, page 37

# keep only published (as in Fig 3 and conclusion)
datpub = dat[pub == TRUE, .(v, mu_scatter, tabs)]

set.seed(1)
small = datpub[sample(nrow(datpub), npubsig), ]  %>% 
  mutate(
      v = factor(
        v, levels = c(TRUE,FALSE)
        , labels = c('True Factor (Non-Null)', 'False Factor (Null)')
      )
    ) 

# plot settings
texty = 180
textsize = 7
xlimnum = c(1.6,10)
ylimnum = c(-0.5,300)
vjustnum = -0.5
hjustnum = 0

color_h = c(MATRED, MATPURPLE, 'black')
p1 = ggplot(small, aes(x=tabs, y=mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 3) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  coord_cartesian(xlim = xlimnum, ylim = ylimnum) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  theme(
    legend.position = c(.75, .15)
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  ) +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)")) +
  geom_vline(xintercept = h_bonf, linetype = 'solid', color = color_h[1]) +
  annotate(geom="text", label="Bonferroni (5%)"
    , x=h_bonf, y=texty, vjust=vjustnum, hjust=hjustnum
    , angle = 90, size = textsize, color = color_h[1]) +
  geom_vline(xintercept = h_by13, linetype = 'dashed', color = color_h[2]) +
  annotate(geom="text", label="BY Theorem 1.3 (1%)"
    , x=h_by13, y=texty, vjust=vjustnum, hjust=hjustnum
    , angle = 90, size = textsize, color = color_h[2])  +
  geom_vline(xintercept = h_by13_5, linetype = 'twodash', color = color_h[3]) +
  annotate(geom="text", label="BY Theorem 1.3 (5%)"
    , x=h_by13_5, y=texty, vjust=vjustnum, hjust=hjustnum
    , angle = 90, size = textsize, color = color_h[3]) 

ggsave('../results/hlz-pub-scatter.pdf', p1, width = 12, height = 8)


# check to console

# HLZ's Figure 3 implied FDR
sum(small$tabs > h_by13_5)
sum(datpub$tabs > h_by13_5)/nrow(datpub)
sum(small$tabs > h_by13_5)/nrow(small)*0.05+(1-sum(small$tabs > h_by13_5)/nrow(small))

# HLZ's SMM FDR overall
sum(small$v == 'False Factor (Null)')
sum(small$v == 'False Factor (Null)')/nrow(small)
sum(datpub$v == FALSE)/nrow(datpub)


sum(small$tabs < h_bonf)
sum(small$tabs < h_by13)
sum(small$tabs < h_by13_5)

sum(small$v == 'False Factor (Null)')

sum(small$v == 'False Factor (Null)')/nrow(small)

