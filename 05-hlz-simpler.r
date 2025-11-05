# ABOUTME: Provides a simplified HLZ simulation for two key comparison figures.
# Inputs:
#   - No external data files; uses random draws plus tidyverse/data.table/ggplot2/extrafont/latex2exp
# Outputs:
#   - results/hlz-simp-standard.pdf
#   - results/hlz-simp-post-truth.pdf
# How to run:
#   Rscript 05-hlz-simpler.r
#   Rscript 05-hlz-simpler.r --vanilla

# Setup -------------------------------------------------------------------
rm(list = ls())

# uncomment for pretty fonts
# install.packages('extrafont')
# extrafont::font_import()

library(here)
library(tidyverse)
library(data.table)
library(ggplot2)
library(extrafont)
library(latex2exp)

here::i_am("05-hlz-simpler.r")

results_dir <- here("results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATPURPLE = rgb(0.5800, 0.4000, 0.7400)

color3point0 = 'darkorchid'
label_3pt0 = '3.0 hurdle'


chen_theme =   theme_minimal() +
  theme(
    text = element_text(family = "Palatino Linotype")
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
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

# true / false labels (for v)
label_true = 'Alt: [Exp Ret]>0'
label_false =  'Null: [Exp Ret]=0'

nport = 1e4
nsim = 1000
n = nport*nsim

set.seed(121)

# hlz baseline
v  = runif(n) > 0.444 # Table 5, rho = 0.2, p_0
se = 1500/sqrt(12)/sqrt(240) # page 29, N = 240 and page 30, annual volatilty of 15% 
mu = rexp(n, 1/55.5); mu[!v] = 0 # Table 5, rho = 0.2, \lambda
theta = mu / se
theta_scatter = theta; theta_scatter[!v] = rnorm(sum(!v), 0, 0.1)
mu_scatter = mu; mu_scatter[!v] = rnorm(sum(!v), 0, 0.05)*se

# Fabian Winkler's fast method for common correlations
rho = 0.2
sigc = sqrt(rho)
sige = sqrt(1-sigc^2)

# simulate common component (in blocks)
simi = matrix(1:nsim, nrow = nsim, ncol = nport) %>% t() %>% as.vector
c = matrix(rnorm(nsim, 0, sigc), nrow = nsim, ncol = nport) %>% t() %>% 
  as.vector
# add idiosyncratic noise
e = rnorm(nport*nsim, 0, sige)
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
  simi, c, Z, theta, v, theta_scatter
) %>% 
  mutate(
    t = theta + Z
    , tabs = abs(t)
    , mu = mu
    , mu_scatter = mu_scatter
  ) %>% 
  mutate(
    v = factor(
      v, levels = c(TRUE,FALSE)
      , labels = c(label_true, label_false)
    )
  ) %>%   
  mutate(
    tselect = tabs
  )

# fit shrinkage
datsum = dat %>% 
  mutate(
    tgroup = ntile(tselect,1000)
  ) %>% 
  group_by(tgroup) %>% 
  summarize(
    tselect_left = min(tselect), tselect_right = max(tselect)
    , Etselect = mean(tselect), Etheta = mean(theta), n = dplyr::n()
    , nfalse = sum(v == label_false) 
  ) 

# fit "cdf"
datsum = datsum %>% 
  arrange(-Etselect) %>% 
  mutate(
    nfalse_cum = cumsum(nfalse)
    , n_cum = cumsum(n)
    , fdr_tselect_left = nfalse_cum / n_cum * 100
  ) %>% 
  arrange(Etselect)

# find hurdles
hurdle_05 = min(datsum$tselect_left[which(datsum$fdr_tselect_left < 5)])
hurdle_01 = min(datsum$tselect_left[which(datsum$fdr_tselect_left < 1)])
hurdle_3pt0 = 3.0

# Plot setup 
# find a subset for plotting
# nplot = 1500 # close to HLZ's estimate of "total factors" (ignoring unidentified scale)
nplot = 1378 #  # Table 5, rho = 0.2, M
set.seed(11)
small = dat[sample(1:n,nplot),]
small %>% filter(tabs > hurdle_01) %>% summarize(sum(v==label_false), n())

# plot sizing (shared by below)
texty = 250
textsize = 7
linesize = 1.1

# Baseline scatter for post-truth charts -----------------------------------

p2 = ggplot(small, aes(x = tselect, y = mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = c(MATBLUE, MATRED)) +
  coord_cartesian(xlim = c(-0.1, 10), ylim = c(-0.5, 300)) +
  scale_x_continuous(breaks = seq(-10, 20, 2)) +
  scale_y_continuous(breaks = seq(0, 500, 50)) +
  chen_theme +
  xlab(TeX('|t-statistic|')) +
  ylab(TeX('Expected Return (bps p.m.)')) +
  theme(legend.position = 'none')

# POST TRUTH --------------------------------------------------------------------

# Start with simplified hurdles chart -------------------------------------------

q1 = p2 +
  chen_theme +
  geom_vline(xintercept = hurdle_3pt0, size = linesize, color = color3point0, linetype = 'longdash') +
  annotate(
    geom = "text",
    label = label_3pt0,
    x = hurdle_3pt0, y = texty, vjust = -1,
    family = "Palatino Linotype", angle = 90, size = textsize, color = color3point0
  ) +
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  ) 

ggsave(file.path(results_dir, "hlz-simp-standard.pdf"), q1, width = 12, height = 8, device = cairo_pdf)  

# Relabel using post-truth -----------------------------------------------------------------

# label observations using post-truth method
lab_hlz_true = '"Alt": Significant'
lab_hlz_false = '"Null": Insignificant'

# post-truth labels
small = small %>% 
  mutate(
    v_hlz_3point0 = case_when(
      tabs > hurdle_3pt0 ~ TRUE,
      tabs <= hurdle_3pt0 ~ FALSE
    ),
    v_hlz_3point0 = factor(
      v_hlz_3point0, levels = c(TRUE, FALSE), labels = c(lab_hlz_true, lab_hlz_false)
    )
  )

# relabel 
q2_3point0 = ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v_hlz_3point0, color = v_hlz_3point0, shape = v_hlz_3point0), size = 2.5) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  theme(
    legend.position = c(.80, .15)
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  ) +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)")) +
  geom_vline(xintercept = hurdle_3pt0, size = linesize, color = color3point0, linetype = 'longdash') +
  annotate(
    geom = "text",
    label = label_3pt0,
    x = hurdle_3pt0, y = texty, vjust = -1,
    family = "Palatino Linotype", angle = 90, size = textsize, color = color3point0
  ) +
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  )
 
ggsave(file.path(results_dir, "hlz-simp-post-truth.pdf"), q2_3point0, width = 12, height = 8, device = cairo_pdf)    


# Some numbers for the paper ---------------------------------------------------
sig_share =  dat %>% filter(tabs>2) %>% summarize(sig_share = mean(tabs>2.27)) %>% 
  pull(sig_share)

sig_share*0.05  + (1-sig_share)

dat[tabs>2, ] %>% summarize(mean(v==label_false), n())

small[tabs>2.95]
small[tabs>2.27]
small[(tabs>2.27) & (v==label_false)]

small %>% 
  summarize(
    sum(tabs>2.95)
    , sum(tabs>2.27)
    , sum((tabs>2.27) & (v==label_false))
    , sum(tabs<3.8 & tabs>2)
    , sum(tabs>2)
    , sum((tabs>2) & (v==label_false))
  ) %>% 
  t()

-1*qnorm(0.05/296/2)

dat %>% 
  summarize(
    sum((tabs<3.8)&(tabs>2))/sum(tabs>2) 
  )

dat %>% 
  summarize(
    sum((tabs<3.0)&(tabs>2))/sum(tabs>2) 
  )
