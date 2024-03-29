# Setup -------------------------------------------------------------------
rm(list = ls())

# uncomment for pretty fonts
# install.packages('extrafont')
# extrafont::font_import()

library(tidyverse)
library(data.table)
library(ggplot2)
library(extrafont)
library(latex2exp)

dir.create('../results/')

MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

chen_theme =   theme_minimal() +
  theme(
    text = element_text(family = "Palatino Linotype")
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
    
    # Font sizes
    # , axis.title.x = element_text(size = 26),
    # axis.title.y = element_text(size = 26),
    # axis.text.x = element_text(size = 22),
    # axis.text.y = element_text(size = 22),
    # legend.text = element_text(size = 18),

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
label_true = 'True: [Exp Ret]>0'
label_false =  'False: [Exp Ret]=0'

nport = 1e4
nsim = 500
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
  Z, theta, v, theta_scatter
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
hurdle_bonf05 = qnorm(1-0.05/300/2) # assumes everything is published, as in HLZ's conclusion text


# Plot setup ----------------------------------------------

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

# Plot slowly ----------------------------------------------------

p0 = ggplot(small %>% 
              filter(v==label_false), aes(x=tselect,y=mu_scatter)) +
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values=c(MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)")) +
  theme(legend.position = 'none')

ggsave('../results/slow-0.pdf', width = 12, height = 8, device = cairo_pdf)

p1 = ggplot(small %>% 
              filter(v==label_false), aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) +
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values=c(MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)")) +
  theme(legend.position = 'none')

ggsave('../results/slow-1.pdf', width = 12, height = 8, device =cairo_pdf)

p2 =  ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) + 
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)")) +
  theme(legend.position = 'none')

ggsave('../results/slow-2.pdf', width = 12, height = 8, device =cairo_pdf)

p3 = p2 + chen_theme

ggsave('../results/slow-3.pdf', width = 12, height = 8, device =cairo_pdf)

# plot FDR 1% line --------------------------------------------------------

p4 = p3 +
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash') +  
  annotate(geom="text", 
           label="Discoveries ->", 
           x=3+1.1, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 0, size = textsize, color = 'black'
  ) 

ggsave('../results/slow-4.pdf', width = 12, height = 8, device =cairo_pdf)

p5 = p4 +  
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) +
  annotate("rect", xmin = 3, xmax = 3.5, 
           ymin = -10, ymax = 5, size=2,
           color='orange', fill=NA) +
  annotate("segment", x=4.2, xend=3.5,
           y=12, yend=5, size=2,
           color='orange',
           arrow=arrow(type = "closed", 
                       length = unit(0.2, "inches"))) 
  annotate("text", x=39/10, y=20, hjust=0,
           color='orange', size=7,
           label='False Discoveries', family = "Palatino Linotype") 

ggsave('../results/slow-5.pdf', width = 12, height = 8, device =cairo_pdf)

# plot FDR 5% line ------------------------------------------------------------
p6 = p5 + 
  geom_vline(xintercept = 2.27, size = linesize, color=MATYELLOW
             , linetype='longdash') 

ggsave('../results/slow-6.pdf', width = 12, height = 8, device =cairo_pdf)

p7 = p6 +
  annotate(geom="text", label="FDR = 5%", 
           x=2.40, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATYELLOW
  ) + 
  annotate("rect", xmin = 2.3, xmax = 3.5, 
           ymin = -10, ymax = 5, size=2,
           color='orange', fill=NA) +
  annotate("segment", x=4.2, xend=3.5,
           y=12, yend=5, size=2,
           color='orange',
           arrow=arrow(type = "closed", 
                       length = unit(0.2, "inches")))+
  annotate("text", x=39/10, y=20, hjust=0,
           color='orange', size=7,
           label='False Discoveries', family = "Palatino Linotype")

ggsave('../results/slow-7.pdf', width = 12, height = 8, device =cairo_pdf)

# plot classical hurdle --------------------------------------------------------------

p8 = p7  + 
  geom_vline(xintercept = 1.96, size = linesize) 

ggsave('../results/slow-8.pdf', width = 12, height = 8, device =cairo_pdf)

p9 = p8 +
  annotate(geom="text", label="FDR = 9%", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  ) + 
  annotate("rect", xmin = 2, xmax = 3.5, 
           ymin = -10, ymax = 5, size=2,
           color='orange', fill=NA) +
  annotate("segment", x=4.2, xend=3.5,
           y=12, yend=5, size=2,
           color='orange',
           arrow=arrow(type = "closed", 
                       length = unit(0.2, "inches")))+
  annotate("text", x=39/10, y=20, hjust=0,
           color='orange', size=7,
           label='False Discoveries', family = "Palatino Linotype")

ggsave('../results/slow-9.pdf', width = 12, height = 8, device =cairo_pdf)

# POST TRUTH --------------------------------------------------------------------

# Start with 3 lines chart -------------------------------------------------------

q1 = p3 +
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash')+  
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) + 
  geom_vline(xintercept = 2.27, size = linesize, color=MATYELLOW
             , linetype='longdash') +
  annotate(geom="text", label="FDR = 5%", 
           x=2.40, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATYELLOW
  ) + 
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="FDR = 9%", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  ) 
  
ggsave('../results/post-truth-1.pdf', width = 12, height = 8, device =cairo_pdf)

# add Bonferroni
qtemp = q1 + geom_vline(xintercept = hurdle_bonf05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=hurdle_bonf05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) 
ggsave('../results/post-truth-2.pdf', width = 12, height = 8, device =cairo_pdf)  

# Relabel using post-truth -----------------------------------------------------------------

# label left of bonf as false
lab_hlz_false = '"False": Insignificant'
lab_hlz_true = '"True": Significant'

# post-truth labels
small = small %>% 
  mutate(
    v_hlz = case_when(
      tabs > hurdle_bonf05 ~ TRUE
      , tabs <= hurdle_bonf05 ~ FALSE
    )
  ) %>% mutate(
    v_hlz = factor(
      v_hlz, levels = c(TRUE,FALSE), labels = c(lab_hlz_true, lab_hlz_false)
    )
  ) %>% 
  mutate(
    v_hlz_alt = case_when(
      tabs > hurdle_05 ~ TRUE
      , tabs <= hurdle_05 ~ FALSE
    )
  ) %>% mutate(
    v_hlz_alt = factor(
      v_hlz_alt, levels = c(TRUE,FALSE), labels = c(lab_hlz_true, lab_hlz_false)
    )
  ) 

# relabel
q2 = ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v_hlz, color = v_hlz, shape = v_hlz), size = 2.5) +
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
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash')+  
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) + 
  geom_vline(xintercept = 2.27, size = linesize, color=MATYELLOW
             , linetype='longdash') +
  annotate(geom="text", label="FDR = 5%", 
           x=2.40, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATYELLOW
  ) + 
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="FDR = 9%", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  )  + 
  geom_vline(xintercept = hurdle_bonf05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=hurdle_bonf05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) 
 
ggsave('../results/post-truth-3.pdf', width = 12, height = 8, device =cairo_pdf)  

# even if exp ret is super high
tempcolor = 'red'
qtemp = q2 + 
  annotate("rect", xmin = 32/10, xmax = 38/10, 
           ymin = 100, ymax = 180, size=2,
           color=tempcolor, fill=NA) +
  annotate("segment", x=0/10, xend=33/10,
           y=180, yend=140, size=2,
           color=tempcolor,
           arrow=arrow(type = "closed", 
                       length = unit(0.2, "inches")))
ggsave('../results/post-truth-4.pdf', width = 12, height = 8, device =cairo_pdf)    

# Relabel for FDR 5% post-truth -------------------------------------------

# alternative relabelling
p5 = ggplot(small, aes(x=tselect,y=mu_scatter))  +
  geom_point(aes(group = v_hlz_alt, color = v_hlz_alt, shape = v_hlz_alt), size = 2.5) +
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
  geom_vline(xintercept = hurdle_bonf05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=hurdle_bonf05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) +
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash')+  
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) + 
  geom_vline(xintercept = 2.27, size = linesize, color=MATYELLOW
             , linetype='longdash') +
  annotate(geom="text", label="FDR = 5%", 
           x=2.40, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATYELLOW
  ) + 
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="FDR = 9%", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  )   + 
  annotate("rect", xmin = 32/10, xmax = 38/10, 
           ymin = 100, ymax = 180, size=2,
           color=tempcolor, fill=NA) +
  annotate("segment", x=0/10, xend=33/10,
           y=180, yend=140, size=2,
           color=tempcolor,
           arrow=arrow(type = "closed", 
                       length = unit(0.2, "inches")))

ggsave('../results/post-truth-5.pdf', width = 12, height = 8, device =cairo_pdf)  


# Some numbers for the paper ---------------------------------------------------
sig_share =  dat %>% filter(tabs>2) %>% summarize(sig_share = mean(tabs>2.27)) %>% 
  pull(sig_share)

sig_share*0.05  + (1-sig_share)

dat[tabs>2, ] %>% summarize(mean(v==label_false), n())

small[tabs>2.95]
small[tabs>2.27]
small[(tabs>2.27) & (v==label_false)]
hurdle_bonf05

small %>% 
  summarize(
    sum(tabs>2.95)
    , sum(tabs>2.27)
    , sum((tabs>2.27) & (v==label_false))
    , sum(tabs<3.8 & tabs>2)
  ) %>% 
  t()

-1*qnorm(0.05/296/2)


dat %>% 
  summarize(
    sum((tabs<3.8)&(tabs>2))/sum(tabs>2) 
  )