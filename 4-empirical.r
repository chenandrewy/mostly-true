# ENVIRONMENT ====

rm(list = ls())
library(tidyverse)
library(data.table)
library(googledrive)
library(readxl)
library(RColorBrewer)
library(lubridate)
library(boot)
library(e1071)
library(ggthemes)
library(gridExtra)
library(latex2exp)
library(Hmisc)

source('0-functions.r')

load('../data/emp_data.Rdata')

tdatpub = emp_sum

# ESTIMATE FDR ====

# settings
tbarlist = quantile(tdatpub$t, seq(0.1,0.9,0.1))

# estimate bias
est_mix = estimate_mixture(tdatpub$t, tgood = 2.6, pnull = 0.9, shape = 2, 1)
est_yz = list(
  tgood = 2.6, Pr_tgood = 0.30, C = 3.3
)
est_exp = estimate_exponential(tdatpub$t, 2.6)


# PLOT ====

# fdr calculations
tbarlist = seq(1,6,0.25) 

fdr = estimate_fdr(
  tdatpub$t, tbarlist = tbarlist, C = 1
) %>% 
  transmute(tbar, dr=dr, fdr_bh = fdrhat) %>% 
  mutate(
    fdryz =fdr_bh*est_yz$C
    , fdrmix =fdr_bh*est_mix$C
    , fdrexp =fdr_bh*est_exp$C
  ) %>% 
  mutate_at(
    .vars = vars(c(-tbar))
    , .funs = ~round(.*100,1)
  ) 

# prep plot
plotme = fdr %>% 
  select(tbar, starts_with('fdr'), -fdr_bh) %>% 
  pivot_longer(-tbar, names_to = 'type', values_to = 'fdr') %>% 
  mutate(
    type = factor(
      type
      , levels = c("fdryz","fdrmix","fdrexp")
      , labels = c('Yan-Zheng Bound','Conservative Mix','Exponential')
    )
  )


# plot
legtitle = 'Publication Bias Adjustment'
ggplot(
  plotme
  , aes(x=tbar, y=fdr, group = type)
) +
  geom_line(aes(linetype=type, color = type), size=2.5) +
  coord_cartesian(ylim = c(0,100)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted")) +
  scale_color_manual(values=c("#619CFF","#00BA38", "#F8766D")) +
  theme_economist_white(gray_bg = FALSE) +
  theme(
    axis.title = element_text(size = 20)
    , axis.text = element_text(size = 14)      
    , legend.title = element_text(size = 14)
    , legend.text = element_text(size = 14)
    , legend.background = element_rect(colour = 'black', fill = 'white', linetype = 'solid')
  )  +
  labs(
    x = TeX('\\bar{t}')
    , y = TeX('FDR Upper Bound for $|t_i|>\\bar{t}$ (\\%)')
    , linetype = legtitle, color = legtitle
  )  +
  theme(
    legend.position = c(70,80)/100
    , legend.key.width = unit(3,'cm')
  ) +
  scale_x_continuous(breaks = 1:6)

pdfscale = 0.7
ggsave('../results/emp-fdr.pdf', width = 10*pdfscale, height = 8*pdfscale)

# TABLE ====

texme = fdr %>% 
  mutate(
    blank = numeric(length(fdr$tbar))*NA
  ) %>% 
  transmute(
    tbar = tbar
    , dr = dr
    , blank
    , fdr_yz
    , fdr_mix
    , fdr_exp
  ) 



# Produces latex table code
temp = latex(
  texme
  , file = '../results/tab-emp.tex'
  , table.env = F
  , first.hline.double = FALSE
  , rowname=NULL
  , na.blank = T
  , already.math.col.names = T
  , cdec = c(2,1, 1,1,1)
)


