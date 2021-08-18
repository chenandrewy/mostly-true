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
library(truncdist)
library(gridExtra)
library(Hmisc)

source('0-functions.r')

load('../data/emp_data.Rdata')

tdatpub = emp_sum

# ESTIMATE FDR ====

# settings
tbarlist = seq(0,6,0.25) 
tbarlist = quantile(tdatpub$t, seq(0.1,0.9,0.1))

# estimate bias
est_mix = estimate_mixture(tdatpub$t, tgood = 2.6, pnull = 0.9, shape = 2, 1)
est_yz = list(
  tgood = 2.6, Pr_tgood = 0.30, C = 3.3
)
est_exp = estimate_exponential(tdatpub$t, 2.6)

# fdr calculations
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


# PLOT ====

ggplot(
  fdr %>% 
    select(tbar, starts_with('fdr')) %>% 
    pivot_longer(-tbar, names_to = 'type', values_to = 'fdr')
  , aes(x=tbar, y=fdr, group = type)
) +
  geom_line(aes(linetype=type, color = type)) +
  coord_cartesian(ylim = c(0,100))


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


