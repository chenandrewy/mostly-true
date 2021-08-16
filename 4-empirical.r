# ==== ENVIRONMENT ====

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

source('0-functions.r')

ret = fread('../data/clean_ret.csv')

# univariate
tdatpub = ret %>% 
  filter(insamp) %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , nmonth = n()
    , t = rbar/vol*sqrt(nmonth)
  ) %>% 
  filter(t>1.96)

## SETTINGS
tbarlist = seq(0,6,0.25) 

tbarlist = quantile(tdatpub$t, seq(0.1,0.9,0.1))

# fdrhat using exp
est_bias = estimate_exponential(tdatpub$t, 2.6)
fdr_exp = estimate_fdr(
  tdatpub$t, tbarlist = tbarlist, C = est_bias$C
) %>% 
  transmute(tbar, dr, fdrhat_exp = fdrhat)

# fdrhat using mix
est_bias = estimate_mixture(tdatpub$t, tgood = 2.6, pnull = 0.9, shape = 2, 1)

est_bias

fdr_mix = estimate_fdr(
  tdatpub$t, tbarlist = tbarlist, C = est_bias$C
) %>% 
  transmute(tbar, fdrhat_mix = fdrhat)

est_all = fdr_exp %>% 
  left_join(fdr_mix, by = 'tbar')


est_all

ggplot(
  est_all %>% 
    select(tbar, starts_with('fdr')) %>% 
    pivot_longer(-tbar, names_to = 'type', values_to = 'fdr')
  , aes(x=tbar, y=fdr, group = type)
) +
  geom_line(aes(linetype=type, color = type)) +
  coord_cartesian(ylim = c(0,1))


