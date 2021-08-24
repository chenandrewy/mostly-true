# 2021 08 estimate bias adjustment 

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
source('0-functions.r')

load('../data/emp_data.Rdata')


qlist = seq(0.1,0.9,0.1)

# SUMMARY STATS ====

## univariate
tab_univar_quantiles = emp_sum %>% select(-signalname) %>% 
  apply(2, quantile, probs=qlist) %>% 
  t()


## correlations

# first do pairwise complete
retwideall = pivot_wider(
  emp_ret
  , c(signalname,date,ret), names_from = signalname, values_from = ret
) %>% 
  select(-date)
corpairwise = cor(retwideall, use = 'pairwise.complete.obs')

# then do balanced matrix with minimum missing
retwidebal = pivot_wider(
  emp_ret
  , c(signalname,date,ret), names_from = signalname, values_from = ret
) %>% 
  select(-date)
corbalanced = cor(retwidebal, use = 'complete.obs')

tab_corr_quantiles = rbind(
  quantile(
    corpairwise[lower.tri(corpairwise)], probs = qlist
  ) 
  , quantile(
    corbalanced[lower.tri(corbalanced)], probs = qlist
  ) 
) 
row.names(tab_corr_quantiles) = c('pairwise','subset complete')


## table prep


# 
finalMatrix <- rbind(tab_univar_quantiles, tab_corr_quantiles)

# Need to make these numbers shorter. Otherwise, they won't fit in one page.
finalMatrix <- round(finalMatrix, 1)
finalMatrix[3,] = as.integer(round(finalMatrix[3,],0))
finalMatrix

# Take transpose
finalMatrix <- t(finalMatrix)
# Rearrange columns
order <- c("t", "rbar", "vol", "nmonth", "pairwise", "subset complete")
finalMatrix <- finalMatrix[,order]
#Take transpose again
finalMatrix <- t(finalMatrix)

# % character in column names will cause issues with latex. Remove them.
colnames(finalMatrix) <- gsub("\\%", "", colnames(finalMatrix))

# Specify rownames
rownames(finalMatrix) <- c(Hmisc::latexTranslate("| t |"),
                           "Mean Return",
                           "Volatility",
                           "Num of Obs",
                           "Pairwise",
                           "Subset Complete")


## table 

# Produces latex table code
capture.output(
  Hmisc::latex(finalMatrix,
               file = "",
               title = '',
               table.env=F,
               cgroup = "Percentile",
               n.cgroup = 9,
               rgroup = c("Section 1", "Section 2"),
               n.rgroup = c(length(row.names(tab_univar_quantiles)),
                            length(row.names(tab_corr_quantiles))),
               cgroupTexCmd = "normalfont",
               rgroupTexCmd = "normalfont",
               col.just = c('c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'), 
               colhead = colnames(finalMatrix),
               math.row.names = F,
               first.hline.double = FALSE,
               insert.bottom = "",
               rdec = c(2, 2, 2, 0, 2, 2)),
  file = "../results/mattab-sum.tex")


# NON-PARAMETRIC FDR ====

# settings
tgood = 2.6
nulldf = 100
C = 6.7

tbarlist = quantile(emp_sum$t, seq(0.1,0.6,0.1))
tbarlist = tbarlist[tbarlist > tgood]

fdr = estimate_fdr(emp_sum$t, tbarlist, C=C, nulldf=nulldf)

tabme = fdr %>% 
  mutate(
    pval = 2*pt(-tbar,nulldf)
    , fdrbh = fdrhat/C
  ) %>% 
  select(
    tbar, pval, dr, fdrbh, fdrhat
  ) %>% 
  mutate(
    tbar = round(tbar,2)
  ) %>% 
  mutate_at(
    .vars = vars(c(-tbar))
    , .funs = ~round(.*100,1)
  ) %>% 
  remove_rownames() %>% 
  t()


# Produces latex table code
temp = Hmisc::latex(
  tabme
  , file = '../results/mattab-nonpar.tex'
  , table.env = F
  , first.hline.double = FALSE
  , title = NULL
  , rowname= c(
    ''
    ,'$p$-value for $\\bar{t}$'
    ,'Share of Selected $|t_i|>\\bar{t}$ (\\%)'
    ,'$\\widehat{\\text{FDR}}_{N,\\text{BHS}}(\\bar{t})$ (\\%)'
    ,'$\\widehat{\\text{FDR}}_{N}(\\bar{t})$ (\\%)'
    )
  , already.math.row.names = T
  , na.blank = T
  , rdec = c(2,1, 1,1,1)
)
