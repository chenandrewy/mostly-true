# 2021 08 Andrew
# Generate data for BH style FDR estimates 
# makes data/clean_ret.csv
# also makes summary stats

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

### USER ENTRY
# root of April 2021 release on Gdrive
pathRelease = 'https://drive.google.com/drive/folders/1I6nMmo8k_zGCcp9tUvmMedKTAkb9734R'
url_prefix = 'https://drive.google.com/uc?export=download&id='


# login to gdrive
# this prompts a login
pathRelease %>% drive_ls()


# ==== DOWNLOAD DATA =====

# download monthly returns
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='Portfolios') %>% drive_ls() %>% 
  filter(name=='Full Sets OP') %>% drive_ls() %>% 
  filter(name=='PredictorPortsFull.csv')

drive_download(target_dribble, path = '../data/PredictorPortsFull.csv', overwrite = T)

# download header info
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='SignalDocumentation.xlsx')

drive_download(target_dribble, path = '../data/SignalDocumentation.xlsx', overwrite = T)

# ==== PROCESS DATA ====
ret = fread('../data/PredictorPortsFull.csv') %>% 
  filter(port=='LS') %>% 
  select(signalname, date, ret)

header = read_excel('../data/SignalDocumentation.xlsx',sheet='BasicInfo') %>% 
  left_join(
    read_excel('../data/SignalDocumentation.xlsx',sheet='AddInfo')
    , by = c('Acronym','Authors')
  ) %>% 
  rename(signalname=Acronym)

ret = ret %>% 
left_join(
  header %>% select(signalname,SampleStartYear,SampleEndYear)
  , by = c('signalname')
) %>% 
  mutate(
    insamp = year(date) >= SampleStartYear & year(date) <= SampleEndYear
  ) %>% 
  select(signalname,date,ret,insamp)


# ==== SAVE TO DISK ====
write_csv(ret, '../data/clean_ret.csv')


# ==== SUMMARY STATS ====

# --- user entry ---
min_nmonth = 200
min_nsignal = 150

# re-read for modularity
ret = fread('../data/clean_ret.csv')

# univariate
signalsum = ret %>% 
  filter(insamp) %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , nmonth = n()
    , t = rbar/vol*sqrt(nmonth)
  )

## correlations

# first do pairwise complete
retwideall = pivot_wider(
  ret
  , c(signalname,date,ret), names_from = signalname, values_from = ret
) %>% 
  select(-date)
corpairwise = cor(retwideall, use = 'pairwise.complete.obs')

# then do balanced matrix with minimum missing
monthsum = ret %>% 
  group_by(date) %>% 
  summarize(nsignal = sum(!is.na(ret))) 

retwide = ret %>%
  left_join(
    signalsum,  by = 'signalname'
  ) %>% 
  left_join(
    monthsum, by = 'date'
  ) %>% 
  filter(
    nmonth >= min_nmonth, nsignal >= min_nsignal, !is.na(ret)
  ) %>% 
  select(signalname, date, ret) %>% 
  pivot_wider(
    c(signalname,date,ret), names_from = signalname, values_from = ret
  ) %>% 
  select(-date) %>% 
  filter(complete.cases(.))

corbalanced = cor(retwide, use = 'complete.obs')

# save balanced matrix 
write_csv(retwide, '../data/balanced_ret.csv')


## tables ====
qlist = seq(0.1,0.9,0.1)

tab_univar_quantiles = signalsum %>% select(-signalname) %>% 
  apply(2, quantile, probs=qlist) %>% 
  t()


tab_corr_quantiles = rbind(
  quantile(
    corpairwise[lower.tri(corpairwise)], probs = qlist
  ) 
  , quantile(
    corbalanced[lower.tri(corbalanced)], probs = qlist
  ) 
  ) 
row.names(tab_corr_quantiles) = c('pairwise','subset complete')


tab_univar_quantiles

tab_corr_quantiles




# ====
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
finalMatrix

# LATEX TABLE ------------------------------------------------------------

# Produces latex table code
capture.output(
  Hmisc::latex(finalMatrix,
               file = "",
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
               insert.bottom = ""),
  file = "table.tex")
