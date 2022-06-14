# 2021 08 Andrew
# Generate data for BH style FDR estimates 
# makes '../data/emp_data.Rdata'

# ENVIRONMENT ====

rm(list = ls())
source('0-functions.r')
detach_all()
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
# root of March 2022 release on Gdrive
pathRelease = 'https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo'


# login to gdrive
# this prompts a login
pathRelease %>% drive_ls()

# parameters for balanced panel
min_nmonth = 200
min_nsignal = 150


# DOWNLOAD DATA =====

# download monthly returns
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='Portfolios') %>% drive_ls() %>% 
  filter(name=='Full Sets OP') %>% drive_ls() %>% 
  filter(name=='PredictorPortsFull.csv')

drive_download(target_dribble, path = '../data/PredictorPortsFull.csv', overwrite = T)


# download header info
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='SignalDoc.csv')

drive_download(target_dribble, path = '../data/SignalDoc.csv', overwrite = T)

# PROCESS DATA ====

## import into R ====

ret0 = fread('../data/PredictorPortsFull.csv') %>% 
  filter(port=='LS') %>% 
  select(signalname, date, ret)

header = fread('../data/SignalDoc.csv') %>% 
  rename(signalname=Acronym)


## construct benchmark data ====
# add sample info
cz_ret = ret0 %>% 
  left_join(
    header %>% select(signalname,SampleStartYear,SampleEndYear)
    , by = c('signalname')
  ) %>% 
  mutate(
    insamp = year(date) >= SampleStartYear & year(date) <= SampleEndYear
  ) %>% 
  select(signalname,date,ret,insamp)

# repeatedly used summary stats (in-sample)
cz_sum = cz_ret %>% 
  filter(insamp) %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , nmonth = n()
    , traw = rbar/vol*sqrt(nmonth)
    , tabs = abs(traw)
  )



# ADD YZ DATA ====
# created 2022 04 to read yz data
# yz data send via email from Sterling Yan in 2019
# this is simpler b/c no Gdrive stuff, no in-sample stuff, no balancing

library(haven) # for read_sas

temp = read_sas('../data_yan_zheng/unzipped/Yan_Zheng_RFS_Data.sas7bdat')

# clean, select ew or vw
yz_ret = temp %>% 
  mutate(
    signalname = paste(transformation, fsvariable, sep = '.')
  ) %>% 
  transmute(
    signalname, date = DATE, ret = 100*ddiff_ew
  )

# find summary stats
yz_sum = yz_ret %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret), vol = sd(ret), nmonth = n(), traw = rbar/vol*sqrt(nmonth), tabs = abs(traw)
  )

# SAVE TO DISK ====


save(
  list = c(
    'cz_ret','cz_sum', 'yz_ret', 'yz_sum'
    )
  , file = '../data/emp_data.Rdata'
)

