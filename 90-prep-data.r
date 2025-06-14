# 2021 08 Andrew
# Generate data for BH style FDR estimates 
# makes '../data/emp_data.Rdata'

# ENVIRONMENT -------------------------------

rm(list = ls())

# Set working directory to unbreakable-bh folder
if (basename(getwd()) != "unbreakable-bh") {
  # Try to find unbreakable-bh directory
  if (dir.exists("unbreakable-bh")) {
    setwd("unbreakable-bh")
  } else if (dir.exists("../unbreakable-bh")) {
    setwd("../unbreakable-bh")  
  } else {
    stop("Please run this script from the unbreakable-bh directory or its parent directory.")
  }
}

source('functions.r')
library(googledrive)

### USER ENTRY
# root of March 2022 release on Gdrive
pathRelease = 'https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo'

# login to gdrive
# this prompts a login
pathRelease %>% drive_ls()

# DOWNLOAD DATA -------------------------------

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

# PROCESS DATA -------------------------------

## import into R -------------------------------

ret0 = fread('../data/PredictorPortsFull.csv') %>% 
  filter(port=='LS') %>% 
  select(signalname, date, ret)

header = fread('../data/SignalDoc.csv') %>% 
  rename(signalname=Acronym)

## construct benchmark data -------------------------------
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
  ) %>% setDT()

# ADD CLZ DATA -------------------------------
# created 2023 12 
# this replaces the yz data

# copy-paste from browser via Chen's website
url_clz = 'https://drive.google.com/drive/folders/16RqeHNyU5gcqjRUvqSeOfQCxu_B2mfcZ'

target_dribble = url_clz %>% drive_ls() %>% 
  filter(name=='DataMinedLongShortReturnsEW.csv')

drive_download(target_dribble, path = '../data/CLZ_raw.csv', overwrite = T)

# clean up
temp0 = fread('../data/CLZ_raw.csv') 

clz_ret = temp0 %>% 
  mutate(date = paste(year, month, '28', sep = '-')
      , date = as.Date(date, format = '%Y-%m-%d')) %>% 
  transmute(signalname = signalid, date, ret) %>% 
  # filter for July 1963 or later
  filter(date >= '1963-07-01')

# repeatedly used summary stats (in-sample)
clz_sum = clz_ret %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , nmonth = n()
    , traw = rbar/vol*sqrt(nmonth)
    , tabs = abs(traw)
  ) %>% setDT()

# ADD CLZ VW DATA -------------------------------
# created 2024 03
# trying to kick it up a notch

# copy-paste from browser via Chen's website
url_clz = 'https://drive.google.com/drive/folders/16RqeHNyU5gcqjRUvqSeOfQCxu_B2mfcZ'

target_dribble = url_clz %>% drive_ls() %>% 
  filter(name=='DataMinedLongShortReturnsVW.csv')

drive_download(target_dribble, path = '../data/CLZvw_raw.csv', overwrite = T)

# clean up
temp0 = fread('../data/CLZvw_raw.csv') 

clzvw_ret = temp0 %>% 
  mutate(date = paste(year, month, '28', sep = '-')
      , date = as.Date(date, format = '%Y-%m-%d')) %>% 
  transmute(signalname = signalid, date, ret) %>% 
  # filter for July 1963 or later
  filter(date >= '1963-07-01')

# repeatedly used summary stats (in-sample)
clzvw_sum = clzvw_ret %>% 
  group_by(signalname) %>% 
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , nmonth = n()
    , traw = rbar/vol*sqrt(nmonth)
    , tabs = abs(traw)
  ) %>% setDT()

# SAVE TO DISK -------------------------------

save(
  list = c(
    'cz_ret','cz_sum', 'clz_ret', 'clz_sum', 'clzvw_ret', 'clzvw_sum'
    )
  , file = '../data/emp_data.Rdata'
)


# CHECK IF YZ DATA IS AVAILABLE ===

if (file.exists('../data_yan_zheng/unzipped/Yan_Zheng_RFS_Data.sas7bdat')){
  # import sas data
  library(haven)
  yzraw = read_sas('../data_yan_zheng/unzipped/Yan_Zheng_RFS_Data.sas7bdat')
  setDT(yzraw)

  # clean
  yzraw[ , signalname := paste0(transformation, '|', fsvariable)][
    , date := DATE][
    , c('DATE', 'transformation','fsvariable') := NULL]
  
  # separate
  yz_ret = yzraw %>% transmute(date,signalname,ret=ddiff_ew)
  yzvw_ret = yzraw %>% transmute(date,signalname,ret=ddiff_vw)

  # summarize
  yz_sum = yz_ret %>% 
    group_by(signalname) %>% 
    summarize(
      rbar = mean(ret)
      , vol = sd(ret)
      , nmonth = n()
      , traw = rbar/vol*sqrt(nmonth)
      , tabs = abs(traw)
    ) %>% setDT()
  yzvw_sum = yzvw_ret %>% 
    group_by(signalname) %>% 
    summarize(
      rbar = mean(ret)
      , vol = sd(ret)
      , nmonth = n()
      , traw = rbar/vol*sqrt(nmonth)
      , tabs = abs(traw)
    ) %>% setDT()
  
  ## SAVE TO DISK -------------------------------
  save(
    list = c(
      'cz_ret','cz_sum'
      , 'clz_ret', 'clz_sum', 'clzvw_ret', 'clzvw_sum'
      , 'yz_ret', 'yz_sum', 'yzvw_ret', 'yzvw_sum'
      )
    , file = '../data/emp_data.Rdata'
  )  

} # end if yz data is availabl
