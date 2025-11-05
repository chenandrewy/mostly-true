# ABOUTME: Downloads and processes predictor return data from Google Drive and Dropbox
# ABOUTME: Creates cleaned datasets and CLZ bootstraps for Chen et al., CLZ, and Yan-Zheng predictors
# Inputs:
#   - config-and-functions.r (project configuration and helpers)
#   - googledrive package (for authentication)
#   - here package (for project-relative paths)
#   - Chen et al. predictors data (downloaded from Google Drive)
#   - CLZ predictor data (downloaded from Google Drive)
#   - Yan-Zheng predictor data (downloaded from Dropbox)
# Outputs:
#   - data/emp_data.Rdata (contains cz_ret, cz_sum, clz_ret, clz_sum, clzvw_ret, clzvw_sum, yz_ret, yz_sum, yzvw_ret, yzvw_sum)
#   - data/bootact.Rdata (bootstrap samples for CLZ returns)
#   - data/PredictorPortsFull.csv (downloaded from Google Drive)
#   - data/SignalDoc.csv (downloaded from Google Drive)
#   - data/CLZ_raw.csv (downloaded from Google Drive)
#   - data/CLZvw_raw.csv (downloaded from Google Drive)
#   - data/Yan_Zheng_RFS_Data.sas7bdat (downloaded from Dropbox)
# How to run:
#   Rscript 01-prep-data.r
#   (will prompt for Google Drive authentication in browser)

# ENVIRONMENT -------------------------------

rm(list = ls())

library(here)
library(googledrive)
library(haven)
here::i_am("01-prep-data.r")

source(here("config-and-functions.r"))

paths <- project_paths()
data_dir <- paths$data

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

drive_download(target_dribble, path = file.path(data_dir, 'PredictorPortsFull.csv'), overwrite = T)

# download header info
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='SignalDoc.csv')

drive_download(target_dribble, path = file.path(data_dir, 'SignalDoc.csv'), overwrite = T)

# PROCESS DATA -------------------------------

## import into R -------------------------------

ret0 = fread(file.path(data_dir, 'PredictorPortsFull.csv')) %>% 
  filter(port=='LS') %>% 
  select(signalname, date, ret)

header = fread(file.path(data_dir, 'SignalDoc.csv')) %>% 
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

drive_download(target_dribble, path = file.path(data_dir, 'CLZ_raw.csv'), overwrite = T)

# clean up
temp0 = fread(file.path(data_dir, 'CLZ_raw.csv')) 

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

drive_download(target_dribble, path = file.path(data_dir, 'CLZvw_raw.csv'), overwrite = T)

# clean up
temp0 = fread(file.path(data_dir, 'CLZvw_raw.csv')) 

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

# ADD YZ DATA -------------------------------

# download from Dropbox
url_yz = 'https://www.dropbox.com/scl/fi/tv7zhkgrb7jbl4mr2ccdn/Yan_Zheng_RFS_Data.sas7bdat?rlkey=zxcrit8ptfe4mdq0q5blbx71c&st=s5ox6lts&dl=1'
download.file(url_yz, destfile = file.path(data_dir, 'Yan_Zheng_RFS_Data.sas7bdat'), mode = 'wb')

# import sas data
yzraw = read_sas(file.path(data_dir, 'Yan_Zheng_RFS_Data.sas7bdat'))
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

# SAVE TO DISK -------------------------------

save(
  list = c(
    'cz_ret','cz_sum'
    , 'clz_ret', 'clz_sum', 'clzvw_ret', 'clzvw_sum'
    , 'yz_ret', 'yz_sum', 'yzvw_ret', 'yzvw_sum'
    )
  , file = file.path(data_dir, 'emp_data.Rdata')
)

# BOOTSTRAP CLZ RETURNS -------------------------------

set_boot = list(
  nboot = 1000,
  min_obs_pct = 50,
  ncore = 8
)

boot_timer_start = Sys.time()

bootact = bootstrap_flex(
  clz_ret,
  nboot = set_boot$nboot,
  min_obs_pct = set_boot$min_obs_pct,
  demean = FALSE,
  ncore = set_boot$ncore
)

boot_timer_end = Sys.time()

print(paste0('CLZ bootstrap finished in minutes: ', difftime(boot_timer_end, boot_timer_start, units = 'mins')))

save(bootact, file = file.path(data_dir, 'bootact.Rdata'))
