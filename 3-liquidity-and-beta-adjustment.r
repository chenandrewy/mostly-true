# 2024 03: to correct Harvey and Liu's (2020) errors in calculating the 
# fraction of true strategies in Yan-Zheng 2017

# Setup -----------------------------------------------------------------------

rm(list = ls())
source('0-functions.r')
load('../data/emp_data.Rdata')


# merge vw and ew data
yz_ret1 = clz_ret %>% mutate(sweight = 'ew') %>%
rbind(
    clzvw_ret %>% mutate(sweight = 'vw') 
) %>% 
    select(sweight, signalname, date, ret) 

## Download FF + Carhart factors ====

# Download FF3
download.file("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_CSV.zip"
    , '../data/deleteme.zip')
unzip('../data/deleteme.zip', exdir = '../data/')

# Download Momentum
download.file("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Momentum_Factor_CSV.zip"
    , '../data/deleteme.zip')
unzip('../data/deleteme.zip', exdir = '../data/')

# read csvs and merge
ff_fac = fread('../data/F-F_Research_Data_Factors.csv') %>% 
    left_join(fread('../data/F-F_Momentum_Factor.csv'), by = 'V1') 
colnames(ff_fac) = c('yearm','mktrf','smb','hml','rf','mom')

# clean
ff_fac = ff_fac %>% 
    mutate(year = yearm %/% 100, month = yearm %% 100) %>% 
    filter(!is.na(mktrf+smb+hml+mom)) %>% 
    select(year, month, mktrf, smb, hml, mom, rf)

# Find alphas (takes a couple minutes) ----------------------------------------

yz_ret2 = copy(yz_ret1)

# merge
yz_ret2[ , c('year','month') := list(year(date), month(date))]
yz_ret2[ff_fac, on = c('year','month')
    , `:=` (mktrf = i.mktrf, smb = i.smb, hml = i.hml, mom = i.mom, rf = i.rf)]
yz_ret2[ , retrf := ret - rf]

# define list of factor models
#   ret is already a long-short so I shouldn't subtract rf (?)
modelstrlist = c(
    'ret ~ 1'
    ,'ret ~ mktrf'
    , 'ret ~ mktrf + smb + hml'
    , 'ret ~ mktrf + smb + hml + mom'
    )
modelnamelist = c('Raw', 'CAPM', 'FF3', 'FF3+Mom')

# find t-stats for each model
regest = list()
print('Finding t-stats for model: ')
for (modelstr in modelstrlist) {
    print(modelstr)
    modelform = as.formula(modelstr)
    regest[[modelstr]] = yz_ret2[, list(
        tstat = summary(lm(modelform, data = .SD))[['coefficients']]['(Intercept)', 't value']
    )
    , by = c('sweight','signalname')] %>% 
        mutate(model = modelnamelist[match(modelstr, modelstrlist)])
}
regest = do.call(rbind, regest)

# Find FDRmax -----------------------------------------------------------------

# find pFmax 

# pmin is what HL call theta
# convert to tmax using
# pmin = Pr(|Z|<tmax) = (Pr(Z \in [-tmax, tmax]))
#   = 2*(Pr(Z < tmax)-0.5)

tmaxlist = c(
    0.5
    , 1.0
    , qnorm(0.4/2+0.5) # pmin = 0.4
    , qnorm(0.6/2+0.5) # pmin = 0.4
    , qnorm(0.8/2+0.5) # pmin = 0.4
) %>% print()

storeydat = list()
for (tmax in tmaxlist){
    storeydat[[as.character(tmax)]] = regest %>% 
        as_tibble() %>% 
        group_by(sweight, model) %>% 
        summarize(Pr_lt_tmax = mean(abs(tstat) < tmax)) %>% 
        transmute(sweight, model, pmin = 2*(pnorm(tmax)-0.5), tmax, Pr_lt_tmax) %>% 
        mutate(pFmax = Pr_lt_tmax/pmin)
}
storeydat = do.call(rbind, storeydat) %>% 
    arrange(sweight, model, tmax)
setDT(storeydat)

# find FDRs
tmin = 2
# F_null = function(tmin) 2*(1-pnorm(tmin)) # better to use the ez 5% for consistency
F_null_ez = 0.05
tabdat = regest[ , .(Pr_gt_tmin = mean(abs(tstat) > tmin)), by = c('sweight','model')] %>% 
    left_join(storeydat, by = c('sweight','model')) %>% 
    mutate(FDRmax = F_null_ez/Pr_gt_tmin*pFmax
        , FDRmaxez = F_null_ez/Pr_gt_tmin) %>% 
    as_tibble() %>% 
    arrange(tmax, sweight, model) 

# export to latex -------------------------------------------------------------

tmaxselect = 0.5

# make pretty
tab = tabdat %>% 
    filter(tmax == tmaxselect) %>% 
    filter(model != 'Raw') %>% 
    select(sweight, model
        , Pr_gt_tmin, FDRmaxez
        , Pr_lt_tmax, pFmax,  FDRmax) %>% 
    pivot_longer(cols = -c(sweight, model)) %>% 
    mutate(value = round(100*value, 1)) %>%
    pivot_wider(names_from=c('sweight','model'), values_from=value) %>% 
    mutate(name = case_when(
        name == 'Pr_lt_tmax' ~ paste0('$\\Pr(|t|\\le', round(tmaxselect,1), ')$ (\\%)')
        , name == 'pFmax' ~ '$\\Pr(F)$ max (\\%)'
        , name == 'Pr_gt_tmin' ~ paste0('$\\Pr(|t|>', tmin, ')$ (\\%)')
        , name == 'FDRmax' ~ '$\\FDRez$ Storey Bound (\\%)'
        , name == 'FDRmaxez' ~ '$\\FDRez$ Easy Bound (\\%)'
    ))

# export to latex
library(kableExtra)
tab %>% 
    kable('latex', booktabs = T, linesep = '', escape = F, digits = 1
    ) %>% 
    cat(file='../results/temp.tex')

# read in temp.tex and modify manually
tex = readLines('../results/temp.tex') 
tex[4] = ' & 1-Factor & 3-Factor & 4-Factor & 1-Factor & 3-Factor & 4-Factor \\\\'
tex = tex %>% append('\\\\ \\hline', after=7)
tex = tex  %>% append(' & \\multicolumn{3}{c}{Equal-Weighted}  & \\multicolumn{3}{c}{Value-Weighted} \\\\', after = 3) 

writeLines(tex, con='../results/vw-ffn-storey.tex')

# export to latex with raw ----------------------------------------------------

tmaxselect = 0.5

# make pretty
tab = tabdat %>% 
    filter(tmax == tmaxselect) %>% 
    mutate(model = if_else(model=='Raw','0_Raw',model)) %>%
    arrange(sweight,model) %>% 
    select(sweight, model
        , Pr_gt_tmin, FDRmaxez
        , Pr_lt_tmax, pFmax,  FDRmax) %>% 
    pivot_longer(cols = -c(sweight, model)) %>% 
    mutate(value = round(100*value, 1)) %>%
    pivot_wider(names_from=c('sweight','model'), values_from=value) %>% 
    mutate(name = case_when(
        name == 'Pr_lt_tmax' ~ paste0('$\\Pr(|t|\\le', round(tmaxselect,1), ')$ ')
        , name == 'pFmax' ~ '$\\Pr(F)$ max '
        , name == 'Pr_gt_tmin' ~ paste0('$\\Pr(|t|>', tmin, ')$ ')
        , name == 'FDRmax' ~ '$\\FDRez$ max Storey'
        , name == 'FDRmaxez' ~ '$\\FDRez$ max Easy'
    ))

# add blank column
tab2 = cbind(tab[ , 1:5]
    , matrix('', nrow(tab), 1)
    , tab[ , 6:9])

# export to latex
library(kableExtra)
tab2 %>% 
    kable('latex', booktabs = T, linesep = '', escape = F, digits = 1
    ) %>% 
    cat(file='../results/temp.tex')

# read in temp.tex and modify manually
tex = readLines('../results/temp.tex') 
tex[4] = ' & Raw & CAPM & FF3 & 4-Fac & & Raw & CAPM & FF3 & 4-Fac \\\\'
tex[5] = '\\cline{2-5} \\cline{7-10}'
tex = tex %>% append('\\\\ \\cline{2-5} \\cline{7-10}', after=7)
tex = tex  %>% append(' & \\multicolumn{4}{c}{Equal-Weighted}  & &  \\multicolumn{4}{c}{Value-Weighted} \\\\', after = 3) 

writeLines(tex, con='../results/vw-ffn-storey-raw.tex')







