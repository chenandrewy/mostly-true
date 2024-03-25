# 2024 03: to correct Harvey and Liu's (2020) errors in calculating the 
# fraction of true strategies in Yan-Zheng 2017

# Setup -----------------------------------------------------------------------

rm(list = ls())
source('0-functions.r')

## Import YZ data ====
# created 2022 04 to read yz data
# yz data send via email from Sterling Yan in 2019
# this is simpler b/c no Gdrive stuff, no in-sample stuff, no balancing

library(haven) # for read_sas
yz_ret0 = read_sas('../data_yan_zheng/unzipped/Yan_Zheng_RFS_Data.sas7bdat')
setDT(yz_ret0)

# clean
yz_ret1 = copy(yz_ret0)
yz_ret1[ , ':=' (date = DATE, signalname=paste(transformation, fsvariable, sep = '.'))] 
yz_ret1 = yz_ret1 %>% transmute(signalname, date, ret = 100*ddiff_ew, sweight = 'ew') %>% 
    rbind(
        yz_ret1 %>% transmute(signalname, date, ret = 100*ddiff_vw, sweight = 'vw')
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

# Find alphas -----------------------------------------------------------------
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

# Find P(F) -----------------------------------------------------------

# pmin is what HL call theta
# convert to tmax using
# pmin = Pr(|Z|<tmax) = (Pr(Z \in [-tmax, tmax]))
#   = 2*(Pr(Z < tmax)-0.5)

tmaxlist = c(
    0.5
    , qnorm(0.4/2+0.5) # pmin = 0.4
    , qnorm(0.6/2+0.5) # pmin = 0.4
    , qnorm(0.8/2+0.5) # pmin = 0.4
) %>% print()


# estimate pFmax
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

storeydat 

# find t-stat hurdles
for (modelcur in modelnamelist){
    if (modelcur==modelnamelist[1]){hurdledat = data.table()}
    for (sweightcur in c('ew','vw')){
        # define emprical CDF
        Fcur = ecdf(abs(regest[sweight == sweightcur & model == model]$tstat))
        tabs_emp_max = max(abs(regest[sweight == sweightcur & model == model]$tstat))
        
        # loop over storey tuning parameters
        for (storeytmax in tmaxlist){
            # get pFmax
            pFmaxcur = storeydat[sweight == sweightcur & model == modelcur & tmax == storeytmax]$pFmax
            
            # solve FDRmaxfun(tmin) = alpha
            alphalist = c(0.01, 0.05, 0.1)
            for (alpha in alphalist){
                solveme = function(tmin) 2*(1-pnorm(tmin))/(1-Fcur(tmin))*pFmaxcur - alpha
                thurdle = uniroot(solveme, c(0.1,tabs_emp_max-0.1))$root
                pct_signif = 100*(1-Fcur(thurdle))
                pct_signif_null = 100*2*(1-pnorm(thurdle))

                # finally, store
                tempdat = data.table(
                    sweight = sweightcur, model = modelcur, tmax = storeytmax
                    , FDRmax = alpha, thurdle = thurdle, pct_signif = pct_signif
                    , pct_signif_null = pct_signif_null
                )
                hurdledat = rbind(hurdledat, tempdat)                
            }            
        }
    }
}


# increase print width
options(width = 120)

# make simpler table
tab1 = storeydat %>% filter(pmin==0.8, model != 'Raw') %>% 
    left_join(
        hurdledat[FDRmax==0.05] %>% transmute(sweight, model, tmax
            , thurdle05 = thurdle
            , pct_signif_null05 = pct_signif_null
            , pct_signif05 = pct_signif)
        , by = c('sweight','model','tmax')
    ) %>% 
    left_join(
        hurdledat[FDRmax==0.1] %>% transmute(sweight, model, tmax
            , thurdle10 = thurdle
            , pct_signif_null10 = pct_signif_null
            , pct_signif10 = pct_signif)
        , by = c('sweight','model','tmax')
    ) %>%     
    pivot_longer(cols=c('tmax','pmin','Pr_lt_tmax','pFmax'
        , ends_with('05'), ends_with('10'))) %>% 
    # round numbers
    mutate(
        value = case_when(
            name %in% c('pmin','Pr_lt_tmax','pFmax') ~ round(100*value, 1)
            , name %in% c('tmax','thurdle05','thurdle10') ~ round(value, 2)
            , grepl('pct_', name)  ~ round(value, 1)
            , TRUE ~ round(value,1)
        )
    ) %>% 
    pivot_wider(names_from=c('sweight','model'), values_from=value)

# format
tmaxnum = storeydat %>% filter(pmin==0.8) %>% pull(tmax) %>% unique()


# using a frame that lists the rows we want
tab2 = tibble(
    name = c('Pr_lt_tmax', 'pFmax'
        , 'blank1'
        , 'FDRmax05', 'thurdle05', 'pct_signif_null05', 'pct_signif05'
        , 'blank2'
        , 'FDRmax10', 'thurdle10', 'pct_signif_null10', 'pct_signif10')
) %>% 
    left_join(tab1, by = 'name') %>% 
    # clean up row names
    mutate(name = case_when(
        name == 'Pr_lt_tmax' ~ paste0('$\\Pr(|t|\\le', round(tmaxnum,1), ')$ (\\%)')
        , name == 'pFmax' ~ '$\\Pr(F)$ max (\\%)'
        , grepl('thurdle', name) ~ '$h^\\ast_{\\text{FDRmax}}$'
        , grepl('pct_signif_null', name) ~ '$\\Pr(|t| >h^\\ast_{\\text{FDRmax}} | F)$ (\\%)'
        , grepl('pct_signif', name) ~ '$\\Pr(|t| >h^\\ast_{\\text{FDRmax}})$ (\\%)'
        , TRUE ~ name
    )) %>% 
    print()

# export to latex
library(kableExtra)
tab2 %>% 
    kable('latex', booktabs = T, linesep = '', escape = F, digits = 1
    ) %>% 
    cat(file='../results/temp.tex')

# read in temp.tex and modify manually
tex = readLines('../results/temp.tex') 
tex[4] = ' & 1-Factor & 3-Factor & 4-Factor & 1-Factor & 3-Factor & 4-Factor \\\\'
tex[8] = '\\\\'
tex[9] = 'FDRmax = 5\\% \\\\ \\hline'
tex[13] = '\\\\'
tex[14] = 'FDRmax = 10\\% \\\\ \\hline'
tex = tex %>% 
    append(' & \\multicolumn{3}{c}{Equal-Weighted}  & \\multicolumn{3}{c}{Value-Weighted} \\\\', after = 3) 

tex
writeLines(tex, con='../results/yz-storey.tex')

pnorm(1.3) - pnorm(-1.3)


# pFmax based on Yan-Zheng's Table 1 -----------------------------------------
# t9010 = c(-3.48, 2.41)
t9010 = c(-1.62,1.58)
diff(pnorm(t9010))
(90-10)/diff(pnorm(t9010))

