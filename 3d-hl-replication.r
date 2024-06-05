# Setup -----------------------------------------------------------------------

rm(list = ls())
source('0-functions.r')
load('../data/emp_data.Rdata')

# merge vw and ew data
ret1 = yz_ret %>% mutate(sweight = 'ew') %>%
rbind(
    yzvw_ret %>% mutate(sweight = 'vw') 
) %>% 
    select(sweight, signalname, date, ret) 
setDT(ret1)

## Download FF + Carhart factors ====

# Download FF3
download.file("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_CSV.zip"
    , '../data/deleteme.zip')
unzip('../data/deleteme.zip', exdir = '../data/')

# Download Momentum
download.file("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Momentum_Factor_CSV.zip"
    , '../data/deleteme.zip')
unzip('../data/deleteme.zip', exdir = '../data/')

## read csvs and merge ====
ff_fac = fread('../data/F-F_Research_Data_Factors.csv') %>% 
    left_join(fread('../data/F-F_Momentum_Factor.csv'), by = 'V1') 
colnames(ff_fac) = c('yearm','mktrf','smb','hml','rf','mom')

# clean
ff_fac = ff_fac %>% 
    mutate(year = yearm %/% 100, month = yearm %% 100) %>% 
    filter(!is.na(mktrf+smb+hml+mom)) %>% 
    select(year, month, mktrf, smb, hml, mom, rf)

# Find alphas (takes a couple minutes) ----------------------------------------

ret2 = copy(ret1)

# merge
ret2[ , c('year','month') := list(year(date), month(date))]
ret2[ff_fac, on = c('year','month')
    , `:=` (mktrf = i.mktrf, smb = i.smb, hml = i.hml, mom = i.mom, rf = i.rf)]
ret2[ , retrf := ret - rf]

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
    regest[[modelstr]] = ret2[, list(
        tstat = summary(lm(modelform, data = .SD))[['coefficients']]['(Intercept)', 't value']
    )
    , by = c('sweight','signalname')] %>% 
        mutate(model = modelnamelist[match(modelstr, modelstrlist)])
}
regest = do.call(rbind, regest)

regest = regest %>% 
    mutate(pval = 2*(1-pnorm(abs(tstat)))) %>% 
    select(sweight, signalname, model, tstat, pval)

# Calculate EZ FDR bound ------------------------------------------------------

regest2 = copy(regest)

# Calculate denominator 
setorder(regest2, sweight, model, pval)
regest2[ , Pr_pval2_lt_pval := (1:.N)/.N, by = c('sweight','model')]

# Calculate EZ Bound
regest2[ , FDR_ez := pval/Pr_pval2_lt_pval]

# BY 1.3 Bound ----------------------------------------------------------------

# find penalty
temp = regest2 %>% 
    group_by(sweight, model) %>%
    summarize(penalty = sum(1/(1:n())), .groups = 'drop') 

# merge on
temp2 = merge(regest2, temp, by = c('sweight','model')) %>% 
    mutate(FDRhat = FDR_ez * penalty) %>% 
    group_by(sweight, model) 

# find bounds for 5 and 10 percent FDRmax
for (FDRmax in c(0.05, 0.10)) {
    if(FDRmax == 0.05) tabout = data.table()
    temp3 = temp2 %>% 
        filter(FDRhat <= FDRmax) %>% 
        filter(row_number() == n()) %>% 
        mutate(pct_signif = 100*Pr_pval2_lt_pval
            , pct_FDRmax = 100*FDRmax) %>% 
        select(pct_FDRmax, sweight, model, pct_signif) %>% 
        setDT()
    tabout = rbind(tabout, temp3)
}

tab_BY13 = tabout
rm(list = ls(pattern = "temp"))

# Hand calculation of Storey --------------------------------------------------

# estimate pF
temp1 =  regest2 %>%
    group_by(sweight, model) %>%
    summarize(pval_gt_theta = mean(pval>0.8)) %>%
    ungroup() %>% 
    mutate(pF = pval_gt_theta/(1-0.8), pT = 1-pF) %>% 
    setDT()

# add pF to regest2 
temp2 = merge(regest2, temp1, by = c('sweight','model')) %>% 
    mutate(FDRhat = FDR_ez * pF) %>%
    group_by(sweight, model)

for (FDRmax in c(0.05, 0.10)) {
    if(FDRmax == 0.05) tabout = data.table()
    temp3 = temp2 %>% 
        filter(FDRhat <= FDRmax) %>% 
        filter(row_number() == n()) %>% 
        mutate(pct_signif = 100*Pr_pval2_lt_pval
            , pct_FDRmax = 100*FDRmax) %>% 
        select(pct_FDRmax, sweight, model, pct_signif) %>% 
        setDT()
    tabout = rbind(tabout, temp3)
}

tab_Storey = tabout %>% 
    left_join(temp1 %>% transmute(sweight, model, pct_alt = 100*pT)
        , by = c('sweight','model')) 
rm(list = ls(pattern = "temp"))

# Export latex ----------------------------------------------------

# rearrange
tab_all = tab_Storey %>% mutate(meth = 'Storey') %>% 
    bind_rows(
        tab_BY13 %>% mutate(meth = 'BY1.3')
    ) 

tab_signif = tab_all %>% 
    filter(model!='Raw') %>% select(-pct_alt) %>% 
    pivot_wider(names_from = c('sweight','model')
        , values_from = c('pct_signif')) %>% 
    arrange(pct_FDRmax, meth)

tab_alt = tab_all %>% 
    filter(model!='Raw', pct_FDRmax==5, meth=='Storey') %>% 
    select(-pct_signif) %>% 
    pivot_wider(names_from = c('sweight','model')
        , values_from = c('pct_alt')) %>% 
    arrange(pct_FDRmax, meth) %>% 
    mutate(pct_FDRmax = 'pct_alt')

tab = tab_signif %>% rbind(tab_alt) 

# add blank rows
tab2 = bind_rows(
    data.table(meth = 'pct_signif_05')
    , tab[1:2,]
    , data.table(meth = 'pct_signif_10')
    , tab[3:4,]
    , data.table(meth = 'pct_alt')
    , tab[5, ]
) %>% 
select(-pct_FDRmax) 

# export to latex
library(kableExtra)
tab2 %>% 
  kable('latex', booktabs = T, linesep = '', escape = F, digits = 1
  ) %>% 
  cat(file='../results/temp.tex')

# read in temp.tex and modify manually
tex = readLines('../results/temp.tex') 
tex[4] = ' & CAPM & FF3 & 4-fac & CAPM & FF3 & 4-fac \\\\'
tex[6] = '\\multicolumn{7}{l}{Percent Significant (FDR $\\le$ 5\\%)} \\\\ \\cline{1-7}'
tex[9] = '\\multicolumn{7}{l}{Percent Significant (FDR $\\le$ 10\\%)} \\\\ \\cline{1-7}'
tex[12] = '\\multicolumn{7}{l}{Minimum Percent \`\`True Strategies\'\' (Non-Null)} \\\\ \\cline{1-7}'

# replace names using grep
tex = tex %>% gsub('BY1.3', 'BY 2001 Thm 1.3', .)  %>% 
    gsub('Storey', 'Storey 2002', .)

# add ew vw header
tex = tex  %>% append(' & \\multicolumn{3}{c}{Equal-Weighted}  &  \\multicolumn{3}{c}{Value-Weighted} \\\\', after = 3) 

writeLines(tex, con='../results/hl-rep.tex')



