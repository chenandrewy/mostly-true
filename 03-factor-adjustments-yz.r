#%% Setup -----------------------------------------------------------------------

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
load('../data/emp_data.Rdata')

# merge vw and ew data
ret1 = yz_ret %>% mutate(sweight = 'ew') %>%
rbind(
    yzvw_ret %>% mutate(sweight = 'vw') 
) %>% 
    select(sweight, signalname, date, ret) 
setDT(ret1)

## Download FF + Carhart factors 

# Download FF3
download.file("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_CSV.zip"
    , '../data/deleteme.zip')
unzip('../data/deleteme.zip', exdir = '../data/')

# Download Momentum
download.file("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Momentum_Factor_CSV.zip"
    , '../data/deleteme.zip')
unzip('../data/deleteme.zip', exdir = '../data/')

## read csvs and merge 
ff_fac = fread('../data/F-F_Research_Data_Factors.csv') %>% 
    left_join(fread('../data/F-F_Momentum_Factor.csv'), by = 'V1') 
colnames(ff_fac) = c('yearm','mktrf','smb','hml','rf','mom')

# clean
ff_fac = ff_fac %>% 
    mutate(year = yearm %/% 100, month = yearm %% 100) %>% 
    filter(!is.na(mktrf+smb+hml+mom)) %>% 
    select(year, month, mktrf, smb, hml, mom, rf)

# Find alphas (takes a couple minutes) ===

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

#%% Calculate Visual Bound ------------------------------------------------------

regest2 = copy(regest)

tab_viz = regest2 %>% 
  group_by(sweight, model) %>% 
  summarize(
    Pr_tgt_2 = mean(abs(tstat) > 2)
    , Pr_tlt_05 = mean(abs(tstat) < 0.5)
    , Pr_tlt_1 = mean(abs(tstat) < 1)
  ) %>% 
  ungroup() %>% 
  mutate(
    FDRmax_ez = 0.05/Pr_tgt_2
    , pFmax = Pr_tlt_05/0.38
    , pFmax_1 = Pr_tlt_1/0.68
    , FDRmax_viz = FDRmax_ez*pFmax
    , FDRmax_viz_1 = FDRmax_ez*pFmax_1
  )

#%% Calculate Storey following HL ------------------------------------------------------

regest2 = copy(regest)

# Calculate denominator 
setorder(regest2, sweight, model, pval)
regest2[ , Pr_pval2_lt_pval := (1:.N)/.N, by = c('sweight','model')]

# Calculate EZ Bound
regest2[ , FDR_ez := pval/Pr_pval2_lt_pval]

# Hand calculation of Storey ===

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

#%% Export latex ----------------------------------------------------

# convert to long and bind
temp1 = tab_viz %>% pivot_longer(cols = -c(sweight, model)) %>% 
  mutate(value = value*100)
temp2 = tab_Storey %>% select(-pct_alt) %>% pivot_wider(
    id_cols = c(sweight, model), names_from = pct_FDRmax, values_from = pct_signif
    , names_prefix = 'pct_signif_'
    ) %>% 
    pivot_longer(cols = -c(sweight, model)) 
tab = temp1 %>% bind_rows(temp2)

# make the nice wide table
stat_select = c('Pr_tgt_2', 'Pr_tlt_1'
    , 'pFmax', 'FDRmax_viz_1'
    ,'pct_signif_5', 'pct_signif_10')
stat_label = c('Share of $|t_i| > 2.0$', 'Share of $|t_i| < 1.0$'
    , '$\\Pr(\\nullt_i)$ Upper Bound', '$\\FDRez$ Upper Bound'
    , 'FDR $\\le$ 5\\%', 'FDR $\\le$ 10\\%')
model_select = c('CAPM', 'FF3', 'FF3+Mom')

tabwide = tab %>% filter(name %in% stat_select
    , model %in% model_select) %>% 
    mutate(name = factor(name, levels = stat_select, labels = stat_label)
        , model = factor(model, levels = model_select)) %>% 
    arrange(sweight, model, name) %>% 
    pivot_wider(names_from = c('sweight','model'), values_from = value) %>% 
    print()

# export to latex
library(kableExtra)
tabwide %>% 
  kable('latex', booktabs = T, linesep = '', escape = F, digits = 1
  ) %>% 
  cat(file='../results/temp.tex')

# read in temp.tex and modify manually
tex = readLines('../results/temp.tex') 

# add rows
tex = tex %>% append('& \\multicolumn{3}{c}{Equal-Weighted}  &  \\multicolumn{3}{c}{Value-Weighted} \\\\', after = 3)
tex = tex %>% append(paste(
    "\\vspace{-1.5ex} \\\\ "
    , "\\multicolumn{7}{l}{Panel (a): Visual Bound on $\\FDRez$}"
    , "\\\\ \\midrule"
    , sep = " "
), after = 6)
tex = tex %>% append(paste(
    "\\\\"
    , "\\multicolumn{7}{l}{Panel (b): Percent Significant using Storey (2002)}"
    , "\\\\ \\midrule"
    , sep = " "
), after = 11)

# edit rows
tex[5] = ' & CAPM & FF3 & 4-fac & CAPM & FF3 & 4-fac \\\\'
tex[6] = '\\cline{2-7}'

writeLines(tex, con='../results/yz-fdr.tex')




