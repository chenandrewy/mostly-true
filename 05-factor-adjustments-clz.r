# ABOUTME: Computes factor-adjusted statistics for CLZ predictors and produces LaTeX output.
# Inputs:
#   - config-and-functions.r (project helper functions)
#   - data/emp_data.Rdata (predictor returns from 01-prep-data.r)
#   - Internet: Ken French factor ZIP archives
# Outputs:
#   - data/F-F_Research_Data_Factors.csv
#   - data/F-F_Momentum_Factor.csv
#   - results/temp.tex
#   - results/vw-ffn-visual-raw.tex
# How to run:
#   Rscript 05-factor-adjustments-clz.r
#   Rscript 05-factor-adjustments-clz.r --vanilla

# Setup -----------------------------------------------------------------------

rm(list = ls())

library(here)
here::i_am("05-factor-adjustments-clz.r")

source(here("config-and-functions.r"))

paths <- project_paths()
data_dir <- paths$data
results_dir <- paths$results

temp_zip_path <- file.path(data_dir, "deleteme.zip")
ff3_csv_path <- file.path(data_dir, "F-F_Research_Data_Factors.csv")
momentum_csv_path <- file.path(data_dir, "F-F_Momentum_Factor.csv")
temp_tex_path <- file.path(results_dir, "temp.tex")
vw_ffn_visual_raw_path <- file.path(results_dir, "vw-ffn-visual-raw.tex")

load(file.path(data_dir, "emp_data.Rdata"))

# merge vw and ew data
ret1 = clz_ret %>% mutate(sweight = 'ew') %>%
rbind(
    clzvw_ret %>% mutate(sweight = 'vw') 
) %>% 
    select(sweight, signalname, date, ret) 

## Download FF + Carhart factors ====

# Download FF3
download.file(
    "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_CSV.zip",
    temp_zip_path
)
unzip(temp_zip_path, exdir = data_dir)

# Download Momentum
download.file(
    "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Momentum_Factor_CSV.zip",
    temp_zip_path
)
unzip(temp_zip_path, exdir = data_dir)

## read csvs and merge ====
ff_fac = fread(ff3_csv_path) %>% 
    left_join(fread(momentum_csv_path), by = 'V1') 
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
# Fnull = function(tmin) 2*(1-pnorm(tmin)) # better to use the ez 5% for consistency
F_null_ez = 0.05
tabdat = regest[ , .(Pr_gt_tmin = mean(abs(tstat) > tmin)), by = c('sweight','model')] %>% 
    left_join(storeydat, by = c('sweight','model')) %>% 
    mutate(FDRmax = F_null_ez/Pr_gt_tmin*pFmax
        , FDRmaxez = F_null_ez/Pr_gt_tmin) %>% 
    as_tibble() %>% 
    arrange(tmax, sweight, model) 

# Export latex ----------------------------------------------------

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
    name == 'Pr_lt_tmax' ~ paste0('Share of $|t_i| <', round(tmaxselect,1), '$ ')
    , name == 'pFmax' ~ '$\\Pr(\\nullt_i)$ Bound'
    , name == 'Pr_gt_tmin' ~ paste0('Share of $|t_i| >', tmin, '$ ')
    , name == 'FDRmax' ~ 'Visual Bound on $\\FDRez$'
    , name == 'FDRmaxez' ~ 'Easy Bound on $\\FDRez$'
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
  cat(file = temp_tex_path)

# read in temp.tex and modify manually
tex = readLines(temp_tex_path) 
tex[4] = ' & Raw & CAPM & FF3 & 4-Fac & & Raw & CAPM & FF3 & 4-Fac \\\\'
tex[5] = '\\cline{2-5} \\cline{7-10}'
tex = tex %>% append('\\\\ \\cline{2-5} \\cline{7-10}', after=7)
tex = tex  %>% append(' & \\multicolumn{4}{c}{Equal-Weighted}  & &  \\multicolumn{4}{c}{Value-Weighted} \\\\', after = 3) 

writeLines(tex, con = vw_ffn_visual_raw_path)




    
