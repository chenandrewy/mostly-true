# 2021 08 Andrew
# frequently used functions for bh with pub bias

# PATHS AND LIBRARIES ====

dir.create('../data/', showWarnings = F)
dir.create('../results/', showWarnings = F)
dir.create('../results/sim-theory-free/', showWarnings = F)
dir.create('../results/sim-extrap-pub/', showWarnings = F)

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(latex2exp)

# RANDOM VARIABLES ====

# DENSITY
dmix = function(t,pnull,shape, scale){
  sigma = 1
  pnull*2/sqrt(2*pi*sigma^2)*exp(-t^2/(2*sigma^2)) + 
    (1-pnull)*dgamma(t,shape,1/scale)
}

# RV
rmix = function(n, pnull, shape, scale){
  inull = runif(n) < pnull
  t = rgamma(n, shape, 1/scale)
  t[inull] = abs(rnorm(sum(inull)))
  output = t
}

simmix = function(nsim,pnull,shape,scale,sigma){
  nnull = sum(runif(nsim) < pnull)
  t = abs(rnorm(nnull,0,sigma))  %>% 
    c(
      rgamma(nsim-nnull,shape,1/scale)
    )
  null = logical(nsim)
  null[1:nnull] = T
  p = 2*pnorm(-t)
  
  simdat = data.frame(
    t=t, null=null, p=p
  )
  
}



# expectation of truncated mix
etmix = function(pnull,shape,scale,tmin,ord=1){
  tempd = function(t){dmix(t,pnull,shape,scale)}
  Fh = integrate(tempd, 0, tmin)$value
  intme = function(t){
    t^ord*tempd(t)/(1-Fh)
  }
  temp = integrate(intme,tmin,Inf)
  temp$value
}



# FDR ESTIMATION  ====

estimate_fdr = function(
  t
  , tbarlist = seq(0,6,0.5)    
  , C = 1
  , nulldf = 200
  , null = numeric(length(t))*NA
){
  
  # initialize
  fdrlist = numeric(length(tbarlist))*NA
  drlist = fdrlist
  fdrhatlist = fdrlist
  
  # estimate for each tbar in tbarlist
  for (ti in 1:length(tbarlist)){
    i = abs(t)>tbarlist[ti]
    
    # discovery rate
    drlist[ti] = sum(i)/length(t)
    
    # find fdrs
    if (sum(i) > 0){
      fdrhatlist[ti] = 2*pt(-tbarlist[ti],nulldf)/drlist[ti]*C
      fdrlist[ti] = mean(null[i])
      fdrhatlist[ti] = min(fdrhatlist[ti], 1)
      
    } else{
      fdrlist[ti] = 0
      fdrhatlist[ti] = 0
    }
  }
  est = data.frame(
    tbar = tbarlist
    , dr = drlist
    , fdr_actual =  fdrlist
    , fdrhat = fdrhatlist 
  )
} # end function

estimate_fdr_parametric = function(
  pnull, shape, scale
  , tbarlist = seq(0,6,0.5)    
  , nulldf = 200
){
  
  # initialize
  drlist = numeric(length(tbarlist))*NA
  fdrhatlist = numeric(length(tbarlist))*NA
  fdrlist = numeric(length(tbarlist))*NA
  dmix1 = function(t){
    dmix(t, pnull, shape, scale)
  }

  # estimate for each tbar in tbarlist
  for (ti in 1:length(tbarlist)){

    # discovery rate
    drlist[ti] = integrate(dmix1,tbarlist[ti],Inf)$value
    
    # find fdrs
    if (drlist[ti] > 0){
      fdrhatlist[ti] = 2*pt(-tbarlist[ti],nulldf)/drlist[ti]
      fdrhatlist[ti] = min(fdrhatlist[ti], 1)
      
    } else{
      fdrhatlist[ti] = 0
    }
  }
  
  est = data.frame(
    tbar = tbarlist
    , dr = drlist
    , fdr_actual =  fdrlist
    , fdrhat = fdrhatlist 
  )
} # end function


estimate_mixture = function(
  t_emp, tgood
  , pnull, shape, sigma
){
  
  # |t| is mix of folded normal and gamma
  # folded normal becomes dirac delta siwh sigma = 0.001
  
  # with diract delta, pnull becomes a free parameter that has
  # absoultey no efefct on shape.  With a standard normal, however,
  # a larger pnull implies a larger scale, I think because 
  # a large pnull increases the densitiy of t-stats near tgood
  # and to keep that density near the data a higher scale is needed
  # to offset it
  
  # note sigma and shape are globals within this function
  
  ## estimate
  ord = 1
  minme = function(scale){
    (
      etmix(pnull,shape,scale,tgood,ord) - mean(t_emp[t_emp>tgood]^ord)
    )^2
  }
  
  est = optimize(minme,c(0.1/shape,6/shape))
  
  # check 
  # etmix(pnull,shape,est$minimum,tgood,ord)
  
  # find Pr(|t|>tgood)
  dmix1 = function(t){dmix(t,pnull,shape,est$minimum)}
  Pr_tgood = integrate(dmix1,tgood,Inf)$value
  
  ## pack results
  fitmix = list(
    tgood = tgood  
    , pnull = pnull
    , shape = shape
    , sigma =  sigma
    , scalehat = est$minimum
    , objhat = est$objective
    , Pr_tgood = Pr_tgood
    , C = 1/Pr_tgood
  )
  
} # end function

estimate_exponential = function(
  t_emp,tgood
){
  # this is super easy and hardly requires its own function
  # but I find it confusing if it doesn't have its own function
  scalehat = mean(t_emp[t_emp>tgood]) - tgood
  Pr_tgood = 1-pexp(tgood, 1/scalehat)
  
  fitexp = list(
    tgood = tgood
    , scalehat = scalehat
    , Pr_tgood = Pr_tgood
    , C = 1/Pr_tgood
  )
}

# SIMULATION ====




average_many_sims = function(
  emat, N, T_, exmethod, expar 
  , pnull, Emualt, tgood, smarg
  , nsim 
  , tgoodhat, pnullhat, shapehat, nulldf
  , tbarlist
){
  simmany = data.frame()
  for (simi in 1:nsim){
    print(paste0('simulation number ', simi))
    
    estat = estatsim(emat, N, T_, exmethod, expar)
    tdat = estat_to_tdat(N, T_, pnull, Emualt, estat)
    tdatpub = tdat_to_tdatpub(tdat, tgood = tgood, smarg = smarg )
    
    # actual fdr
    fdr_actual = estimate_fdr(
      tdat$t, tbarlist = tbarlist, null = tdat$null
    ) %>% 
      transmute(tbar, dr_actual = dr, fdr_actual)
    
    # fdrhat using exp
    bias_exp = estimate_mixture(tdatpub$t, tgood = tgoodhat, pnull = 0, shape = 1, 1)
    fdr_exp = estimate_fdr_parametric(
      pnull =  0, shape = 1, bias_exp$scalehat
      , tbarlist = tbarlist,nulldf = nulldf
    ) %>%
      transmute(tbar, fdrhat_exp = fdrhat)
    
    # fdrhat using mix
    bias_mix = estimate_mixture(tdatpub$t, tgood = tgoodhat, pnull = pnullhat, shape = shapehat, 1)
    
    fdr_mix = estimate_fdr_parametric(
      pnullhat, shapehat, bias_mix$scalehat
      , tbarlist = tbarlist,nulldf = nulldf
    ) %>%
      transmute(tbar, fdrhat_mix = fdrhat)
    
    # bind fdrs together
    est_all = fdr_actual %>% 
      left_join(fdr_exp, by = 'tbar') %>% 
      left_join(fdr_mix, by = 'tbar')
    
    est_all$simi = simi
    
    simmany = rbind(simmany, est_all)  
    
  } # end for simi
  
  # average across simulations
  manysum = simmany %>% 
    group_by(tbar) %>% 
    summarise_at(
      vars(2:(dim(simmany)[2])-2), mean, na.rm=F
    ) %>% 
    mutate(
      Ndisc = dr_actual*N
    )
  
} # end average_many_sims






# ==== PROB THEORY STUFF ====
# conditional expectation of folded norm
fe_normf = function(x){abs(x)*dnorm(x)} # function to integrate to get expectation
e_tnormf = function(xmin){ # integrate expecation
  integrate(f, xmin, Inf)
}





# TOOLS ====

detach_all = function(){
  invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
}

histcomp = function(
  dat1,dat2
  ,label1 = '1',label2 = '2'
  ,edge=seq(-3,15,0.5),tgood=-Inf
){
  t = dat1
  t = t[t>min(edge) & t<max(edge)]
  hemp = hist(t,edge)
  
  t = dat2
  t = t[t>min(edge) & t<max(edge)]
  hsim = hist(t,edge)
  
  plotme = rbind(
    data.frame(
      t = hemp$mids
      , f = hemp$density / sum(dat1>tgood) * length(dat1)
      , group = label1
    )
    , data.frame(
      t = hsim$mids
      , f = hsim$density / sum(dat2>tgood) * length(dat2)
      , group = label2
    )
  )
  
  p1 = ggplot(plotme, aes(x=t, y=f, fill=group)) +
    geom_bar(stat='identity', position='identity',alpha=0.6, show.legend = F) 
  
  p2 = ggplot(
    plotme %>% filter(t>tgood)
    , aes(x=t, y=f, fill=group)
  ) +
    geom_bar(stat='identity', position='identity',alpha=0.6, show.legend = F) 
  
  
  grid.arrange(p1, p2, nrow=1)  
} # end function



# AESTHETICS ====

library(latex2exp)
library(extrafont)


MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

NICEBLUE = "#619CFF"
NICEGREEN = "#00BA38"
NICERED = "#F8766D"

