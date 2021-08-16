# 2021 08 Andrew
# frequently used functions for bh with pub bias

# ==== FDR ESTIMATION  ====

run_bh_plus = function(t,C=1,qlist=c(0.01, 0.05, 0.10, 0.20)){
  
  Ns = length(t)
  psort = sort(2*pnorm(-t))
  fdrhat = psort/(1:Ns/Ns)*C
  
  istar = integer(length(qlist))
  for (qi in 1:length(qlist)){
    istar[qi] = sum(fdrhat < qlist[qi])
  }
  
  dt_fdr = data.table(
    psort = psort
    , tsort = sort(t, decreasing = T)
    , fdrhat = fdrhat
  )
  
  dt_hurdle = data.frame(
    fdr_ub = qlist
    , phurdle = psort[istar]  
    , thurdle = -1*qnorm(psort[istar]/2)
  )
  
  bh = list(
    fdr = dt_fdr
    , hurdle = dt_hurdle
  )
  
} # end function

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


# ==== PLOT TWO HISTGORAMS ====

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
      , f = hemp$density / sum(dat1>tgood) * length(t_emp)
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


# ==== PROB THEORY STUFF ====
# conditional expectation of folded norm
fe_normf = function(x){abs(x)*dnorm(x)} # function to integrate to get expectation
e_tnormf = function(xmin){ # integrate expecation
  integrate(f, xmin, Inf)
}


# ==== MODELS OF |t| ====

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


## setup
# declare mixture
dmix = function(t,pnull,scale){
  pnull*2/sqrt(2*pi*sigma^2)*exp(-t^2/(2*sigma^2)) + 
    (1-pnull)*dgamma(t,shape,1/scale)
}

# expectation of truncated mix
etmix = function(pnull,scale,tmin,ord=1){
  tempd = function(t){dmix(t,pnull,scale)}
  Fh = integrate(tempd, 0, tmin)$value
  intme = function(t){
    t^ord*tempd(t)/(1-Fh)
  }
  temp = integrate(intme,tmin,Inf)
  temp$value
}

# # debug: check expectation
# scale = 4
# etmix(pnull,scale,tgood)
# 
# nsim = 1e5
# tsim = rmix(nsim,pnull,scale)
# mean(tsim[tsim>tgood])

## estimate
ord = 1
minme = function(scale){
  (
    etmix(pnull,scale,tgood,ord) - mean(t_emp[t_emp>tgood]^ord)
  )^2
}

est = optimize(minme,c(0.1/shape,6/shape))

# find Pr(|t|>tgood)
tempd = function(t){dmix(t,pnull,est$minimum)}
Pr_tgood = integrate(tempd,tgood,Inf)$value

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
  Pr_tgood = pexp(tgood, 1/scalehat)
  
  fitexp = list(
    tgood = tgood
    , scalehat = scalehat
    , Pr_tgood = Pr_tgood
    , C = 1/Pr_tgood
  )
}