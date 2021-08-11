# 2021 08 Andrew
# frequently used functions for bh with pub bias

# ==== BENJAMINI HOCHBERG PLUS ====

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




