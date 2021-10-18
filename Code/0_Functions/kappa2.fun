kappa2<-function(xtab)
{
  #computes Cohen's kappa and variance estimates, 95% CI
  #from a symmetric agreement table xtab
  #return value is a list with elements kap and vark
  #See Bishop et al. Disc. Mult. Analysis pp. 395-397 1975
  
  #Originally coded by C. Randin, UNIL
  #Adjusted by N.E.Zimmermann, WSL

  kap <- 0
  if(nrow(xtab) != ncol(xtab)) {
    #cat("\n freq. table not symmetric.\n")
    #return(kap, vark, totn)
  }
  # compute obs props:
  totn <- sum(xtab)
  # calc marginals, expected diag props.
  obs <- xtab/totn 
  pi <- apply(xtab, 1, sum)/totn
  pj <- apply(xtab, 2, sum)/totn
  exp <- outer(pi, pj)
  theta.2 <- sum(diag(exp))
  theta.1 <- sum(diag(obs))
  kap <- (theta.1 - theta.2)/(1 - theta.2)
  return(kap)
}
