## Pleiotropy for quantitative trait
## Author: DJ Schaid
## Date: 7/27/2016

pleio.q.fit <- function(y, geno){

  ## handle missing data
  
  nonmiss.geno <- !is.na(geno)
  nonmiss.y <- apply(!is.na(y), 1, all)
  nonmiss <- nonmiss.y & nonmiss.geno
  y <- y[nonmiss,]
  geno <- geno[nonmiss]
  
  n.subj <- nrow(y)
  n.traits <- ncol(y)

  
  ## center y and geno
  y <- center(y)
  geno <- geno - mean(geno)
  
  ## fit model without intercept, since x & y both centered
  fit <-  lm(y ~ -1 + geno) 
  beta.ols <- as.vector(fit$coeff)
  var.res <- var(fit$residuals)

  y.vec <- as.vector(y)
  
  ## cholesky of var.res
  lmat <- t(chol(var.res))
  linv <- solve(lmat)
  
  ident.n <- diag(n.subj)
  
  ## transform y and x by decorrelating
  
  x <- kronecker(linv, geno)
  y <- kronecker(linv, ident.n) %*% y.vec
  
  xx.inv <- solve(t(x) %*% x)

  obj.fit <- list(x=x,
                  xx.inv = xx.inv,
                  beta.ols = beta.ols,
                  n.traits = n.traits)

  return(obj.fit)
}
