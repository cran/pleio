
## Pleiotropy for quantitative trait, sequential test of 0, then 1, etc.
## Author: DJ Schaid
## Date: 7/27/2016
pleio.q.sequential <- function(obj.fit, pval.threshold){
  pval <- pval.threshold / 2
  n.traits <- obj.fit$n.traits
 
  count <- 0
  save <- NULL
  while(pval < pval.threshold & count < n.traits){
    save <- pleio.q.test(obj.fit, count.nonzero.beta=count)
    pval <- save$pval
    index.beta <- save$index.nonzero.beta
    count <- count + 1
  }

  ## decrement count to account for "+1" in above loop, in case
  ## pval > pval.threshold when count === 0
  
  count <- count - 1
  
  return(list(pval=pval, index.beta=index.beta))
}
