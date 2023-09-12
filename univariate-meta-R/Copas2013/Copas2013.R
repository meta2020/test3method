##
## Copas 2013 Log-likelihood 
##
##----------------------------

library(metafor)
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

y <- dat$yi
v <- dat$vi

Copas2013  <- function(y, v, p){ ## y is t-type statistic
               
  # n  <- length(y)
  b.range <- 10/(min(y)-max(y)) 

## OPTIMIZATION USING POPULATION MODEL
  
  llk.p <- function(par) {
    
    u <- par[1]
    t2 <- par[2]
    
    ll <- -0.5 * ( sum( log((v + t2)), na.rm = TRUE ) + sum( (y-u)^2 / (v + t2), na.rm = TRUE )) 
    
    -ll
  }
  
  int <- try(nlminb(start = c(0, 0.1), objective = llk.p, 
                       lower=c(-Inf, .Machine$double.eps), upper = c(Inf, Inf)), 
                silent=T)$par
  
  int <- round(int,2)
  
  
## OPTIMIZATION USING OBSERVED MODEL 
  llk.o <- function(par){
    
    u <- par[1]
    t2 <- par[2]
    b <- par[3]

    t <- y / sqrt(v)

    ## ESTIMATE ALPHA    
    
    f.a <- function(a) mean(1 / pnorm((a + b * (u/sqrt(v))) / sqrt( 1 + b^2 * (1 + t2 /v)) ), na.rm = TRUE) - 1/p
    
    a.p <- suppressWarnings(pracma::newton(f.a, 0)$root)

    ## FINAL LOG-LIKELIHOOD
    
    ll <- -0.5 * ( sum( log( (v + t2)), na.rm = TRUE ) + sum( (y-u)^2 / (v + t2), na.rm = TRUE)) + 
      sum( log(pnorm(a.p + b * t)), na.rm = TRUE) - 
      sum( log(pnorm((a.p + b * (u/sqrt(v))) / sqrt( 1 + b^2 * (1 + t2 /v)) )), na.rm = TRUE)
    
    -ll
  
  }

result <- try(nlminb(c(int, 0),  llk.o, lower=c(-Inf, 0 , -1), upper = c(Inf,5,1)), silent=T)

result$initial <- int

result

}


Copas2013(y, v, 1)



       
       
