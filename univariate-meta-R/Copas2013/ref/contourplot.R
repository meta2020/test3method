      Copas.plot  <- function(y,x,start){

       beta.range <- 10/(min(y)-max(y)) 

 alpha.all = seq(-0.25,2,0.01)
  beta.all = seq(beta.range,0,0.01)

   loglik2 = function(param){
    theta <- param[ 1 ]
      tau <- param[ 2 ]

    ll <- -0.5*sum(log(1+tau^2*x^2))-
           0.5*sum((y-theta*x)^2/(1+tau^2*x^2))+
               sum(log(pnorm(alpha+beta*y)))-
               sum(log(pnorm((alpha+beta*theta*x)/sqrt(1+beta^2*(1+tau^2*x^2)))))
    ll <- -ll
  
return(ll)
  }

mu = matrix(NA,length(alpha.all),length(beta.all))
for (i in seq(along=alpha.all)){
 for(j in seq(along=beta.all)){
  alpha = alpha.all[i]  
  beta  = beta.all[j]
  mu[i,j] = nlminb(start,loglik2,lower=c(-Inf,0),upper=c(Inf,0))$par[1]
 }
}
 contour(alpha.all,beta.all,mu,xlab=~alpha,ylab=~beta)
}
