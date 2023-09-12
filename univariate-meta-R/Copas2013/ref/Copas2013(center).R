Copas2013.center  <- function(y,x,p,start){
  
               n  <- length(y)
             y_m  <- mean(y)
       beta.range <- 10/(min(y)-max(y)) 

#Estimation with the reduced model
   loglik = function(param){
     
    theta <- param[ 1 ]
      tau <- param[ 2 ]
     beta <- param[ 3 ]

#Estimate alpha            
     func = function(alpha){
            sum(pnorm((alpha+beta*theta*x)/sqrt(1+beta^2*(1+tau^2*x^2)))^-1)-n/p}
    alpha = newton(func,0)$root

#likelihood of the reduced model
    ll <- -0.5*sum(log(1+tau^2*x^2))-
           0.5*sum((y-theta*x)^2/(1+tau^2*x^2))+
               sum(log(pnorm(alpha+beta*(y-y_m))))-
               sum(log(pnorm((alpha+beta*theta*x)/sqrt(1+beta^2*(1+tau^2*x^2)))))
    ll <- -ll
  
return(ll)
  }

result <- try(nlminb(start,objective=loglik,lower=c(-Inf,0,beta.range),upper=c(Inf,Inf,0)),silent=T)
if (class(result)!="try-error"){

#extract the logliklihood at MLE, note: here the function returns the -loglik
   MLE <- result$objective

#calculate confidence interval with the likelihood ratio test
            LR <- function(param){
              
            lr <- (MLE+2-loglik(param))^2
       return(lr)
            }
            
lb  <- try(nlminb(start=result$par,objective=LR,lower=c(-Inf,0,beta.range),upper=c(result$par[1],Inf,0))$par[1],silent=T)
ub  <- try(nlminb(start=result$par,objective=LR,lower=c(result$par[1],0,beta.range),upper=c(Inf,Inf,0))$par[1],silent=T)
lb  <- ifelse (class(lb)=="try-error",NA,lb)
ub  <- ifelse (class(ub)=="try-error",NA,ub)

est <- cbind(result$par[1],lb,ub)
 } else {est <- c(NA,NA,NA)}

return(est)
}

