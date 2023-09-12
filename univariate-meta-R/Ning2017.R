################################
## COPAS-LIKE SELECTION MODEL ##
################################
# Obtains maximum likelihood estimates based on the Copas-like 
# selection model in Ning et al. (2017)

# INPUTS:
# y = given vector of observed treatment effects
# s = given vector of within-study standard errors
# init = optional initialization values for (theta, tau, rho, gamma0, gamma1) 
#        If the user does not provide, these are estimated from the data.
# maxit = maximum number of iterations in the optimization

# OUTPUTS:
# theta.hat = estimate for theta
# tau.hat = estimate for tau
# rho.hat = estimate for rho
# gamma0.hat = estimate for gamma0
# gamma1.hat = estimate for gamma1
# H = estimate of Hessian matrix for (theta.hat, tau.hat, rho.hat, gamma0.hat, gamma1.hat)
#     Square root of diagonal entries can be used to estimate standard errors
# conv = convergence. "1"=converged, "0"=failed to converge

CopasLikeSelection <- function(y, s, init = NULL, tol=1e-20, maxit=1000){
  
  # Check that y and s are of the same length
  if(length(y) != length(s))
      stop("Please ensure that 'y' and 's' are of the same length.")
  # Check that all values in s are > 0
  if(any(s <= 0))
      stop("Please ensure that all observed standard errors 's' are greater than zero.")
  
  # Bind 'y' and 's'
  data = data.frame(cbind(y,s))
  max.s = max(s)
  
  # If no initial values were specified
  if(is.null(init)){
  
    ### for initial estimates of (1) theta (2) tau (3) rho (4) gamma0 and (5) gamma1
    init = rep(0,5)
    
    # initial estimate of theta
    init[1]=mean(data[,1])
    # moment estimator for tau
    init[2]=sqrt(abs(sd(data[,1])^2-(sd(data[,2])^2+mean(data[,2])^2)))
    
    rho.test=seq(-0.98,0.98,by=0.05)
    
    test = vector(mode = "list", length = length(rho.test))
    for(k in 1:length(test)){
      test[[k]]=optimal.gammas(data,theta=init[1],
                               tau=init[2],rho=rho.test[k])
    }
    # Find the rho that maximizes the log-likelihood
    loglik.max.index = 1
    loglik.max = as.double(test[[1]][2])
    for(k in 2:length(test)){
      if(as.double(test[[k]][2]) > loglik.max){
        loglik.max = as.double(test[[k]][2])
        loglik.max.index = k
      }
    }  
    # initial value for rho
    init[3] = rho.test[loglik.max.index]
    # initial value for gamma0
    init[4] = as.double(test[[loglik.max.index]][1])
    # initial value for gamma1
    init[5] = as.double(test[[loglik.max.index]][2])
  }
  
  # If initial values were specified but not of correct length
  if(!is.null(init)){
    if(length(init) != 5)
      stop("Please enter a vector of length 5 for initial values of (theta, tau, rho, gamma0, gamma1) in this order.")
    
    if(init[3]>=1){
      init[3] = 0.98
    } else if(init[3]<=-1){
      init[3] = -0.98
    }
  }

  par.new = init # Initialize values
  counter = 0    # start counter
  
  while (counter <= maxit) {
  
    counter = counter + 1
    par.old = par.new
    
    # E-step for EM algorithm
    myM= E.m(par.old,data)  

    # M-step for EM algorithm
    output = stats::optim(par.old,E.loglik,data=data,m=myM,method="L-BFGS-B",lower=c(-Inf,1e-3,-0.99,-Inf,0), 
                   upper=c(Inf,2,0.99,Inf,Inf), control = list(maxit = 1000),hessian=TRUE)
    
    
    par.new = output$par
    
    #cat("iter = ", counter, "P = ", par.new, "\n")
    
    if (max(abs(par.new-par.old)) < tol){  ###if convergence achieved
      
      # cat("\nSuccessfully Converged\n")
      # cat("iter = ", counter, "par = ", par.new, "\n")
      
      theta.hat = par.new[1]
      tau.hat = par.new[2]
      rho.hat = par.new[3]
      gamma0.hat = par.new[4]
      gamma1.hat = par.new[5]

      # Extract standard errors for theta and tau
      H = solve(output$hessian)
      
      return(list(theta.hat = theta.hat, 
                  tau.hat = tau.hat,
                  rho.hat = rho.hat,
                  gamma0.hat = gamma0.hat,
                  gamma1.hat = gamma1.hat,
                  H = H,
                  conv=1))
    
    }
  }
  
  # print("Convergence Failed")
  theta.hat = par.new[1]
  tau.hat = par.new[2]
  rho.hat = par.new[3]
  gamma0.hat = par.new[4]
  gamma1.hat = par.new[5]
  H = solve(output$hessian)

  return(list(theta.hat = theta.hat, 
              tau.hat = tau.hat,
              rho.hat = rho.hat,
              gamma0.hat = gamma0.hat,
              gamma1.hat = gamma1.hat,
              H = H,
              conv=0)) 
}  





##
## HELPER FUNCTIONS
##

######################
## HELPER FUNCTIONS ##
######################

# for the E-step in the EM algorithm for the Copas-like selection model
E.m=function(para,data){ # used for E step in EM 
  theta=para[1] ;   tau=para[2] ;   rho=para[3]
  gamma0=para[4];    gamma1=para[5]
  theta.hat=data[,1];s=data[,2]
  temp1=gamma0+gamma1/s
  temp2=rho*s*(theta.hat-theta)/(tau^2+s^2)
  denom = sqrt(1-rho^2*s^2/(tau^2+s^2))
  #calculating P[Z>0|theta] = \Phi[v_i]
  p1=pnorm((temp1+temp2)/denom)
  p2=pnorm(-(temp1+temp2)/denom)
  #calculating E[m]=(1-p)/p (expected # of unpublished studies)
  return(p2/p1)  
}

# the log likelihood function for the M-step in the EM algorithm
# for the Copas-like selection model
E.loglik=function(para,data,m){ 
  theta=para[1] ;   tau=para[2] ;   rho=para[3]
  gamma0=para[4];    gamma1=para[5]
  theta.hat=data[,1]; s=data[,2]
  junk1=gamma0 + gamma1/s + rho*s*(theta.hat-theta)/(tau^2+s^2)
  junk2=sqrt(1-rho^2*s^2/(tau^2+s^2))
  v=junk1/junk2
  ###bounding V?
  v[v < -37] <- -37
  v[v > 37] <- 37
  p1=pnorm(v)
  p2=pnorm(-v)
  ###full log likelihood 
  result=log(p1)+m*log(p2)-(m+1)/2*log(tau^2+s^2)-(m+1)/2*(theta.hat-theta)^2/(tau^2+s^2)
  return(-sum(result)) ###often prefer to minimize neg. log likelihood instead of max
}


# the log likelihood function for the standard random effects meta-analysis
RE.loglik=function(para,data){ 
  theta=para[1]    
  tau.sq=para[2]^2 
  y.obs = data[,1]
  s.obs.sq = data[,2]^2
  
  ## full log likelihood 
  result=-log(s.obs.sq+tau.sq)-(y.obs-theta)^2/(s.obs.sq+tau.sq) 
  return(-sum(result)) ###often prefer to minimize neg. log likelihood instead of max
}


## Estimate rho given other parameters (theta,tau,gamma0,gamma1)
optimal.gammas=function(data,theta,tau,rho){
  
  # Initial (gamma0,gamma1)
  gamma.params = c(-1,0.1)
  # Maximize w.r.t. (gamma0,gamma1)
  output = stats::optim(gamma.params,loglik.Copas,data=data,theta=theta,tau=tau,
                        rho=rho,method="L-BFGS-B",lower=c(-2,0.01), 
                        upper=c(1.5,0.9), control = list(maxit = 1000),hessian=TRUE)
  
  return(list(gamma0.hat=output$par[1],
              gamma1.hat=output$par[2],
              loglik=(-1)*output$val))
}

## Log-likelihood of original Copas selection model
loglik.Copas=function(gamma.params,theta,tau,rho,data){
  
  theta.hat=data[,1]; s=data[,2];
  gamma0=gamma.params[1]
  gamma1=gamma.params[2]
  
  ## Compute v
  junk1 = gamma0+(gamma1/s)+(rho*s*(theta.hat-theta))/(tau^2+s^2)
  junk2 = sqrt(1-(rho^2*s^2)/(tau^2+s^2))
  ## Bound v
  v=junk1/junk2
  v[v < -37] <- -37
  v[v > 37] <- 37
  
  term1 = -0.5*log(tau^2+s^2)
  term2 = -((theta.hat-theta)^2)/(2*(tau^2+s^2))
  term3 = -log(pnorm(gamma0+gamma1/s))
  term4 = log(pnorm(v))
  
  # Minimize the negative of the sum
  return(-sum(term1+term2+term3+term4))
}