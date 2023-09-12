gradient.function  <- function(param){
           
            mu     <- param[ 1 ]
           tau     <- param[ 2 ]
           rho     <- param[ 3 ]
        gamma0     <- param[ 4 ]
        gamma1     <- param[ 5 ] 
   
         ni.pub    <- ni[zi==1]
       ni.unpub    <- ni[zi==0]    
         ui.pub    <- gamma0 + gamma1 * sqrt ( ni.pub )
       ui.unpub    <- gamma0 + gamma1 * sqrt ( ni.unpub )
  ui.unpub[ui.unpub< -37] <- -37

     rho_tilde     <- si*rho/sqrt(tau^2+si^2)            
           vi      <- ( ui.pub + rho_tilde * ( yi - mu ) / sqrt ( tau^2 + si^2)) / sqrt ( 1 - rho_tilde^2 )
      vi[vi < -37] <- -37
        lambda     <- function(x){dnorm(x)/pnorm(x)}
          
        iterm1     <- tau^2 + si^2
        iterm2     <- yi - mu
#------------------------derivation of random effect model ------------------------#
  diff.tau         <- 2*tau*(-0.5*iterm1^-1+iterm2^2/(2*iterm1^2))
  diff.mu          <- iterm2/iterm1
#------------------------------derivation of vi for published studies--------------#
        vi.top     <- ui.pub+rho*si*(yi-mu)/iterm1
     vi.bottom     <- sqrt(1-rho_tilde^2)

    diff.vi.mu     <- vi.bottom^-1*(-rho_tilde/sqrt(iterm1))

diff.vi.top.tau    <- -2*tau*rho*si*iterm2/iterm1^2
diff.vi.bottom.tau <- (vi.bottom^-1)*rho_tilde*(tau*rho*si/iterm1^1.5)
   diff.vi.tau     <- (1-rho_tilde^2)^-1*(diff.vi.top.tau*vi.bottom-diff.vi.bottom.tau*vi.top)

diff.vi.top.rho    <- si*iterm2/iterm1
diff.vi.bottom.rho <- -(vi.bottom^-1)*rho*si^2/iterm1
    diff.vi.rho    <- vi.bottom^-2*(diff.vi.top.rho*vi.bottom-diff.vi.bottom.rho*vi.top)

       grad.mu     <- sum ( diff.mu + lambda(vi)*diff.vi.mu)
       
       grad.tau    <- sum ( diff.tau + lambda(vi)*diff.vi.tau)
       
       grad.rho    <- sum ( lambda(vi)*diff.vi.rho)   
       
    grad.gamma0    <- sum( lambda(vi)*vi.bottom^-1)-sum(dnorm(ui.unpub)/(1-pnorm(ui.unpub)))
    grad.gamma1    <- sum( sqrt(ni.pub)*lambda(vi)*vi.bottom^-1)-sum(sqrt(ni.unpub)*dnorm(ui.unpub)/(1-pnorm(ui.unpub)))

    gradient       <- c(-grad.mu,-grad.tau,-grad.rho,-grad.gamma0,-grad.gamma1)
    gradient
}

hessian  <- function(param){
           
            mu     <- param[ 1 ]
           tau     <- param[ 2 ]
           rho     <- param[ 3 ]
        gamma0     <- param[ 4 ]
        gamma1     <- param[ 5 ] 

   
         ni.pub    <- ni[zi==1]
       ni.unpub    <- ni[zi==0]    
         ui.pub    <- gamma0 + gamma1 * sqrt ( ni.pub )
       ui.unpub    <- gamma0 + gamma1 * sqrt ( ni.unpub )
  ui.unpub[ui.unpub< -37] <- -37

     rho_tilde     <- si*rho/sqrt(tau^2+si^2)            
           vi      <- ( ui.pub + rho_tilde * ( yi - mu ) / sqrt ( tau^2 + si^2)) / sqrt ( 1 - rho_tilde^2 )
      vi[vi < -37] <- -37
        lambda     <- function(x){dnorm(x)/pnorm(x)}
          
        iterm1     <- tau^2 + si^2
        iterm2     <- yi - mu
#------------------------derivation of random effect model ------------------------#
       diff2.mu    <- -1/iterm1
       diff.mu.tau <- -iterm2*2*tau/iterm1^2
       diff2.tau   <- -2*tau^2/iterm1^2-(iterm1-iterm2^2)/iterm1^2+4*tau^2*(iterm1-iterm2^2)/iterm1^3
#------------------------------derivation of vi for published studies--------------#
        vi.top     <- ui.pub+rho*si*(yi-mu)/iterm1
     vi.bottom     <- sqrt(1-rho_tilde^2)

    diff.vi.mu     <- vi.bottom^-1*(-rho_tilde/sqrt(iterm1))
   diff2.vi.mu     <- diff.vi.mu^2*lambda(vi)*(-vi-lambda(vi))
      hes.mu2      <- sum(diff2.mu+diff2.vi.mu)

diff.vi.top.tau    <- -2*tau*rho*si*iterm2/iterm1^2
diff.vi.bottom.tau <- (vi.bottom^-1)*rho_tilde*(tau*rho*si/iterm1^1.5)
   diff.vi.tau     <- (1-rho_tilde^2)^-1*(diff.vi.top.tau*vi.bottom-diff.vi.bottom.tau*vi.top)
diff.vi.mu.tau     <- (1-rho_tilde^2)^-1*(2*tau*rho*si/iterm1^2*vi.bottom+diff.vi.bottom.tau*(rho*si/iterm1))

diff.vi.top.rho    <- si*iterm2/iterm1
diff.vi.bottom.rho <- -(vi.bottom^-1)*rho*si^2/iterm1
    diff.vi.rho    <- vi.bottom^-2*(diff.vi.top.rho*vi.bottom-diff.vi.bottom.rho*vi.top)
diff.vi.mu.rho     <- (1-rho_tilde^2)^-1*(rho*si/iterm1*diff.vi.bottom.rho-vi.bottom*si/iterm1)


      hes.mu.tau   <- hes.tau.mu <- sum(-diff.vi.mu*diff.vi.tau*(vi*lambda(vi)+lambda(vi)^2)+
                                    lambda(vi)*diff.vi.mu.tau+diff.mu.tau)
     
      hes.mu.rho   <- hes.rho.mu <- sum(-diff.vi.mu*diff.vi.rho*(vi*lambda(vi)+lambda(vi)^2)+
                                    lambda(vi)*diff.vi.mu.rho)

      hes.mu.r0    <- hes.r0.mu  <- sum(-diff.vi.mu*vi.bottom^-1*(vi*lambda(vi)+lambda(vi)^2))

      hes.mu.r1    <- hes.r1.mu  <- sum(-diff.vi.mu*sqrt(ni.pub)/vi.bottom*(vi*lambda(vi)+lambda(vi)^2))

diff.vi.top.tau    <- -2*tau*rho*si*iterm2/iterm1^2 
diff.vi.bottom.tau <- (vi.bottom^-1)*rho_tilde*(tau*rho*si/iterm1^1.5)
   diff.vi.tau     <- (1-rho_tilde^2)^-1*(diff.vi.top.tau*vi.bottom-diff.vi.bottom.tau*vi.top)
   diff.vi.tau.p1 <- -2*rho*si*(rho*si*ui.pub+2*iterm2)*tau^2/(iterm1^3*vi.bottom^3)
diff.vi.tau.p2 <- -rho*si*((rho*si*ui.pub+2*iterm2)*tau^2+rho*si^3*ui.pub+(2-rho^2)*si^2*
                   yi+(mu*rho^2-2*mu)*si^2)/(iterm1^3*vi.bottom^3)
diff.vi.tau.p3 <- 6*rho*si*tau^2*((rho*si*ui.pub+2*iterm2)*tau^2+rho*si^3*ui.pub+(2-rho^2)*si^2*
                   yi+(mu*rho^2-2*mu)*si^2)/(iterm1^4*vi.bottom^3)
diff.vi.tau.p4 <- 3*rho^3*si^3*tau^2*((rho*si*ui.pub+2*iterm2)*tau^2+rho*si^3*ui.pub+(2-rho^2)*si^2*
                   yi+(mu*rho^2-2*mu)*si^2)/(iterm1^5*vi.bottom^5)
diff.vi.tau.com<- diff.vi.tau.p1+diff.vi.tau.p2+diff.vi.tau.p3+diff.vi.tau.p4
 diff2.vi.tau     <- diff.vi.tau^2*lambda(vi)*(-vi-lambda(vi))+lambda(vi)*(diff.vi.tau.com)

diff.vi.tau.rho.p1 <- -3*rho^3*si^4*tau*vi.top/(iterm1^3*vi.bottom^5)
diff.vi.tau.rho.p2 <- -3*rho^2*si^3*tau*iterm2/(iterm1^3*vi.bottom^3)
diff.vi.tau.rho.p3 <- -2*rho*si^2*tau*vi.top/(iterm1^2*vi.bottom^3)
diff.vi.tau.rho.p4 <- -2*si*tau*iterm2/(iterm1^2*vi.bottom)

diff.vi.tau.rho    <- diff.vi.tau.rho.p1+diff.vi.tau.rho.p2+
                      diff.vi.tau.rho.p3+diff.vi.tau.rho.p4

      hes.tau2     <- sum(diff2.tau+diff2.vi.tau)

   hes.tau.rho     <- hes.rho.tau <- sum(-diff.vi.tau*diff.vi.rho*(vi*lambda(vi)+lambda(vi)^2)+
                                    lambda(vi)*diff.vi.tau.rho)

   hes.tau.r0      <- hes.r0.tau  <- sum(-diff.vi.tau*vi.bottom^-1*(vi*lambda(vi)+lambda(vi)^2)+
                                    lambda(vi)*(-diff.vi.bottom.tau/vi.bottom^2))

   hes.tau.r1      <- hes.r1.tau  <- sum(-diff.vi.tau*sqrt(ni.pub)*vi.bottom^-1*(vi*lambda(vi)+lambda(vi)^2)+
                                    lambda(vi)*sqrt(ni.pub)*(-diff.vi.bottom.tau/vi.bottom^2))
 diff.vi.top.rho   <- si*iterm2/iterm1
diff.vi.bottom.rho <- -(vi.bottom^-1)*rho*si^2/iterm1
    diff.vi.rho    <- vi.bottom^-2*(diff.vi.top.rho*vi.bottom-diff.vi.bottom.rho*vi.top)
    diff.vi.rho2   <- si^2*ui.pub/(iterm1*vi.bottom^3)+3*rho*si^3*(si*ui.pub*rho+iterm2)/(iterm1^2*vi.bottom^5)

     hes.rho2      <- sum(diff.vi.rho^2*lambda(vi)*(-vi-lambda(vi))+lambda(vi)*diff.vi.rho2)
   hes.rho.r0 <- hes.r0.rho <- sum(-diff.vi.rho*vi.bottom^-1*(vi*lambda(vi)+lambda(vi)^2)+
                                   lambda(vi)*si^2*rho/(iterm1*vi.bottom^3))
   hes.rho.r1 <- hes.r1.rho <- sum(-diff.vi.rho*sqrt(ni.pub)*vi.bottom^-1*(vi*lambda(vi)+lambda(vi)^2)+
                                   lambda(vi)*sqrt(ni.pub)*si^2*rho/(iterm1*vi.bottom^3))

  diff.unp.r02     <- ui.unpub*dnorm(ui.unpub)/(1-pnorm(ui.unpub))-dnorm(ui.unpub)^2/(1-pnorm(ui.unpub))^2

       hes.r02     <- sum(vi.bottom^-2*lambda(vi)*(-vi-lambda(vi)))+sum(diff.unp.r02)
       hes.r12     <- sum(ni.pub*vi.bottom^-2*lambda(vi)*(-vi-lambda(vi)))+sum(diff.unp.r02*ni.unpub)

   hes.r0.r1 <- hes.r1.r0 <- sum(sqrt(ni.pub)*vi.bottom^-2*lambda(vi)*(-vi-lambda(vi)))+
                             sum(diff.unp.r02*sqrt(ni.unpub))


   hes.data <- c(hes.mu2,hes.mu.tau,hes.mu.rho,hes.mu.r0,hes.mu.r1,
                hes.tau.mu,hes.tau2,hes.tau.rho,hes.tau.r0,hes.tau.r1,
                hes.rho.mu,hes.rho.tau,hes.rho2,hes.rho.r0,hes.rho.r1,
                hes.r0.mu,hes.r0.tau,hes.r0.rho,hes.r02,hes.r0.r1,
                hes.r1.mu,hes.r1.tau,hes.r1.rho,hes.r1.r0,hes.r12)
   hessian  <- matrix(-hes.data,ncol=5,nrow=5)
   return(hessian)
}

    loglik    <- function(param){
        mu    <- param[1]
       tau    <- param[2]
       rho    <- param[3]
    gamma0    <- param[4]
    gamma1    <- param[5]

    ni.pub    <- ni[zi==1]
  ni.unpub    <- ni[zi==0]    
    ui.pub    <- gamma0 + gamma1 * sqrt ( ni.pub )
  ui.unpub    <- gamma0 + gamma1 * sqrt ( ni.unpub )

 rho_tilde    <- si * rho / sqrt ( tau^2 + si^2 )

        vi    <- ( ui.pub + rho_tilde * ( yi - mu ) / sqrt ( tau^2 +  si ^2 ) ) / sqrt ( 1 - rho_tilde^2 )
 vi[vi < -37] <- -37
        
        ll    <- sum ( -0.5 * log ( tau^2 + si^2 ) -( yi - mu ) ^2 / ( 2 * ( tau^2 +  si^2 ) ) + log ( pnorm ( vi )))+ 
                 sum ( log ( 1 - pnorm ( ui.unpub)))
        ll    <- -ll
               return(ll)
    }

