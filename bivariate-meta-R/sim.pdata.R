##******************************************************************************
##
## CREATE POPULATIN DATA
##
##******************************************************************************

# par <- c(20, 0, 1.73, 0.5, 0.5, -0.15, 0.5, 0.5); DT  <- sim.pdata(par)

.sim.pdata <- function(par){
  
  
  # par <- set[[1]][1,]
  
  S   <- par[1]
  mu1 <- par[2]
  mu2 <- par[3]
  
  tau1.sq <- par[4]
  tau2.sq <- par[5]
  tau12   <- par[6]
  
  c1.sq <- par[7]
  c2.sq <- par[8]
  
  c1 <- sqrt(c1.sq)
  c2 <- sqrt(c2.sq)
  
  rho1 <- par[9]
  rho2 <- par[10]
  
  v1 <- rnorm(S, abs(mu1)/2, 0.5)^2
  v2 <- rnorm(S, abs(mu2)/2, 0.5)^2
  s1 <- sqrt(v1)
  s2 <- sqrt(v2)
  
  mu <- c(mu1, mu2)
  
  ## Sigma+Omega
  
  SO <- lapply(1:S, function(s){
    
    matrix(c(v1[s]+tau1.sq, tau12, tau12, v2[s]+tau2.sq),2,2)
    
  })
  
  ## CHECK PD
  
  check.PD <- sapply(1:S, function(s){
    
    eigen(SO[[s]])$values
    
  })
  
  if (any(as.vector(check.PD)<= 0)) stop("VAR-COV MATRIX (S+O) is NOT PD")
  
  ## y FROM N(mu, SO)
  
  Y <- t(sapply(1:S, function(s) mvtnorm::rmvnorm(1, mu, SO[[s]])))
  
  ## SENS AND SPEC
  
  X <- plogis(Y)
  
  ## FINAL DATAFRMAE WITH NAME (y1, y2, v1, v2)
  
  ## AUGMENTED T-LNDOR
  ldor.t <- (c1*Y[,1]+c2*Y[,2])/sqrt(c1.sq*v1+c2.sq*v2)
  
  ## MU_IC,SIGMA_IC
  mu.c.obs <- sapply(1:S, function(s){
    
    Phi.12 <- matrix(c(rho1*s1[s], rho2*s2[s]),1,2)
    Phi.22 <- matrix(c(v1[s]+tau1.sq, tau12, tau12, v2[s]+tau2.sq),2,2)
    D <- matrix(Y[s,] - mu, 2, 1)
    Phi.12 %*% solve(Phi.22) %*% D
    
  })
  
  sigma.c <- sapply(1:S, function(s){
    
    Phi.12 <- matrix(c(rho1*s1[s], rho2*s2[s]),1,2)
    Phi.22 <- matrix(c(v1[s]+tau1.sq, tau12, tau12, v2[s]+tau2.sq),2,2)
    1- Phi.12 %*% solve(Phi.22) %*% t(Phi.12)
    
  })
  
  DT <- cbind.data.frame(Y, v1, v2, X, ldor.t, mu.c.obs, sigma.c)
  
  colnames(DT) <- c("y1", "y2", "v1", "v2", "se", "sp", "t.clnDOR", "mu.c.obs", "sigma.c")
  
  DT
}