
library("mixmeta")

## 2. load all functions
rm(list = ls())
files.sources <- list.files("bivariate-meta-R/")
sapply(paste0("bivariate-meta-R/", files.sources), source)

## load data
test.data <- read.csv("influenza.csv")
# test.data <- read.csv("mental.csv")


## contiuity correction: single = adding 0.5 to study with 0; all = adding 0.5 to all the studies
co.data <- correction(test.data, type = "single")
# co.data <- correction(test.data, type = "all")

## data after logit transformation
test.data2 <- logit.data(co.data)


# y1 <- test.data2$y1
# y2 <- test.data2$y2
# v1 <- test.data2$v1
# v2 <- test.data2$v2
# sigma1 <- sqrt(test.data2$v1)
# sigma2 <- sqrt(test.data2$v2)
# n <- nrow(test.data)

## 4 methods

## standard model, MLE estimate, p=1
fit0 <- dtametasa.fc(test.data2, p=1, sauc.type = "sroc")
init <- fit0$par[1:5]
init
print(fit0)


## likelihood bansed method, fix c contrast, c=(0.5, 0.5)
fit3 <- dtametasa.fc(test.data2, p=0.8, sauc.type = "sroc", c1.square = 0.5)
fit3$par
print(fit3)


## Heckman-type conditional likelihood
fit1 <- Piao2019(test.data2, init = c(init, c(0.1,0.1,-0.1,0.1,0.1)))
fit1$par
plogis(fit1$par[1:2])
## change another initial values 
fit1 <- Piao2019(test.data2, init = c(init, c(0.3, 0.4, -0.4, 0.5, 0.7)))
fit1$par
plogis(fit1$par[1:2])

## Empirical full likelihood
fit2 <- Li2021(test.data2, eta1.tilde = fit1$par[8:10], eta0 = c(init, c(0.1,0.1,-0.1,0.1,0.1)))
fit2$par 
fit2$N
plogis(fit2$par[1:2])


fit2 <- Li2021(test.data2, eta1.tilde = c(-0.4, 0.48, 0.66), eta0 = c(init, c(0.1,0.1,-0.1,0.1,0.1)))
fit2$par 
fit2$N
plogis(fit2$par[1:2])

