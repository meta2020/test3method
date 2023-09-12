rm(list = ls())

## 1. load package
library("mixmeta")
library("mada")

## Zhou et al. (2023) can be installed, or load at Step 2.
# devtools::install_github("meta2020/dtametasa")
# library(dtametasa)
# help(dtametasa.fc)

## 2. load all functions

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
## fit0 <- dtametasa.fc(test.data2, p=1, sauc.type = "hroc")
fit0 <- dtametasa.fc(test.data2, p=1, sauc.type = "sroc")
init <- fit0$par[1:5]
init
## see detailed estimates
# print(fit0)

## Same with function in mada package
fit.mada <- reitsma(data=test.data, correction.control = "single", method = "ml")
fit.mada

## likelihood bansed method, fix c contrast, c=(0.5, 0.5)
fit3 <- dtametasa.fc(test.data2, p=0.8, sauc.type = "sroc", c1.square = 0.5)
fit3$par
# print(fit3)

fit3.1 <- dtametasa.rc(test.data2, p=0.8, sauc.type = "sroc")
fit3.1$par
# print(fit3.1)

## Heckman-type conditional likelihood
fit1 <- Piao2019(test.data2, init = c(init, c(0.1,0.1,-0.1,0.1,0.1)))
fit1$par
# plogis(fit1$par[1:2])

## change other initial values 
fit1 <- Piao2019(test.data2, init = c(init, c(0.3, 0.4, -0.4, 0.5, 0.7)))
fit1$par
# plogis(fit1$par[1:2])


## Empirical full likelihood
fit2 <- Li2021(test.data2, eta1.tilde = fit1$par[8:10], eta0 = c(init, c(0.1,0.1,-0.1,0.1,0.1)))
fit2$par 
fit2$N
# plogis(fit2$par[1:2])

## change other initial values 
fit2.1 <- Li2021(test.data2, eta1.tilde = c(-0.4, 0.48, 0.66), eta0 = rep(0, 10))
fit2.1$par 
fit2.1$N
# plogis(fit2$par[1:2])



