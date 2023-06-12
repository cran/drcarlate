## ---- include = FALSE---------------------------------------------------------
library(drcarlate)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# Parameter dgpflag declares three ways to generate random data.
# FuncDGP will generate random data according to the method specified in (i), (ii) or (iii) when dgpflag = 1, 2 or 3 respectively.

# Parameter rndflag declares four ways to randomly assign treatment effects.
# rndflag = 1 - SRS; rndflag = 2 - WEI; rndflag = 3 - BCD; rndflag = 4 - SBR
# Note that CovadpRnd is built into FuncDGP, so it is not necessary to use CovadpRnd alone in the actual operation of generating random data

# Let's take dgpflag=1 and rndflag=2 for example
random_dgp <- FuncDGP(dgptype = 1, rndflag = 2, n = 100, g = 4, pi = c(0.5, 0.5, 0.5, 0.5))
# We can see that the return value of FuncDGP is a list of nine matrices, We can easily extract what we need from it
Y <- random_dgp$Y
X <- random_dgp$X
S <- random_dgp$S
A <- random_dgp$A
Y1 <- random_dgp$Y1
Y0 <- random_dgp$Y0
D1 <- random_dgp$D1
D0 <- random_dgp$D0
D <- random_dgp$D

## -----------------------------------------------------------------------------
# compute estimated LATE
tauhat <- tau(muY1 = Y1, muY0 = Y0, muD1 = D1, muD0 = D0, A = A, S = S, Y = Y, D = D)

#compute estimated treatment assignment probabilities
pihat(A = A, S = S)

# compute estimated standard deviation
stanE(muY1 = Y1, muY0 = Y0, muD1 = D1, muD0 = D0, A = A, S = S, Y = Y, D = D, tauhat = tauhat)

## ---- message=FALSE,results='hide'--------------------------------------------
# let's take dgpflag = 1 for example.
true_value <- TrueValue(dgptype = 1, vIdx = 1:4, n = 100, g = 4, pi = c(0.5, 0.5, 0.5, 0.5))
true_tau <- true_value$tau

## -----------------------------------------------------------------------------
# SRS - WEI - BCD - SBR
true_tau

## -----------------------------------------------------------------------------
# remember that we set dgpflag = 1 and rndflag = 2 before
LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = 4, modelflag = 1, iridge = 0.001)

## -----------------------------------------------------------------------------
# set random seed
set.seed(1)

## ---- results='hide'----------------------------------------------------------
# get true tau
true_tau <- TrueValue(dgptype = 1, vIdx = 1:4, n = 1000, g = 4, pi = c(0.5, 0.5, 0.5, 0.5))

## ---- warning=FALSE-----------------------------------------------------------
# see the output: size, iPert = 0
Output(ii = 1, tau = tauhat[1], dgptype = 1, rndflag = 1, n = 1000, g = 4, 
       pi = c(0.5, 0.5, 0.5, 0.5),
       iPert = 0, iq = 0.05, iridge = 0.001)

# see the output: power, iPert = 1
Output(ii = 1, tau = tauhat[1], dgptype = 1, rndflag = 1, n = 1000, g = 4, 
       pi = c(0.5, 0.5, 0.5, 0.5),
       iPert = 1, iq = 0.05, iridge = 0.001)


## -----------------------------------------------------------------------------
# For example, if we wanted to get the results shown in Panel A in Table 1, we could run the following two functions separately.

# First of all, There are four strata data (g = 4), the probability of data being treated at each strata is equal to 0.5 (pi = c(0.5, 0.5, 0.5, 0.5)), the random data generation process follows the data generation process 1 (DGP = 1), the total sample size is 200 (n = 200), run 10,000 Monte Carlo simulations (iMonte = 10000), the confidence level for hypothesis testing is 5%, and finally，get the size of the Monte Carlo simulation (iPert = 0).
# Time Consumeing. This command will take about 1.4 hours to execute, so do not run it unless necessary.

# JLTZ(iMonte = 10000, dgptype = 1, n = 200, g = 4, pi = c(0.5, 0.5, 0.5, 0.5), iPert = 0, iq = 0.05, iridge = 0.001)

# Second of all, There are four strata data (g = 4), the probability of data being treated at each strata is equal to 0.5 (pi = c(0.5, 0.5, 0.5, 0.5)), the random data generation process follows the digital generation process 1 (DGP = 1), the total sample size was 200 (n = 200), run 10,000 Monte Carlo simulations (iMonte = 10000), the confidence level for hypothesis testing is 5%, and finally，get the power of the Monte Carlo simulation (iPert = 1).
# Time Consumeing. This command will take about 1.4 hours to execute, so do not run it unless necessary.

# JLTZ(iMonte = 10000, dgptype = 1, n = 200, g = 4, pi = c(0.5, 0.5, 0.5, 0.5), iPert = 1, iq = 0.05, iridge = 0.001)

## -----------------------------------------------------------------------------
# With the above Settings in mind, we get the following JLTZ function.

# JLTZ(iMonte = 10000, dgptype = 1, n = 200, g = 4, pi = c(0.5,0.5,0.5,0.5), iPert = 0,iq = 0.05, iridge = 0.001)

# Since it takes about 1.4 hours to run this function, we'll give you the results directly that users can check them out themselves.

#        vProb_d1   vProb_d2   vProb_d3   vProb_d4
# [1,] 0.03139898 0.03548546 0.03071509 0.03139898
# [2,] 0.04356403 0.04189256 0.03951368 0.04356403
# [3,] 0.04307085 0.04238541 0.03919373 0.04307085
# [4,] 0.05622226 0.04632824 0.05327148 0.05622226
# [5,] 0.10833470 0.08854937 0.09486482 0.10833470
# [6,]        NaN        NaN        NaN        NaN
# [7,] 0.03435805 0.03334976 0.03455447 0.03435805
# [8,] 0.04586553 0.05076392 0.05087186 0.04586553



## -----------------------------------------------------------------------------

# Set up ------------------------------------------------------------------
library(drcarlate)
library(pracma)

# Load data ---------------------------------------------------------------
# data_for_final.csv add covariates b_total_income log_b_total_income b_exp_total_30days
# edu_index b_asset_index b_health_index to data_for_tab4.csv

data_table <- drcarlate::data_table

data1 <- as.matrix(data_table)

colnames(data1) <- NULL

# create S for wave 1: number of stratnum, total 41 stratnum
n <- size(data1,1)
S <- zeros(n,1)
for (i in 1:41) {
  S[data1[, size(data1,2)-i+1] == 1] <- 41 - i + 1
}

data1 <- cbind(data1, S)

# create Y, X, S, D: number of outcomes considered is 9:
# save the results for NA, TSLS, L, NL, F
vtauhat <- NaN * ones(8,5)
vsighat <- NaN * ones(8,5)
n_vec <- NaN * ones(8,1)

iridge = 0.01 #tuning parameter for in ridge regression

for (i in 1:8) {
  if (i <= 4) {
  # data_used: 1st col- Y, 2nd col- A, 3rd col- D, 4~(end-1)th col- X,
  # last col is strata number
    data_used <- cbind(data1[, 4+(i-1)*3], data1[, 2:3], data1[, 5:6],
                       data1[, size(data1,2)-1-41], data1[, size(data1,2)])
  } else {
    data_used <- cbind(data1[, 4+(i-1)*3], data1[,2:3], data1[, (4+(i-1)*3+1):(4+(i-1)*3+2)],
                       data1[, size(data1,2)-1-41], data1[, size(data1,2)])
  }

  # delete the outcome variables with NaN
  data_used <- data_used[!is.nan(data_used[,1]), ]

  # check whether a strata has obs less than 10
  for (s in 1:41) {
    if (length(data_used[data_used[, size(data_used,2)] == s,1]) < 10) {
      # delete strata has obs less than 10
     data_used[data_used[, size(data_used,2)] ==s, ] <- matrix()
     data_used <- data_used[!is.na(data_used[,1]),]
    }
  }

  # update n
  n <- size(data_used,1)
  n_vec[i] <- n

  # find the unique value for stratnum
  stratnum <-  sort(base::unique(data_used[, size(data_used,2)]))

  # create A
  A <-  matrix(data_used[, 2])

  # create D
  D <- matrix(data_used[, 3])

  # create Y
  Y <- matrix(data_used[, 1])

  # create S
  S <- matrix(data_used[, size(data_used,2)])

  # create X
  # use baseline total income
  X <- data_used[,4:6]

  # make a index
  print(stringr::str_c("Now i equals to ", i, " !"))

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  %% No adjustment (z) NA
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  muY0_z <- zeros(n,1)
  muY1_z <- zeros(n,1)
  muD0_z <- zeros(n,1)
  muD1_z <- zeros(n,1)

  vtauhat[i,1] <- tau(muY1 = muY1_z, muY0 = muY0_z, muD1 = muD1_z, muD0 = muD0_z,
                      A = A, S = S, Y = Y, D = D, stratnum = stratnum)
  vsighat[i,1] <- stanE(muY1 = muY1_z, muY0 = muY0_z, muD1 = muD1_z, muD0 = muD0_z,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[i,1], stratnum = stratnum)/sqrt(n)

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  %% TSLS
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mS_iv <- zeros(n, length(stratnum))

  for (j in 1:length(stratnum)) {
    s <- stratnum[j]
    mS_iv[S==s, j] <- 1
  }

  mX_iv <- cbind(D, X, mS_iv)
  mZ_iv <- cbind(A, X, mS_iv)
  K <- size(mX_iv,2)

  vPara_iv <- inv(t(mZ_iv) %*% mX_iv) %*% t(mZ_iv) %*% Y
  mE_iv2 <- diag(((Y - mX_iv %*% vPara_iv)^2)[,])
  m0_iv <- inv(t(mZ_iv) %*% mX_iv/n) %*% (t(mZ_iv) %*% mE_iv2 %*% mZ_iv/n) %*% inv(t(mX_iv) %*% mZ_iv/n)
  vtauhat[i,2] <- vPara_iv[1]
  vsighat[i,2] <- sqrt(m0_iv[1,1])/sqrt(n)

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  %% Linear+linear model (L)
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  muY0_a <- NaN*ones(n,1)
  muY1_a <- NaN*ones(n,1)
  muD0_a <- NaN*ones(n,1)
  muD1_a <- NaN*ones(n,1)
  for (j in 1:length(stratnum)) {
    s = stratnum[j]

    result <- LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = s, modelflag = 1, iridge = iridge)

    theta_0s_a <- result[["theta_0s"]]
    theta_1s_a <- result[["theta_1s"]]
    beta_0s_a <- result[["beta_0s"]]
    beta_1s_a <- result[["beta_1s"]]

    muY0_a[S==s] <- X[S==s,] %*% theta_0s_a
    muY1_a[S==s] <- X[S==s,] %*% theta_1s_a
    muD0_a[S==s] <- X[S==s,] %*% beta_0s_a
    muD1_a[S==s] <- X[S==s,] %*% beta_1s_a
  }

  vtauhat[i,3] <- tau(muY1 = muY1_a, muY0 = muY0_a, muD1 = muD1_a, muD0 = muD0_a,
                      A = A, S = S, Y = Y, D = D, stratnum = stratnum)
  vsighat[i,3] <- stanE(muY1 = muY1_a, muY0 = muY0_a, muD1 = muD1_a, muD0 = muD0_a,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[i,3], stratnum = stratnum)/sqrt(n)
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  %% Linear+logistic model (NL)
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  muY0_b <- NaN*ones(n,1)
  muY1_b <- NaN*ones(n,1)
  muD0_b <- NaN*ones(n,1)
  muD1_b <- NaN*ones(n,1)

  for (j in 1:length(stratnum)) {
    s = stratnum[j]

    result <- LinearLogit(Y = Y, D = D, A = A, X = X, S = S, s = s, modelflag = 2, iridge = iridge)

    theta_0s_b <- result[["theta_0s"]]
    theta_1s_b <- result[["theta_1s"]]
    beta_0s_b <- result[["beta_0s"]]
    beta_1s_b <- result[["beta_1s"]]

    muY0_b[S==s] <- X[S==s,] %*% matrix(theta_0s_b[2:length(theta_0s_b)]) + theta_0s_b[1]
    muY1_b[S==s] <- X[S==s,] %*% matrix(theta_1s_b[2:length(theta_0s_b)]) + theta_1s_b[1]
    muD0_b[S==s] <- LogisticReg(x = (X[S==s,] %*% matrix(beta_0s_b[-1]) + beta_0s_b[1]))
    muD1_b[S==s] <- LogisticReg(x = (X[S==s,] %*% matrix(beta_1s_b[-1]) + beta_1s_b[1]))
  }

  vtauhat[i,4] <- tau(muY1 = muY1_b, muY0 = muY0_b, muD1 = muD1_b, muD0 = muD0_b,
                      A = A, S = S, Y = Y, D = D, stratnum = stratnum)
  vsighat[i,4] <- stanE(muY1 = muY1_b, muY0 = muY0_b, muD1 = muD1_b, muD0 = muD0_b,
                        A = A, S = S, Y = Y, D = D,tauhat = vtauhat[i,4], stratnum = stratnum)/sqrt(n)

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  %% Further Efficiency Improvement
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  X_c <- cbind(X, muD1_b, muD0_b)
  muY0_c <- NaN*ones(n,1)
  muY1_c <- NaN*ones(n,1)
  muD0_c <- NaN*ones(n,1)
  muD1_c <- NaN*ones(n,1)

  for (j in 1:length(stratnum)) {
    s = stratnum[j]

    result <- LinearLogit(Y = Y, D = D, A = A, X = X_c, S = S, s = s, modelflag = 1, iridge = iridge)

    theta_0s_c <- result[["theta_0s"]]
    theta_1s_c <- result[["theta_1s"]]
    beta_0s_c <- result[["beta_0s"]]
    beta_1s_c <- result[["beta_1s"]]

    muY0_c[S==s] <- X_c[S==s,] %*% theta_0s_c
    muY1_c[S==s] <- X_c[S==s,] %*% theta_1s_c
    muD0_c[S==s] <- X_c[S==s,] %*% beta_0s_c
    muD1_c[S==s] <- X_c[S==s,] %*%beta_1s_c
  }

  vtauhat[i,5] <- tau(muY1 = muY1_c, muY0 = muY0_c, muD1 = muD1_c, muD0 = muD0_c,
                      A = A, S = S, Y = Y, D = D, stratnum = stratnum)
  vsighat[i,5] <- stanE(muY1 = muY1_c, muY0 = muY0_c, muD1 = muD1_c, muD0 = muD0_c,
                        A = A, S = S, Y = Y, D = D, tauhat = vtauhat[i,5], stratnum = stratnum)/sqrt(n)
  }

# show the LATE estimates and standard errors in Table 5 of Jiang et al.(2022). 
# LATE Estimates
vtauhat

# Standard Errors
vsighat

