
######################################  
######  load and prepare data  #######
######################################                                                                                

# load data
setwd("~/Documents/Uni/Master StatWi/Seeschwalben")
dat <- read.csv("slick2-h1-h2.csv")
colnames(dat)[1] <- "ID"

# load packages
library(HiddenMarkov)
library(moveHMM)

# overview
head(dat)
table(dat$ID) 

# for modelling, values exactly equal to 1 need to be changed by adding small quantity
which(dat$tortuosity==1)
ind <- which(dat$tortuosity==1)
dat$tortuosity[ind] <- rep(1.01,length(ind)) 

# prepare data
data <- prepData(dat,type="UTM")
data$log.tort <- log(data$tortuosity)


#######################################
## fit model without slick covariate ##
#######################################

# define function that calculates log-likelihood of 3-state HMM
mllk <- function(theta.star,x,N){
  Gamma <- diag(N)  
  Gamma[!Gamma] <- exp(theta.star[1:6])
  Gamma <- Gamma/rowSums(Gamma) 
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N)) 
  mu <- cumsum(exp(theta.star[7:9]))  
  sigma <- cumsum(exp(theta.star[10:12])) 
  allprobs <- matrix(1,dim(x)[1],N)
  ind <- which(!is.na(x$log.tort)) 
  allprobs[ind,] <- cbind(dgamma(x$log.tort[ind],shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]),
                          dgamma(x$log.tort[ind],shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]),
                          dgamma(x$log.tort[ind],shape=mu[3]^2/sigma[3]^2,scale=sigma[3]^2/mu[3]))
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:dim(x)[1]){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

# numerical likelihood maximisation 
theta.star <- c(rep(-2,6),log(c(0.01,0.02,0.03)),log(c(0.005,0.01,0.02)))
N=3
mod <- nlm(mllk,theta.star,x=data,N=N,print.level=2,iterlim=500,stepmax=10)
theta.star.mle <- mod$estimate

# backtransform working parameters from likelihood optimisation to actual scale of interest
Gamma <- diag(N)
Gamma[!Gamma] <- exp(theta.star.mle[1:6])
Gamma <- Gamma/rowSums(Gamma)
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
mu <- cumsum(exp(theta.star.mle[7:9]))
sigma <- cumsum(exp(theta.star.mle[10:12]))

# display parameter estimates
round(Gamma,4)
mu
sigma
delta


####################################################
## fit model including slick covariate (as dummy) ##
####################################################

# define function that calculates log-likelihood of 3-state HMM with slick covariate
mllk.slick <- function(theta.star,x){
  mu.tort <- cumsum(exp(theta.star[1:3]))
  sigma.tort <- cumsum(exp(theta.star[4:6]))
  delta <- c(1,exp(theta.star[7:8])/sum(c(1,exp(theta.star[7:8]))))
  beta <- cbind(theta.star[9:14],theta.star[15:20])
  allprobs <- matrix(1,dim(x)[1],3)
  ind <- which(!is.na(x$log.tort)) # indices of non-missing step lengths
  allprobs[ind,] <- cbind(dgamma(x$log.tort[ind],shape=mu.tort[1]^2/sigma.tort[1]^2,scale=sigma.tort[1]^2/mu.tort[1]),
                          dgamma(x$log.tort[ind],shape=mu.tort[2]^2/sigma.tort[2]^2,scale=sigma.tort[2]^2/mu.tort[2]),
                          dgamma(x$log.tort[ind],shape=mu.tort[3]^2/sigma.tort[3]^2,scale=sigma.tort[3]^2/mu.tort[3]))
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:dim(x)[1]){
    eta <- beta[,1]+beta[,2]*x[t,"slick"]
    Gamma <- diag(3)
    Gamma[!Gamma] <- exp(eta)
    Gamma <- Gamma/rowSums(Gamma)
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

# numerical likelihood maximisation 
N=3
delta0 <- c(0,0)
beta0 <- rep(-2,N*(N-1))
beta1 <- rep(0,N*(N-1))
theta.star <- c(log(mu),log(sigma),delta0,beta0,beta1)
mod_slick <- nlm(mllk.slick,theta.star,x=data,print.level=2,iterlim=500,stepmax=10)

# backtransform working parameters from likelihood optimisation to actual scale of interest
mu_slick <- cumsum(exp(mod_slick$estimate[1:3]))
sigma_slick <- cumsum(exp(mod_slick$estimate[4:6]))
beta_slick <- cbind(mod_slick$estimate[9:14],mod_slick$estimate[15:20])
delta_slick <- c(1,exp(mod_slick$estimate[7:8])/sum(c(1,exp(mod_slick$estimate[7:8]))))

# display parameter estimates
mu_slick
sigma_slick
beta_slick
delta_slick


######################
## model comparison ##
######################

# AIC
aic <- numeric(2)
aic[1] <- 2*mod$minimum + 2*length(mod$estimate)
aic[2] <- 2*mod_slick$minimum + 2*length(mod_slick$estimate)
aic
which.min(aic) 

# BIC
bic <- numeric(2)
T <- dim(data)[1]
bic[1] <- 2*mod$minimum + log(T)*length(mod$estimate)
bic[2] <- 2*mod_slick$minimum + log(T)*length(mod_slick$estimate)
bic
which.min(bic)

# LRT (H0: both models equally good, H1: the more complex model is much better)
# test statistic:
(-1)*2*mod_slick$minimum-(-1)*2*mod$minimum 
# chi-squared distribution, critical value:
qchisq(0.95,6)
# 146.2831 > 12.59159 -> reject H0


###############################################################
## main results with respect to covariate slick (best model) ##
###############################################################

# TPMs with and without slick
eta <- beta_slick[,1]+beta_slick[,2]*0
Gamma_slick0 <- diag(3)
Gamma_slick0[!Gamma_slick0] <- exp(eta) 
Gamma_slick0 <- Gamma_slick0/rowSums(Gamma_slick0)
round(Gamma_slick0,3)
eta <- beta_slick[,1]+beta_slick[,2]*1
Gamma_slick1 <- diag(3)
Gamma_slick1[!Gamma_slick1] <- exp(eta) 
Gamma_slick1 <- Gamma_slick1/rowSums(Gamma_slick1)
round(Gamma_slick1,3)

# stationary distribution with and without slick
delta0<-solve(t(diag(3)-Gamma_slick0+1),c(1,1,1))
round(delta0,3)
delta1<-solve(t(diag(3)-Gamma_slick1+1),c(1,1,1))
round(delta1,3)


#######################################
## state decoding (under best model) ##
#######################################

# function for state decoding using Viterbi algorithm
viterbi <- function(x,slick,mu,sigma,delta,beta){
  n <- length(x)
  allprobs <- matrix(1,n,3)
  ind <- which(!is.na(x))
  allprobs[ind,] <- cbind(dgamma(x[ind], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]),
                          dgamma(x[ind], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]), 
                          dgamma(x[ind], shape=mu[3]^2/sigma[3]^2,scale=sigma[3]^2/mu[3]))
  xi <- matrix(0,n,3)
  foo <- delta*allprobs[1,]
  xi[1,] <- foo/sum(foo)
  for (t in 2:n) {
    eta <- beta[,1]+beta[,2]*slick[t]
    Gamma <- diag(3)
    Gamma[!Gamma] <- exp(eta)
    Gamma <- Gamma/rowSums(Gamma)
    foo <- apply(xi[t-1,]*Gamma,2,max)*allprobs[t,]
    xi[t,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for (t in (n-1):1){
    eta <- beta[,1]+beta[,2]*slick[t]
    Gamma <- diag(3)
    Gamma[!Gamma] <- exp(eta)
    Gamma <- Gamma/rowSums(Gamma)
    iv[t] <- which.max(Gamma[,iv[t+1]]*xi[t,])
  }
  iv
}

HMMstate<-NULL
for(i in unique(data$ID)){
  x = data$log.tort[which(data$ID==i)]
  slick = data$slick[which(data$ID==i)]
  vitstates <- viterbi(x,slick,mu_slick,sigma_slick,delta_slick,beta_slick)
  HMMstate <- c(HMMstate,vitstates)
}
length(HMMstate)
data$HMMstate <- HMMstate
delta_vit<-prop.table(table(HMMstate))
write.csv(data,"slick2-h1-h2_with_states.csv",row.names=FALSE)


###############################################################################
## plot state-dependent distributions for best model (including slick dummy) ##
###############################################################################

# visualisation of estimated state-dependent distributions
hist(data$log.tort, probability = TRUE, breaks = 1500, col = "light grey", xlab = "", ylab = "",
     main = "", xlim = c(0,0.06), ylim = c(0,130)) 
title(main = "", cex.main = 2 , cex.lab = 1.3, xlab="log(tort)",ylab="Density")
colors <- c("#E69F00", "#009E73", "#56B4E9", "black")
curve(delta_vit[1] * dgamma(x, shape = mu_slick[1]^2/sigma_slick[1]^2,scale = sigma_slick[1]^2/mu_slick[1]), add = TRUE, col = colors[1], lwd = 2 , n=1001) 
curve(delta_vit[2] * dgamma(x, shape = mu_slick[2]^2/sigma_slick[2]^2,scale = sigma_slick[2]^2/mu_slick[2]), add = TRUE, col = colors[2], lwd = 2 , n=1001)
curve(delta_vit[3] * dgamma(x, shape = mu_slick[3]^2/sigma_slick[3]^2,scale = sigma_slick[3]^2/mu_slick[3]), add = TRUE, col = colors[3], lwd = 2 , n=1001)
curve(delta_vit[1] * dgamma(x, shape = mu_slick[1]^2/sigma_slick[1]^2,scale = sigma_slick[1]^2/mu_slick[1])
      + delta_vit[2] * dgamma(x, shape = mu_slick[2]^2/sigma_slick[2]^2,scale = sigma_slick[2]^2/mu_slick[2])
      + delta_vit[3] * dgamma(x, shape = mu_slick[3]^2/sigma_slick[3]^2,scale = sigma_slick[3]^2/mu_slick[3]),
      add = TRUE, col = colors[4], lwd = 2, lty = 2 , n=1001)
legend("topright", inset = c(0,0), legend = c("state 1: flying", "state 2: searching", "state 3: hovering/foraging", "marginal distribution"),pch = 16, col = colors, bty="n")



