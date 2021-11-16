#### Team-Seeschwalbe

## Daten einlesen
dat <- read.csv("slick2-h1-h2.csv")
## Paket laden
library(moveHMM)
## Daten anzeigen
View(dat)
head(dat)

## Anzahl der Beobachtungen pro Tern
table(dat$ID)

## neue Spalte im Datensatz erstellen mit der Frequency
data2 <- dat$ID
x <- data.frame(data2, freq=ave(seq_along(data2), data2, FUN=length))
dat$freq <- x$freq
View(dat)

which(data_frame$Freq >400) #welche Terns haben mehr als 400 Beobachtungen
sum(data_frame$Freq >400) #wie viele Terns haben mehr als 400 Beobachtungen

## Datensatz nur mit Terns, die >400 bzw. >700 Beobachtungspunkte besitzen
dat_400 <- subset(dat, dat$freq > 400)
dat_700 <- subset(dat, dat$freq > 700)

## Daten aufbereiten
data <- prepData(dat_700, type="UTM")

## 2 State HMM fitten
?fitHMM
# mit turning anlge
mod <- fitHMM(data,nbStates=2, stepPar0 = c(0.1,0.25, 0.05, 0.05), anglePar0 = c(0,0,1,1),verbose=1, angleDist="vm",stationary=T)
plot(mod)
mod
# ohne turning angle (läuft viel schneller durch!!)
mod2 <- fitHMM(data,nbStates=2, stepPar0 = c(0.1,0.25, 0.05, 0.05), angleDist="none", verbose=1,stationary=T)
plot(mod2)
mod2


### Modell manuell fitten
mllk <- function(theta.star, x, N){
  Gamma <- diag(N) #identity matrix 
  Gamma[!Gamma] <- exp(theta.star[1:2]) #turn zeros in matrix (off-diagonal entries) into exponentials of theta-stars
  Gamma <- Gamma/rowSums(Gamma) #devide by row sums so that entries of a row sum to 1 (transformation is inverse logit link)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N)) #stationary distribution erstellen
  mu.step <- exp(theta.star[3:4]) #Mittelwert
  sigma <- exp(theta.star[5:6]) #Varianz
  allprobs <- matrix(1,dim(x)[1],N)
  ind <- which(!is.na(x$step)) # indices of non-missing step lengths
    allprobs[ind,] <- cbind(dgamma(x$step[ind],shape=mu.step[1]^2/sigma[1]^2,scale=sigma[1]^2/mu.step[1]),
                            dgamma(x$step[ind],shape=mu.step[2]^2/sigma[2]^2,scale=sigma[2]^2/mu.step[2]))
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

N=2

# finding starting values
hist(data$step, breaks=20, main = "Histogram of the step lengths", xlab = "step lengths", xlim=c(0,0.45))
# starting values
theta.star <- c(rep(-2,(N-1)*N),log(c(0.15,0.25)),log(c(0.05,0.03)))
# numerical maximization 
mod.manuell <- nlm(mllk,theta.star,x=data,N=N,print.level=1,iterlim=500)
theta.star.mle = mod.manuell$estimate

# back-transformation
Gamma <- diag(N)
Gamma[!Gamma] <- exp(theta.star.mle[1:((N-1)*N)])
Gamma <- Gamma/rowSums(Gamma)
delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
mu.step <- exp(theta.star.mle[(N-1)*N+1:N])
sigma <- exp(theta.star.mle[(N-1)*N+(N+1):(2*N)])

round(Gamma, 4) #Matrix auf 4 Nachkommastellen gerundet
mu.step
sigma
delta

## grafische Darstellung
# State dependent distributions
hist(data$step,probability=TRUE,breaks=100,col="light grey",xlab="step length",main="State-dependent distributions", xlim= c(0,0.45))
curve(delta[1]*dgamma(x,shape=mu.step[1]^2/sigma[1]^2,scale=sigma[1]^2/mu.step[1]),add=TRUE,col=2, lwd=2)
curve(delta[2]*dgamma(x,shape=mu.step[2]^2/sigma[2]^2,scale=sigma[2]^2/mu.step[2]),add=TRUE,col=3, lwd=2)
#marginal distribution
curve(delta[1]*dgamma(x,shape=mu.step[1]^2/sigma[1]^2,scale=sigma[1]^2/mu.step[1])+delta[2]*dgamma(x,shape=mu.step[2]^2/sigma[2]^2,scale=sigma[2]^2/mu.step[2]),add=TRUE,col="black", lwd=2)

#### with covariates ####
### incorporate covariate slick
mllk.slick <- function(theta.star,x){
  mu.step <- exp(theta.star[1:2])
  sigma <- exp(theta.star[3:4])
  delta <- c(plogis(theta.star[5]),1-plogis(theta.star[5]))
  beta <- cbind(theta.star[6:7],theta.star[8:9])
  allprobs <- matrix(1,dim(x)[1],2)
  ind <- which(!is.na(x$step)) # indices of non-missing step lengths
  allprobs[ind,] <- cbind(dgamma(x$step[ind],shape=mu.step[1]^2/sigma[1]^2,scale=sigma[1]^2/mu.step[1]),
                          dgamma(x$step[ind],shape=mu.step[2]^2/sigma[2]^2,scale=sigma[2]^2/mu.step[2]))
  #for (j in 1:2){
  #  allprobs[,j] <- dgamma(x$step,shape=mu.step[j]^2/sigma[j]^2,scale=sigma[j]^2/mu.step[j])
  #}
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:dim(x)[1]){
    eta <- beta[,1]+beta[,2]*x[t,"slick"]
    Gamma <- diag(2)
    # tpm has to be updated after every new estimation of the beta values
    Gamma[!Gamma] <- exp(eta) # assign values to the off-diagonal entries of the Gamma-matrix
    Gamma <- Gamma/rowSums(Gamma)
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
return(-l)
}

N=2

#starting values
delta0 <- 0.5
# choose beta0 such that the diagonal elements of the resulting Gamma matrix will be highest (similar to starting values for gamma_ij)
beta0 <- rep(-2,N*(N-1))
# here my initial guess is that our covariate has no effect on the state process and hence I set all beta1 values to zero
beta1 <- rep(0,N*(N-1))

theta.star <- c(log(mu.step),log(sigma),qlogis(delta0), beta0, beta1)
mod_slick <- nlm(mllk.slick,theta.star,x=data,print.level=2)

# natural parameters
mu.step_slick <- exp(mod_slick$estimate[1:2])
sigma_slick <- exp(mod_slick$estimate[3:4])
mu.step_slick
sigma_slick

### Modellüberprüfung
colours <- c("#E69F00", "#56B4E9", "#009E73", "black")
curve(delta3[1]*dgamma(x,shape=mu.step3[1]^2/sigma.step3[1]^2,scale=sigma.step3[1]^2/mu.step3[1])
      +delta3[2]*dgamma(x,shape=mu.step3[2]^2/sigma.step3[2]^2,scale=sigma.step3[2]^2/mu.step3[2])
      +delta3[3]*dgamma(x,shape=mu.step3[3]^2/sigma.step3[3]^2,scale=sigma.step3[3]^2/mu.step3[3]),
      add=TRUE,col="black", lwd=2)
legend("topright", inset = c(0,0), legend= c("State 1","State 2","State 3","Sum of densities"), pch = 16, col=colours)
      
      
