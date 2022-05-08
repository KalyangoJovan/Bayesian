
library(rjags)
library(tidyverse)
library(mvtnorm)
library(invgamma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

#setwd("C:/Users/wambu/Desktop/Exam Bayes2021")
houses<-read.table("house2.txt",header=TRUE)

starting_values <- list(
  chain1 = list(beta = c(0, 4, 1, 6,4), sigma2 = 1),
  chain2 = list(beta = c(1, -1.5, .5, 2.5, 1.87), sigma2 = 2.5) 
)

# 2
#start outer function
gibbs_spl <- function(formula, n_chain = 2, n_iter = 10000, n_burnin = 1000, 
                      starting_values, data){
  
  # start inner function
  gibbs_spl_chain <- function(formula, n_iter = 10000, n_burnin = 1000, 
                              starting_values, data){
    # extract outcome 
    y <- model.frame(formula, data)[, 1]
    # extract model matrix 
    terms <- terms.formula(formula)
    X <- model.matrix(terms, data)
    # pre-calculations
    XX_t <- t(X)%*%X
    XX_t_i <- solve(XX_t)# predictors/ regressors 
    beta_ml <- solve(XX_t, t(X) %*%y)## betas 
    n <- nrow(data)
    # memory for output
    parameters <-  matrix(NA, nrow = (n_iter + n_burnin), ncol =  ncol(X)+1)
    #  starting values
    beta <- starting_values$beta
    sigma2  <- starting_values$sigma2
    # sampling from conditional posterior
    for (i in 1:(n_iter+n_burnin)){
      beta <- MASS::mvrnorm(n = 1, beta_ml, sigma2 * XX_t_i)
      sigma2 <- 1/ rgamma(1, shape = n/2, rate = t(y-X%*%beta)%*%(y-X%*%beta) *.5)
      ## save output
      parameters[i, 1:ncol(X)] <- beta
      parameters[i, ncol(X)+1] <- sqrt(sigma2) # transform variance into sd to resemble freq. mlr
    }
    # recode burnin samples to NA & remove
    parameters[1:n_burnin, ] <- NA
    parameters <- na.omit(as.data.frame(parameters))
    colnames(parameters) <-  c(sapply(0:(ncol(X)-1), function(x)(paste("b", as.character(x), sep = ""))), "sigma")
    return(parameters)
    
  }
  #  outer function again
  output_spl <- list()
  out<- 
    for (i in 1:n_chain){      
      # apply sampling for each chain
      output_spl[[i]] <- gibbs_spl_chain(formula = formula, n_iter = n_iter, n_burnin = n_burnin, starting_values = starting_values[[i]], data = houses)
      out <- data.frame()
    }
  # the final output
  for(i in 1:length(output_spl)){
    output_spl[[i]]$chain <- i
  }
  for(i in 1:length(output_spl)){
    out <- rbind(out, output_spl[[i]])
  }
  # return final output
  return(out)
}


m_gibbs   <- gibbs_spl(formula = Price ~ Taxes+Beds+Baths+Size,
                       n_chain = 2,
                       n_iter = 10000, 
                       n_burnin = 1000, 
                       starting_values, 
                       data = houses)
View(m_gibbs)

#3
#how to run a Metroplios_ Hastings  from Gibbs 
## load the data set we have to use in this case( only Size as my predictor)
#setwd("C:/Users/wambu/Desktop/Exam Bayes2021")
house1 <- read.table("house2.txt",header=TRUE)# am interested to perform a simple linear regreesion equivalent
house1
# the model to compare with 
model1 <- lm(house1$Price~., data = house1)
summary(model1)

trueA <- -79.24
trueB <- 1.69
trueSd <- sqrt(var(house1$Price)) 
sampleSize <- 27
# create independent x-values 
x <- (-(sampleSize)/2):((sampleSize)/2)
#x <- house1$Size
y<-  trueA*x+ trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)
#y <- house1$Price
## Likelihood 
likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)    
}

#  plot the likelihood profile of the slope a
slopevalues <- function(x){return(likelihood(c(x, trueB, trueSd)))}
slopelikelihoods <- lapply(seq(-500, 500, by=.05), slopevalues )
plot (seq(-100, 100, by=.05), slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")
# Prior distribution
prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min=0, max=20, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}
#The product of prior and likelihood is the actual quantity the MCMC will be working on. 
#This function is called the posterior (or to be exact, it's called the posterior after it's 
#normalized, which the MCMC).#Again, here we work with the sum because we work with logarithms.
## posterior we use sum instead of product, 
posterior <- function(param){
  return (likelihood(param) + prior(param))
}
# Metropolis algorithm 

proposalfunction <- function(param){
  return(rnorm(3,mean = param, sd= c(0.001,0.005,0.003)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

yo<-run_metropolis_MCMC(startvalue = c(0,4,20),100000)
yo
startvalue = c(0,4,20)
chain = run_metropolis_MCMC(startvalue, 100000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
acceptance


# summary
par(mfrow = c(2,3))
hist(chain[-(1:burnIn),1],nclass=40, main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]), col="red")
abline(v = trueA, col="red" )
hist(chain[-(1:burnIn),2],nclass=40, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]),col="red")
abline(v = trueB, col="red" )
hist(chain[-(1:burnIn),3],nclass=40, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]),col="red")
abline(v = trueSd, col="red" )
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a" )
abline(h = trueA, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b")
abline(h = trueB, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd")
abline(h = trueSd, col="red" )


#4 Assess convergence of the model
# Trace plot
par(mfrow = c(2, 3))
for(i in 1:6){
  plot(m_gibbs[m_gibbs$chain == 1, i], type = "l",
       xlab = colnames(m_gibbs[-ncol(m_gibbs)])[i],
       ylab = "value"
  )
  lines(m_gibbs[m_gibbs$chain == 2, i], col = "blue", 
        ylab = "value")
  
}

# Desity plot 
par(mfrow = c(2, 3))
for(i in 1:6){
  plot(density(m_gibbs[, i]), main = "", xlab = colnames(m_gibbs[-6])[i])
  abline(v = c(quantile(m_gibbs[,i], .025), quantile(m_gibbs[,i], .95)), col = "blue")
}

# Autocorrelation function 
acfm <- function(model, n_lags = 1000){
  
  n_par <- ncol(model) - 1
  n_chain <- length(unique(model$chain))
  n_lags 
  
  autocors <- list()
  length(autocors) <- n_chain
  
  for (i in 1:n_chain){
    
    autocors[[i]] <- data.frame()
    
    
    for (j in 1:n_par){
      
      par <- model[model$chain == i,j]
      
      for (k in 1:n_lags){
        
        par_t0 <- numeric(length(par)-k)
        par_t1 <- numeric(length(par)-k)
        par_t0 <- par[1:(length(par)-k)]
        par_t1 <- par[(k+1):length(par)]
        
        #return(c(length(par_t0), length(par_t1)))
        
        autocors[[i]][j, k] <- cor(par_t0, par_t1)
        
      }
    }
    autocors[[i]] <- cbind(1, autocors[[i]])
    colnames(autocors[[i]]) <- sapply((0:n_lags+1), function(x)(paste("lag", as.character(x), sep = " ")))
    
  }
  
  return(autocors)
}

# test it on my model
acfm(m_gibbs)
# plot it 
acfm_chain1 <- t(acfm(m_gibbs)[[1]])
colnames(acfm_chain1) <- c("b0", "b1", "b2", "b3","b4", "sigma")
autocor_chain2 <- acfm(m_gibbs)[[2]]

par(mfrow = c(3, 2))
for (i in 1:6){
  plot(acfm_chain1[,i],  ylab = paste( "autocorrelation",  colnames(m_gibbs[-ncol(m_gibbs)])[i]), type = "l", xlab = "lag")
  
}


## gelman function 
gelman_rubin <- function(model){
  m <- length(unique(model$chain))
  n <- nrow(model[model$chain == 1,])
  
  # write  function to compute R 
  gelman_rubin_inner <- function(para, model, m, n){
    var_within <- 0
    sums <- 0
    
    for(i in 1:m){
      x_bar <- mean(model[model$chain == i, colnames(model) == para])
      var_within <- var_within + var(model[model$chain == i, colnames(model) == para])
      x_grand <- mean(model[, para])
      sums <- sums + (x_bar - x_grand)^2
    }
    var_between <- n/(m-1) * sums
    var_total <- (n-1)/n * var_within + (1/n)*var_between
    R <- sqrt(var_total/var_within)
    return(R)
  }
  
  # apply function to all parameters in model
  parameter_names <- colnames(model)[-ncol(model)]
  output <- numeric(length(parameter_names))
  names(output) <- colnames(model)[-ncol(model)]
  
  for (i in 1:length(parameter_names)){
    
    output[i] <- gelman_rubin_inner(parameter_names[i], model, m, n)
    
  }
  
  return(output)
  
}



gelman_rubin_stat <- gelman_rubin(m_gibbs)
gelman_rubin_stat


#The Gelman and Rubin statistics of all parameters were very close to 1 or even equal to 1 f
#or example, b2 and b3. This suggests that there were no issues with convergence,hence, t
#he parameters are independent. .

mc_error <- function(model){
  n_par <- ncol(model)-1
  mcerror <- numeric(n_par)
  pars <- colnames(model[-ncol(model)])
  for (i in 1:ncol(model[, -ncol(model)])){
    mcerror[i] <- sd(model[, pars[i]])/ sqrt(nrow(model*2))
    names(mcerror) <- colnames(model[-ncol(model)])
  }
  return(mcerror)
}

#5. Check a model assumption of linear regression by means of a posterior predictive p-value

# step 2 
theta <- gibbs_spl(formula = Price ~Taxes + Beds + Baths+ Size,
                   n_chain = 2,
                   n_iter = 10000, 
                   n_burnin = 1000, 
                   starting_values, 
                   data = houses)[-7]


# step 3 
# extract outcome 
y <- model.frame(formula = Price~Taxes+Beds+Baths+Size, houses)[, 1]
# extract model matrix 
terms <- terms.formula(Price~Taxes+Beds+Baths+Size)
X <- model.matrix(terms, houses)
# simulate data
output <- list()
for (i in 1:10000){
  
  y_hat  <- X%*%as.matrix(t(theta[i, -6]))
  output[[i]] <- cbind(y_hat = y_hat, X)
  colnames(output[[i]])[1] <- "y_hat"
  
}
#step 4
e <- list()
t <- numeric(10000)
for (i in 1:10000){
  
  e[[i]] <- as.vector(output[[i]][, "y_hat"]) - y
  
}

D <- numeric(10000)
for (i in 1:10000){
  
  D[i] <- abs(mean(e[[i]])- median(e[[i]]))
  
}
# step 5 
mean_beta <- sapply(theta[-ncol(theta)], mean)
y_hat_obs <-  X%*%as.matrix((mean_beta))
e_obs <- y_hat_obs - y
D_obs <- abs(mean(e_obs) - median(e_obs))

p <- mean(ifelse(D > D_obs, 1, 0))
p

hist(D, main = "")
abline(v = D_obs, col = "blue")
#6. Obtain parameter estimates, credible intervals and interpretation
summary_pars <-
  m_gibbs%>% 
  select(b0:sigma) %>% 
  summarise(EAP = sapply(., mean),
            SD = sapply(.,sd)
  )
rownames(summary_pars) <- c("b0", "b1", "b2", "b3", "b4" ,"sigma")
knitr::kable(t(cbind(round(summary_pars, 3), `Gelman Rubin` = round(gelman_rubin_stat, 4), 
                `MC error` = mc_err)), caption = "Parameter Estimates")

## CIs
summary_quantiles <-
  m_gibbs%>% 
  select(b0:sigma) %>% 
  summarise(`2.5%` = sapply(., quantile, 0.025),
            `25%` = sapply(., quantile, 0.25),
            `50%` = sapply(., quantile, .5),
            `75%` = sapply(., quantile, .75),
            `97.5%` = sapply(., quantile, .975))
rownames(summary_quantiles) <- c("b0", "b1", "b2", "b3", "b4", "sigma")
knitr::kable(round(summary_quantiles, 3), caption = "Posterior Quantiles")

# 7 Compare multiple models by means of DIC and Bayes Factor.

## DIC
dic <- function(formula, data, n_iter, starting_values){
  
  m <- gibbs_spl(formula = formula, n_chain = 2,   n_iter = n_iter, n_burnin = 1000, starting_values = starting_values, data = data)
  samp <- m[, -ncol(m)]
  
  # construct likelihood given mean theta
  # extract model matrix 
  terms <- terms.formula(formula)
  X <- model.matrix(terms, data)
  beta <- sapply(samp[, 1:(ncol(samp)-1)], mean)
  # extract outcome
  y <- model.frame(formula, data)[,1]
  # compute predicted values
  y_hat <- as.matrix(X %*% beta)
  # d_hat
  d_hat <-  -2*sum(dnorm(x = y, mean = y_hat, sd = mean(samp$sigma), log = TRUE))
  # dbar
  l <- numeric(nrow(samp))
  
  for (i in 1:nrow(samp)){
    
    y_hat <-   X%*%t(as.matrix(samp[i,-ncol(samp)]))
    l[i] <- -2*sum(dnorm(y, mean = y_hat,
                         sd = (samp$sigma[i]), log = T))
    
  }
  
  # compute dbar, pd * dic
  d_bar <- mean(l)
  pd <- d_bar - d_hat
  dic <- d_hat + 2*pd
  # return output
  return(round(data.frame(mean_deviance = d_hat, 
                          complexity = pd,
                          dic = dic), 3))
}

# Model comparison 
dic(formula = Price~Taxes+Beds+Baths+Size, 
    data = houses, 
    n_iter = 10000,
    starting_values)

dic(formula =Price~Taxes+Beds+Baths, 
    data = houses, 
    n_iter = 10000,
    starting_values)  # starting values here?!

dic(formula =Price~Taxes+Beds, 
    data = houses, 
    n_iter = 10000,
    starting_values)

dic(formula =Price~Taxes+Beds+Size, 
    data = houses, 
    n_iter = 10000,
    starting_values)

dic(formula =Price~Taxes+Beds*Size, 
    data = houses, 
    n_iter = 10000,
    starting_values)

#The dic function, i have built suggests that a model without baths has the lowest value($DIC =289.384 $), this implies that a 
#combination of three parameters produce the best model.
#Bayes Factor
BF_function <- function(model, data){ 
  #linear model and parameters 
  l_model <- lm(model,data)
  Sigma <- vcov(l_model)[-1,-1]
  Maxlik<-coef(l_model)[-1]
  j <- length(Maxlik)
  # bayes factor 
  Posterior<-pmvnorm(lower=rep(0,j), upper=rep(Inf,j), mean=Maxlik, sigma=Sigma)
  Prior<-pmvnorm( lower=rep(0,j),upper=rep(Inf,j),mean=rep(0,4),sigma=Sigma/(4/nrow(data)))
  Bayesfactor<-Posterior[1]/Prior[1]
  return(Bayesfactor)
}
BF_function(model = Price~Taxes + Beds +Baths + Size , data = houses)
# comparison in bayes

BF_function1 <- function(model, data){ 
  #linear model and parameters 
  l_model <- lm(model,data)
  Sigma <- vcov(l_model)[-2,-2]
  Maxlik<-coef(l_model)[-2]
  j <- length(Maxlik)
  # bayes factor 
  Posterior<-pmvnorm(lower=rep(0,j), upper=rep(Inf,j), mean=Maxlik, sigma=Sigma)
  Prior<-pmvnorm( lower=rep(0,j),upper=rep(Inf,j),mean=rep(0,3),sigma=Sigma/(3/nrow(data)))
  Bayesfactor<-Posterior[1]/Prior[1]
  return(Bayesfactor)
}
BF_function1(model = Price~ Beds +Baths + Size , data = houses)
# comparison in bayes


BF_function2<- function(model, data){ 
  #linear model and parameters 
  l_model <- lm(model,data)
  Sigma <- vcov(l_model)[-3,-3]
  Maxlik<-coef(l_model)[-3]
  j <- length(Maxlik)
  # bayes factor 
  Posterior<-pmvnorm(lower=rep(0,j), upper=rep(Inf,j), mean=Maxlik, sigma=Sigma)
  Prior<-pmvnorm( lower=rep(0,j),upper=rep(Inf,j),mean=rep(0,2),sigma=Sigma/(2/nrow(data)))
  Bayesfactor<-Posterior[1]/Prior[1]
  return(Bayesfactor)
}
BF_function2(model = Price~Baths + Size , data = houses)


BF_function3<- function(model, data){ 
  #linear model and parameters 
  l_model <- lm(model,data)
  Sigma <- vcov(l_model)[-4,-4]
  Maxlik<-coef(l_model)[-4]
  j <- length(Maxlik)
  # bayes factor 
  Posterior<-pmvnorm(lower=rep(0,j), upper=rep(Inf,j), mean=Maxlik, sigma=Sigma)
  Prior<-pmvnorm( lower=rep(0,j),upper=rep(Inf,j),mean=rep(0,1),sigma=Sigma/(1/nrow(data)))
  Bayesfactor<-Posterior[1]/Prior[1]
  return(Bayesfactor)
}
BF_function3(model = Price~Size , data = houses)





