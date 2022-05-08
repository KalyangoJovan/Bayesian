# Bayesian statistics - Exercise 4
library(rjags)

# C. Assessing convergence
# -----------------------------------------------------------------------------
# 1. Load the data into R
source("Exercise 4 - Data.txt")

# 2. Create a model object with the jags.model() function
# Specify the model, data and number of chains
model <- jags.model(file = "Exercise 4 - Model B.txt", data = dat, n.chains = 2)

# 3. Use the update function to run the Markov Chain for a burn-in period 
# of 1000 iterations	
update(object = model, n.iter = 1000)

# 4. Use the coda.samples function to get samples from the posterior 
# distribution
parameters <- c("a", "b1", "odds.ratio")   # parameters you want to model
samples <- coda.samples(model = model, variable.names = parameters, n.iter = 10000)

# History plot, autocorrelation plot
library(mcmcplots)
mcmcplot(mcmcout = samples)

# Gelman-Rubin diagnostic 
gelman.plot(samples)


# D. Interpreting results
# -----------------------------------------------------------------------------
# 1. Load the data into R
source("Exercise 4 - Data.txt")

# 2. Create a model object with the jags.model() function
# Specify the model, data and number of chains
# Model changed to model C, with centered predictor variable
model1 <- jags.model(file = "Exercise 4 - Model C.txt", data = dat, n.chains = 2)

# 3. Use the update function to run the Markov Chain for a burn-in period 
# of 1000 iterations	
update(object = model1, n.iter = 1000)

# 4. Use the coda.samples function to get samples from the posterior 
# distribution
parameters <- c("a", "b1", "odds.ratio")   # parameters you want to model
samples <- coda.samples(model = model1, variable.names = parameters, n.iter = 10000)

summary(samples)

# E. Model comparison using the DIC.
# ----------------------------------------------------------------------------

dic.model1 <- dic.samples(model1, 1000, "pD") # Deviance Information Criterion
dic.model1

model2 <- jags.model(file = "Exercise 4 - Model E 2.txt", data = dat, n.chains = 2)
update(object = model2, n.iter = 1000)
parameters <- c("a", "b1","b2", "odds.ratio_gender", "odds.ratio_age")   # parameters you want to model
samples <- coda.samples(model = model2, variable.names = parameters, n.iter = 1000)
dic.model2 <- dic.samples(model2, 1000, "pD") # Deviance Information Criterion
dic.model2
summary(samples)

model3 <- jags.model(file = "Exercise 4 - Model E 3.txt", data = dat, n.chains = 2)
update(object = model3, n.iter = 1000)
parameters <- c("a", "b1","b2", "odds.ratio_gender", "odds.ratio_problems")   # parameters you want to model
samples <- coda.samples(model = model3, variable.names = parameters, n.iter = 1000)
dic.model3 <- dic.samples(model3, 1000, "pD") # Deviance Information Criterion
dic.model3 
summary(samples)

model4 <- jags.model(file = "Exercise 4 - Model E 4.txt", data = dat, n.chains = 2)
update(object = model4, n.iter = 1000)
parameters <- c("a", "b1","b2", "odds.ratio_gender", "odds.ratio_no.diagnosis")   # parameters you want to model
samples <- coda.samples(model = model4, variable.names = parameters, n.iter = 1000)
dic.model4 <- dic.samples(model4, 1000, "pD") # Deviance Information Criterion
dic.model4 
summary(samples)


