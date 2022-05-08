# Bayesian statistics - Exercise 2
library(rjags)

# C : RUNNING THE ANALYSIS 
# -----------------------------------------------------------------------------
# 1. Load the data into R
source("Exercise 2 - Data.txt")

# 2. Create a model object with the jags.model() function
# Specify the model, data and number of chains
model <- jags.model(file = "Exercise 2 - Model_specified.txt", data = dat, n.chains = 2)

# 3. Use the update function to run the Markov Chain for a burn-in period 
# of 1000 iterations	
update(object = model, n.iter = 1000)

# 4. Use the coda.samples function to get samples from the posterior 
# distribution
parameters <- c("theta.PE", "theta.PC", "RR")   # parameters you want to model
samples <- coda.samples(model = model, variable.names = parameters, n.iter = 10000)

# 5. Get an informative summary of the samples
summary(samples)

# D : DISCREPANCY MEASURE FOR THE EMPIRICAL DATA
# -----------------------------------------------------------------------------
# Discrepancy measure for the PE condition
proportion.half1.PE <- mean(dat$LD.PE[1:70])
proportion.half2.PE <- mean(dat$LD.PE[71:141])
diff.PE <- proportion.half1.PE - proportion.half2.PE

# Discrepancy measure for the PC condition
proportion.half1.PC <- mean(dat$LD.PC[1:71])
proportion.half2.PC <- mean(dat$LD.PC[72:143])
diff.PC <- proportion.half1.PC - proportion.half2.PC

# E : DISCREPANCY MEASURES FOR MODEL-PREDICTED (REPLICATED) DATA SETS
# -----------------------------------------------------------------------------
# 1. Extract the parameter estimates
theta.PE.chain1 <- samples[[1]][,"theta.PE"]
theta.PE.chain2 <- samples[[2]][,"theta.PE"]

theta.PC.chain1 <- samples[[1]][,"theta.PC"]
theta.PC.chain2 <- samples[[2]][,"theta.PC"]

# 2. Sample replicated datasets
# Load LaplacesDemon library for rbern() function
library(LaplacesDemon)

# Storage room (each row is a replicated dataset) :
replicated.PE.chain1 <- array(data = NA, dim = c(length(theta.PE.chain1), dat$n.PE))
replicated.PE.chain2 <- array(data = NA, dim = c(length(theta.PE.chain2), dat$n.PE))
replicated.PC.chain1 <- array(data = NA, dim = c(length(theta.PC.chain1), dat$n.PC))
replicated.PC.chain2 <- array(data = NA, dim = c(length(theta.PC.chain2), dat$n.PC))

# Sample replicated datasets
# For each parameter estimate...
for(t in 1:length(theta.PE.chain1)) {
	# ... sample a replicated dataset by sampling n times from the Bernoulli distribution 
	# with as probability of success on each trial the parameter estimate
	replicated.PE.chain1[t,] <- rbern(n = dat$n.PE, prob = theta.PE.chain1[t])
	replicated.PE.chain2[t,] <- rbern(n = dat$n.PE, prob = theta.PE.chain2[t])
	
	replicated.PC.chain1[t,] <- rbern(n = dat$n.PC, prob = theta.PC.chain1[t])
	replicated.PC.chain2[t,] <- rbern(n = dat$n.PC, prob = theta.PC.chain2[t])
}

# 3. Calculate the discrepancy measures for the replicated datasets
# Chain 1 
predicted.proportion.half1.PE.chain1 <- apply(replicated.PE.chain1, 1, function(x) mean(x[1:70]))
predicted.proportion.half2.PE.chain1 <- apply(replicated.PE.chain1, 1, function(x) mean(x[71:141]))
predicted.diff.PE.chain1 <- predicted.proportion.half1.PE.chain1 - predicted.proportion.half2.PE.chain1

predicted.proportion.half1.PC.chain1 <- apply(replicated.PC.chain1, 1, function(x) mean(x[1:71]))
predicted.proportion.half2.PC.chain1 <- apply(replicated.PC.chain1, 1, function(x) mean(x[72:143]))
predicted.diff.PC.chain1 <- predicted.proportion.half1.PC.chain1 - predicted.proportion.half2.PC.chain1

# Chain 2
predicted.proportion.half1.PE.chain2 <- apply(replicated.PE.chain2, 1, function(x) mean(x[1:70]))
predicted.proportion.half2.PE.chain2 <- apply(replicated.PE.chain2, 1, function(x) mean(x[71:141]))
predicted.diff.PE.chain2 <- predicted.proportion.half1.PE.chain2 - predicted.proportion.half2.PE.chain2

predicted.proportion.half1.PC.chain2 <- apply(replicated.PC.chain2, 1, function(x) mean(x[1:71]))
predicted.proportion.half2.PC.chain2 <- apply(replicated.PC.chain2, 1, function(x) mean(x[72:143]))
predicted.diff.PC.chain2 <- predicted.proportion.half1.PC.chain2 - predicted.proportion.half2.PC.chain2

# F : CODE TO OBTAIN POSTERIOR PREDICTIVE P-VALUES
# -----------------------------------------------------------------------------
# Create a 0/1 variable indicating whether the replicated data have more 
# extreme discrepancy measures than our emperical data

# Chain 1 
abs.diff.PE <- abs(diff.PE)
abs.predicted.diff.PE.chain1 <- abs(predicted.diff.PE.chain1)
ppp.PE.chain1 <- ifelse((abs.predicted.diff.PE.chain1 - abs.diff.PE) >= 0, 1, 0)

abs.diff.PC <- abs(diff.PC)
abs.predicted.diff.PC.chain1 <- abs(predicted.diff.PC.chain1)
ppp.PC.chain1 <- ifelse((abs.predicted.diff.PC.chain1 - abs.diff.PC) >= 0, 1, 0)

# Chain 2
abs.predicted.diff.PE.chain2 <- abs(predicted.diff.PE.chain2)
ppp.PE.chain2 <- ifelse((abs.predicted.diff.PE.chain2 - abs.diff.PE) >= 0, 1, 0)

abs.predicted.diff.PC.chain2 <- abs(predicted.diff.PC.chain2)
ppp.PC.chain2 <- ifelse((abs.predicted.diff.PC.chain2 - abs.diff.PC) >= 0, 1, 0)

# Posterior predictive p-values
mean(ppp.PC.chain1)
mean(ppp.PE.chain1)

mean(ppp.PC.chain2)
mean(ppp.PE.chain2)

# OPTIONAL: E-F, POSTERIOR PREDICTIVE P-VALUES FOR BOTH CHAINS WITHIN A FOR LOOP
# -----------------------------------------------------------------------------
# Storage room
ppp.PE <- array(data = NA, dim = c(10000, 2))
ppp.PC <- array(data = NA, dim = c(10000, 2))

# Repeat for both chains 
for(c in 1:2) {
	# Extract parameter estimates
	theta.PC <- samples[[c]][,"theta.PC"]  
	theta.PE <- samples[[c]][,"theta.PE"]
	
	# Sample replicated datasets
	replicated.PE <- array(data = NA, dim = c(length(theta.PE), dat$n.PE))
	replicated.PC <- array(data = NA, dim = c(length(theta.PC), dat$n.PC))
	
	for(t in 1:length(theta.PE)) {
		replicated.PE[t,] <- rbern(n = dat$n.PE, prob = theta.PE[t])
		replicated.PC[t,] <- rbern(n = dat$n.PC, prob = theta.PC[t])
	}
	
	# Discrepancy measure for the replicated datasets
	predicted.proportion.half1.PE <- apply(replicated.PE, 1, function(x) mean(x[1:70]))
	predicted.proportion.half2.PE <- apply(replicated.PE, 1, function(x) mean(x[71:141]))
	predicted.diff.PE <- predicted.proportion.half1.PE - predicted.proportion.half2.PE
	
	predicted.proportion.half1.PC <- apply(replicated.PC, 1, function(x) mean(x[1:71]))
	predicted.proportion.half2.PC <- apply(replicated.PC, 1, function(x) mean(x[72:143]))
	predicted.diff.PC <- predicted.proportion.half1.PC - predicted.proportion.half2.PC
	
	# Create a 0/1 variable indicating whether the replicated data have more 
	# extreme discrepancy measures than our emperical data
	abs.diff.PE <- abs(diff.PE)
	abs.predicted.diff.PE <- abs(predicted.diff.PE)
	ppp.PE[,c] <- ifelse((abs.predicted.diff.PE - abs.diff.PE) >= 0, 1, 0)
	
	abs.diff.PC <- abs(diff.PC)
	abs.predicted.diff.PC <- abs(predicted.diff.PC)
	ppp.PC[,c] <- ifelse((abs.predicted.diff.PC - abs.diff.PC) >= 0, 1, 0)
	
}

# PPP per chain
apply(ppp.PE, 2, mean)
apply(ppp.PC, 2, mean)








