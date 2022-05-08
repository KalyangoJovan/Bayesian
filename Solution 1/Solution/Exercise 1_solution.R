# Bayesian statistics - Exercise 1


# A: INPUT YOUR DATA
# -----------------------------------------------------------------------------
# 1. Load the data into R

# specify the data
dat <- list("y.PE" = 58, "n.PE" = 141, "y.PC" = 40, "n.PC" = 143)

# Or source the data : 
source("Exercise 1 - Data.txt")

# B: SPECIFY YOUR MODEL
# -----------------------------------------------------------------------------
# See file Model - 1 solution C



# C: OBTAIN INITIAL VALUES
# -----------------------------------------------------------------------------
# this step can be skipped in this exercise

# D: OBTAIN SAMPLES FROM THE POSTERIOR DISTRIBUTION OF SAMPLES
# -----------------------------------------------------------------------------
# 1. load the rjags library
library(rjags)   # library used to access JAGS with R

# 2. model, data and chains specification : 
model.def <- jags.model(file = "Exercise 1 - Model solution C.txt", 
												data = dat, n.chains = 2)

# 3. burn-in period : 
update(object = model.def, n.iter = 1000) # burn-in period

# 4. obtain samples from the posterior distribution of the parameters and monitor these : 
parameters <- c("theta.PE", "theta.PC", "RR")   # parameters you want to model
res <- coda.samples(model = model.def, variable.names = parameters, n.iter = 10000)

# E: INSPECTING CONVERGENCE
# -----------------------------------------------------------------------------

#library(mcmcplots)
#mcmcplot(res)

# F: SUBSTANTIVE INTERPRETATION
# -----------------------------------------------------------------------------

# obtain summary statistics
summary(res)

# J: ASSESING PRIOR INFLUENCE
# -----------------------------------------------------------------------------
# See file Model - 1 solution I for Ronald prior

# 1. model, data and chains specification : 
model.def <- jags.model(file = "Exercise 1 - Model solution I.txt", 
                        data = dat, n.chains = 2)

# 2. burn-in period : 
update(object = model.def, n.iter = 1000) # burn-in period

# 3. obtain samples from the posterior distribution of the parameters and monitor these : 
parameters <- c("theta.PE", "theta.PC", "RR")   # parameters you want to model
res <- coda.samples(model = model.def, variable.names = parameters, n.iter = 10000)

# 4.obtain summary statistics
summary(res)

# K: ASSESING PRIOR INFLUENCE, PART 2
# -----------------------------------------------------------------------------
# See file Model - 1 solution K for Jessica prior


# 1. model, data and chains specification : 
model.def <- jags.model(file = "Exercise 1 - Model solution K.txt", 
                        data = dat, n.chains = 2)

# 2. burn-in period : 
update(object = model.def, n.iter = 1000) # burn-in period

# 3. obtain samples from the posterior distribution of the parameters and monitor these : 
parameters <- c("theta.PE", "theta.PC", "RR")   # parameters you want to model
res <- coda.samples(model = model.def, variable.names = parameters, n.iter = 10000)

# 4.obtain summary statistics
summary(res)

