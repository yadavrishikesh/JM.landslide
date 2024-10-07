## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(JM.landslide)

## -----------------------------------------------------------------------------
# Load necessary libraries
library(JM.landslide)
library(spdep)
library(INLA)

# Load Wenchuan landslides data
data("Wenchuan_info_used_for_simulation")

## -----------------------------------------------------------------------------
# Create adjacency graph file
nb2INLA("adjgraph-sim.txt", poly2nb(shp_selected, queen = FALSE, row.names = shp_selected$SU_ID))
adjacensy <- inla.read.graph(filename = "adjgraph-sim.txt")

## -----------------------------------------------------------------------------
# Create the precision matrix
N <- adjacensy$n
diag.Q <- diag(adjacensy$nnbs, N)
A.Q <- matrix(0, nrow = N, ncol = N)
for (i in 1:N) {
  A.Q[i, adjacensy$nbs[[i]]] <- 1
}
Q <- diag.Q - A.Q  # Precision matrix

## -----------------------------------------------------------------------------
# Set fixed parameters
kappa_w1 <- 10
kappa_w2 <- 5
kappa_eta <- 10
kappa_mu <- 5
intercept1 <- 2
intercept2 <- 4
beta <- 1
other.hyper <- list(
  kappa_w1 = kappa_w1,
  kappa_w2 = kappa_w2,
  kappa_eta = kappa_eta,
  kappa_mu = kappa_mu,
  intercept1 = intercept1,
  intercept2 = intercept2,
  beta = beta
)

## -----------------------------------------------------------------------------
# Simulate covariates
Z1 <- mvtnorm::rmvnorm(N, mean = rep(0, 3))
Z2 <- mvtnorm::rmvnorm(N, mean = rep(0, 2))
hyper.mu<- c(20,0.1)

set.seed(1)
beta1 <- runif(ncol(Z1), -1, 1)
beta2 <- runif(ncol(Z2), -1, 1)

## -----------------------------------------------------------------------------
model_type <- "jSp"  # Fixed Effects

# Simulate data using the function
set.seed(1)
Sim_mark_data <- sim_mark_function(
  Q = Q,
  hyper.mu =  hyper.mu,
  other.hyper = other.hyper,
  beta1 = beta1, beta2 = beta2,
  Z1 = as.matrix(Z1), Z2 = as.matrix(Z2),
  mark_dist = "eGPD",
  model_type = model_type
)

## ----fig.width = 8, fig.height = 12-------------------------------------------
# Access results
Y <- Sim_mark_data$Y
A <- Sim_mark_data$A
mu <- Sim_mark_data$mu
eta <- Sim_mark_data$eta
W1 <- Sim_mark_data$W1
W2 <- Sim_mark_data$W2
par(mfrow=c(3,2))
hist(Y)
hist(A)
hist(eta)
hist(mu)
hist(W1)
hist(W2)

## -----------------------------------------------------------------------------
# Define the model type
model_type <- "jSp"
outputs <- mcmc_sampler(
  Y = Sim_mark_data$Y,               # Simulated landslide counts
  A = Sim_mark_data$A,               # Simulated landslide sizes
  Z1 = Z1,             # Covariates for counts
  Z2 = Z2,             # Covariates for sizes
  CV = "WS",          # Type of cross-validation (within sample)
  mark_dist = "eGPD", # Distribution of the mark process
  thr.family = "gamma", # Family of the threshold model
  model_type = model_type, # Model type used in the simulation
  adjacensy = adjacensy, # Adjacency structure for the spatial model
  no.rm.obs = 200,
  N.MCMC = 2000,
  samples.store = 250,
  print.result = FALSE
)

## -----------------------------------------------------------------------------
outputs$JM.info$summry_hyper

## -----------------------------------------------------------------------------
outputs$JM.info$summ_fixed_effects

