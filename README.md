
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JM.landslide: Joint Modeling of landsldie Counst and Sizes

<!-- badges: start -->
<!-- badges: end -->

The goal of `JM.landslide` package is to demonstrate how to jointly
model landslide counts and sizes using a sub-asymptotic distribution for
the landslide sizes, following the joint modeling approach of [Yadav et
al. (2024)](https://arxiv.org/abs/2404.09156) and [Yadav et
al. (2023)](https://doi.org/10.1093/jrsssc/qlad077). The model considers
various covariates and spatial dependencies, allowing for an effective
analysis of landslide counts and sizes. Briefly, the model is defined as
follows:

$$
Y_i \mid \boldsymbol \eta \sim \text{Poisson}[e_i \exp\{\eta_i\}], \quad i=1,\ldots,n;
$$

$$
\boldsymbol\eta = \gamma_1 \boldsymbol{1}_n + \boldsymbol \beta_1 \mathbf{Z}_{1} + \boldsymbol W_{1} + \boldsymbol \varepsilon_{\eta};
$$

$$
A_j \mid \boldsymbol \mu, \boldsymbol \Theta_A \sim F_A[\cdot;\exp\{\mu_j\},\boldsymbol \Theta_{A}],\quad j:Y_j>0;
$$

$$
\boldsymbol\mu = \gamma_2 \boldsymbol{1}_n + \boldsymbol \beta_2 \mathbf{Z}_{2} + \beta \boldsymbol W_1 + \boldsymbol W_{2} + \boldsymbol \varepsilon_{\mu}.
$$

In this model:

- $\boldsymbol \eta$ represents the linear predictor for landslide
  counts.
- $A_j$ denotes the size of landslides.
- $e_i$ are offset terms, which can be set to one or the specific slope
  unit area.

Further details about the model can be found in the Yadav et al. (2024)
and Yadav et al. (2023)

## Installation

You can install the development version of JM.landslide from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("yadavrishikesh/JM.landslide")
#pak::pak("yadavrishikesh/JM.landslide")
```

### Synthetic data Simulation

#### 2. Loading Required Libraries and Data

``` r
# Load necessary libraries
library(JM.landslide)
library(spdep)
library(INLA)

# Load Wenchuan landslides data
data("Wenchuan_info_used_for_simulation")
```

- **Purpose**: This section loads the necessary R packages and the
  dataset used for simulation.
- **Libraries**:
  - `JM.landslide`: Contains functions for joint modeling of landslide
    data.
  - `spdep`: Provides functions for spatial dependence analysis.
  - `INLA`: Used for approximate Bayesian inference for latent models.
- **Data**: The `Wenchuan_info_used_for_simulation` dataset contains the
  necessary covariates for the simulation.

#### 3. Creating the Adjacency Graph

``` r
# Create adjacency graph file
nb2INLA("adjgraph-sim.txt", poly2nb(shp_selected, queen = FALSE, row.names = shp_selected$SU_ID))
adjacensy <- inla.read.graph(filename = "adjgraph-sim.txt")
```

- **Purpose**: This section generates an adjacency graph based on the
  spatial structure of the data.
- **Functionality**:
  - `nb2INLA`: Converts neighborhood structures into a format compatible
    with INLA.
  - `poly2nb`: Generates a list of neighbors for each spatial polygon.
  - `inla.read.graph`: Reads the adjacency graph file into R for further
    analysis.

#### 4. Creating the Precision Matrix

``` r
# Create the precision matrix
N <- adjacensy$n
diag.Q <- diag(adjacensy$nnbs, N)
A.Q <- matrix(0, nrow = N, ncol = N)
for (i in 1:N) {
  A.Q[i, adjacensy$nbs[[i]]] <- 1
}
Q <- diag.Q - A.Q  # Precision matrix
```

- **Purpose**: This section constructs the precision matrix `Q`, which
  is crucial for defining the spatial dependency structure of the model.
- **Explanation**:
  - `diag.Q`: Creates a diagonal matrix where the diagonal elements
    represent the number of neighbors for each node.
  - `A.Q`: Initializes an adjacency matrix where connections between
    nodes are indicated.
  - The final precision matrix `Q` is derived by subtracting the
    adjacency matrix from the diagonal matrix.

#### 5. Setting Fixed Parameters

``` r
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
```

- **Purpose**: This section initializes hyperparameters that will be
  used in the model.
- **Hyperparameters**:
  - `kappa_w1`, `kappa_w2`, `kappa_eta`, `kappa_mu`: Precision
    parameters for the spatial random effects.
  - `intercept1`, `intercept2`, `beta`: Intercept and regression
    coefficients for the model.

#### 6. Simulating Covariates

``` r
# Simulate covariates
Z1 <- mvtnorm::rmvnorm(N, mean = rep(0, 3))
Z2 <- mvtnorm::rmvnorm(N, mean = rep(0, 2))
hyper.mu<- c(20,0.1)

set.seed(1)
beta1 <- runif(ncol(Z1), -1, 1)
beta2 <- runif(ncol(Z2), -1, 1)
```

- **Purpose**: Simulates covariate matrices `Z1` and `Z2` which contain
  fixed effects for the landslide counts and sizes.
- **Explanation**:
  - `rmvnorm`: Generates multivariate normal random variables.
  - `runif`: Generates uniform random variables for the regression
    coefficients.

#### 7. Simulating Data

``` r
model_type <- "FE"  # Fixed Effects

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
```

- **Purpose**: This section simulates the landslide count and size data
  using the specified model.
- **Functionality**:
  - `sim_mark_function`: A custom function that simulates the data based
    on the defined model, precision matrix, hyperparameters, and
    covariates.

#### 8. Accessing Results

``` r
# Access results
Y <- Sim_mark_data$Y
A <- Sim_mark_data$A
mu <- Sim_mark_data$mu
eta <- Sim_mark_data$eta
W1 <- Sim_mark_data$W1
W2 <- Sim_mark_data$W2
```

- **Purpose**: This section extracts the simulated results from the
  output of the simulation function.
- **Explanation**:
  - `Y`: Simulated landslide counts.
  - `A`: Simulated landslide sizes.
  - `mu`, `eta`: Latent variables representing the predictors for counts
    and sizes.
  - `W1`, `W2`: Spatially structured latent random effects.

## Fitting the Model

After simulating the data, we can proceed to fit the model using the
MCMC sampler provided in the `JM.landslide` package. This allows us to
analyze the simulated landslide counts and sizes, estimating the model
parameters based on the generated data.

#### Code to Fit the Model

``` r
# Define the model type
model_type <- "jSp"

# Load the simulated data
load(paste0("SimulatedData_", model_type, ".RData"))

# Load the JM.landslide library
library(JM.landslide)

# Access the help documentation for the mcmc_sampler function
?JM.landslide::mcmc_sampler

# Fit the model using the MCMC sampler
outputs <- mcmc_sampler(
  Y = Y,               # Simulated landslide counts
  A = A,               # Simulated landslide sizes
  Z1 = Z1,             # Covariates for counts
  Z2 = Z2,             # Covariates for sizes
  CV = "WS",          # Type of cross-validation (within sample)
  mark_dist = "eGPD", # Distribution of the mark process
  thr.family = "gamma", # Family of the threshold model
  model_type = model_type, # Model type used in the simulation
  adjacensy = adjacensy # Adjacency structure for the spatial model
)
```

#### Explanation

- **Defining the Model Type**: The `model_type <- "jSp"` line specifies
  that we are working with a joint model with spatial components.

- **Loading the Simulated Data**: The
  `load(paste0("SimulatedData_", model_type, ".RData"))` command loads
  the simulated data generated earlier. This allows us to use the
  simulated counts and sizes in the model fitting process.

- **Loading Required Libraries**: The `library(JM.landslide)` statement
  loads the `JM.landslide` package, which contains the functions needed
  for fitting the model.

- **Model Fitting**: The `mcmc_sampler` function is then called to fit
  the model:

  - **Parameters**:
    - `Y` and `A`: These are the simulated landslide counts and sizes,
      respectively.
    - `Z1` and `Z2`: Covariate matrices used to model the effects on
      counts and sizes.
    - `CV`: Specifies the type of cross-validation method used, here set
      to “WS” for within-sample.
    - `mark_dist`: Indicates the distribution for the mark process, set
      to “eGPD” (extended Generalized Pareto Distribution).
    - `thr.family`: Specifies the family of the threshold model, in this
      case, “gamma”.
    - `model_type`: Indicates the type of model being used, which aligns
      with the previously defined `model_type`.
    - `adjacensy`: Provides the adjacency structure for the spatial
      dependencies in the model.

### 10. Conclusion

With the model now fitted, we can proceed to analyze the results and
assess the performance of the simulated data under the chosen model. The
MCMC sampler estimates the parameters, allowing us to evaluate the
impact of various covariates on landslide occurrences and sizes.

In this vignette, we have demonstrated the complete process from
simulating landslide data to fitting a statistical model using the
`JM.landslide` package. This structured approach allows researchers to
analyze spatial data effectively, accounting for dependencies and
covariates.

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

## References

1.  Yadav, R., Lombardo, L., & Huser, R. (2024). *Statistics of extremes
    for natural hazards: landslides and earthquakes*. arXiv preprint
    [arXiv:2404.09156](https://arxiv.org/abs/2404.09156).
2.  Yadav, R., Huser, R., Opitz, T., & Lombardo, L. (2024) Rishikesh
    Yadav, Raphael Huser, Thomas Opitz, and Luigi Lombardo (2023).
    *Joint modeling of landslide counts and sizes using spatial marked
    point processes with sub-asymptotic mark distributions.*
    `Journal of the Royal Statistical Society Series C (JRSSC): Applied Statistics, qlad077`.
    <https://doi.org/10.1093/jrsssc/qlad077>