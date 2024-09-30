#' Update the latent parameters in parallel (matrix forms) using the Random Walk Metropolis (RWM) in parallel
#'
#' This function performs updates for latent parameters in a parallel manner using the
#' Random Walk Metropolis (RWM) algorithm. It updates the latent parameters based on
#' random walk proposals in a Gibbs sampling step to achieve convergence.
#'
#' @param par_cur A matrix of dimension ns x nt containing the current parameter values.
#' @param par_name Name of the latent parameter vectors being updated here.
#' @param loglik_data Log-likelihood at the data level, which depends on the underlying
#'   latent parameters.
#' @param loglik_latent Log-likelihood at the latent level (may be considered as a prior).
#' @param ns Number of spatial locations.
#' @param nt Number of time points.
#' @param var_prop Tuning parameter (variance) of the normal random walk proposals.
#' @param ... Additional arguments to be passed to the log-likelihood functions.
#'
#' @return A list containing the updated parameter matrix and the number of times a
#'   parameter is accepted.
#'
#' @examples
#' # Example usage of the latent_MH_update_parallel function:
#' par_cur <- matrix(rnorm(5 * 4), nrow = 5, ncol = 4)
#' par_name <- "latent_parameter"
#' ns <- 5
#' nt <- 4
#' var_prop <- 0.1
#' acc_par <- matrix(runif(ns * nt), nrow = ns, ncol = nt)
#' loglik_data <- function(latent_parameter) -0.5 * (latent_parameter - 1)^2
#' loglik_latent <- function(latent_parameter) dnorm(latent_parameter, mean = 1, sd = 0.1, log = TRUE)
#' updated_params <- latent_MH_update_parallel(par_cur, par_name, loglik_data, loglik_latent, ns, nt, var_prop, acc_par)
#'
#' @seealso
#' Other functions used in the Gibbs sampling process.

latent_MH_update_parallel <-
  function(par_cur,
           par_name,
           loglik_data,
           loglik_latent,
           ns,
           nt,
           var_prop,
           ...) {
    #browser()
    ellipsis <- list(...)
    args_loglik <- formalArgs(loglik_data)
    args_prior<- formalArgs(loglik_latent)
    stopifnot(
      is.finite(par_cur),
      isTRUE(all(args_loglik == formalArgs(loglik_data))),
      isTRUE(all(args_prior == formalArgs(loglik_latent))),
      dim(var_prop) == c(ns, nt)
    )
    
    ## current parameter value
    loglik_args <- ellipsis[names(ellipsis) %in% args_loglik]
    prior_args<- ellipsis[names(ellipsis) %in% args_prior]
    
    loglik_args[[par_name]] <- par_cur
    prior_args[[par_name]] <- par_cur
    
    loglik_curr <-
      do.call(what = loglik_data, args = loglik_args) +
      do.call(what = loglik_latent, args = prior_args)
    ## proposing the parameter
    par_prop <- matrix(rnorm(n = ns * nt,
                             mean = par_cur,
                             sd = sqrt(var_prop)), nrow = ns, ncol=nt)
    
    loglik_args[[par_name]] <- par_prop
    prior_args[[par_name]] <- par_prop
    
    loglik_prop <-
      do.call(what = loglik_data, args = loglik_args) +
      do.call(what = loglik_latent, args = prior_args)
    
    log_MH_ratio <-
      loglik_prop - loglik_curr
    
    par_new<- ifelse((matrix(log(runif(ns * nt)), nrow = ns, ncol = nt) < as.matrix(log_MH_ratio)) & (!is.na(as.matrix(log_MH_ratio))),
                     par_prop, par_cur)
    
    
    # par_new<- ifelse((matrix(log(runif(ns * nt)), nrow = ns, ncol = nt) < as.matrix(log_MH_ratio)) ,
    #                    par_prop, par_cur)
    
    return(par_new)
  }

#' Update the latent parameters in parallel (matrix forms) using the Random Walk Metropolis (RWM) in parallel
#'
#' This function performs updates for latent parameters in a parallel manner using the
#' Random Walk Metropolis (RWM) algorithm. It updates the latent parameters based on
#' random walk proposals in a Gibbs sampling step to achieve convergence.
#'
#' @param par_cur A matrix of dimension ns x nt containing the current parameter values.
#' @param par_name Name of the latent parameter vectors being updated here.
#' @param loglik_data Log-likelihood at the data level, which depends on the underlying
#'   latent parameters.
#' @param loglik_latent Log-likelihood at the latent level (may be considered as a prior).
#' @param ns Number of spatial locations.
#' @param nt Number of time points.
#' @param var_prop Tuning parameter (variance) of the normal random walk proposals.
#' @param ... Additional arguments to be passed to the log-likelihood functions.
#'
#' @return A list containing the updated parameter matrix and the number of times a
#'   parameter is accepted.
#'
#' @examples
#' # Example usage of the latent_MH_update_parallel function:
#' par_cur <- matrix(rnorm(5 * 4), nrow = 5, ncol = 4)
#' par_name <- "latent_parameter"
#' ns <- 5
#' nt <- 4
#' var_prop <- 0.1
#' acc_par <- matrix(runif(ns * nt), nrow = ns, ncol = nt)
#' loglik_data <- function(latent_parameter) -0.5 * (latent_parameter - 1)^2
#' loglik_latent <- function(latent_parameter) dnorm(latent_parameter, mean = 1, sd = 0.1, log = TRUE)
#' updated_params <- latent_MH_update_parallel(par_cur, par_name, loglik_data, loglik_latent, ns, nt, var_prop, acc_par)
#'
#' @seealso
#' Other functions used in the Gibbs sampling process.

latent_MH_update_parallel_constrain <-
  function(par_cur,
           par_name,
           loglik_data,
           loglik_latent,
           ns,
           nt,
           var_prop,
           transform,
           lb, 
           ub,
           ...) {
    #browser()
    ellipsis <- list(...)
    args_loglik <- formalArgs(loglik_data)
    args_prior<- formalArgs(loglik_latent)
    stopifnot(
      is.finite(par_cur),
      isTRUE(all(args_loglik == formalArgs(loglik_data))),
      isTRUE(all(args_prior == formalArgs(loglik_latent))),
      dim(var_prop) == c(ns, nt)
    )
    
    ## current parameter value
    loglik_args <- ellipsis[names(ellipsis) %in% args_loglik]
    prior_args<- ellipsis[names(ellipsis) %in% args_prior]
    
    loglik_args[[par_name]] <- par_cur
    prior_args[[par_name]] <- par_cur
    
    loglik_curr <-
      do.call(what = loglik_data, args = loglik_args) +
      do.call(what = loglik_latent, args = prior_args)
    
    # transforming the parameters to constrain scales  
    trpar_curr <- ifelse(transform, transfo(par = par_cur, lb = lb, ub = ub), par_cur)
    # Sample a proposal from the normal: proposing the parameter 
    trpar_prop<- matrix(rnorm(n = ns * nt,
                              mean = trpar_curr,
                              sd = sqrt(var_prop)), nrow = ns, ncol=nt)
    # Compute the difference of log-Jacobin
    # adj <- ifelse(transform, jac_inv_transfo(tpar = trpar_prop, lb = lb, ub = ub, log=TRUE) - 
    #                 jac_inv_transfo(tpar = trpar_curr, lb = lb, ub = ub, log=TRUE), 0)
    adj <- jac_inv_transfo(tpar = trpar_prop, lb = lb, ub = ub, log=TRUE) - 
      jac_inv_transfo(tpar = trpar_curr, lb = lb, ub = ub, log=TRUE)
    ### transforming back to original scale
    par_prop <- ifelse(transform, inv_transfo(tpar = trpar_prop, lb = lb, ub = ub), trpar_prop)
    
    loglik_args[[par_name]] <- par_prop
    prior_args[[par_name]] <- par_prop
    
    loglik_prop <-
      do.call(what = loglik_data, args = loglik_args) +
      do.call(what = loglik_latent, args = prior_args)
    
    log_MH_ratio <-
      loglik_prop - loglik_curr + adj
    
    par_new<- ifelse((matrix(log(runif(ns * nt)), nrow = ns, ncol = nt) < as.matrix(log_MH_ratio)) & (!is.na(as.matrix(log_MH_ratio))),
                     par_prop, par_cur)
    # par_new<- ifelse((matrix(log(runif(ns * nt)), nrow = ns, ncol = nt) < as.matrix(log_MH_ratio)) ,
    #                    par_prop, par_cur)
    
    return(par_new)
  }



#' Univariate (block) updates based on random walk proposals
#'
#' @param par_curr current parameter values
#' @param par_name name of the parameter of interest. While writing the log-likelihood and log-prior the parameter of interest should have same parameter values 
#' @param loglik log-likelihood function 
#' @param logprior log-prior function 
#' @param var_markdist variance of the random walk proposals
#' @param lb lower support of the parameter of interest 
#' @param ub upper support of the parameter of interest 
#' @param transform logical, if TRUE transformation are applied are to convert the current parameter space to whole real line and then performed Gaussian proposals. If FALSE then a truncated normal proposals are used in the (lb, ub)
#' @param ... additional parameter 
#'
#' @return return an updated values of parameter of interest
block_rw_update <-
  function(par_curr,
           par_name,
           loglik,
           logprior,
           var_markdist,
           lb = -Inf,
           ub = Inf,
           transform = FALSE,
           ...){
    #browser()
    ellipsis <- list(...)
    args_loglik <- formalArgs(loglik)
    args_logprior<- formalArgs(logprior)
    
    stopifnot(is.finite(par_curr))
    # if(lb == -Inf & ub == Inf){
    #   transform <- FALSE
    # }
    # Parameter value to be returned
    par_new <- par_curr # if all things fail
    # Copy arguments into a list to call function
    loglik_args <- ellipsis[names(ellipsis) %in% args_loglik]
    logprior_args<- ellipsis[names(ellipsis) %in% args_logprior]
    # Override value of parameter
    loglik_args[[par_name]] <- par_curr
    logprior_args[[par_name]]<- par_curr
    
    logpost_curr <-
      do.call(what = loglik,
              args = loglik_args) +
      do.call(what = logprior,
              args = logprior_args)
    if(!transform){
      par_prop <- rtnorm(n = length(par_curr),
                         a = lb,
                         b = ub,
                         mean = par_curr,
                         sd = sqrt(var_markdist))
      adj <- sum(dtnorm(par_curr,
                        a = lb,
                        b = ub,
                        mean = par_prop,
                        sd = sqrt(var_markdist),
                        log = TRUE)) -
        sum(dtnorm(par_prop,
                   a = lb,
                   b = ub,
                   mean = par_curr,
                   sd = sqrt(var_markdist),
                   log = TRUE))
    } else{ # Transformation is TRUE
      # Y = f(X)
      # Transform parameter to the real line
      trpar_curr <- transfo(par = par_curr, lb = lb, ub = ub)
      # Sample a proposal from the normal
      trpar_prop <- rnorm(n = length(par_curr), mean = trpar_curr, sd = sqrt(var_markdist))
      # Compute the difference of log-Jacobin
      adj <- sum(jac_inv_transfo(tpar = trpar_prop, lb = lb, ub = ub, log=TRUE)) -
        sum(jac_inv_transfo(tpar = trpar_curr, lb = lb, ub = ub, log=TRUE))
      ### tranforming back to origial scale
      par_prop <- inv_transfo(tpar = trpar_prop, lb = lb, ub = ub)
    }
    # Reverse move
    loglik_args[[par_name]] <- par_prop
    logprior_args[[par_name]] <- par_prop
    
    logpost_prop <-
      do.call(what = loglik,
              args = loglik_args) +
      do.call(what = logprior,
              args = logprior_args)
    log_MH_ratio <-
      logpost_prop - logpost_curr + adj ### because adjusted in the Jacobean
    
    if((log_MH_ratio > log(runif(1))) & (!is.na(log_MH_ratio))){
      par_new <- par_prop
    }
    
    # if((log_MH_ratio > log(runif(1)))){
    #   par_new <- par_prop
    # } 
    return(par_new)
  }










