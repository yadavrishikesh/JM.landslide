#' #' @index_MCMC_iter: index of MCMC iterations
#' #' @sigma2_adapt: the tuning parameter in the MH and MALA algorithm
#' #' @target_accept: target acceptance probability, which is 0.224 for random wala MH and 0.57 for MALA\
#' #' @rate_adapt: acceptance rate in MCMC iterations calculated at every fixed number of iterations
#' #' @burn_in1: the first burning period in MCMC algorithm
#' #' @burn_in2: the second burning period in MCMC algorithm
#' #' @adapt_seq: the sequence, at which we update the @sigma2_adapt
#' #' @adpat: the number of MCMC iterations after which we update the @sigma2_adapt
#' #' @adapt_param: the scale parameter in adaptive algorithm
#' #' @lower.acc: The lower acceptance probability for RWM and MALA which is 0.15, and 0.50, respectively
#' #' @upper.acc: The lower acceptance probability for RWM and MALA which is 0.30, and 0.65, respectively
#' adpative_function<- function(index_MCMC_iter, sigma2_adapt, target_accept, rate_adapt, burn_in1, burn_in2, adapt, adpat_param, adapt_seq, lower.acc, upper.acc){
#'   if(index_MCMC_iter < burn_in1 + burn_in2){
#'     if(index_MCMC_iter %in% adapt_seq){
#'       if (index_MCMC_iter< burn_in1 ){
#'         sigma2_adapt<- exp(((rate_adapt/adapt)-target_accept)/adpat_param) * sigma2_adapt
#'       } else {
#'         sigma2_adapt<- ifelse((((rate_adapt/adapt)>upper.acc) | ((rate_adapt/adapt)<lower.acc)), 
#'                               exp(((rate_adapt/adapt)-target_accept)/adpat_param) * sigma2_adapt, sigma2_adapt)
#'       }
#'     }
#'   } else {
#'     sigma2_adapt<- sigma2_adapt
#'   }
#'   return(sigma2_adapt)
#' }





#' Adaptive function for updating the variance components of MH and MALA algorithm
#'
#' @param index_MCMC_iter index of MCMC iterations
#' @param sigma2_adapt the tuning parameter in the MH and MALA algorithm
#' @param target_accept target acceptance probability, which is 0.224 for random wala MH and 0.57 for MALA
#' @param rate_adapt acceptance rate in MCMC iterations calculated at every fixed number of iterations
#' @param burn_in1 the first burning period in MCMC algorithm
#' @param burn_in2 the second burning period in MCMC algorithm
#' @param adapt the sequence, at which we update the @sigma2_adapt
#' @param adpat_param the number of MCMC iterations after which we update the @sigma2_adapt
#' @param adapt_seq the scale parameter in adaptive algorithm
#' @param lower.acc The lower acceptance probability for RWM and MALA which is 0.15, and 0.50, respectively
#' @param upper.acc The lower acceptance probability for RWM and MALA which is 0.30, and 0.65, respectively

adpative_function <- function(index_MCMC_iter,
                              sigma2_adapt,
                              target_accept,
                              rate_adapt,
                              burn_in1,
                              burn_in2,
                              adapt,
                              adpat_param,
                              adapt_seq,
                              lower.acc,
                              upper.acc) {
  if (index_MCMC_iter < burn_in1 + burn_in2) {
    if (index_MCMC_iter %in% adapt_seq) {
      if (index_MCMC_iter < burn_in1) {
        sigma2_adapt <- exp(((rate_adapt / adapt) - target_accept) / adpat_param) * sigma2_adapt
      } else {
        sigma2_adapt <- ifelse((((
          rate_adapt / adapt
        ) > upper.acc) | ((
          rate_adapt / adapt
        ) < lower.acc)), exp(((rate_adapt / adapt) - target_accept
        ) / adpat_param) * sigma2_adapt, sigma2_adapt)
      }
    }
  } else {
    sigma2_adapt <- sigma2_adapt
  }
  return(sigma2_adapt)
}


