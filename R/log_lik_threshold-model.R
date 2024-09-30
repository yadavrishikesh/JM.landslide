#' log-likelihood of the mark-distribution process a the data level; i.e., the likelihood of density A | theta_A, mu
#'
#' @param cur_par  hyper-parameters vector of the continuous mark distributions 
#' @param mu log-median of the latent size process. More explicitly exp(mu) is the median of the size distribution A
#' @param A the size data 
#' @param sum_dens whether to calculate the log-density individually or the sum of log-likelihood 
#' @param ind_zeros_counts 
#' @param thr.family 
#'
#' @return 
log_lik_thr<- function(
    cur_par,
    mu, 
    A, 
    thr.family,
    ind_zeros_counts,
    sum_dens=TRUE) {
  #browser()
  if(thr.family=="gamma"){
    k<- cur_par[1]
    if (sum_dens==TRUE){
      log.post_A<- sum(ifelse(ind_zeros_counts, 0, dgamma(x = A, shape = k, scale = exp(mu)/k , log = TRUE)))
    } else {
      log.post_A<- ifelse(ind_zeros_counts, 0, dgamma(x = A, shape = k, scale = exp(mu)/k , log = TRUE))
    }
  } 
  if(thr.family=="logNormal"){
    k<- cur_par[1]
    if (sum_dens==TRUE){
      log.post_A<- sum(ifelse(ind_zeros_counts, 0, dnorm(log(A), mean=mu, sd=sqrt(1/k), log = TRUE)))
    } else {
      log.post_A<- ifelse(ind_zeros_counts, 0, dnorm(log(A), mean=mu, sd=sqrt(1/k), log = TRUE))
    }
  }
  return(log.post_A)
}

