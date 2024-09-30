#' log-likelihood of the mark-distribution process a the data level; i.e., the likelihood of density A | theta_A, mu
#'
#' @param cur_par  hyper-parameters vector of the continuous mark distributions 
#' @param mu log-median of the latent size process. More explicitly exp(mu) is the median of the size distribution A
#' @param A the size data 
#' @param sum_dens whether to calculate the log-density individually or the sum of log-likelihood 
#' @param ind_zeros_counts 
#'
#' @return 
log_lik_indicator_model<- function(mu, A) {
  log.post_A<-  A * log(logistic(mu)) + (1-A) * log(1- logistic(mu))  # dbinom(A, size=1, prob = exp(mu)/(1+exp(mu)), log = TRUE)
return(log.post_A)
}



#' Title
#'
#' @param mu_mean 
#' @param kappa.mu 
#' @param ind_zero 
#' @param ind.NA 
#'
#' @return
#' @export
#'
#' @examples
update_mu<- function(mu_mean, kappa.mu, ind_zero, ind.NA, CV){
  #browser()
  sim<- rep(NA, length(A))
  sim[ind_zero]<- truncnorm::rtruncnorm(n = 1, a=-Inf, b=0, mean = mu_mean[ind_zero], sd=sqrt(1/kappa.mu))
  sim[!ind_zero]<- truncnorm::rtruncnorm(n = 1, a=0, b=Inf, mean = mu_mean[!ind_zero], sd=sqrt(1/kappa.mu))
  if(CV=="OOS"){
    sim[ind.NA]<- truncnorm::rtruncnorm(n =1,  a=-Inf, b=Inf, mean = mu_mean[ind.NA], sd=sqrt(1/kappa.mu))
  }
  return(sim)
}
