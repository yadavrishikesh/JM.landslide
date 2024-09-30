#' log-likelihood of the mark-distribution process at the data level; i.e., the likelihood of density A | theta_A, mu
#'
#' @param cur_par  hyper-parameters vector of the continuous mark distributions 
#' @param mu log-median of the latent size process. More explicitly exp(mu) is the median of the size distribution A
#' @param A the size data 
#' @param mark_dist choice of mark_dist of distributions for the size process.
#' @param sum_dens whether to calculate the log-density individually or the sum of log-likelihood 
#'
#' @return 
log_lik_shape_mixture<- function(
    cur_par,
    mu, 
    A, 
    mark_dist, 
    threshold,
    sum_dens=TRUE
) {
  ## browser()
  n1<- length(mu)
  if(mark_dist=="eGPD"){
    k<- cur_par[1]
    xi<- cur_par[2]
    sigma<- exp(mu) / evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
    
    if (sum_dens==TRUE){
      log.post_A<- rep(0, n1)
      ind.postv<- A>0
      log.post_A[ind.postv] <- dEGPD1(x = A[ind.postv], k = k, xi = xi, sigma = sigma[ind.postv], log = TRUE)
      log.post_A<- sum(log.post_A)
      
    } else {
      log.post_A<- rep(0, n1)
      ind.postv<- A>0
      log.post_A[ind.postv] <- dEGPD1(x = A[ind.postv], k = k, xi = xi, sigma = sigma[ind.postv], log = TRUE)
      
    }
  } else if (mark_dist=="bGPD"){
    k<- cur_par[1]
    if (sum_dens==TRUE){
      
      log.post_A<- rep(0, n1)
      ind<- (A < threshold) & (A!=0)
      log.post_A[ind]<- dbeta((A[ind]/threshold[ind]), shape1 = k, shape2 = (threshold[ind] / exp(mu[ind]) - 1) * k, log = TRUE) - log(threshold[ind])
      log.post_A<- sum(log.post_A)
      
    } else {
      log.post_A<- rep(0, n1)
      ind<- (A < threshold) & (A!=0)
      log.post_A[ind]<- dbeta((A[ind]/threshold[ind]), shape1 = k, shape2 = (threshold[ind] / exp(mu[ind]) - 1) * k, log = TRUE) - log(threshold[ind])
    }
  } else if(mark_dist=="tgGPD"){
    k<- cur_par[1]
    if (sum_dens==TRUE){
      
      log.post_A<- rep(0, n1)
      ind<- (A < threshold) & (A!=0)
      log.post_A[ind]<-  dgamma(A[ind], shape = k, scale = exp(mu[ind])/k, log = TRUE) - pgamma(threshold[ind], shape = k, scale =  exp(mu[ind])/k, log.p = TRUE)
      log.post_A<- sum(log.post_A)
    } else {
      
      log.post_A<- rep(0, n1)
      ind<- (A < threshold) & (A!=0)
      log.post_A[ind]<-  dgamma(A[ind], shape = k, scale = exp(mu[ind])/k, log = TRUE) - 
        pgamma(threshold[ind], shape = k, scale =  exp(mu[ind])/k, log.p = TRUE)
      
    }
  }
  return(log.post_A)
}





#' Title
#'
#' @param thr.exceed.ind 
#' @param thr.prob 
#' @param cur_par 
#' @param mu 
#' @param A 
#' @param mark_dist 
#' @param ind_zeros_counts 
#' @param threshold 
#' @param sum_dens 
#'
#' @return
#' @export
#'
#' @examples
log_lik_GPD_param<- function(
    cur_par,
    A, 
    threshold
) {
  sigma.GP<- cur_par[1]
  xi<- cur_par[2]
  thr.exceed.ind<- A > threshold
  #log.post_A<-  thr.prob[thr.exceed.ind] * evd::dgpd(x=A[thr.exceed.ind], loc = threshold[thr.exceed.ind], scale = sigma.GP, shape = xi, log = TRUE)
  log.post_A<- sum(evd::dgpd(x=A[thr.exceed.ind], 
                             loc = threshold[thr.exceed.ind],
                             scale = sigma.GP, shape = xi, log = TRUE))
  return(log.post_A)
}



