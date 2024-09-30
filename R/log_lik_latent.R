#####################################################################
######### Log-likelihood for latent mu in case of indicator model ######################
#####################################################################


#' Title
#'
#' @param A 
#' @param mu 
#' @param intercept2 
#' @param W1 
#' @param W2 
#' @param beta 
#' @param beta2 
#' @param kappa_mu 
#' @param log.hyper.mu 
#' @param Z2 
#' @param family 
#'
#' @return
#' @export
#'
#' @examples
log_lik_latent_mu_thr_ind<- function(mu, intercept2, W2, beta2,  kappa_mu, Z2){
  #browser()
  #log.post.process_mu<- -0.5*kappa_mu*sum(mu^2) + kappa_mu * sum(mu *(intercept2+Z2%*%beta2+beta*(A2%*%W1)+A2%*%W2)) 
  loglik<- -0.5*kappa_mu* mu^2 + kappa_mu * mu *(intercept2 + Z2 %*% beta2 + W2)
  return(loglik)
  
}

#####################################################################
######### Log-likelihood for latent mu and etas for joint model ######################
#####################################################################


#' Title
#'
#' @param A 
#' @param mu 
#' @param intercept2 
#' @param W1 
#' @param W2 
#' @param beta 
#' @param beta2 
#' @param kappa_mu 
#' @param log.hyper.mu 
#' @param Z2 
#' @param family 
#'
#' @return
#' @export
#'
#' @examples
log_lik_latent_mu_JM<- function(mu, intercept2, W1, W2, beta, beta2,  kappa_mu, Z2){
  #browser()
  #log.post.process_mu<- -0.5*kappa_mu*sum(mu^2) + kappa_mu * sum(mu *(intercept2+Z2%*%beta2+beta*(A2%*%W1)+A2%*%W2)) 
  loglik<- -0.5*kappa_mu* mu^2 + kappa_mu * mu *(intercept2 + Z2 %*% beta2 + beta* W1 + W2)
  return(loglik)
  
}



#' Title
#'
#' @param Y 
#' @param eta 
#' @param intercept1 
#' @param beta1 
#' @param kappa_eta 
#' @param W1 
#' @param Z1 
#' @param A1 
#'
#' @return
#' @export
#'
#' @examples
log_lik_latent_eta_JM<- function(eta, intercept1, beta1, kappa_eta, W1, Z1){
  loglik<- -0.5 * kappa_eta * eta^2 + kappa_eta * eta  * (intercept1 + Z1%*%beta1 + W1)
  return(loglik)
}

#' Title
#'
#' @param Y 
#' @param eta 
#' @param intercept1 
#' @param beta1 
#' @param kappa_eta 
#' @param W1 
#' @param Z1 
#' @param A1 
#'
#' @return
#' @export
#'
#' @examples
log_lik_eta_JM<- function(Y, eta){
  loglik<-  Y*eta - exp(eta)
  return(loglik)
}