###############################################################
################ For threshold #####################
###############################################################

#' Title
#'
#' @param node_index 
#' @param W2 
#' @param mu_latent 
#' @param kappa_w2 
#' @param kappa_mu 
#' @param nbd_info 
#' @param no_of_nbd 
#'
#' @return
#' @export
#'
#' @examples
W2_sim_Gibbs_componentwise_thr<- function(node_index, W2, mean_mu_latent, kappa_w2, kappa_mu, nbd_info, no_of_nbd){
  mu_W2_index<- sum(W2[nbd_info[[node_index]]])/no_of_nbd[node_index]
  var_w2<- 1/(no_of_nbd[node_index] * kappa_w2)
  var_mu<- 1/kappa_mu
  
  var_sim<- (var_mu * var_w2 /(var_mu + var_w2))
  mean_sim<- var_sim * (mu_W2_index / var_w2 + mean_mu_latent/var_mu)
  sim<- rnorm(1, mean= mean_sim, sd= sqrt(var_sim))
  return(sim)
}





###############################################################
################ For Indictaor models #####################
###############################################################

#' Title
#'
#' @param node_index 
#' @param W2 
#' @param mu_latent 
#' @param kappa_w2 
#' @param kappa_mu 
#' @param nbd_info 
#' @param no_of_nbd 
#'
#' @return
#' @export
#'
#' @examples
W2_sim_Gibbs_componentwise_thr_ind<- function(node_index, W2, mean_mu_latent, kappa_w2, kappa_mu, nbd_info, no_of_nbd){
  mu_W2_index<- sum(W2[nbd_info[[node_index]]])/no_of_nbd[node_index]
  var_w2<- 1/(no_of_nbd[node_index] * kappa_w2)
  var_mu<- 1/kappa_mu
  
  var_sim<- (var_mu * var_w2 /(var_mu + var_w2))
  mean_sim<- var_sim * (mu_W2_index / var_w2 + mean_mu_latent/var_mu)
  sim<- rnorm(1, mean= mean_sim, sd= sqrt(var_sim))
  return(sim)
}





###############################################################
################ For Joint Models ###########################
###############################################################

#' Title
#'
#' @param node_index 
#' @param W1 
#' @param mean_eta_latent 
#' @param mean_mu_latent 
#' @param kappa_w1 
#' @param kappa_eta 
#' @param kappa_mu 
#' @param beta 
#' @param nbd_info 
#' @param no_of_nbd 
#'
#' @return
#' @export
#'
#' @examples
W1_sim_Gibbs_componentwise_JM<- function(node_index, W1, mean_eta_latent, mean_mu_latent, kappa_w1, kappa_eta, kappa_mu, beta, nbd_info, no_of_nbd){
  #browser()
  mu_W1_index<- sum(W1[nbd_info[[node_index]]])/no_of_nbd[node_index]
  var_w1<- 1/(no_of_nbd[node_index] * kappa_w1)
  var_eta<- 1/kappa_eta
  var_mu<- 1/kappa_mu
  
  var_sim<- (1/ (beta^2 / var_mu +  1 / var_eta + 1 / var_w1)) 
  mean_sim<- var_sim * (mu_W1_index / var_w1 +  mean_eta_latent / var_eta + beta * mean_mu_latent/var_mu )
  sim<- rnorm(1, mean= mean_sim, sd= sqrt(var_sim))
  return(sim)
}



#' Title
#'
#' @param node_index 
#' @param W2 
#' @param mu_latent 
#' @param kappa_w2 
#' @param kappa_mu 
#' @param nbd_info 
#' @param no_of_nbd 
#'
#' @return
#' @export
#'
#' @examples
W2_sim_Gibbs_componentwise_JM<- function(node_index, W2, mean_mu_latent, kappa_w2, kappa_mu, nbd_info, no_of_nbd){
  mu_W2_index<- sum(W2[nbd_info[[node_index]]])/no_of_nbd[node_index]
  var_w2<- 1/(no_of_nbd[node_index] * kappa_w2)
  var_mu<- 1/kappa_mu
  
  var_sim<- (var_mu * var_w2 /(var_mu + var_w2))
  mean_sim<- var_sim * (mu_W2_index / var_w2 + mean_mu_latent/var_mu)
  sim<- rnorm(1, mean= mean_sim, sd= sqrt(var_sim))
  return(sim)
}


