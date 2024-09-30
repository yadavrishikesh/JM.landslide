######################################################################################
###################### updates for threshold model  #############################
######################################################################################

#' Function to simulate the kappa_w2 for threshold model
#'
#' @param W2 
#' @param node1 
#' @param node2 
#'
#' @return
#' @export
#'
#' @examples
kappa_w2_sim_thr<-function(W2, node1, node2, hyper_fixed){
  N<-length(W2)
  sim_kappa_w2<-rgamma(1, shape=hyper_fixed$kappa_w2[1]+0.5*(N-1), 
                       rate= hyper_fixed$kappa_w2[2] + 0.5*sum((W2[node1]-W2[node2])^2))
  return(sim_kappa_w2)
}

#' Function to simulate the kappa_mu for threshold model
#'
#' @param mu 
#' @param intercept2 
#' @param beta2 
#' @param beta 
#' @param W2 
#' @param Z2 
#'
#' @return
#' @export
#'
#' @examples
kappa_mu_sim_thr<-function(mu, intercept2, beta2, W2, Z2, hyper_fixed){
  n2<-length(mu)
  sim_kappa_mu<- rgamma(1, shape=hyper_fixed$kappa_mu[1] + 0.5 * n2, 
                        rate= hyper_fixed$kappa_mu[2] + 0.5 * (sum((mu-intercept2- Z2%*%beta2 -W2)^2)))
  return(sim_kappa_mu)
}


#' Function to simulate the intercept2 for threshold model
#'
#' @param mu 
#' @param beta2 
#' @param beta 
#' @param W2 
#' @param kappa_mu 
#' @param Z2 
#'
#' @return
#' @export
#'
#' @examples
intercept2_sim_thr<-function(mu, beta2, beta,W2, kappa_mu, Z2, hyper_fixed){
  n2<-length(mu)
  # prec.intercept2<-hyper_fixed[11]+n2*kappa_mu
  # mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2-beta*(A2%*%W1)-A2%*%W2)/prec.intercept2
  prec.intercept2<-hyper_fixed$intercept2+n2*kappa_mu
  mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2- W2)/prec.intercept2
  sim_intercept2<-rnorm(n=1, mean=mean.intercept2, sd=sqrt(1/prec.intercept2))
  return(sim_intercept2)
}


#' Function to simulate beta2 for threshold model
#'
#' @param mu 
#' @param intercept2 
#' @param beta 
#' @param kappa_mu 
#' @param W2 
#' @param Z2 
#' @param Z2.crossprd 
#'
#' @return
#' @export
#'
#' @examples
beta2_sim_thr<-function(mu, intercept2, beta, kappa_mu,  W2, Z2, Z2.crossprd, hyper_fixed){
  q<-ncol(Z2)
  latent.cov.inv<- hyper_fixed$beta2 * diag(1,q) + kappa_mu * Z2.crossprd 
  latent.mean.part<- kappa_mu * (t(Z2)%*%(mu-intercept2 - W2))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(q)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
}

######################################################################################
###################### updates for indicator models  #############################
######################################################################################


#' Function to simulate the kappa_w2 for threshold indicator models
#'
#' @param W2 
#' @param node1 
#' @param node2 
#'
#' @return
#' @export
#'
#' @examples
kappa_w2_sim_thr_ind<-function(W2, node1, node2, hyper_fixed){
  N<-length(W2)
  sim_kappa_w2<-rgamma(1, shape=hyper_fixed$kappa_w2[1] + 0.5*(N-1), 
                       rate= hyper_fixed$kappa_w2[2] + 0.5*sum((W2[node1]-W2[node2])^2))
  return(sim_kappa_w2)
}


#' Function to simulate the kappa_mu for threshold indicator models
#'
#' @param mu 
#' @param intercept2 
#' @param beta2 
#' @param W2 
#' @param Z2 
#'
#' @return
#' @export
#'
#' @examples
kappa_mu_sim_thr_ind<-function(mu, intercept2, beta2, W2, Z2, hyper_fixed){
  n2<-length(mu)
  sim_kappa_mu<- rgamma(1, shape=hyper_fixed$kappa_mu[1] + 0.5 * n2, 
                        rate= hyper_fixed$kappa_mu[2] + 0.5 * (sum((mu-intercept2- Z2%*%beta2 -W2)^2)))
  return(sim_kappa_mu)
}



#' Function to simulate the intercept2  for threshold indicator models
#'
#' @param mu 
#' @param beta2 
#' @param beta 
#' @param W2 
#' @param kappa_mu 
#' @param Z2 
#'
#' @return
#' @export
#'
#' @examples
intercept2_sim_thr_ind<-function(mu, beta2, W2, kappa_mu, Z2, hyper_fixed){
  n2<-length(mu)
  prec.intercept2<-hyper_fixed$intercept2+n2*kappa_mu
  mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2- W2)/prec.intercept2
  sim_intercept2<-rnorm(n=1, mean=mean.intercept2, sd=sqrt(1/prec.intercept2))
  return(sim_intercept2)
}


#' Function to simulate beta2  for threshold indicator models
#'
#' @param mu 
#' @param intercept2 
#' @param beta 
#' @param kappa_mu 
#' @param W2 
#' @param Z2 
#' @param Z2.crossprd 
#'
#' @return
#' @export
#'
#' @examples
beta2_sim_thr_ind<-function(mu, intercept2, beta, kappa_mu,  W2, Z2, Z2.crossprd, hyper_fixed){
  q<-ncol(Z2)
  latent.cov.inv<- hyper_fixed$beta2 * diag(1,q) + kappa_mu * Z2.crossprd 
  latent.mean.part<- kappa_mu * (t(Z2)%*%(mu-intercept2 - W2))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(q)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
  
}


######################################################################################
###################### updates for joint models #############################
######################################################################################

#' Function to simulate the kappa_eta for joint model
#'
#' @param eta 
#' @param intercept1 
#' @param beta1 
#' @param W1 
#' @param Z1 
#'
#' @return
#' @export
#'
#' @examples
kappa_eta_sim_JM<-function(eta, intercept1, beta1, W1, Z1, hyper_fixed){
  n1<-length(eta)
  sim_kappa_eta<-rgamma(1, shape=hyper_fixed$kappa_eta[1]+0.5*n1, 
                        rate= hyper_fixed$kappa_eta[2] + 0.5*(sum((eta-intercept1-Z1%*%beta1- W1)^2)))
  
  return(sim_kappa_eta)
}


#' Function to simulate the kappa_w1 for joint model
#'
#' @param W1 
#' @param node1 
#' @param node2 
#'
#' @return
#' @export
#'
#' @examples
kappa_w1_sim_JM<-function(W1, node1, node2, hyper_fixed){
 # browser()
  N<-length(W1)
  sim_kappa_w1<- rgamma(1, shape=hyper_fixed$kappa_w1[1]+0.5*(N-1), 
                       rate= hyper_fixed$kappa_w1[2] + 0.5*sum((W1[node1]-W1[node2])^2))
  return(sim_kappa_w1)
}


#' Function to simulate the kappa_w2 for joint model
#'
#' @param W2 
#' @param node1 
#' @param node2 
#'
#' @return
#' @export
#'
#' @examples
kappa_w2_sim_JM<-function(W2, node1, node2, hyper_fixed){
  N<-length(W2)
  sim_kappa_w2<-rgamma(1, shape=hyper_fixed$kappa_w2[1]+0.5*(N-1), 
                       rate= hyper_fixed$kappa_w2[2] + 0.5*sum((W2[node1]-W2[node2])^2))
  return(sim_kappa_w2)
}


#' Function to simulate the kappa_mu for joint model
#'
#' @param mu 
#' @param intercept2 
#' @param beta2 
#' @param beta 
#' @param W1 
#' @param W2 
#' @param Z2 
#'
#' @return
#' @export
#'
#' @examples
kappa_mu_sim_JM<-function(mu, intercept2, beta2, beta, W1, W2, Z2, hyper_fixed){
  n2<-length(mu)
  sim_kappa_mu<-rgamma(1, shape=hyper_fixed$kappa_mu[1] + 0.5 * n2, 
                       rate= hyper_fixed$kappa_mu[2] + 0.5 * (sum((mu-intercept2- Z2%*%beta2 - beta* W1-W2)^2)))
  return(sim_kappa_mu)
}



#' Function to simulate the beta for joint model
#'
#' @param mu 
#' @param intercept2 
#' @param kappa_mu 
#' @param beta2 
#' @param W1 
#' @param W2 
#' @param Z2 
#'
#' @return
#' @export
#'
#' @examples
beta_sim_JM<-function(mu, intercept2, kappa_mu, beta2, W1, W2, Z2, hyper_fixed){
  #prec_beta<-kappa_mu * (t(W1) %*% t.A2.A2 %*% W1) + hyper_fixed[9]
  #mean_beta<-kappa_mu*(t(W1)%*% t(A2)%*% (mu-intercept2-Z2%*%beta2- W2))/prec_beta  ## may be we can simplify it later
  prec_beta<-kappa_mu * sum(W1^2) + hyper_fixed$beta
  mean_beta<-kappa_mu*(sum(W1 * (mu-intercept2-Z2%*%beta2- W2)))/prec_beta  ## may be we can simplify it later
  sim_beta<-rnorm(n=1, mean=mean_beta, sd=sqrt(1/prec_beta))
  return(sim_beta)
}


#' Function to simulate the intercept1 for joint model
#'
#' @param eta 
#' @param beta1 
#' @param W1 
#' @param kappa_eta 
#' @param Z1 
#'
#' @return
#' @export
#'
#' @examples
intercept1_sim_JM<-function(eta, beta1, W1, kappa_eta, Z1, hyper_fixed){
  n1<-length(eta)
  prec.intercept1<- hyper_fixed$intercept1+n1*kappa_eta
  mean.intercept1<- kappa_eta*sum((eta-Z1%*%beta1- W1))/prec.intercept1
  sim_intercept1<-rnorm(n=1, mean=mean.intercept1, sd=sqrt(1/prec.intercept1))
  return(sim_intercept1)
}


#' function to simulate the intercept2 for joint model
#'
#' @param mu 
#' @param beta2 
#' @param beta 
#' @param W1 
#' @param W2 
#' @param kappa_mu 
#' @param Z2 
#'
#' @return
#' @export
#'
#' @examples
intercept2_sim_JM<-function(mu, beta2, beta, W1, W2, kappa_mu, Z2, hyper_fixed){
  n2<-length(mu)
  # prec.intercept2<-hyper_fixed[11]+n2*kappa_mu
  # mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2-beta*(A2%*%W1)-A2%*%W2)/prec.intercept2
  prec.intercept2<-hyper_fixed$intercept1+n2*kappa_mu
  mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2-beta*W1- W2)/prec.intercept2
  sim_intercept2<-rnorm(n=1, mean=mean.intercept2, sd=sqrt(1/prec.intercept2))
  return(sim_intercept2)
}


#' function to simulate beta1 for joint model
#'
#' @param eta 
#' @param intercept1 
#' @param W1 
#' @param kappa_eta 
#' @param Z1 
#'
#' @return
#' @export
#'
#' @examples
beta1_sim_JM<-function(eta, intercept1, W1, kappa_eta, Z1, hyper_fixed){
  p<-ncol(Z1)
  latent.cov.inv<- hyper_fixed$beta1 * diag(1,p) + kappa_eta * (t(Z1)%*%Z1)
  latent.mean.part<- kappa_eta * (t(Z1)%*%(eta-intercept1- W1))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(p)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
  
}

#' Function to simulate beta2 for joint model
#'
#' @param mu 
#' @param intercept2 
#' @param beta 
#' @param kappa_mu 
#' @param W1 
#' @param W2 
#' @param Z2 
#'
#' @return
#' @export
#'
#' @examples
beta2_sim_JM<-function(mu, intercept2, beta, kappa_mu, W1, W2, Z2, hyper_fixed){
  q<-ncol(Z2)
  latent.cov.inv<- hyper_fixed$beta2 * diag(1,q) + kappa_mu * (t(Z2)%*%Z2)
  latent.mean.part<- kappa_mu * (t(Z2)%*%(mu-intercept2 - beta * W1 - W2))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(q)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
  
}















