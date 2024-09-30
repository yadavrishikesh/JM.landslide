######################################################################################
##########  Initial values for threshold model ##########################################
######################################################################################
init_fun_log.hyper.mu_thr_model<-function(thr.family, seed){
  set.seed(seed)
  if(thr.family=="gamma"){
    log.hyper.mu<-c(runif(1, 0, 1))
    log.hyper.mu.name<- c(expression(k))
  }
  if(thr.family=="logNormal"){
    log.hyper.mu<-c(runif(1, 0, 1))
    log.hyper.mu.name<- c(expression(k))
  }
  return(list(log.hyper.mu.init=log.hyper.mu, log.hyper.mu.name=log.hyper.mu.name))
}

init_fun_all_other_param_thr_model<-function(Z2, A, thr.family, seed, simulation){
  set.seed(seed)
  if(thr.family=="gamma"){
    kappa_w2_init<- runif(1, 1, 5)
    kappa_mu_init<- runif(1, 1, 5)
    intercept2.init<- runif(1, -1, 1)
    
    n1<- nrow(Z2)
    
    beta2.init<- rep(runif(1, -1, 1), ncol(Z2))
    mu_init<- rep(runif(1, -0.1, 0.1), n1) 
    w2_init<- rep(runif(1, -1, 1), times=n1)
    init1<-c( kappa_w2_init, kappa_mu_init, intercept2.init,
              beta2.init, w2_init, mu_init)
    
    model.param.name<- c(expression(kappa[w]), expression(kappa[mu]),
                         "intercept-e", paste0("beta2_",1:ncol(Z2)), 
                         expression(mu[1]), expression(mu[2]), 
                         expression(W2[1]), expression(W2[2]))
  } 
  if(thr.family=="logNormal"){
    kappa_w2_init<- runif(1, 1, 5)
    kappa_mu_init<- runif(1, 1, 5)
    intercept2.init<- runif(1, -1, 1)
    
    n1<- nrow(Z2)
    
    beta2.init<- rep(runif(1, -1, 1), ncol(Z2))
    mu_init<- rep(runif(1, -0.1, 0.1), n1) 
    w2_init<- rep(runif(1, -1, 1), times=n1)
    init1<-c( kappa_w2_init, kappa_mu_init, intercept2.init,
              beta2.init, w2_init, mu_init)
    
    model.param.name<- c(expression(kappa[w]), expression(kappa[mu]),
                         "intercept-e", paste0("beta2_",1:ncol(Z2)), 
                         expression(mu[1]), expression(mu[2]), 
                         expression(W2[1]), expression(W2[2]))
  } 
  return(list(init.all.other.param=init1, model.param.name.all.other.param=model.param.name))
}


######################################################################################
##########  Initial values for indicator model  ##########################################
######################################################################################

init_fun_all_other_param_indicator_model<-function(Z2, A, seed){
  set.seed(seed)
  kappa_w2_init<- runif(1, 1, 5)
  kappa_mu_init<- runif(1, 1, 5)
  intercept2.init<- runif(1, -1, 1)
  
  n1<- nrow(Z2)
  
  beta2.init<- rep(runif(1, -1, 1), ncol(Z2))
  mu_init<- rep(runif(1, -0.1, 0.1), n1) 
  w2_init<- rep(runif(1, -1, 1), times=n1)
  init1<-c( kappa_w2_init, kappa_mu_init, intercept2.init,
            beta2.init, w2_init, mu_init)
  
  model.param.name<- c(expression(kappa[w]), expression(kappa[mu]),
                       "intercept-e", paste0("beta2_",1:ncol(Z2)), 
                       expression(mu[1]), expression(mu[2]), 
                       expression(W2[1]), expression(W2[2]))
  
  return(list(init.all.other.param=init1, 
              model.param.name.all.other.param=model.param.name))
}



######################################################################################
##########  Initial values for joint model  ##########################################
######################################################################################


#' Initial values and parameter names for hyperparameter in the bulk
#'
#' @param mark_dist 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
init_fun_hyper.mu_JM<-function(mark_dist, seed){
  set.seed(seed)
  if (mark_dist=="eGPD"){
    log.hyper.mu<- list("k" =runif(1, 0, 10), "xi"=runif(1, 0, 0.2))
    log.hyper.mu.name<- c(expression(k), expression(xi))
  } else if(mark_dist=="bGPD"){
    log.hyper.mu<- list("k" =runif(1, 0, 10))
    log.hyper.mu.name<- c(expression(k)) 
  } else if(mark_dist=="tgGPD"){
    log.hyper.mu<-  list("k" =runif(1, 0, 10))
    log.hyper.mu.name<- c(expression(k)) 
  }
  return(list(log.hyper.mu.init=log.hyper.mu, 
              log.hyper.mu.name=log.hyper.mu.name))
}


#' Initial values and parameter names for hyperparameter for GP
#'
#' @param mark_dist 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
init_fun_hyper.GP_JM<-function(mark_dist, seed){
  set.seed(seed)
  if (mark_dist=="eGPD"){
    log.hyper.mu<- NULL
    log.hyper.mu.name<- NULL
  } else if(mark_dist=="bGPD"){
    log.hyper.mu<- list("sigma.GP" =runif(1, 0, 10), "xi" = runif(1, 0, 0.2))
    log.hyper.mu.name<- c(expression(sigma[GP]), expression(xi)) 
  } else if(mark_dist=="tgGPD"){
    log.hyper.mu<-  list("sigma.GP" =runif(1, 0, 10), "xi" = runif(1, 0, 0.2))
    log.hyper.mu.name<- c(expression(sigma[GP]), expression(xi)) 
  }
  return(list(log.hyper.mu.init=log.hyper.mu, 
              log.hyper.mu.name=log.hyper.mu.name))
}

#' nitial values of all the other parameters, that does not depends on the choice of mark distributions 
#'
#' @param Z1 
#' @param Z2 
#' @param A 
#' @param Y 
#' @param seed 
#' @param simulation 
#' @param threshold 
#' @param thr.acces.ind 
#'
#' @return
#' @export
#'
#' @examples
init_fun_all_other_param_JM<-function(Z1, Z2, A, Y, seed, simulation, threshold, mark_dist, model_type, thr.acces.ind){
   #browser()
  set.seed(seed)
  kappa_w1_init<- runif(1, 1, 5)
  kappa_w2_init<- runif(1, 1, 5)
  kappa_eta_init<- runif(1, 1, 5)
  kappa_mu_init<- runif(1, 1, 5)
  intercept1.init<- runif(1, -1, 1)
  intercept2.init<- runif(1, -1, 1)
  beta.init<- runif(1, -1, 1)
  
  n.Q<- length(Y)
  
  beta1.init<-rep(runif(1, -1, 1), ncol(Z1))
  beta2.init<-rep(runif(1, -1, 1), ncol(Z2))
  eta_init<-rep(runif(1, -1, 1), times=length(Y))
  W1_init<-rep(runif(1, -1, 1), times=n.Q)
  if(mark_dist=="eGPD"){
    mu_init<- rep(runif(1, -0.1, 0.1), nrow(Z1))
  } else {
    mu_init<- ifelse(thr.acces.ind, runif(1, -0.1, 0.1) + log(threshold), runif(1, -0.1, 0.1))
  }
  w2_init<-rep(runif(1, -1, 1), times=n.Q)
  init<- list("kappa_w1_init" = kappa_w1_init, "kappa_w2_init"=kappa_w2_init, "kappa_eta_init"=kappa_eta_init, 
              "kappa_mu_init"=kappa_mu_init, "intercept1.init"=intercept1.init, "intercept2.init"=intercept2.init,
              "beta.init"=beta.init, "beta1.init"=beta1.init, "beta2.init"=beta2.init, "eta_init"=eta_init, 
              "W1_init"=W1_init, "mu_init"=mu_init,  "w2_init"=w2_init)
  if(model_type=="FE"){
    model.param.name<-c(expression(kappa[eta]), expression(kappa[mu]),
                        "intercept1","intercept2", paste0("beta1_",1:ncol(Z1)), paste0("beta2_",1:ncol(Z2)), 
                        expression(eta[1]), expression(eta[2]), 
                        expression(mu[1]))
  } else {
    model.param.name<-c(expression(kappa[w1]), expression(kappa[w2]), expression(kappa[eta]), expression(kappa[mu]),
                        "intercept1","intercept2", expression(beta), paste0("beta1_",1:ncol(Z1)), paste0("beta2_",1:ncol(Z2)), 
                        expression(eta[1]), expression(eta[2]), expression(W1[1]), expression(mu[1]), expression(W2[1]))
  }
  return(list(init.all.other.param=init, 
              model.param.name.all.other.param=model.param.name))
}






