#############################################
########## Main MCMC function  ##############
#############################################
mcmc_sampler_threhshold_model<-function(N.MCMC, 
                                  A,
                                  ind.NA,
                                  Z2,
                                  thin,
                                  adapt,
                                  burn_in1, 
                                  burn_in2, 
                                  tun.mu, 
                                  tun.hyper.mu,
                                  hyper_fixed,
                                  print.result,
                                  traceplot, 
                                  CV, 
                                  true.values, 
                                  simulation, 
                                  nbd_info, 
                                  no_of_nbd, 
                                  node.set,
                                  ind_zeros_counts, 
                                  q.probs, 
                                  thr.family, 
                                  hyper.mu_adapt_seq2,
                                  mu_adapt_seq2,
                                  eta_adapt_seq2, 
                                  init.seed)
{ 
 
  
  n2<- sum(!ind_zeros_counts)
  q<-ncol(Z2)
  n1<- nrow(Z2)
  n.Q<- n1
  
  Z2.crossprd<- (t(Z2)%*%Z2)
  ####### Imputation in case of within sample diagnostics
  imputed.A.WSD<-rep(0, sum(!ind_zeros_counts)) # Imputed post. mean of Y in case of within sample diganostics=Impute.A.WSD/ N-(burn_in1+burn_in2)
  imputed.A.WSD.samples<-rep(0, sum(!ind_zeros_counts)) # Used for std of A in case of within sample diganostics
  
  ####### posterior summary of latent parameters
  post.sum.mean.mu <- rep(0, times=n1)
  post.sum.squre.mu<- rep(0, times=n1)
  post.sum.mean.w2<- rep(0, times=n1)
  post.sum.squre.w2<- rep(0, times=n1)
  
  imputed.A.WSD<- rep(0, times=n2)
  imputed.A.squre<- rep(0, times=n2)
  
  ### storing the probability of interest
  post.mean.quantile<-  matrix(0, nrow = n2, ncol=length(q.probs))
  post.squre.quantile<-  matrix(0, nrow = n2, ncol=length(q.probs))
  
  ## Storing the tuning parameteres
  tun.hyper.mu<- rep(tun.hyper.mu, 1) ### rep(tun.hyper.mu, 2) if the kappa and xi are updated separately
  sigma.matrix<-matrix(nrow=floor(N.MCMC/adapt),ncol=2,byrow=TRUE) ## storing the adaptive scale parameters
  sigma.matrix[1,]<- c(tun.mu, tun.hyper.mu)
  
  init.hyper.mu_and_name<- init_fun_log.hyper.mu_thr_model(thr.family=thr.family, seed=init.seed)
  cur.samples.log.hyper.mu<- init.hyper.mu_and_name$log.hyper.mu.init #current samples for gamma1
  init.all.other.param_and.name<- init_fun_all_other_param_thr_model(Z2=Z2, A=A, seed=init.seed, thr.family=thr.family,
                                                           simulation=simulation)
  
  init.other.param<- init.all.other.param_and.name$init.all.other.param
  cur.samples.kappa_w2<-init.other.param[1] #current samples for \kappa_w2
  cur.samples.kappa_mu<-init.other.param[2] #current samples for \kappa_w2
  cur.samples.intercept2<-init.other.param[3]
  cur.samples.beta2<-init.other.param[(3+1):(4-1+q)]    #current samples for beta2
  cur.samples.mu<-init.other.param[((4+q)):((4-1+q+n1))] #current samples for mu
  cur.samples.w2<-init.other.param[((4+q+n1)):((4-1+q+n1+n.Q))]  #current samples for w2
  cur.samples.w2<- cur.samples.w2-mean(cur.samples.w2) ## making sure that sum of W2 is zero
  
  ### saving the Chain for all the hyperparameters and some latent parmeters
  samples <- matrix(nrow=floor(N.MCMC/thin),ncol=length(cur.samples.log.hyper.mu)+3+q+4,byrow=TRUE) # storing the samples
  samples[1,]<- c(cur.samples.log.hyper.mu, init.other.param[c(1:3, 4:(3+q), (4+q):(3+q+2), (4+q+n1):(3+q+n1+2))]) ## saving only few samples
  
  tun.mu<- rep(tun.mu, n1)
  
  j<-1
  l<-1
  m<-1
  k<-1
  for (i in 1:(N.MCMC-1)){
    # if((i%%(adapt))-1==0 & i< (burn_in1+burn_in2+2)){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
    #   rate.latent<-0
    #   rate.hyper<-0
    # }
    
    if(((i%%(adapt))-1==0) & (i< (burn_in1+burn_in2+2))){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
      rate.mu<- rep(0, n1)
      rate.hyper.mu<- rep(0,1)
    }
    
    if((i%%thin)-1==0){  ## by is the thinning steps and adapt is the number of iterations after which i update the variance of MALA and random walk algorithms
      if(simulation==TRUE){
        par(mfrow=c(4,4),oma=c(0,0,2,0),mar=c(4,5,1,1))
      } else {
        par(mfrow=c(5,5),oma=c(0,0,2,0),mar=c(4,5,1,1))
      }
      if(print.result==TRUE){
        if(i< (burn_in1+burn_in2+2)){
          cat(paste0(" THRESHOLD MODEL:", "\n", 
                      " Iteration: ",i, "\n",
                      " Accep rate mu[1] = ", round(rate.mu[1]/((i%%(adapt))+1), digits = 3),
                      " | sigma mu[1] = ", round(tun.mu[1], digits = 5), "\n",
                     " Accep rate mu[2] = ", round(rate.mu[2]/((i%%(adapt))+1), digits = 3),
                     " | sigma mu[2] = ", round(tun.mu[2], digits = 5), "\n",
                      " Accep rate hyper mark=", round(rate.hyper.mu/((i%%(adapt))+1), digits = 3),
                      " | sigma hyper mark=", round(tun.hyper.mu, digits = 5), "\n",
                     "--------------------------------------------------------------------------", "\n",
                     sep=""))
        } else{
          
          cat(paste0(" THRESHOLD MODEL:", "\n", 
                     " Iteration: ",i, "\n",
                     " Accep rate mu[1] = ", round(rate.mu[1]/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma mu[1] = ", round(tun.mu[1], digits = 5), "\n",
                     " Accep rate mu[2] = ", round(rate.mu[2]/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma mu[2] = ", round(tun.mu[2], digits = 5), "\n",
                     " Accep rate hyper mark=", round(rate.hyper.mu/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma hyper mark=", round(tun.hyper.mu, digits = 5), "\n",
                     "--------------------------------------------------------------------------", "\n",
                     sep="")) 
        }
      }
      
      if(traceplot){
        model.param.name<-c(init.hyper.mu_and_name$log.hyper.mu.name, init.all.other.param_and.name$model.param.name.all.other.param)
        ## Plotting the traceplots
        for (ll  in 1: length(model.param.name)) {
          if(simulation==TRUE){
            plot(thin*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=model.param.name[ll]) # Plot for alpha.tilde
            abline(h=true.values[ll], col=2)
          }  else {
            plot(thin*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=model.param.name[ll]) # Plot for alpha.tilde
            
          }
        }

      }
      if((i%%adapt)-1==0){
        #plot(sigma.matrix[1:k,1],sigma.matrix[1:k,2],xlab="sigma.latent",ylab="sigma.hyper2") #plot for the scale parameters chosen adaptively
        k<-k+1
      }
      l<-l+1
    }
    ################################ imputations ############################################################
    
    if(CV=="OOS"){
      if(thr.family=="gamma"){
        imputed_A<- rgamma(n = sum(!ind_zeros_counts), scale = exp(cur.samples.mu[!ind_zeros_counts])/ k,
                           shape = k)
      } else if(thr.family=="logNormal"){
        imputed_A<- exp(rnorm(n = sum(!ind_zeros_counts), mean = cur.samples.mu[!ind_zeros_counts],
                              sd =sqrt(1/k)))
      }
      A.all<-  rep(0, length=n1)
      A.all[!ind_zeros_counts]<- imputed_A
      A[ind.NA]<- A.all[ind.NA]
    } 
    
    #### saving the results from imputed one after the burning samples 
    if(i>burn_in1+burn_in2){
      k<- cur.samples.log.hyper.mu[1]
      if(thr.family=="gamma"){
        imputed_A<- rgamma(n = sum(!ind_zeros_counts), scale = exp(cur.samples.mu[!ind_zeros_counts])/ k,
                           shape = k)
      } else if(thr.family=="logNormal"){
        imputed_A<- exp(rnorm(n = sum(!ind_zeros_counts), mean = cur.samples.mu[!ind_zeros_counts],
                              sd =sqrt(1/k)))
      }
      
      
      if(thr.family=="gamma"){
        imputed.A.WSD<- imputed.A.WSD + imputed_A
        imputed.A.squre<- imputed.A.squre + imputed_A^2
        post.sum.mean.mu<- post.sum.mean.mu + exp(cur.samples.mu) ### estimated mean 
        post.sum.squre.mu<- post.sum.squre.mu + exp(cur.samples.mu)^2 ### estimated mean counts 
        post.sum.mean.w2<- post.sum.mean.w2 + cur.samples.w2
        post.sum.squre.w2 <- post.sum.squre.w2 + cur.samples.w2^2
        qunatile_est<- matrix(NA, nrow= sum(!ind_zeros_counts), ncol=length(q.probs))
        
        for (pp in 1:length(q.probs)) {
          qunatile_est[,pp]<-  qgamma(q.probs[pp], scale = exp(cur.samples.mu[!ind_zeros_counts])/ k,
                                      shape = k)
        }
      }
      if(thr.family=="logNormal"){
        post.sum.mean.mu<- post.sum.mean.mu + cur.samples.mu ### estimated mean 
        post.sum.squre.mu<- post.sum.squre.mu + cur.samples.mu^2 ### estimated mean counts 
        imputed.A.WSD<- imputed.A.WSD + imputed_A
        imputed.A.squre<- imputed.A.squre + imputed_A^2
        post.sum.mean.w2<- post.sum.mean.w2 + cur.samples.w2
        post.sum.squre.w2 <- post.sum.squre.w2 + cur.samples.w2^2
        qunatile_est<- matrix(NA, nrow= sum(!ind_zeros_counts), ncol=length(q.probs))
        
        for (pp in 1:length(q.probs)) {
          qunatile_est[,pp]<-  exp(qnorm(q.probs[pp], mean = cur.samples.mu[!ind_zeros_counts],
                                         sd = sqrt(1/k)))
        } 
      }
      
      post.mean.quantile<- post.mean.quantile + qunatile_est
      post.squre.quantile<- post.squre.quantile + qunatile_est^2
      
    }
    
   # browser()
    ############ updating all the model parameters #######
    ###### Proposing new parameters for kappa_w2 (Gibbs steps) #####
    cur.samples.kappa_w2<- kappa_w2_sim_thr(W2=cur.samples.w2, node1 = node.set[,1], node2 = node.set[,2], hyper_fixed=hyper_fixed)  # Proposed hyperparameters using uniform random walks
    
    #cur.samples.kappa_a<- kappa_a  # Proposed hyperparameters using uniform random walks
    #### Proposing new parameters for kappa_a (Gibbs steps) #####
    cur.samples.kappa_mu<- kappa_mu_sim_thr(mu=cur.samples.mu, intercept2 = cur.samples.intercept2, beta2 = cur.samples.beta2,
                                        W2=cur.samples.w2, Z2=Z2, hyper_fixed) # Proposed hyperparameters using uniform random walks
    
    #cur.samples.kappa_mu<- true.values[3]  # Proposed hyperparameters using uniform random walks
    
    #### Proposing new parameters for hyper.mu MH steps #####
    prop_hyper_marks<- block_rw_update(par_curr=cur.samples.log.hyper.mu,
                                       par_name="cur_par",
                                       loglik=log_lik_thr,
                                       logprior=gamma_prior_markdist_thr,
                                       var_markdist=tun.hyper.mu,
                                       lb = 0,
                                       ub = Inf,
                                       transform = TRUE,
                                       mu=cur.samples.mu,
                                       A=A,
                                       thr.family=thr.family,
                                       ind_zeros_counts=ind_zeros_counts,
                                       sum_dens=TRUE,
                                       hyper.mu_fixed=c(1/10,1/100))
    
    rate.hyper.mu<- ifelse(any(prop_hyper_marks==cur.samples.log.hyper.mu), rate.hyper.mu , rate.hyper.mu + 1)
    tun.hyper.mu<- adpative_function(index_MCMC_iter=i, sigma2_adapt=tun.hyper.mu, target_accept=0.40,
                                        rate_adapt=rate.hyper.mu, burn_in1=burn_in1, burn_in2=burn_in2, adapt_seq=hyper.mu_adapt_seq2,
                                        adapt=adapt, adpat_param=1, lower.acc=0.35, upper.acc=0.55)
    cur.samples.log.hyper.mu<- prop_hyper_marks
    #cur.samples.log.hyper.mu[1]<- 0.5
    #cur.samples.log.hyper.mu<- true.values[1:2]
    #cur.samples.log.hyper.mu<- true.values[1]
    
    ###### Proposing new parameters for intercept2 (Gibbs steps) ######
    cur.samples.intercept2<- intercept2_sim_thr(mu=cur.samples.mu, beta2=cur.samples.beta2,
                                            W2=cur.samples.w2, kappa_mu = cur.samples.kappa_mu, Z2=Z2, hyper_fixed=hyper_fixed)
    
    ###### updating beta2 using Gibbs #########
    cur.samples.beta2<- beta2_sim_thr(mu=cur.samples.mu, intercept2=cur.samples.intercept2, kappa_mu=cur.samples.kappa_mu, 
                                  W2=cur.samples.w2, Z2=Z2, Z2.crossprd=Z2.crossprd, hyper_fixed=hyper_fixed)
    
    ###### updating mu using MALA #########
    prop_mu<- latent_MH_update_parallel(par_cur = cur.samples.mu,
                                        par_name = "mu",
                                        loglik_data = log_lik_thr,
                                        loglik_latent = log_lik_latent_mu_thr_ind,
                                        var_prop = tun.mu,
                                        A = A,
                                        cur_par = cur.samples.log.hyper.mu,
                                        intercept2 = cur.samples.intercept2,
                                        W2=cur.samples.w2,
                                        beta2=cur.samples.beta2,
                                        kappa_mu = cur.samples.kappa_mu,
                                        Z2=Z2,
                                        ns=n1,
                                        nt =1,
                                        thr.family=thr.family,
                                        sum_dens=FALSE,
                                        ind_zeros_counts = ind_zeros_counts
    )
    
    rate.mu<- ifelse(cur.samples.mu==prop_mu, rate.mu, rate.mu+1)
    tun.mu<-adpative_function(index_MCMC_iter=i, sigma2_adapt=tun.mu, target_accept=0.40,
                                 rate_adapt=rate.mu, burn_in1=burn_in1, burn_in2=burn_in2,
                                 adapt=adapt, adpat_param=1, adapt_seq=mu_adapt_seq2, lower.acc=0.30, upper.acc=0.50)
    
    cur.samples.mu<- prop_mu
    #cur.samples.mu<- mu
    ###### updating W2 using Gibbs #########
    mean_mu_latent_vec<- cur.samples.mu - cur.samples.intercept2 -  Z2 %*% cur.samples.beta2 
    for (ind in 1:n1) {
      cur.samples.w2[ind]<- W2_sim_Gibbs_componentwise_thr(node_index = ind,  
                                                       W2 = cur.samples.w2, 
                                                       mean_mu_latent = mean_mu_latent_vec[ind] ,
                                                       kappa_w2 = cur.samples.kappa_w2, 
                                                       kappa_mu = cur.samples.kappa_mu , 
                                                       nbd_info = nbd_info, 
                                                       no_of_nbd = no_of_nbd)
    }
    #cur.samples.w2<- W2
    cur.samples.w2<- cur.samples.w2- mean(cur.samples.w2) ## making sure that sum of W2 is zero
    
    ### Saving the samples after thinning the samples at every by iterations .  
    if((i%%thin)-1==0){
      samples[j,]<- c(cur.samples.log.hyper.mu,
                      cur.samples.kappa_w2,
                      cur.samples.kappa_mu,
                      cur.samples.intercept2,
                      cur.samples.beta2,
                      cur.samples.mu[1:2],
                      cur.samples.w2[1:2])
      
      j=j+1
    }
    ### storing the adaptive tuning parameters
    if((i%%adapt)-1==0){ # to save allexp the scale parameter of the MALA
      sigma.matrix[m,]<- c(tun.mu[1], tun.hyper.mu)
      m=m+1
    }
  }
  
  
  return(list("samples"=samples, ### saving the sample for the hyperparameters to see the traceplots
              "imputed.A.WSD"=imputed.A.WSD,
              "imputed.A.squre"=imputed.A.squre,
              "post.sum.mean.mu" = post.sum.mean.mu,  
              "post.sum.squre.mu" = post.sum.squre.mu, 
              "post.mean.quantile"= post.mean.quantile,
              "post.squre.quantile" = post.squre.quantile,
              "post.sum.mean.w2" =post.sum.mean.w2,
              "post.sum.squre.w2"= post.sum.squre.w2,
              "tuning_param_x_hyper"=sigma.matrix ## tuning parameters values
  )) ## acceptance rate for eta
}

