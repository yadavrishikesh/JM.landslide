#############################################
########## Main MCMC function  ##############
#############################################
mcmc_sampler_joint_model_jSp<-function(N.MCMC, 
                                       Y,
                                       ind_NA_Y, 
                                       A, 
                                       ind_NA_A, 
                                       Z1, 
                                       Z2, 
                                       model_type,
                                       Sim_data,
                                       thin, 
                                       adapt, 
                                       burn_in1, 
                                       burn_in2, 
                                       tun.eta, 
                                       tun.mu,
                                       tun.hyper.mu, 
                                       tun.hyper.GP,
                                       mark_dist,
                                       hyper_fixed,
                                       print.result=TRUE, 
                                       traceplot=FALSE, 
                                       model.base,
                                       CV,
                                       true.values=NULL, 
                                       simulation, 
                                       nbd_info, 
                                       no_of_nbd, 
                                       node.set,
                                       ind_zeros_counts, 
                                       threshold, 
                                       thr.acces.ind, 
                                       thr.prob, 
                                       q.probs, 
                                       hyper.mu_adapt_seq2,
                                       mu_adapt_seq2,
                                       eta_adapt_seq2, 
                                       samples.store,
                                       init.seed)
{ 
  n1<-length(Y)
  n2<- length(A)
  n.Q<-length(Y)
  p<-ncol(Z1)
  q<-ncol(Z2)
  
  ####### Imputation in case of within sample diagnostics
  imputed.Y.WSD<-rep(0, n1)  # Imputed post. mean of Y in case of within sample diganostics=Impute.Y.WSD / N-(burn_in1+burn_in2)
  imputed.Y2.WSD<-rep(0, n1)  # Imputed post. mean of Y in case of within sample diganostics=Impute.Y.WSD / N-(burn_in1+burn_in2)
  data.mean.A<- mean(A, na.rm=TRUE)
  data.mean.Y<- mean(Y, na.rm=TRUE)
  imputed.A.WSD<- rep(0, sum(!ind_zeros_counts)) # Imputed post. mean of Y in case of within sample diganostics=Impute.A.WSD/ N-(burn_in1+burn_in2)
  imputed.A2.WSD<-rep(0, sum(!ind_zeros_counts)) # Imputed post. mean of Y in case of within sample diganostics=Impute.A.WSD/ N-(burn_in1+burn_in2)
  
  no.samples<- samples.store
  samples.save<-  floor(seq(from=burn_in1+burn_in2+1, to=N.MCMC, length.out=no.samples))
  samples.save[no.samples]<- samples.save[no.samples]-1 ### last samples 
  
  imputed.Y.WSD.samples<- array(NA, dim = c(no.samples, n1))
  imputed.A.WSD.samples<- array(NA, dim = c(no.samples,  sum(!ind_zeros_counts)))
  ####### Imputation in case of out of sample diagnostics
  imputed.Y.OSD<- rep(0, times=sum(ind_NA_Y))  # Imputed  summation values for Y
  imputed.A.OSD<-rep(0, times=sum(ind_NA_A))  # Imputed  summation values for A
  imputed.Y.OSD.samples<- array(NA, dim = c(no.samples, sum(ind_NA_Y)))
  imputed.A.OSD.samples<- array(NA, dim = c(no.samples,  sum(ind_NA_Y)))
  
  ####### posterior summary of latent parameters
  post.sum.mean.mu <- rep(0, times=n2)
  post.sum.squre.mu<- rep(0, times=n2)
  post.sum.mean.eta<- rep(0, times=n2)
  post.sum.squre.eta<- rep(0, times=n2)
  post.sum.mean.w1<- rep(0, times=n2)
  post.sum.squre.w1<- rep(0, times=n2)
  post.sum.mean.w2<- rep(0, times=n2)
  post.sum.squre.w2<- rep(0, times=n2)
  
  ### storing the probability of interest
  post.mean.condprob<-  matrix(0, nrow = n1, ncol=length(q.probs))
  post.mean.uncondprob<-  matrix(0, nrow = n1, ncol=length(q.probs))
  post.squre.condprob<-  matrix(0, nrow = n1, ncol=length(q.probs))
  post.squre.uncondprob<-  matrix(0, nrow = n1, ncol=length(q.probs))
  
  
  ## Storing the tuning parameteres
  tun.hyper.mu<- rep(tun.hyper.mu, 1) ### rep(tun.hyper.mu, 2) if the kappa and xi are updated separately
  sigma.matrix<-matrix(nrow=floor(N.MCMC/adapt),ncol=4,byrow=TRUE) ## storing the adaptive scale parameters
  sigma.matrix[1,]<-c(tun.eta, tun.mu, tun.hyper.mu, tun.hyper.GP)
  
  ### initial values 
  cur.samples.log.hyper.mu<- as.numeric(init_fun_hyper.mu_JM(mark_dist=mark_dist, seed=init.seed)$log.hyper.mu.init) #current samples for gamma1
  cur.samples.hyper.GP<- as.numeric(init_fun_hyper.GP_JM(mark_dist=mark_dist, seed=init.seed)$log.hyper.mu.init)

  init.all.other.param_and.name<- init_fun_all_other_param_JM(Z1=Z1, Z2=Z2, A=A, Y=Y, seed=init.seed,
                                                           simulation=simulation, threshold = threshold, mark_dist=mark_dist,
                                                           model_type=model_type, thr.acces.ind=thr.acces.ind)
  cur.samples.kappa_w1<- init.all.other.param_and.name$init.all.other.param$kappa_w1_init
  cur.samples.kappa_w2<- init.all.other.param_and.name$init.all.other.param$kappa_w2_init
  cur.samples.kappa_eta<- init.all.other.param_and.name$init.all.other.param$kappa_eta_init
  cur.samples.kappa_mu<- init.all.other.param_and.name$init.all.other.param$kappa_mu_init
  cur.samples.intercept1<- init.all.other.param_and.name$init.all.other.param$intercept1.init
  cur.samples.intercept2<- init.all.other.param_and.name$init.all.other.param$intercept2.init
  cur.samples.beta<- init.all.other.param_and.name$init.all.other.param$beta.init
  cur.samples.beta1<- init.all.other.param_and.name$init.all.other.param$beta1.init
  cur.samples.beta2<- init.all.other.param_and.name$init.all.other.param$beta2.init  ##urrent samples for beta2
  cur.samples.eta<- init.all.other.param_and.name$init.all.other.param$eta_init    #current samples for eta
  cur.samples.w1<- init.all.other.param_and.name$init.all.other.param$W1_init   #current samples for w1
  cur.samples.mu<- init.all.other.param_and.name$init.all.other.param$mu_init #current samples for mu
  cur.samples.w2<- init.all.other.param_and.name$init.all.other.param$w2_init  #current samples for w2
  
  ### saving the Chain for all the hyperparameters and some latent parmeters
  samples <- matrix(nrow=floor(N.MCMC/thin),ncol=length(cur.samples.log.hyper.mu) + length(cur.samples.hyper.GP) + 7+p+q+5,byrow=TRUE) # storing the samples
  samples[1,]<- c(cur.samples.log.hyper.mu, cur.samples.hyper.GP, cur.samples.kappa_w1, cur.samples.kappa_w2, 
                  cur.samples.kappa_eta, cur.samples.kappa_mu, cur.samples.intercept1, cur.samples.intercept2, cur.samples.beta,
                  cur.samples.beta1, cur.samples.beta2, cur.samples.eta[1:2], cur.samples.w1[1], cur.samples.mu[1], cur.samples.w2[1])
  
  tun.mu<- rep(tun.mu, n1)
  tun.eta<- rep(tun.eta, n2)
  
  
  j<-1
  l<-1
  m<-1
  k<-1
  ls<- 1
  for (i in 1:(N.MCMC-1)){
    if(((i%%(adapt))-1==0) & (i< (burn_in1+burn_in2+2))){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
      rate.eta<- rep(0, n1)
      rate.mu<- rep(0, n1)
      rate.hyper.mu<- rep(0,1)
      rate.hyper.GP<- rep(0,1)
    }
    
    if((i%%thin)-1==0){  ## thin is the thinning steps and adapt is the number of iterations after which i update the variance of MALA and random walk algorithms
      par(mfrow=c(7,7),oma=c(0,0,2,0),mar=c(4,5,1,1))
      if(print.result==TRUE){
        if(i< (burn_in1+burn_in2+2)){
          cat(paste0(" JOINT SPATIAL MODEL:", "\n", 
                     " Iteration: ",i, "\n",
                     " Accep rate eta = ",round(rate.eta[1]/((i%%(adapt))+1), digits = 3),
                     " | sigma eta = ",round(tun.eta[1], digits = 5), "\n",
                     " Accep rate mu = ", round(rate.mu[1]/((i%%(adapt))+1), digits = 3),
                     " | sigma mu = ", round(tun.mu[1], digits = 5), "\n",
                     " Accep rate hyper mark = ",round(rate.hyper.mu/((i%%(adapt))+1), digits = 3),
                     " | sigma hyper mark = ", round(tun.hyper.mu, digits = 5), "\n",
                     " Accep rate hyper GP = ", round(rate.hyper.GP/((i%%(adapt))+1), digits = 3),
                     " | sigma hyper GP = ", round(tun.hyper.GP, digits = 3), "\n",
                     "--------------------------------------------------------------------------", "\n",
                     sep=""))
        } else{
          cat(paste0(" JOINT SPATIAL MODEL:", "\n", 
                     " Iteration: ",i, "\n",
                     " Accep rate eta = ",round(rate.eta[1]/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma eta = ",round(tun.eta[1], digits = 5), "\n",
                     " Accep rate mu = ", round(rate.mu[1]/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma mu = ", round(tun.mu[1], digits = 5), "\n",
                     " Accep rate hyper mark = ",round(rate.hyper.mu/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma hyper mark = ", round(tun.hyper.mu, digits = 5), "\n",
                     " Accep rate hyper GP = ", round(rate.hyper.GP/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma hyper GP = ", round(tun.hyper.GP, digits = 3), "\n",
                     "--------------------------------------------------------------------------", "\n",
                     sep=""))
        }
      }
      
      if(traceplot){
        model.param.name<- c(init_fun_hyper.mu_JM(mark_dist=mark_dist, seed=init.seed)$log.hyper.mu.name,
                             init_fun_hyper.GP_JM(mark_dist=mark_dist, seed=init.seed)$log.hyper.mu.name,
                             init.all.other.param_and.name$model.param.name.all.other.param)
        ## Plotting the traceplots
        
        if(simulation==TRUE){
          par(mfrow=c(5,4))
        }
        for (ll  in 1: length(model.param.name)) {
          if(simulation==TRUE){
            plot(thin*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=model.param.name[ll]) # Plot for alpha.tilde
            abline(h=true.values[ll], col=2)
          }  else {
            plot(thin*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=model.param.name[ll]) # Plot for alpha.tilde
            
          }
        }
        #hist(cur.samples.mu)
        
      }
      
      if((i%%adapt)-1==0){
        #plot(sigma.matrix[1:k,1],sigma.matrix[1:k,2],xlab="sigma.latent",ylab="sigma.hyper2") #plot for the scale parameters chosen adaptively
        k<-k+1
      }
      
      l<-l+1
    }
    
    
    ################################ imputations ############################################################
   #browser()
    
    if(CV=="OOS"){
      # browser()
      ### Imputation of Y
      imputed_Y<- impute.NA.Y(ind_NA_Y = ind_NA_Y, eta = cur.samples.eta, CV=CV)
      Y[ind_NA_Y]<- imputed_Y   #replacing the NA by the imputed values
      
      ### Imputation of A
      imputed_A<- impute.NA.A(CV=CV, ind_NA_A = ind_NA_A, ind_zeros_counts=ind_zeros_counts, mu=cur.samples.mu, thr.prob=thr.prob,
                              cur_par = c(cur.samples.log.hyper.mu, cur.samples.hyper.GP), mark_dist = mark_dist, threshold=threshold 
      )
      # A.all<-  rep(0, length=n1)
      # A.all[!ind_zeros_counts]<- imputed_A
      # A[ind_NA_A]<- A.all[ind_NA_A]
      A[ind_NA_A]<- as.numeric(imputed_A$imputed_NA_A)
      
      if(i>burn_in1+burn_in2){ ## storing the samples: Posterior predictive mean and standard variability
        imputed.Y.OSD<-imputed.Y.OSD + imputed_Y / data.mean.Y
        imputed.A.OSD<-imputed.A.OSD + sqrt(imputed_A$imputed_NA_A)
      }
      
      if(i == samples.save[ls]){
        imputed.Y.OSD.samples[ls,]<- imputed_Y
        imputed.A.OSD.samples[ls,]<- imputed_A$imputed_NA_A
        ls<- ls + 1
      }
    }
    
    
    
    if(CV=="WS"){
      if(i>burn_in1+burn_in2){
        imputed_Y<- impute.NA.Y(ind_NA_Y = ind_NA_Y, eta = cur.samples.eta, CV=CV)
        imputed_A<- impute.NA.A(CV=CV, ind_NA_A = ind_NA_A, ind_zeros_counts=ind_zeros_counts, mu=cur.samples.mu, thr.prob=thr.prob,
                                cur_par = c(cur.samples.log.hyper.mu,cur.samples.hyper.GP), mark_dist = mark_dist, threshold=threshold)
        
        imputed.Y.WSD<- imputed.Y.WSD + imputed_Y/data.mean.Y
        imputed.Y2.WSD<- imputed.Y2.WSD + (imputed_Y^2)/(data.mean.Y^2)
        
        imputed.A.WSD<-imputed.A.WSD + sqrt(imputed_A$imputed_NA_A)/data.mean.A
        imputed.A2.WSD<-imputed.A2.WSD + (sqrt(imputed_A$imputed_NA_A)^2)/(data.mean.A^2)
      }
      
      if(i == samples.save[ls]){
        imputed.Y.WSD.samples[ls,]<- imputed_Y
        imputed.A.WSD.samples[ls,]<- imputed_A$imputed_NA_A
        ls<- ls + 1
      }
      
    }
    
    
    ########################## saving the posterior mean and standard deviations of some interesting quantity ########################
    if(i>burn_in1+burn_in2){
      ### posteriors summary of latent mean vectors mu and eta
      post.sum.mean.mu <- post.sum.mean.mu + cur.samples.mu
      post.sum.squre.mu<- post.sum.squre.mu + cur.samples.mu^2
      post.sum.mean.eta<- post.sum.mean.eta + cur.samples.eta
      post.sum.squre.eta<- post.sum.squre.eta + cur.samples.eta^2
      
      ### posteriors summary of ICAR vectors W1 and W2
      post.sum.mean.w1<- post.sum.mean.w1 + cur.samples.w1
      post.sum.squre.w1<- post.sum.squre.w1 + cur.samples.w1^2
      post.sum.mean.w2<- post.sum.mean.w2 + cur.samples.w2
      post.sum.squre.w2<- post.sum.squre.w2 + cur.samples.w2^2
      
      #### probability and conditional probability estimations 
      cond.probs.est<- uncond.probs.est<- matrix(NA, nrow=n1, ncol=length(q.probs))
      if(mark_dist=="eGPD"){
        k<- cur.samples.log.hyper.mu[1]
        xi<- cur.samples.log.hyper.mu[2]
        sigma<- exp(cur.samples.mu) / evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
        for (pp in 1:length(q.probs)) {
          cond.probs.est[,pp]<-  1- pEGPD1(x = q.probs[pp], k = k, xi =xi, sigma = sigma, log = FALSE)
          uncond.probs.est[,pp]<- cond.probs.est[,pp] * (1-exp(-exp(cur.samples.eta)))
        }
        post.mean.condprob<- post.mean.condprob + cond.probs.est
        post.mean.uncondprob<- post.mean.uncondprob + uncond.probs.est
        post.squre.condprob<- post.squre.condprob + cond.probs.est^2
        post.squre.uncondprob<- post.squre.uncondprob + uncond.probs.est^2
      } 
      # else if(mark_dist=="bGPD"){
      #   k<- cur.samples.log.hyper.mu[1]
      #   sigma.GP<- cur.samples.hyper.GP[1]
      #   xi<- cur.samples.hyper.GP[2]
      #   sigma<- ifelse(A > threshold, sigma.GP, (threshold / exp(cur.samples.mu) - 1) * k)
      #   #thr.probs.r<- mean(A > threshold)
      #   for (pp in 1:length(q.probs)) {
      #     cond.probs.est[,pp]<-  1- cond.pbGPD(x= q.probs[pp], mu=cur.samples.mu, u.thr = threshold, k = k, xi = xi,
      #                                     sigma = sigma, sigma.GP=sigma.GP, thr.prob = thr.prob) 
      #     uncond.probs.est[,pp]<- cond.probs.est[,pp] * (1-exp(-exp(cur.samples.eta)))
      #   }
      # } else { ## for tgGPD
      #   k<- cur.samples.log.hyper.mu[1]
      #   sigma.GP<- cur.samples.hyper.GP[1]
      #   xi<- cur.samples.hyper.GP[2]
      #   sigma<- ifelse(A > threshold, sigma.GP, exp(cur.samples.mu)/k)
      #   #thr.probs.r<- mean(A > threshold)
      #   for (pp in 1:length(q.probs)) {
      #     cond.probs.est[,pp]<-  1 - cond.ptgGPD(x= q.probs[pp], mu=cur.samples.mu, u.thr = threshold, k = k, xi = xi,
      #                                          sigma = sigma, sigma.GP=sigma.GP, thr.prob = thr.prob)  
      #     uncond.probs.est[,pp]<- cond.probs.est[,pp] * (1-exp(-exp(cur.samples.eta)))
      #   } 
      # }
      
      # post.mean.condprob<- post.mean.condprob + cond.probs.est
      # post.mean.uncondprob<- post.mean.uncondprob + uncond.probs.est
      # post.squre.condprob<- post.squre.condprob + cond.probs.est^2
      # post.squre.uncondprob<- post.squre.uncondprob + uncond.probs.est^2
      
    }
    #browser()
    ###################### updating all the model parameters ##################
   # browser()
    #### Proposing new parameters for kappa_eta (Gibbs steps) #####
    cur.samples.kappa_eta<- kappa_eta_sim_JM(eta = cur.samples.eta, intercept1 = cur.samples.intercept1, 
                                         beta1 = cur.samples.beta1, W1=cur.samples.w1, Z1=Z1, hyper_fixed=hyper_fixed)
    #### Proposing new parameters for kappa_w1 (Gibbs steps) #####
    cur.samples.kappa_w1<- kappa_w1_sim_JM(W1=cur.samples.w1, node1 = node.set[,1], node2 = node.set[,2], hyper_fixed=hyper_fixed)  
    #cur.samples.kappa_w1<- kappa_w1
    ###### Proposing new parameters for kappa_w2 (Gibbs steps) #####
    cur.samples.kappa_w2<- kappa_w2_sim_JM(W2=cur.samples.w2, node1 = node.set[,1], node2 = node.set[,2], hyper_fixed=hyper_fixed)  # Proposed hyperparameters using uniform random walks
    
    #cur.samples.kappa_a<- kappa_a  # Proposed hyperparameters using uniform random walks
    #### Proposing new parameters for kappa_a (Gibbs steps) #####
    cur.samples.kappa_mu<- kappa_mu_sim_JM(mu=cur.samples.mu, intercept2 = cur.samples.intercept2, beta2 = cur.samples.beta2,
                                        beta=cur.samples.beta, W1=cur.samples.w1, W2=cur.samples.w2, Z2=Z2, hyper_fixed=hyper_fixed) # Proposed hyperparameters using uniform random walks
    #cur.samples.kappa_mu<- 3  # Proposed hyperparameters using uniform random walks
    
    #### Proposing new parameters for hyper.mu MH steps #####
    prop_hyper_marks<- block_rw_update(par_curr=cur.samples.log.hyper.mu,
                                       par_name="cur_par",
                                       loglik=log_lik_shape_mixture,
                                       logprior=gamma_prior_markdist_JM,
                                       var_markdist=tun.hyper.mu,
                                       lb = lb_ub_markdist(mark_dist = mark_dist)$lb,
                                       ub = lb_ub_markdist(mark_dist = mark_dist)$ub,
                                       transform = TRUE,
                                       mu=cur.samples.mu,
                                       A=A,
                                       mark_dist=mark_dist,
                                       ind_zeros_counts=ind_zeros_counts,
                                       sum_dens=TRUE,
                                       threshold=threshold,
                                       hyper.mu_fixed=c(1/10,1/100))

    rate.hyper.mu<- ifelse(any(prop_hyper_marks==cur.samples.log.hyper.mu), rate.hyper.mu , rate.hyper.mu + 1)
    tun.hyper.mu<- adpative_function(index_MCMC_iter=i, sigma2_adapt=tun.hyper.mu, target_accept=0.40,
                                     rate_adapt=rate.hyper.mu, burn_in1=burn_in1, burn_in2=burn_in2, adapt_seq=hyper.mu_adapt_seq2,
                                     adapt=adapt, adpat_param=1, lower.acc=0.35, upper.acc=0.55)
    cur.samples.log.hyper.mu<- prop_hyper_marks
    
    #cur.samples.log.hyper.mu<- true.values[1:2]
    #cur.samples.log.hyper.mu<- c(1200, 0.01)
    #  print(cur.samples.log.hyper.mu[2])
     #print(cur.samples.log.hyper.mu)
    
    #### Proposing new parameters for hyper.GP MH steps #####
    if(mark_dist=="bGPD" | mark_dist=="tgGPD"){
      prop_hyper_GP<- block_rw_update(par_curr=cur.samples.hyper.GP,
                                      par_name="cur_par",
                                      loglik=log_lik_GPD_param,
                                      logprior=gamma_prior_GPD_param,
                                      var_markdist=tun.hyper.GP,
                                      lb = lb_ub_GP(mark_dist = mark_dist)$lb,
                                      ub = lb_ub_GP(mark_dist = mark_dist)$ub,
                                      transform = TRUE,
                                      A=A,
                                      threshold=threshold,
                                      hyper.mu_fixed=c(1/10,1/100))
      
      rate.hyper.GP<- ifelse(any(prop_hyper_GP==cur.samples.hyper.GP), rate.hyper.GP , rate.hyper.GP + 1)
      tun.hyper.GP<- adpative_function(index_MCMC_iter=i, sigma2_adapt=tun.hyper.GP, target_accept=0.40,
                                       rate_adapt=rate.hyper.GP, burn_in1=burn_in1, burn_in2=burn_in2,
                                       adapt_seq=hyper.mu_adapt_seq2,
                                       adapt=adapt, adpat_param=1, lower.acc=0.35, upper.acc=0.55)
      cur.samples.hyper.GP<- prop_hyper_GP
    } else {
      cur.samples.hyper.GP<- NULL
    }
    
    ###### Proposing new parameters for beta (Gibbs steps) ######
    if(model.base==TRUE){
      cur.samples.beta=0
    } else {
      cur.samples.beta<- beta_sim_JM(mu=cur.samples.mu, intercept2 = cur.samples.intercept2, kappa_mu = cur.samples.kappa_mu,
                                  beta2=cur.samples.beta2, W1=cur.samples.w1, W2=cur.samples.w2, Z2=Z2, hyper_fixed=hyper_fixed)
    }
    #cur.samples.log.hyper.mu[1]<-log(5)
    #cur.samples.beta<-beta
    ###### Proposing new parameters for intercept1 (Gibbs steps) ######
    cur.samples.intercept1<- intercept1_sim_JM(eta = cur.samples.eta, beta1 = cur.samples.beta1,
                                            W1=cur.samples.w1, kappa_eta = cur.samples.kappa_eta, Z1=Z1, hyper_fixed=hyper_fixed)
    # alternative way to estimate intercept2
    #cur.samples.intercept1<-mean(cur.samples.eta-Z1%*%cur.samples.beta1-A1%*%cur.samples.w1)
    
    ###### Proposing new parameters for intercept2 (Gibbs steps) ######
    cur.samples.intercept2<- intercept2_sim_JM(mu=cur.samples.mu, beta2=cur.samples.beta2, beta=cur.samples.beta,
                                            W1=cur.samples.w1, W2=cur.samples.w2, kappa_mu = cur.samples.kappa_mu, Z2=Z2, hyper_fixed=hyper_fixed)
    ## alternative way to estimate intercept2
    #cur.samples.intercept2<-mean(cur.samples.mu-Z2%*%cur.samples.beta2-cur.samples.beta*(A2%*%cur.samples.w1)-A2%*%cur.samples.w2)
    
    ###### Proposing new parameters for beta1 (Gibbs steps) ######
    cur.samples.beta1<- beta1_sim_JM(eta = cur.samples.eta, intercept1 = cur.samples.intercept1, W1=cur.samples.w1,
                                  kappa_eta = cur.samples.kappa_eta, Z1=Z1, hyper_fixed=hyper_fixed)
    
    #print(cur.samples.beta1)
    ###### updating beta2 using Gibbs #########
    cur.samples.beta2<- beta2_sim_JM(mu=cur.samples.mu, intercept2=cur.samples.intercept2, beta=cur.samples.beta, kappa_mu=cur.samples.kappa_mu, 
                                 W1=cur.samples.w1, W2=cur.samples.w2, Z2=Z2, hyper_fixed=hyper_fixed)
    
    ###### updating mu using MALA #########
    #browser()
    if(mark_dist=="eGPD"){
      prop_mu<- latent_MH_update_parallel(par_cur = cur.samples.mu,
                                          par_name = "mu",
                                          loglik_data = log_lik_shape_mixture,
                                          loglik_latent = log_lik_latent_mu_JM,
                                          var_prop = tun.mu,
                                          A = A,
                                          cur_par = cur.samples.log.hyper.mu,
                                          intercept2 = cur.samples.intercept2,
                                          W1=cur.samples.w1,
                                          W2=cur.samples.w2,
                                          beta=cur.samples.beta,
                                          beta2=cur.samples.beta2,
                                          kappa_mu = cur.samples.kappa_mu,
                                          Z2=Z2,
                                          ns=n2,
                                          nt =1,
                                          sum_dens=FALSE,
                                          ind_zeros_counts = ind_zeros_counts,
                                          threshold=threshold,
                                          mark_dist=mark_dist)
    } else{
      ### lower and upper limits of mus
      lbs<- rep(-Inf, n1)
      ubs<- rep(Inf,n1)
      inds_less_thr<- (A <threshold) & (A!=0)
      ubs[inds_less_thr]<- log(threshold[inds_less_thr])

      prop_mu<- latent_MH_update_parallel_constrain(par_cur = cur.samples.mu,
                                                    par_name = "mu",
                                                    loglik_data = log_lik_shape_mixture,
                                                    loglik_latent = log_lik_latent_mu_JM,
                                                    var_prop = tun.mu,
                                                    A = A,
                                                    cur_par = cur.samples.log.hyper.mu,
                                                    intercept2 = cur.samples.intercept2,
                                                    W1=cur.samples.w1,
                                                    W2=cur.samples.w2,
                                                    beta=cur.samples.beta,
                                                    beta2=cur.samples.beta2,
                                                    kappa_mu = cur.samples.kappa_mu,
                                                    Z2=Z2,
                                                    ns=n2,
                                                    transform = ifelse(inds_less_thr, TRUE, FALSE),
                                                    lb = lbs,
                                                    ub = ubs,
                                                    nt =1,
                                                    sum_dens=FALSE,
                                                    ind_zeros_counts = ind_zeros_counts,
                                                    threshold=threshold,
                                                    mark_dist=mark_dist)
    }

    rate.mu<- ifelse(cur.samples.mu==prop_mu, rate.mu, rate.mu+1)
    tun.mu<-adpative_function(index_MCMC_iter=i, sigma2_adapt=tun.mu, target_accept=0.40,
                              rate_adapt=rate.mu, burn_in1=burn_in1, burn_in2=burn_in2,
                              adapt=adapt, adpat_param=1, adapt_seq=mu_adapt_seq2, lower.acc=0.30, upper.acc=0.50)

    cur.samples.mu<- prop_mu

    ###### updating eta using MALA #########
    #cur.samples.mu<- mu
    ###### updating eta using MALA #########
    prop_eta<- latent_MH_update_parallel(par_cur = cur.samples.eta,
                                         par_name = "eta",
                                         loglik_data = log_lik_eta_JM,
                                         loglik_latent = log_lik_latent_eta_JM,
                                         ns = n1,
                                         nt = 1,
                                         var_prop = tun.eta,
                                         Y = Y,
                                         intercept1 = cur.samples.intercept1,
                                         beta1 = cur.samples.beta1,
                                         kappa_eta = cur.samples.kappa_eta,
                                         W1=cur.samples.w1,
                                         Z1=Z1)

    rate.eta<- ifelse(cur.samples.eta==prop_eta, rate.eta, rate.eta+1)
    tun.eta<-adpative_function(index_MCMC_iter=i, sigma2_adapt=tun.eta, target_accept=0.40,
                               rate_adapt=rate.eta, burn_in1=burn_in1, burn_in2=burn_in2, adapt_seq=eta_adapt_seq2,
                               adapt=adapt, adpat_param=1, lower.acc=0.30, upper.acc=0.50)
    cur.samples.eta<- prop_eta

    
    ###### updating W1 using Gibbs #########
    mean_eta_latent_vec<- cur.samples.eta - cur.samples.intercept1  - Z1 %*% cur.samples.beta1
    mean_mu_latent_vec<- cur.samples.mu - cur.samples.intercept2 -  Z2 %*% cur.samples.beta2 - cur.samples.w2
    
    # ### updating W1s in a for loop sequentially  using Gibbs sampler steps
    for (ind in 1:n1) {
      cur.samples.w1[ind]<- W1_sim_Gibbs_componentwise_JM(node_index = ind,
                                                       W1 = cur.samples.w1,
                                                       mean_eta_latent = mean_eta_latent_vec[ind],
                                                       mean_mu_latent =  mean_mu_latent_vec[ind],
                                                       kappa_w1 = cur.samples.kappa_w1,
                                                       kappa_eta = cur.samples.kappa_eta,
                                                       kappa_mu = cur.samples.kappa_mu,
                                                       beta = cur.samples.beta,
                                                       nbd_info = nbd_info,
                                                       no_of_nbd = no_of_nbd)

    }
    ### updating W1s in parallel  using Gibbs sampler steps
   # cur.samples.w1<- Sim_mark_data$W1
   cur.samples.w1<- cur.samples.w1 - mean(cur.samples.w1) ## making sure that sum of W1 is zero
    
    #print(range(cur.samples.w1))
    ###### updating W2 using Gibbs #########
    mean_mu_latent_vec<- cur.samples.mu - cur.samples.intercept2 -  Z2 %*% cur.samples.beta2 - cur.samples.beta * cur.samples.w1
    for (ind in 1:n1) {
      cur.samples.w2[ind]<- W2_sim_Gibbs_componentwise_JM(node_index = ind,
                                                       W2 = cur.samples.w2,
                                                       mean_mu_latent = mean_mu_latent_vec[ind] ,
                                                       kappa_w2 = cur.samples.kappa_w2,
                                                       kappa_mu = cur.samples.kappa_mu ,
                                                       nbd_info = nbd_info,
                                                       no_of_nbd = no_of_nbd)
    }
    #cur.samples.w2<- Sim_mark_data$W2
    cur.samples.w2<- cur.samples.w2- mean(cur.samples.w2) ## making sure that sum of W2 is zero
    ### Saving the samples after thinning the samples at every thin iterations .  
    if((i%%thin)-1==0){
      # print(cur.samples.log.hyper.mu)
      samples[j,]<- c(cur.samples.log.hyper.mu,
                      cur.samples.hyper.GP,
                      cur.samples.kappa_w1, 
                      cur.samples.kappa_w2,
                      cur.samples.kappa_eta,
                      cur.samples.kappa_mu,
                      cur.samples.intercept1,
                      cur.samples.intercept2,
                      cur.samples.beta,
                      cur.samples.beta1,
                      cur.samples.beta2,
                      cur.samples.eta[1:2],
                      cur.samples.w1[1], 
                      cur.samples.mu[1],
                      cur.samples.w2[1])
      
      j=j+1
    }
    ### storing the adaptive tuning parameters
    if((i%%adapt)-1==0){ # to save allexp the scale parameter of the MALA
      sigma.matrix[m,]<- c(tun.eta[1], tun.mu[1], tun.hyper.mu, tun.hyper.GP)
      m=m+1
    }
  }
# browser()
  
  ######### saving results for qqplots ################
  # WS
  if(CV=="WS"){
    true_sorted.Y <- sort(Y, decreasing = FALSE)
    estimated_samples_sorted.Y <- t(apply(imputed.Y.WSD.samples, 1, sort, decreasing = FALSE))
    estimated_mean.Y <- apply(estimated_samples_sorted.Y, 2, FUN = mean, na.rm=TRUE)
    lower_ci.Y <- apply(estimated_samples_sorted.Y, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
    upper_ci.Y <- apply(estimated_samples_sorted.Y, 2, function(x) quantile(x, 0.975 ,na.rm=TRUE))
    
    true_sorted.A <- sort(A, decreasing = FALSE)
    estimated_samples_sorted.A <- t(apply(imputed.A.WSD.samples, 1, sort, decreasing = FALSE))
    estimated_mean.A <- apply(estimated_samples_sorted.A, MARGIN = 2, FUN = mean, na.rm=TRUE)
    lower_ci.A <- apply(estimated_samples_sorted.A, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
    upper_ci.A <- apply(estimated_samples_sorted.A, 2, function(x) quantile(x, 0.975, na.rm=TRUE))
    
    WS_qqplots<- list(true.Y = true_sorted.Y, est.Y=estimated_mean.Y, lci.Y=lower_ci.Y, uci.Y=upper_ci.Y, 
                      true.A=true_sorted.A, est.A=estimated_mean.A, lci.A=lower_ci.A, uci.A = upper_ci.A)
    
  } else{
    WS_qqplots <- NULL 
  }
  
  # OOS
  if(CV=="OOS"){
    true.Y <- Y[ind_NA_Y]
    true_sorted.Y<- sort(true.Y, decreasing = FALSE)
    
    estimated_samples_sorted.Y <- t(apply(imputed.Y.OSD.samples, 1, sort, decreasing = FALSE))
    estimated_mean.Y <- apply(estimated_samples_sorted.Y, 2, FUN = mean, na.rm=TRUE)
    lower_ci.Y <- apply(estimated_samples_sorted.Y, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
    upper_ci.Y <- apply(estimated_samples_sorted.Y, 2, function(x) quantile(x, 0.975, na.rm=TRUE))
    
    true.A <- A[ind_NA_A]
    true_sorted.A<- sort(true.A, decreasing = FALSE)
    estimated_samples_sorted.A <- t(apply(imputed.A.OSD.samples, 1, sort, decreasing = FALSE))
    estimated_mean.A <- apply(estimated_samples_sorted.A, 2, FUN = mean, na.rm=TRUE)
    lower_ci.A <- apply(estimated_samples_sorted.A, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
    upper_ci.A <- apply(estimated_samples_sorted.A, 2, function(x) quantile(x, 0.975, na.rm=TRUE))
    
    OOS_qqplots<- list(true.Y = true_sorted.Y, est.Y=estimated_mean.Y, lci.Y=lower_ci.Y, uci.Y=upper_ci.Y, 
                       true.A=true_sorted.A, est.A=estimated_mean.A, lci.A=lower_ci.A, uci.A = upper_ci.A)
    
  } else{
    OOS_qqplots <- NULL 
  }
  
  
  ########## estimates and their standard errors ##########
  #### WS
  if(CV=="WS"){
    WS_with_CIs = list(post.mean.Y =  apply(imputed.Y.WSD.samples, MARGIN = 2, FUN = mean, na.rm=TRUE),
                       post.sd.Y= apply(imputed.Y.WSD.samples, MARGIN = 2, FUN = sd, na.rm=TRUE),
                       post.lci.Y= apply(imputed.Y.WSD.samples, MARGIN = 2, FUN = quantile, probs=0.025, na.rm=TRUE),
                       post.uci.Y=apply(imputed.Y.WSD.samples, MARGIN = 2, FUN = quantile, probs=0.975, na.rm=TRUE),
                       post.mean.A= apply(imputed.A.WSD.samples, MARGIN = 2, FUN = mean, na.rm=TRUE),
                       post.sd.A= apply(imputed.A.WSD.samples, MARGIN = 2, FUN = sd, na.rm=TRUE),
                       post.lci.A= apply(imputed.A.WSD.samples, MARGIN = 2, FUN = quantile, probs=0.025, na.rm=TRUE),
                       post.uci.A=apply(imputed.A.WSD.samples, MARGIN = 2, FUN = quantile, probs=0.975, na.rm=TRUE)
    )
  } else {
    WS_with_CIs<- NULL  
  }
  ###  OOS
  if(CV=="OOS"){
    OOS_with_CIs = list(post.mean.Y =  apply(imputed.Y.OSD.samples, MARGIN = 2, FUN = mean, na.rm=TRUE),
                        post.sd.Y= apply(imputed.Y.OSD.samples, MARGIN = 2, FUN = sd, na.rm=TRUE),
                        post.lci.Y= apply(imputed.Y.OSD.samples, MARGIN = 2, FUN = quantile, probs=0.025, na.rm=TRUE),
                        post.uci.Y=apply(imputed.Y.OSD.samples, MARGIN = 2, FUN = quantile, probs=0.975, na.rm=TRUE),
                        post.mean.A= apply(imputed.A.OSD.samples, MARGIN = 2, FUN = mean, na.rm=TRUE),
                        post.sd.A= apply(imputed.A.OSD.samples, MARGIN = 2, FUN = sd, na.rm=TRUE),
                        post.lci.A= apply(imputed.A.OSD.samples, MARGIN = 2, FUN = quantile, probs=0.025, na.rm=TRUE),
                        post.uci.A=apply(imputed.A.OSD.samples, MARGIN = 2, FUN = quantile, probs=0.975, na.rm=TRUE)
    )
  } else {
    OOS_with_CIs<- NULL  
  }
  
  
  
  return(list("samples"=samples, ### saving the sample for the hyperparameters to see the traceplots
              "OOS_with_CIs" =  OOS_with_CIs,
              "WS_with_CIs" =  WS_with_CIs,
              "OOS_qqplots" =  OOS_qqplots,
              "WS_qqplots" = WS_qqplots,
              "imputed.Y.WSD"=imputed.Y.WSD, 
              "imputed.A.WSD"=imputed.A.WSD,
              "imputed.Y.OSD"=imputed.Y.OSD, 
              "imputed.A.OSD"=imputed.A.OSD,
              "post.sum.mean.mu" = post.sum.mean.mu,  
              "post.sum.squre.mu" = post.sum.squre.mu, 
              "post.sum.mean.eta" = post.sum.mean.eta,
              "post.sum.squre.eta"= post.sum.squre.eta,
              "post.mean.condprob"= post.mean.condprob,
              "post.mean.uncondprob" = post.mean.uncondprob,
              "post.squre.condprob" = post.squre.condprob,
              "post.squre.uncondprob"=  post.squre.uncondprob,
              "post.sum.mean.w1"= post.sum.mean.w1, 
              "post.sum.squre.w1" = post.sum.squre.w1,
              "post.sum.mean.w2" =post.sum.mean.w2,
              "post.sum.squre.w2"= post.sum.squre.w2,
              "tuning_param_x_hyper"=sigma.matrix, ## tuning parameters values
              "Acc.rate eta"=rate.eta/(N.MCMC-(burn_in1+burn_in2)))) ## acceptance rate for eta
}

