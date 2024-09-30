

####################################################
########## for threshold model #########################
####################################################


#' prior distributions for the mark distributions
#'
#' @param cur_par 
#' @param mark_dist 
#' @param hyper.mu_fixed 
#'
#' @return
#' @export
#'
#' @examples
gamma_prior_markdist_thr<- function(cur_par,
                                thr.family,
                                hyper.mu_fixed){
  
  if(thr.family=="gamma"){
    k<- cur_par[1]
    logpriordens <- dgamma(k, shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE) 
  }
  if(thr.family=="logNormal"){
    k<- cur_par[1]
    logpriordens <- dgamma(k, shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE) 
  }
  
  return(logpriordens)
}






#####################################################
########## for joint model #########################
####################################################

#' prior distributions for the mark distributions
#'
#' @param cur_par 
#' @param mark_dist 
#' @param hyper.mu_fixed 
#'
#' @return
#' @export
#'
#' @examples
gamma_prior_markdist_JM<- function(cur_par,
                                mark_dist,
                                hyper.mu_fixed){
  
  k<- cur_par[1]
  xi<- cur_par[2]
  if(mark_dist=="eGPD"){ 
    logpriordens <- dgamma(k, shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE) +
      dgamma(xi, shape = 1, rate =15, log = TRUE)
  } else if (mark_dist=="bGPD"){
    logpriordens<- dgamma(k, shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)  
    
  } else if (mark_dist=="tgGPD"){
    logpriordens<- dgamma(k, shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)
  }
  
  return(logpriordens)
}




#' Title
#'
#' @param cur_par 
#' @param mark_dist 
#' @param hyper.mu_fixed 
#'
#' @return
#' @export
#'
#' @examples
gamma_prior_GPD_param<- function(cur_par,
                                 hyper.mu_fixed){
  sigma.GP<- cur_par[1]
  xi<- cur_par[2]
  logpriordens<- dgamma(xi, shape = 1, rate =15, log = TRUE) + 
    dgamma(sigma.GP, shape = hyper.mu_fixed[1], rate = hyper.mu_fixed[2], log = TRUE)
  return(logpriordens)
}







#' Title
#'
#' @param mark_dist 
#'
#' @return
#' @export
#'
#' @examples
lb_ub_markdist<- function(mark_dist){
  if(mark_dist=="eGPD"){ lb=c(0,0); ub=c(Inf, 1)
  } else if (mark_dist=="bGPD"){lb=c(0); ub=c(Inf)
  } else if (mark_dist=="tgGPD"){lb=c(0); ub=c(Inf)
  }
  return(list("ub"=ub, "lb"=lb))
}





#' Title
#'
#' @param mark_dist 
#'
#' @return
#' @export
#'
#' @examples
lb_ub_GP<- function(mark_dist){
  if (mark_dist=="bGPD"){lb=c(0,0); ub=c(Inf, 1)
  } else if (mark_dist=="tgGPD"){lb=c(0,0); ub=c(Inf, 1)
  }
  return(list("ub"=ub, "lb"=lb))
}
