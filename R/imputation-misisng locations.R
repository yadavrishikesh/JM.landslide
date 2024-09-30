#' Imputing NA for the landslides counts
#'
#' @param ind_NA_Y 
#' @param eta 
#' @param CV 
#'
#' @return
#' @export
#'
#' @examples
impute.NA.Y<-function(ind_NA_Y, eta, CV){
  n1<-length(eta)
  if(CV=="OOS"){
    imputed_NA_Y<- rpois(sum(ind_NA_Y), lambda = exp(eta)[ind_NA_Y])
  } else if(CV=="WS"){
    imputed_NA_Y<- rpois(n1, lambda = exp(eta))
  }
  return(imputed_NA_Y)
}


#' Imputing NA for landslides sizes in case of with sample (WS)
#'
#' @param ind_zeros_counts 
#' @param mu 
#' @param cur_par 
#' @param mark_dist 
#' @param threshold 
#' @param thr.prob 
#' @param thr.exceed.ind 
#'
#' @return
#' @export
#'
#' @examples
impute.NA.A<-function(CV,
                      ind_NA_A,
                      ind_zeros_counts,
                      mu, 
                      thr.prob,
                      cur_par,  
                      mark_dist, 
                      threshold,
                      thr.exceed.ind=NULL
){
  #browser()
  n1<- length(mu)
  if (mark_dist=="eGPD"){
    k<- cur_par[1]
    xi<- cur_par[2]
    if(CV=="WS"){
      mu<- mu[!ind_zeros_counts]
      sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
      imputed_NA_A<- rEGPD1(n=length(mu), k=k, xi=xi, sigma=sigma)
    } else{
      mu<- mu[ind_NA_A]
      sigma<- exp(mu)/evd::qgpd(0.5^(1/k), scale = 1, shape = xi)
      imputed_NA_A<- rEGPD1(n=length(mu), k=k, xi=xi, sigma=sigma) 
    }
    
  } else if (mark_dist=="bGPD"){ ### I assume parsimonious gamma distribution whenever we have missing observations
    k<- cur_par[1]
    sigma.GP<- cur_par[2]
    xi<- cur_par[3]
    
    if(CV=="WS"){
      mu<- mu[!ind_zeros_counts]
      threshold<- threshold[!ind_zeros_counts]
      thr.exceed.ind<- exp(mu) > threshold
      sigma<- ifelse(thr.exceed.ind, sigma.GP, (threshold / exp(mu) - 1) * k)
      
      imputed_NA_A<-  ifelse(thr.exceed.ind, 
                             threshold + evd::rgpd(n=1, loc = 0, scale = sigma.GP, shape = xi), 
                             threshold * rbeta(n=1, shape1 = k, shape2 = sigma)
      ) 
    } else {
      mu<- mu[ind_NA_A]
      threshold<- threshold[ind_NA_A]
      thr.exceed.ind<- exp(mu) > threshold
      sigma<- ifelse(thr.exceed.ind, sigma.GP, (threshold / exp(mu) - 1) * k)
      
      imputed_NA_A<-  ifelse(thr.exceed.ind, 
                             threshold + evd::rgpd(n=1, loc = 0, scale = sigma.GP, shape = xi), 
                             threshold * rbeta(n=1, shape1 = k, shape2 = sigma))
    }
    
  } else if (family=="tgGPD"){
    k<- cur_par[1]
    sigma.GP<- cur_par[2]
    xi<- cur_par[3]
    
    
    if(CV=="WS"){
      mu<- mu[!ind_zeros_counts]
      threshold<- threshold[!ind_zeros_counts]
      thr.exceed.ind<- exp(mu) > threshold
      sigma<- ifelse(thr.exceed.ind,  sigma.GP, exp(mu)/k)
      
      imputed_NA_A<-  ifelse(thr.exceed.ind, 
                             threshold + evd::rgpd(n=1, loc = 0, scale = sigma.GP, shape = xi), 
                             rgammat(n=1, upper.bound  = threshold,  shape = k, scale = sigma))
      
    } else {
      mu<- mu[ind_NA_A]
      threshold<- threshold[ind_NA_A]
      thr.exceed.ind<- exp(mu) > threshold
      sigma<- ifelse(thr.exceed.ind,  sigma.GP, exp(mu)/k)
      
      imputed_NA_A<- ifelse(thr.exceed.ind, 
                            threshold + evd::rgpd(n=1, loc = 0, scale = sigma.GP, shape = xi), 
                            rgammat(n=1, upper.bound  = threshold,  shape = k, scale = sigma))
    }
  }
  return(list ("thr.exceed.ind"=thr.exceed.ind, 
               "imputed_NA_A"=imputed_NA_A))
}



