#' cdf of truncated Gamma distribution function
#'
#' @param x 
#' @param range 
#' @param shape 
#' @param scale 
#' @param log.p 
#'
#' @return
#' @export
#'
#' @examples
pgammat<- function(x, upper.bound, shape, scale, log.p=TRUE){
  F.a <- 0
  F.b <- pgamma(upper.bound, shape = shape, scale = scale)
  if(log.p){
    cdf<- pgamma(x, shape = shape, scale = scale, log.p = TRUE)-log(F.b - F.a)
  } else{
    cdf<- pgamma(x, shape = shape, scale = scale)/(F.b - F.a)
  }
  return(cdf)
}



#' density of truncated gamma distribution 
#'
#' @param x 
#' @param range 
#' @param shape 
#' @param scale 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
dgammat<- function(x, upper.bound, shape, scale, log=TRUE){
  F.a <- 0
  F.b <- pgamma(upper.bound, shape = shape, scale = scale)
  if(log){density<- dgamma(x, shape = shape, scale = scale, log = TRUE) - log((F.b - F.a))
  } else{
    density<- dgamma(x, shape = shape, scale = scale) / (F.b - F.a)
  }
  return(density)
}



#' simulation from truncated Gamma distribution function  
#'
#' @param n 
#' @param range 
#' @param shape 
#' @param scale 
#'
#' @return
#' @export
#'
#' @examples
rgammat <- function(n, upper.bound, shape, scale = 1) {
  F.a <- 0
  F.b <- pgamma(upper.bound, shape = shape, scale = scale)
  u <- runif(n, min = F.a, max = F.b)
  x<-qgamma(u, shape = shape, scale = scale)
  return(x)
  
}



#' density function of beta distribution with support in (0, c)
#'
#' @param x 
#' @param u.thr
#' @param shape1 
#' @param shape2 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
dbetat<- function(x, u.thr, shape1, shape2, log=TRUE){
  # browser()
  if(log){
    dens<- dbeta(x/u.thr, shape1 = shape1, shape2 = shape2, log=TRUE) - log(u.thr)
  } else {
    dens<- dbeta(x/u.thr, shape1 = shape1, shape2 = shape2)/u.thr
  }
  return(dens)
}


#' distribution function of beta distribution with support in (0, c)
#'
#' @param x 
#' @param c 
#' @param shape1 
#' @param shape2 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
pbetat<- function(x, u.thr, shape1, shape2, log.p=TRUE){
  if(log.p){
    dens<- pbeta(x/u.thr, shape1 = shape1, shape2 = shape2, log.p=TRUE)
  } else {
    dens<- pbeta(x/u.thr, shape1 = shape1, shape2 = shape2)
  }
  return(dens)
}


#' simulate from the beta distribution with support in (0, c)
#'
#' @param n 
#' @param c 
#' @param shape1 
#' @param shape2 
#'
#' @return
#' @export
#'
#' @examples
rbetat<- function(n, u.thr, shape1, shape2){
  x<- u.thr * rbeta(n=n, shape1 = shape1, shape2 = shape2)
  return(x)
}


#' density function of extended GPD of type 1
#'
#' @param x 
#' @param k 
#' @param xi 
#' @param sigma 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
dEGPD1<-function(x, k, xi, sigma, log=TRUE){
 # browser()
  if(log){
    dens<- log(k) + (evd::dgpd(x=x, loc=0, scale=sigma, shape=xi, log = TRUE)) + 
      ((k-1)* log(evd::pgpd(q=x, loc=0, scale=sigma, shape=xi, lower.tail=TRUE)))
  } else{
    dens<- k * (evd::dgpd(x=x, loc=0, scale=sigma, shape=xi)) * (evd::pgpd(q=x, loc=0, scale=sigma, shape=xi, lower.tail=TRUE))^(k-1)
  }
  return(dens)
}


#' distribution function of extended GPD of type 1
#'
#' @param x 
#' @param k 
#' @param xi 
#' @param sigma 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
pEGPD1<-function(x, k, xi, sigma, log=TRUE){
  if(log){
    cdf<- (k)* log(evd::pgpd(q=x, loc=0, scale=sigma, shape=xi, lower.tail=TRUE))
  } else{
    cdf<-  (evd::pgpd(q=x, loc=0, scale=sigma, shape=xi, lower.tail=TRUE))^(k) 
  }
  return(cdf)
}


### simulating random number from extended GPD of type 1
#' Title
#'
#' @param n 
#' @param k 
#' @param xi 
#' @param sigma 
#'
#' @return
#' @export
#'
#' @examples
rEGPD1<-function(n, k, xi, sigma){
  X<- (sigma/xi) * (((1-(runif(n=n))^(1/k))^(-xi))-1)
  return(X)
}




#' density function of mixture of truncated gamma and GPD distribution 
#'
#' @param x 
#' @param range 
#' @param k 
#' @param xi 
#' @param sigma 
#' @param thr.prob 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
dtgGPD<- function(x, u.thr, k, xi, sigma, sigma.GP, thr.prob, log=TRUE){
  if(log){
    dens<- rep(0, length(x))
    ind<- (x < u.thr) & (x >0)
    dens[ind]<- log(1 - thr.prob) + dgammat(x = x[ind], range = c(0,u.thr), shape = k, scale = sigma[ind], log=TRUE)
    ind.e<- x>u.thr
    dens[ind.e]<- log(thr.prob) + evd::dgpd(x=x[ind.e], loc= u.thr, scale = sigma.GP, shape=xi, log=TRUE)
  } else{
    dens<- rep(0, length(x))
    ind<- (x < u.thr) & (x >0)
    dens[ind]<- log(1 - thr.prob) * dgammat(x = x[ind], range = c(0,u.thr), shape = k, scale = sigma[ind])
    ind.e<- x>u.thr
    dens[ind.e]<- (thr.prob) * evd::dgpd(x=x[ind.e], loc= u.thr, scale = sigma.GP, shape=xi)
  }
  return(dens)
}


#' distribution function of mixture of truncated gamma and GPD distribution 
#'
#' @param x 
#' @param u.thr 
#' @param k 
#' @param xi 
#' @param sigma 
#' @param thr.prob 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
ptgGPD<- function(x, u.thr, k, xi, sigma, sigma.GP, thr.prob, log.p=TRUE){
  if(log.p){
    cdf<- rep(0, length(x))
    ind<- (x < u.thr) & (x >0)
    cdf[ind]<- log(1 - thr.prob) + pgammat(x = x[ind], range = c(0,u.thr), shape = k, scale = sigma[ind], log.p=TRUE)
    ind.e<- x>u.thr
    cdf[ind.e]<- log((1 - thr.prob) + thr.prob  * evd::pgpd(x=x[ind.e], loc= u.thr, shape=xi, scale = sigma.GP))
  } else{
    cdf<- rep(0, length(x))
    ind<- (x < u.thr) & (x >0)
    cdf[ind]<- (1 - thr.prob) * pgammat(x = x[ind], range = c(0,u.thr), shape = k, scale = sigma[ind])
    ind.e<- x>u.thr
    cdf[ind.e]<- (1 - thr.prob)  + (thr.prob) * evd::pgpd(x=x[ind.e], loc= u.thr, scale = sigma.GP, shape=xi)
  }
  return(cdf)
}


#' Title
#'
#' @param x 
#' @param mu 
#' @param u.thr 
#' @param k 
#' @param xi 
#' @param sigma 
#' @param thr.prob 
#'
#' @return
#' @export
#'
#' @examples
cond.ptgGPD<- function(x, mu, u.thr, k, xi, sigma, sigma.GP, thr.prob){
  #browser()
  cdf<- rep(NA, length(mu))
  if(x>u.thr){
    ind<- exp(as.numeric(mu)) < u.thr
    cdf[ind]<- 1
    cdf[!ind]<-  1 - thr.prob + thr.prob * evd::pgpd(x=x, loc= u.thr, scale = sigma.GP, shape=xi)
  } else{
    ind<- exp(as.numeric(mu)) > u.thr
    cdf[ind]<- 0 
    cdf[!ind]<-  (1 - thr.prob) * pgammat(x = x, range = c(0,u.thr), shape = k, scale = sigma[!ind]) 
  }
  return(cdf)
  
}


#' simulating from mixture of truncated gamma and GPD distribution 
#'
#' @param u.thr 
#' @param k 
#' @param xi 
#' @param n 
#' @param ind 
#' @param sigma 
#'
#' @return
#' @export
#'
#' @examples
rtgGPD<- function(n, u.thr, k, xi, sigma, sigma.GP, ind){
  n1<- sum(ind)
  x<- ifelse(ind, evd::rgpd(n=n1, loc=u.thr, scale = sigma.GP, shape = xi), 
             rgammat(n=n-n1, range = c(0, u.thr), shape = k, scale = sigma))
  return(x)
}


#' Density function of mixture of beta in (0, u.thr)  and GPD
#' 
#' @param x a vector of observations  
#' @param u.thr 
#' @param k  shape1 parameter for the  beta distribution 
#' @param xi shape parameter of GP 
#' @param sigma a vectors of scale parameters with same length as the length of x 
#' @param thr.prob 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
dbGPD<- function(x, u.thr, k, xi, sigma, sigma.GP, thr.prob, log=TRUE){
  #browser()
  if(log){
    dens<- rep(0, length(x))
    ind<- (x < u.thr) & (x >0)
    dens[ind]<- log(1 - thr.prob) + dbetat(x = x[ind], u.thr =u.thr, shape1 = k, shape2 = sigma[ind], log=TRUE)
    ind.e<- x>u.thr
    dens[ind.e]<- log(thr.prob) + evd::dgpd(x=x[ind.e], loc= u.thr, scale = sigma.GP, shape=xi, log=TRUE)
  } else{
    dens<- rep(0, length(x))
    ind<- (x < u.thr) & (x >0)
    dens[ind]<- (1 - thr.prob) * dbetat(x = x[ind], u.thr =u.thr, shape1 = k, shape2 = sigma[ind])
    ind.e<- x>u.thr
    dens[ind.e]<- (thr.prob) * evd::dgpd(x=x[ind.e], loc= u.thr, scale = sigma.GP, shape=xi)
  }
  return(dens)
  
}

#' distribution function of mixture of beta in (0, u.thr)  and GPD
#' 
#' @param x 
#' @param u.thr 
#' @param k 
#' @param xi 
#' @param sigma 
#' @param thr.prob 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
pbGPD<- function(x, u.thr, k, xi, sigma, sigma.GP, thr.prob, log.p=TRUE){
  if(log.p){
    cdf<- rep(0, length(x))
    ind<- (x < u.thr) & (x >0)
    cdf[ind]<- log(1 - thr.prob) + pbetat(x = x[ind], u.thr =u.thr, shape1 = k, shape2 = sigma[ind], log.p=TRUE)
    ind.e<- x>u.thr
    cdf[ind.e]<- log((1 - thr.prob) + thr.prob * evd::pgpd(x=x[ind.e], loc= u.thr, shape=xi, scale = sigma.GP, lower.tail = TRUE))
  } else{
    cdf<- rep(0, length(x))
    ind<- (x < u.thr) & (x >0)
    cdf[ind]<- (1 - thr.prob) * dbetat(x = x[ind], u.thr =u.thr, shape1 = k, shape2 = sigma[ind])
    ind.e<- x>u.thr
    cdf[ind.e]<-  (1 - thr.prob) + (thr.prob) * evd::pgpd(x=x[ind.e], loc= u.thr, shape=xi, scale = sigma.GP)
  }
  return(cdf)
  
}


#' Title
#'
#' @param x  x is high quantile, greater than u.thr
#' @param mu 
#' @param u.thr 
#' @param k 
#' @param xi 
#' @param sigma 
#' @param thr.prob 
#' @param log.p 
#'
#' @return
#' @export
#'
#' @examples
cond.pbGPD<- function(x, mu, u.thr, k, xi, sigma, sigma.GP, thr.prob){
  #browser()
  cdf<- rep(NA, length(mu))
  if(x>u.thr){
    ind<- exp(as.numeric(mu)) < u.thr
    cdf[ind]<- 1
    cdf[!ind]<-  1 - thr.prob + thr.prob * evd::pgpd(q=x, loc= u.thr, shape=xi, scale =sigma.GP , lower.tail = TRUE)
  } else{
    ind<- exp(as.numeric(mu)) > u.thr
    cdf[ind]<- 0 
    cdf[!ind]<-  (1 - thr.prob) * pbetat(x = x, u.thr =u.thr, shape1 = k, shape2 = sigma[!ind], log.p=TRUE) 
  }
  
  
  cdf<- ifelse(x>u.thr, ifelse(exp(as.numeric(mu)) < u.thr, 1, 1 - thr.prob + thr.prob * evd::pgpd(q=x, loc= u.thr, shape=xi, scale =sigma.GP , lower.tail = TRUE)),
               ifelse(exp(as.numeric(mu)) > u.thr, 0, (1 - thr.prob) * pbetat(x = x, u.thr =u.thr, shape1 = k, shape2 = sigma[!ind], log.p=TRUE) ))
  
  return(cdf)
}



#' simulating from mixture of beta in (0, u.thr)  and GPD
#'
#' @param u.thr 
#' @param k 
#' @param xi 
#' @param n 
#' @param ind 
#' @param sigma 
#'
#' @return
#' @export
#'
#' @examples
rbGPD<- function(n, u.thr, k, xi, sigma, sigma.GP, ind){
  n1<- sum(ind)
  x<- ifelse(ind, evd::rgpd(n=n1, loc=u.thr, scale = sigma.GP, shape = xi), 
             rbetat(n=n - n1, u.thr = u.thr, shape1 = k, shape2 = sigma))
  return(x)
}


#' Univariate truncated normal sampler
#'
#' This implementations works for a < 37
#' or b > 37
#' @param n [integer] sample size
#' @param mean [numeric] scalar location parameter
#' @param sd [numeric] scalar scale parameter
#' @param a [numeric] lower bound
#' @param b [numeric] upper bound
#' @export
rtnorm <- function(n = 1,
                   mean = 0,
                   sd = 1,
                   a = -Inf,
                   b = Inf
){
  stopifnot(length(a) == 1L,
            length(b) == 1L,
            length(mean) == 1L,
            length(sd) == 1L,
            isTRUE(a < b),
            isTRUE(sd > 0))
  a_std <- (a - mean) / sd
  b_std <- (b - mean) / sd
  if(b_std < -37 | a_std > 37){
    warning("Interval requested is beyond numerical tolerance.\nUse \"TruncatedNormal\" package \"rtnorm\" instead for rare events simulation.")
    return(TruncatedNormal::rtnorm(n = n, mu = mean, sd = sd, lb = a, ub = b))
  }
  if(a_std < 0){
    Fa <- pnorm(a_std)
    Fb <- pnorm(b_std)
    mean + sd * qnorm(Fa + runif(n) * (Fb - Fa))
  } else {
    Fa <- pnorm(a_std, lower.tail = FALSE)
    Fb <- pnorm(b_std, lower.tail = FALSE)
    mean + sd * (-qnorm(Fa - (Fa - Fb) * runif(n)))
  }
}

#' Density of a univariate truncated Normal
#' @inheritParams rtnorm
#' @param x [numeric] observation vector
#' @param log [logical] if \code{TRUE}, return the log density
#' @export
dtnorm <- function(x,
                   mean = 0,
                   sd = 1,
                   a = -Inf,
                   b = Inf,
                   log = FALSE){
  stopifnot(length(a) == 1L,
            length(b) == 1L,
            length(mean) == 1L,
            length(sd) == 1L,
            isTRUE(a < b),
            isTRUE(sd > 0))
  dens <- dnorm(x, mean = mean, sd = sd, log = TRUE)
  # if(is.finite(a) & is.finite(b)){
  dens <- dens  - TruncatedNormal::lnNpr(a = (a-mean)/sd, b = (b-mean)/sd)
  # log(
  #   pnorm(q = b, mean = mean, sd = sd) -
  #     pnorm(q = a, mean = mean, sd = sd))
  # } else if(is.infinite(a) & is.finite(b)){
  #   dens <- dens - pnorm(q = b, mean = mean, sd = sd, log.p = TRUE)
  # } else if(is.finite(a) & is.infinite(b)){
  #   dens <- dens - pnorm(q = a, mean = mean, sd = sd,
  #                        log.p = TRUE, lower.tail = FALSE)
  # }
  return(dens)
}

#' Transform parameter to unconstrained scale 
#'
#' @param par scalar parameter
#' @param lb lower bound
#' @param ub upper bound
#' @reference Section 56 of the Stan Reference manual version 2.9 at \url{https://github.com/stan-dev/stan/releases/download/v2.9.0/stan-reference-2.9.0.pdf}
#' @description 
transfo <- Vectorize(function(par, lb, ub){
  stopifnot(length(par) == 1L,
            length(lb) == 1L,
            length(ub) == 1L,
            isTRUE(lb < ub))
  if(lb == -Inf & ub == Inf){
    return(par)
  } else if(lb > -Inf & ub == Inf){
    return(log(par - lb))
  } else if(lb == -Inf & ub < Inf){
    return(log(ub - par))
  } else if(lb > -Inf & ub < Inf){
    return(qlogis((par - lb) / (ub - lb)))
  }
}, vectorize.args=c("par", "lb", "ub"))

#' Jacobean of the transformations made in \transfo
#'
#' @param tpar the transformed parameter
#' @param lb lower bound of parameter
#' @param ub upper bound of the parameter
#' @param log log=TRUE return log-Jacobean
#'
#' @return
#' @export
#'
#' @examples
jac_inv_transfo <- Vectorize(function(tpar, lb, ub, log = FALSE){
  if(lb == -Inf & ub == Inf){
    ljac <- 0  ## log(tpar) should be zero, why there is log(tpar) we need to check
  }
  if(lb > -Inf & ub == Inf){
    ljac <- tpar
  } else if(lb == -Inf & ub < Inf){
    ljac <- tpar
  } else if(lb > -Inf & ub < Inf){
    ljac <- log(ub - lb) + plogis(tpar, log.p = TRUE) + plogis(tpar, log.p = TRUE, lower.tail = FALSE)
  }
  if(log){
    return(ljac)
  } else{
    return(exp(ljac))
  }
}, vectorize.args=c("tpar", "lb", "ub"))


#' Gradient of the Jacobean 
#' 
#' @param tpar the transformed parameter
#' @param lb lower bound of parameter
#' @param ub upper bound of the parameter
#'
#' @return
#' @export
#'
#' @examples
dlogjac_inv_transfo <- Vectorize(function(tpar, lb, ub){
  if(lb == -Inf & ub == Inf){
    return(0)
  }
  if(lb > -Inf & ub == Inf){
    return(1)
  } else if(lb == -Inf & ub < Inf){
    return(1)
  } else{
    -1 + 2*plogis(-tpar)
  }
}, vectorize.args=c("tpar", "lb", "ub"))


#' Transform the parameters back to original scales
#'
#' @param tpar underlying parameter
#' @param lb  lower support of the parameters
#' @param ub upper support of the parameters
#'
#' @return transformed parameter

inv_transfo <- Vectorize(function(tpar, lb, ub){
  stopifnot(length(tpar) == 1L,
            length(lb) == 1L,
            length(ub) == 1L,
            isTRUE(lb < ub))
  if(lb == -Inf & ub == Inf){
    return(tpar)
  } else if(lb > -Inf & ub == Inf){
    return(exp(tpar) + lb)
  } else if(lb == -Inf & ub < Inf){
    return(ub - exp(tpar))
  } else{
    return(lb + (ub - lb)*plogis(tpar))
  }
}, vectorize.args=c("tpar", "lb", "ub"))



#' Truncated Eigen Square Root (or inverse square root) 
#'
#' Designed for Covaraice or Precision matricies - throws error if 
#' there is any large negative Eigenvalues. 
#'
#' @param Sigma p x p matrix to compute truncated matrix square root
#' @param inv (boolean) whether to compute inverse of square root
#' @return p x k matrix (where k is rank of Sigma) 
trunc_eigen_sqrt <- function(Sigma, inv){
  es <- eigen(Sigma)
  
  # Small negative eigenvalues can occur just due to numerical error and should
  # be set back to zero. 
  es$values[(es$values < 1e-12) & (es$values > -1e-12)] <- 0
  
  # If any eigen values are large negative throw error (something wrong)
  if (any(es$values < -1e-12)) stop("Non-trivial negative eigenvalues present")
  
  # calculate square root and reveal rank (k)
  k <- sum(es$values > 0)
  if (!inv){
    L <- es$vectors %*% diag(sqrt(es$values))  
  } else if (inv) {
    L <- es$vectors %*% diag(1/sqrt(es$values))  
  }
  return(L[,1:k, drop=F])
}

#' Covariance parameterization (Eigen decomposition)
#' @param n number of samples to draw
#' @param mu p-vector mean
#' @param Sigma covariance matrix (p x p)
#' @return matrix of dimension p x n of samples
rMVNormC_eigen <- function(n, mu, Sigma){
  p <- length(mu)
  L <- trunc_eigen_sqrt(Sigma, inv=FALSE)
  k <- ncol(L)
  Z <- matrix(rnorm(k*n), k, n)
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
}

#' Precision parameterization (Eigen decomposition)
#' @param n number of samples to draw
#' @param mu p-vector mean
#' @param Sigma precision matrix (p x p)
#' @return matrix of dimension p x n of samples
rMVNormP_eigen <- function(n, mu, Sigma){
  p <- length(mu)
  L <- trunc_eigen_sqrt(Sigma, inv=TRUE) # only difference vs. cov parameterization, this is done in terms of covariance matrix because eigen values of a matrix and its inverse is the same
  k <- ncol(L)
  Z <- matrix(rnorm(k*n), k, n)
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
}
# n <- 10000
# mu <- 1:3
# Omega <- MASS::ginv(Sigma) # use pseudoinverse 
# x1 <- rMVNormC_eigen(n, mu, Sigma)
# x2 <- rMVNormP_eigen(n, mu, Omega)
# 
# # Create function that tests for equality with high tolerance due to 
# # random number generation
# weak_equal <- function(x, y) all.equal(x, y, tolerance=.05)
# 
# # check row means match up with mu and agree with eachother
# weak_equal(rowMeans(x1), mu) & weak_equal(rowMeans(x2), mu)
# 
# # check empirical covariances agree
# weak_equal(var(t(x1)), var(t(x2)))



