#' @title PDF of Generalized Beta Prime distribution
#' @description This function computes the value of the
#' probability density function of a generalized
#' beta prime (GBP) distribution
#' @param x vector of quantiles
#' @param alpha shape parameter
#' @param beta shape parameter
#' @param c power parameter
#' @param d scale parameter
#' @return value of pdf at x
#' @export
dgbetapr <- function(x, alpha, beta, c, d){
  return(
    c*(x/d)^(alpha*c-1) * (1+(x/d)^c)^(-alpha-beta) / (d*beta(alpha,beta))
  )
}



#' @title CDF of Generalized Beta Prime distribution
#' @description This function computes the value of the
#' distribution function of a generalized
#' beta prime (GBP) distribution
#' @param x vector of quantiles
#' @param alpha shape parameter
#' @param beta shape parameter
#' @param c power parameter
#' @param d scale parameter
#' @return value of cdf at x
#' @export
pgbetapr <- function(x, alpha, beta, c, d){

  u = (x/d)^c
  v = u/(1+u)

  return(
    stats::pbeta(v,alpha,beta)
  )
}

#' @title Quantile Function of Generalized Beta Prime distribution
#' @description This function computes the quantile of a generalized
#' beta prime (GBP) distribution
#' @param p vector of probabilities
#' @param alpha shape parameter
#' @param beta shape parameter
#' @param c power parameter
#' @param d scale parameter
#' @return quantile corresponding to p
#' @export

qgbetapr <- function(p, alpha, beta, c, d){

  apply.fun <- function(p){
    dist <- function(logw){
      return(
        (r2d2glmm::pgbetapr(exp(logw), alpha, beta, c, d) - p)^2
      )
    }

    out <- stats::optimize(dist, lower=log(1/1000), upper=log(1000))$min
    return(exp(out))
  }

  return(sapply(p, apply.fun))

}


#' @title Random Number Generator of Generalized Beta Prime distribution
#' @description This function generates a random number from a generalized
#' beta prime (GBP) distribution
#' @param n number of observations
#' @param alpha shape parameter
#' @param beta shape parameter
#' @param c power parameter
#' @param d scale parameter
#' @return quantile corresponding to p
#' @export

rgbetapr <- function(n, alpha, beta, c, d){
  V <- stats::rbeta(n, alpha, beta)
  return(
    d*(V/(1-V))^(1/c)
  )
}

## Functions for using R2D2 with LOGISTIC REGRESSION

#' @title Expit
#' @description This function computes the logistic function
#' @param x value of x
#' @return inverse logit of x
#' @export

expit <- function(x){
  return(
    1/(1+exp(-x))
  )
}

#' @title Logit
#' @description This function computes the logistic function
#' @param x value of x
#' @return inverse logit of x
#' @export

logit <- function(x){
  return(
    log(x/(1-x))
  )
}


#' @title R2 for a given W
#' @description This function computes the value of R2 for
#' a given W, beta0 and link function
#' @param W value of W
#' @param b0 value of beta0
#' @param link link function used to compute R2. Default options are "Gaussian", "Logistic",
#'  "Poisson", "Poisson_offsets", "ZIP", "NegBinom" and "Weibull". Also can have user-inputted link
#'  function using "Arbitrary".
#' @param theta dispersion parameter (if necessary)
#' @param mu if link="Arbitrary", then this is the mean function of the linear predictor.
#' May depend on theta
#' @param sigma2 if link="Arbitrary", then this is the variance function of the linear predictor.
#' May depend on theta
#' @param disp logical. TRUE if arbitrary link function depends on a dispersion parameter
#' @return value of R2
#' @export

W_to_R2 <- function(W, b0=0, link=c("Gaussian", "Logistic", "Poisson", "Poisson_offsets",
                                  "ZIP", "NegBinom", "Weibull", "Arbitrary"), theta=1,
                    mu, sigma2, disp=FALSE){

  if(link=="Gaussian"){
    return(W/(W+theta))

  }else if(link=="Logistic"){
    K=1000
    i.seq <- seq(1,K-1,1)/K

    apply.fun <- function(x){
      M = sum((r2d2glmm::expit(stats::qnorm(i.seq, b0, sqrt(x))))^2)/(K-1)-(sum(r2d2glmm::expit(stats::qnorm(i.seq, b0, sqrt(x))))/(K-1))^2
      V = sum(r2d2glmm::expit(stats::qnorm(i.seq, b0, sqrt(x)))*(1-r2d2glmm::expit(stats::qnorm(i.seq, b0, sqrt(x)))))/(K-1)

      return(M/(M+V))
    }
    return(sapply(W, apply.fun))

  }else if(link=="Poisson"){
    R2 <- (exp(W)-1)/(exp(W)-1+exp(-b0-0.5*W))
    return(R2)

  }else if(link=="Poisson_offsets"){
    R2 <- (theta*exp(W)-1)/(theta*exp(W)-1+(theta)^(-1/2)*exp(-b0-0.5*W))
    R2min = (theta-1) / (theta-1+theta^(-1/2)*exp(-b0))
    return((R2-R2min)/(1-R2min))

  }else if(link=="ZIP"){
    R2 <- (1-theta)*(exp(W)-1)/((1-theta)*(exp(W)-1)+exp(-b0-0.5*W)+theta*exp(W)) /(1-theta)
    return(R2)

  }else if(link=="NegBinom"){
    R2 <- (exp(W)-1)/(exp(W)-1+theta*exp(-b0-0.5*W))
    return(R2)

  }else if (link=="Weibull"){
    
    r = gamma(1+2/theta)/gamma(1+1/theta)^2
    
    R2 <- (exp(W)-1)*(2+r) / (exp(W)*(2+r)-1) 
    
    
  }else if(link=="Arbitrary"){
    K=1000
    i.seq <- seq(1,K-1,1)/K

    if(disp){
      apply.fun <- function(x){
        input <- stats::qnorm(i.seq, b0, sqrt(x))

        M = sum((mu(input, theta=theta))^2)/(K-1)-(sum(mu(input, theta=theta))/(K-1))^2
        V = sum(sigma2(input, theta=theta))/(K-1)

        return(M/(M+V))
      }
      return(sapply(W, apply.fun))
    }else{
      apply.fun <- function(x){
        input <- stats::qnorm(i.seq, b0, sqrt(x))

        M = sum((mu(input))^2)/(K-1)-(sum(mu(input))/(K-1))^2
        V = sum(sigma2(input))/(K-1)

        return(M/(M+V))
      }
      return(sapply(W, apply.fun))
    }
  }

}

#' @title W for a given R2
#' @description This function computes the value of W for
#' a given R2, b0 and link function
#' @param R2 value of R2
#' @param b0 value of beta0
#' @param link link function used to compute R2. Default options are "Gaussian", "Logistic",
#'  "Poisson", "Poisson_offsets", "ZIP", "NegBinom" and "Weibull". Also can have user-inputted link
#'  function using "Arbitrary".
#' @param theta dispersion parameter (if necessary)
#' @param mu if link="Arbitrary", then this is the mean function of the linear predictor.
#' May depend on theta
#' @param sigma2 if link="Arbitrary", then this is the variance function of the linear predictor.
#' May depend on theta
#' @param disp logical. TRUE if arbitrary link function depends on a dispersion parameter
#' @return value of W
#' @export

R2_to_W <- function(R2, b0=0, link=c("Gaussian", "Logistic", "Poisson", "Poisson_offsets",
                                         "ZIP", "NegBinom", "Weibull", "Arbitrary"), theta=1,
                          mu, sigma2,disp=FALSE){

  apply.fun <- function(R2, bb){
    dist <- function(logw,R2){
      abs(r2d2glmm::W_to_R2(W=exp(logw), b0=bb, link=link,theta=theta,
                            mu=mu,sigma2=sigma2,disp=disp)-R2)
    }
    out <- stats::optimize(dist,
                    lower=log(1/100000),
                    upper=log(100000),
                    maximum=FALSE,
                    R2=R2)$min
    return(exp(out))}

  return(
    sapply(R2, apply.fun, bb=b0)
  )

}




#' @title PDF of W
#' @description This function computes the value of the
#' probability density function of W induced by a
#' Beta(a,b) prior on R2 for any link function
#' @param w vector of quantiles
#' @param a value of a for prior R2~Beta(a,b)
#' @param b value of b for prior R2~Beta(a,b)
#' @param b0 value of beta0
#' @param link link function used to compute R2. Default options are "Gaussian", "Logistic",
#' "Poisson", "Poisson_offsets", "ZIP", "NegBinom" and "Weibull". Also can have user-inputted link
#' function using "Arbitrary".
#' @param theta dispersion parameter (if necessary)
#' @param mu if link="Arbitrary", then this is the mean function of the linear predictor.
#' May depend on theta
#' @param sigma2 if link="Arbitrary", then this is the variance function of the linear predictor.
#' May depend on theta
#' @param disp logical. TRUE if arbitrary link function depends on a dispersion parameter
#' @return value of PDF at W
#' @export

dw <- function(w, a=1, b=1, b0=0, link=c("Gaussian", "Logistic", "Poisson", "Poisson_offsets",
                                               "ZIP", "NegBinom", "Weibull", "Arbitrary"),
               theta=1,mu, sigma2, disp=FALSE){

  if(link=="Gaussian"){
    return(r2d2glmm::dgbetapr(w,a,b,1,theta))

  }else if (link=="Logistic"){
    delta=0.001
    return(
      (r2d2glmm::pw(w=w+delta, a=a, b=b, b0=b0, link=link,theta=theta,mu=mu,sigma2=sigma2,disp=disp) -
         r2d2glmm::pw(w=w-delta, a=a, b=b, b0=b0, link=link,theta=theta,mu=mu,sigma2=sigma2,disp=disp)) / (2*delta)
      )

  }else if(link=="Poisson"){
    return(
      (exp(w)-1)^(a-1) * exp(-b*(b0+w/2)) * (3*exp(w)-1) / (beta(a,b) * 2 * (exp(w)-1+exp(-b0-w/2))^(a+b))
    )

  }else if(link=="Poisson_offsets"){
    return(
      (theta^(a/2) * exp(a*b0) * (1+exp(b0)*(theta-1)*sqrt(theta))^b *
         (1-theta+exp(w/2)*(theta*exp(w)-1))^(a-1)*(3*theta*exp(3*w/2)-exp(w/2))) /
        (
          2*beta(a,b)*(1+sqrt(theta)*exp(b0+w/2)*(theta*exp(w)-1))^(a+b)
        )
    )

  }else if(link=="ZIP"){
    return(
      (exp(w)-1)^(a-1) * exp(-b*(b0+w/2))*(1+theta*exp(b0+3*w/2))^b * (3*exp(w)-1+2*theta*exp(b0+3*w/2)) /
        (beta(a,b) * 2 * (exp(w)-1+exp(-b0-w/2)+theta)^(a+b) * (1+theta*exp(b0+3*w/2)))
    )

  }else if(link=="NegBinom"){
    return(
      theta^b*(exp(w)-1)^(a-1) * exp(-b*(b0+w/2)) * (3*exp(w)-1) / (beta(a,b) * 2 * (exp(w)-1+theta*exp(-b0-w/2))^(a+b))
    )

  }else if(link=="Weibull"){
    
    r = 1/(2+gamma(1+2/theta)/gamma(1+1/theta)^2)
    
    return(
      (1-r)^b*exp(w)*(exp(w)-1)^(a-1) / (beta(a,b)*(exp(w)-r)^(a+b))
    )
    
    
  }else if(link=="Arbitrary"){
    delta=0.001
    return(
      (r2d2glmm::pw(w=w+delta, a=a, b=b, b0=b0, link=link,theta=theta,mu=mu,sigma2=sigma2,disp=disp) -
         r2d2glmm::pw(w=w-delta, a=a, b=b, b0=b0, link=link,theta=theta,mu=mu,sigma2=sigma2,disp=disp)) / (2*delta)
    )
  }
}


#' @title CDF of W
#' @description This function computes the value of the
#' distribution function of W induced by a
#' Beta(a,b) prior on R2 for arbitrary link function
#' @param w vector of quantiles
#' @param a value of a for prior R2~Beta(a,b)
#' @param b value of b for prior R2~Beta(a,b)
#' @param b0 value of beta0
#' @param link link function used to compute R2. Default options are "Gaussian", "Logistic",
#' "Poisson", "Poisson_offsets", "ZIP", "NegBinom" and "Weibull". Also can have user-inputted link
#' function using "Arbitrary".
#' @param theta dispersion parameter (if necessary)
#' @param mu if link="Arbitrary", then this is the mean function of the linear predictor.
#' May depend on theta
#' @param sigma2 if link="Arbitrary", then this is the variance function of the linear predictor.
#' May depend on theta
#' @param disp logical. TRUE if arbitrary link function depends on a dispersion parameter
#' @return value of CDF at W
#' @export

pw <- function(w, a=1, b=1, b0=0, link=c("Gaussian", "Logistic", "Poisson", "Poisson_offsets",
                                         "ZIP", "NegBinom", "Weibull", "Arbitrary"),
               theta=1, mu, sigma2, disp=FALSE){

  R2 <- r2d2glmm::W_to_R2(W=w, b0=b0,link=link,theta=theta,mu=mu,sigma2=sigma2,disp=disp)

  return(
    stats::pbeta(R2, a, b)
  )
}



#' @title Quantile Function of W
#' @description This function computes quantile of W induced by a
#' Beta(a,b) prior on R2 for any link function
#' @param p vector of probabilities
#' @param a value of a for prior R2~Beta(a,b)
#' @param b value of b for prior R2~Beta(a,b)
#' @param b0 value of beta0
#' @param link link function used to compute R2. Default options are "Gaussian", "Logistic",
#' "Poisson", "Poisson_offsets", "ZIP", "NegBinom" and "Weibull". Also can have user-inputted link
#' function using "Arbitrary".
#' @param theta dispersion parameter (if necessary)
#' @param mu if link="Arbitrary", then this is the mean function of the linear predictor.
#' May depend on theta
#' @param sigma2 if link="Arbitrary", then this is the variance function of the linear predictor.
#' May depend on theta
#' @param disp logical. TRUE if arbitrary link function depends on a dispersion parameter
#' @return quantile of W at p
#' @export

qw <- function(p, a=1, b=1, b0=0, link=c("Gaussian", "Logistic", "Poisson", "Poisson_offsets",
                                         "ZIP", "NegBinom", "Weibull", "Arbitrary"),
               theta=1, mu, sigma2, disp=FALSE){

  apply.fun <- function(p){
    dist <- function(logw){
      return(
        (r2d2glmm::pw(w=exp(logw), a=a, b=b, b0=b0,link=link,
                      theta=theta,mu=mu,sigma2=sigma2,disp=disp) - p)^2
      )
    }

    out <- stats::optimize(dist, lower=log(1/100000), upper=log(100000))$min

    return(exp(out))
  }

  return(
    sapply(p, apply.fun)
  )

}



#' @title GBP Approximation
#' @description This function finds the closest generalized beta prime (GBP)
#' distribution to the true pdf of W as measured by the
#' Pearson chi-squared divergence
#' @param a value of a for prior R2~Beta(a,b)
#' @param b value of b for prior R2~Beta(a,b)
#' @param b0 value of beta0
#' @param link link function used to compute R2. Default options are "Gaussian", "Logistic",
#' "Poisson", "Poisson_offsets", "ZIP", "NegBinom" and "Weibull". Also can have user-inputted link
#' function using "Arbitrary".
#' @param theta dispersion parameter (if necessary)
#' @param mu if link="Arbitrary", then this is the mean function of the linear predictor.
#' May depend on theta
#' @param sigma2 if link="Arbitrary", then this is the variance function of the linear predictor.
#' May depend on theta
#' @param disp logical. TRUE if arbitrary link function depends on a dispersion parameter
#' @param lambda tuning parameter
#' @return parameters (alpha,beta,c,d) for GBP
#' @export

WGBP <- function(a=1, b=1, lambda=0.25, b0=0,
                 link=c("Logistic", "Poisson", "Poisson_offsets",
                                                           "ZIP", "NegBinom", "Weibull", "Arbitrary"),
                 theta=1, mu, sigma2, disp=FALSE){

  tau <- seq(0.01, 0.99, length=100)

  w <- r2d2glmm::qw(p=tau, a=a, b=b, b0=b0,link=link,theta=theta,mu=mu,sigma2=sigma2,disp=disp)
  py <- r2d2glmm::dw(w=w, a=a, b=b, b0=b0,link=link,theta=theta,mu=mu,sigma2=sigma2,disp=disp)
  
  w <- w[!is.na(py)]
  py <- py[!is.na(py)]

  dist <- function(logparam){# Distance between W and W'
    px <- r2d2glmm::dgbetapr(w, exp(logparam[1]),exp(logparam[2]),exp(logparam[3]), exp(logparam[4]))[!is.na(py)]

    return(sum((1-px/py)^2) + lambda *( (logparam[1]-log(a))^2 + (logparam[2]-log(b))^2 +
                                          (logparam[3]-log(1))^2 + (logparam[4]-log(1))^2  ))

  }

  return(
    exp(stats::optim(c(log(1), log(1),  log(1), log(1)), dist)$par)
  )
}





