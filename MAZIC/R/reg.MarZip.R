#' Log Likelihood for Marginalized ZIP
#'
#' @description The function returns the log likelihood of a Marginalized ZIP
#'
#' @param A vector of response, a matrix of covariates, names of the parameter, 
#' and tuning parameter for Ridge penalty term (default is 0) 
#'
#' @return Log Likelihood of Marginalized ZIP
#' @export

#### Regression LLK Function for MZIP ####

reg.marzip = function(y, X, par, lambda = 0) {
  
  #### Extract parameters ####
  
  logitpars = par[grep('zero', names(par))]   
  poispars  = par[grep('mean', names(par))]     
  
  #### Regressing Logit Link Part ####
  
  Xlogit  = X
  LPlogit = Xlogit %*% logitpars
  logi0   = plogis(LPlogit)  
  
  #### Regressing Log Link Part adjusted for Marginalized Model ####
  
  Xpois  = X
  mupois = exp(Xpois %*% poispars - log(1 - logi0)) ## - log(1 - logi0) is the adjustment ##
  
  #### Log Likelihoods (LLK) for Logit and Log Part ####
  
  logliklogit = log(logi0 + exp(log(1 - logi0) - mupois))
  loglikpois  = log(1 - logi0) + dpois(y, lambda = mupois, log = TRUE)
  
  #### Estimating Joint LLK for Zero and Count Parts ####
  
  y.zero = y == 0  
  y.count = y > 0
  loglik = sum(logliklogit[y.zero]) + sum(loglikpois[y.count])
  return(-loglik + lambda * (par %*% par))
}