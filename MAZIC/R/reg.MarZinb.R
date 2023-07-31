#' Log Likelihood for Marginalized ZINB
#'
#' @description The function returns the log likelihood of a Marginalized ZIP
#'
#' @param A vector of response, a matrix of covariates, names of the parameter, 
#' and tuning parameter for Ridge penalty term (default is 0) 
#'
#' @return Log Likelihood of Marginalized ZINB
#' @export

#### Regression LLK Function for MZINB ####

reg.marzinb = function(y, X, par, lambda = 0){
  
  #### Extract parameters ####
  
  logitpars = par[grep('zero', names(par))]   
  mupars = par[grep('mean', names(par))]
  theta = exp(par[grep('theta', names(par))])
  
  #### Regressing Logit Link Part ####
  
  Xlogit  = X
  LPlogit = Xlogit %*% logitpars
  logi0   = plogis(LPlogit)  
  
  #### Regressing Log Link Part adjusted for Marginalized Model ####
  
  Xnb  = X
  munb = exp(Xnb %*% mupars - log(1 - logi0)) ## - log(1 - logi0) is the adjustment ##
  
  #### Log Likelihoods (LLK) for Logit and Log Part ####
  
  d = (1/(1 + (theta * munb)))^(1/theta)
  logliklogit = log(logi0 + ((1 - logi0) *  d))
  logliknb  = log(1 - logi0) + dnbinom(y, mu = munb, size = theta, log = TRUE)
  
  #### Estimating Joint LLK for Zero and Count Parts ####
  
  y.zero = y == 0  
  y.count = y > 0
  loglik = sum(logliklogit[y.zero]) + sum(logliknb[y.count])
  return(-loglik + lambda * (par %*% par))
}