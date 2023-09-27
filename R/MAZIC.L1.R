#' Lasso estimates for MAZIC
#'
#' @description The function returns the lasso estimates of each covariate
#'
#' @param The dataframe, the formula of the model, option for using starting values, optimization method,
#' names of the distribution, and verbose option
#'
#' @return The lasso estimates
#' @export

MAZIC.L1 = function(data, formula, startVals = TRUE, opt.method = c('Nelder-Mead', 'BFGS'),
                    mod = c("marzip","marzinb","marzigp"), verbose = TRUE)
{
  
  #### Obtain response variable and design matrix ####
  
  regData = orgMetadata(data = data, formula = formula)
  y = regData$y
  covs = regData$design
  if(verbose){
    message('Snapshot of the Design Matrix of this Data')
    print(head(covs))
  }
  
  #### Getting Start values ####
  starts = getStartVals(y = y, 
                        covs = covs, 
                        startVals = startVals, 
                        mar.mod = mod)
  
  estimates = tryCatch({MAZIC.MLE(y = y, x = covs, startVals = starts, opt.method = opt.method,
                                  mod = mod, verbose = TRUE)},
                       error = function(e){
                         return(NULL)
                       })
  if(is.null(estimates)){
    estimates = MAZIC.ridge(y = y, x = covs, startVals = starts, opt.method = opt.method,
                            mod = mod, lambda = 0.001, verbose = TRUE)
  }
  
  beta = estimates$beta.estimate
  vcov = estimates$vcov
  #### Removing Theta for Downstream Calculations ####
  
  nX = (2 * ncol(data[, covnames]) + 2)
  theta.ind = nX + 1
  vcov = vcov[-theta.ind, -theta.ind, drop = FALSE]
  
  #### Removing Intercepts for Downstream Calculations ####
  
  beta = beta[!names(beta) %in% c("zero.(Intercept)", "mean.(Intercept)")]
  p = length(covnames)
  sigma.11 = vcov[c(1,p+2),c(1,p+2)]
  sigma.12 = vcov[c(1,p+2),-c(1,p+2)]
  sigma.22 = vcov[-c(1,p+2),-c(1,p+2)]
  Sigma = try(sigma.22-t(sigma.12)%*%solve(sigma.11)%*%sigma.12,silent=T)
  #### Converting Sigma to PD in case it not PD ####
  
  e = eigen(Sigma)
  if(round(det(Sigma),4)>0) 
  {
    cov.star=e$vectors%*%diag(1/sqrt(abs(e$values)))%*%t(e$vectors)
  } else {
    cov.star=e$vectors%*%diag(1/sqrt(makePD(Sigma)))%*%t(e$vectors)
  }
  cov.star = as.matrix(sqrtm(cov.star))
  cov.star = scale(cov.star)
  y.star = as.vector(cov.star%*%beta)
  y.star = y.star - mean(y.star)
  colnames(cov.star) = colnames(Sigma)
  names(y.star) = colnames(Sigma)
  
  if(verbose){
    message('Snapshot of Y.star for LSA')
    print(y.star)
    message('Snapshot of scaled X.star for LSA')
    print(cov.star[1:5, 1:5])
  }
  
  #### Fitting Penalty Model ####
  
  fitlass = cv.glmnet(cov.star, y.star)
  lam = fitlass$lambda.min
  print(lam)
  fit = glmnet(cov.star, y.star, lambda = lam)
  #### Collecting Penalized Results ####
  
  coef = as.matrix(coef(fit, s = lam))
  colnames(coef) = 'Penalized.Estimates'
  print(coef)
  tLL = fit$nulldev - deviance(fit)
  k = fit$df
  n = fit$nobs
  aic = -2*tLL+2*k
  bic = log(n)*k - tLL
  
  # tLL = -deviance(fit) # 2*log-likelihood
  # k = dim(model.matrix(fit))[2]
  # n = nobs(fit)
  # aicc = -tLL+2*k+2*k*(k+1)/(n-k-1)
  # aic = -tLL+2*k
  
  # bic = log(n)*k - tLL
  print('Done BIC')
  #### Giving Back Intercept Terms ####
  
  coef = as.data.frame(coef)
  coef = coef[!rownames(coef) %in% '(Intercept)', , drop = FALSE]
  tmp = sigma.12 %*% solve(sigma.22) %*% (beta - as.vector(coef$Penalized.Estimates))
  coef.zero = coef[grep('zero', rownames(coef)), , drop = FALSE]
  coef.mean = coef[grep('mean', rownames(coef)), , drop = FALSE]
  Penalized.Coeff.zero = c(tmp[1, 1], coef.zero$Penalized.Estimates)
  names(Penalized.Coeff.zero)[2:length(Penalized.Coeff.zero)] = rownames(coef.zero)
  Penalized.Coeff.mean = c(tmp[2, 1], coef.mean$Penalized.Estimates)
  names(Penalized.Coeff.mean)[2:length(Penalized.Coeff.mean)] = rownames(coef.mean)
  All.Penalized.Coeffs = list(zero = Penalized.Coeff.zero, mean = Penalized.Coeff.mean)

  print('Done Calculation')
  if(verbose){
    message('Snapshot of Penalized Coefficients (using glmnet) for variables in zero part')
    print(Penalized.Coeff.zero)
    message('Snapshot of Penalized Coefficients (using glmnet) for variables in mean part')
    print(Penalized.Coeff.mean)
    message('Log-Likelihood of the model')
    print(tLL)
    message('AIC of the model')
    print(aic)
    message('BIC of the model')
    print(bic)
  }
  
  #### Giving Back Intercept Terms ####
  
  fintab = list(Penalized.Coeff.estimates = All.Penalized.Coeffs, aic = aic, bic = bic, 
                loglik = tLL)
  return(fintab)
}
