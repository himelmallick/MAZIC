
library(haven)
library(pscl)
library(zic)
library(MASS)
library(VGAM)
library(expm)
library(grpreg)
library(stringr)

#' PD approximation of a matrix.
#'
#' @description The function uses eigen values and eigen vectors to recreate a pd matrix.
#'
#' @param mat A matrix input.
#'
#' @return The function returns a modified matrix which is pd.

makePD = function(mat){
  N = nrow(mat)
  HC = mat
  D = eigen(mat)
  E = D$values
  U = D$vectors
  v = as.numeric(E < 0)
  idx = which(v==1)
  m = sum(v) # number of negative values
  if(m > 0){
    S = sum(v*E)*2
    W = (S*S*100)+1
    P = min(abs(E[idx])) # smallest positive value
    for(i in idx){
      C = E[i]
      E[i] = P * (S-C)*(S-C)/W
    }
  }
  return(E)
}

#### Organize Design Matrix ####

orgMetadata = function(data, formula){
  
  #### Extracting variables from formula ####
  
  vars = all.vars(formula)
  yvar = all.vars(formula)[1]
  desvars = all.vars(formula)[2:length(vars)]
  matchDot = match('.', desvars)
  matchDot = ifelse(is.na(matchDot), FALSE, TRUE)
  
  #### Getting response variable values ####
  
  y = as.numeric(as.matrix(data[, yvar]))
  if(sum(y < 0) > 0){
    stop('There are negative numbers in the response/outcome variable. The 
         response/outcome variable can only have positive integers')
  }
  
  #### Getting design matrix ####
  
  if(matchDot){
    desvars = colnames(data[, !colnames(data) %in% yvar])
    covs = data[, desvars]
    design = as.formula(paste('~', paste(desvars, collapse = "+"), sep=''))
    designMat = model.matrix(design, data = covs)
  }else if(!matchDot){
    covs = data[, desvars, drop = FALSE]
    design = as.formula(paste('~', paste(desvars, collapse = "+"), sep=''))
    designMat = model.matrix(design, data = covs)
  }else{
    stop('Please specify correct independent variable names in the formula argument')
  }
  orgData = list(y = y, design = designMat)
  return(orgData)
}

#### Get Start Values for Optim ####

getStartVals = function(y, covs, startVals, mar.mod){
  
  #### Poisson Starting Values ####
  
  if(mar.mod == 'marzip'){
    if(!startVals){
      starts = c(logit = rep(0, ncol(covs)), 
                 log = rep(0, ncol(covs)))
      names(starts) = c(paste0('zero.', colnames(covs)), 
                        paste0('mean.', colnames(covs)))
    }else{
      dmat = as.data.frame(cbind(y, covs))
      dmat = dmat[, -2]
      form = as.formula(paste('y', '~', paste(names(dmat)[2:ncol(dmat)], 
                                              collapse = "+"), sep = ''))
      fit.start = glm(formula = form, data = dmat, family = 'poisson')
      starts = c(logit = coef(fit.start), 
                 log = coef(fit.start))
      names(starts) = c(paste0('zero.', names(coef(fit.start))), 
                        paste0('mean.', names(coef(fit.start))))
    }
  }else if(mar.mod == 'marzinb'){
    if(!startVals){
      starts = c(logit = rep(0, ncol(covs)), 
                 log = rep(0, ncol(covs)),
                 theta = 0)
      names(starts) = c(paste0('zero.', colnames(covs)), 
                        paste0('mean.', colnames(covs)),
                        'theta')
    }else{
      dmat = as.data.frame(cbind(y, covs))
      dmat = dmat[, -2]
      form = as.formula(paste('y', '~', paste(names(dmat)[2:ncol(dmat)], 
                                              collapse = "+"), sep = ''))
      fit.start = glm.nb(formula = form, data = dmat)
      starts = c(logit = coef(fit.start), 
                 log = coef(fit.start),
                 theta = fit.start$theta)
      names(starts) = c(paste0('zero.', names(coef(fit.start))), 
                        paste0('mean.', names(coef(fit.start))),
                        'theta')
    }
  }else if(mar.mod == 'marzigp'){
    if(!startVals){
      starts = c(logit = rep(0, ncol(covs)), 
                 log = rep(0, ncol(covs)),
                 theta = 0)
      names(starts) = c(paste0('zero.', colnames(covs)), 
                        paste0('mean.', colnames(covs)),
                        'theta')
    }else{
      dmat = as.data.frame(cbind(y, covs))
      dmat = dmat[, -2]
      form = as.formula(paste('y', '~', paste(names(dmat)[2:ncol(dmat)], 
                                              collapse = "+"), sep = ''))
      fit.start = vglm(formula = form, genpoisson2, data = dmat, trace = FALSE)
      fit.start = coef(fit.start, matrix = TRUE)
      starts = c(logit = fit.start[, 1], 
                 log = fit.start[, 1],
                 theta = fit.start[, 2][1])
      names(starts) = c(paste0('zero.', names(fit.start[, 1])), 
                        paste0('mean.', names(fit.start[, 1])),
                        'theta')
    }
  }else{
    stop("The specified marginalized model is not supported")
  }
  return(starts)
}

#### Regression LLK Function for MZIP ####

reg.marzigp = function(y, X, par, lambda = 0) {
  
  #### Extract parameters ####
  
  logitpars = par[grep('zero', names(par))]   
  mupars = par[grep('mean', names(par))]
  theta = exp(par[grep('theta', names(par))])
  
  #### Regressing Logit Link Part ####
  
  Xlogit  = X
  LPlogit = Xlogit %*% logitpars
  logi0   = plogis(LPlogit)  
  
  #### Regressing Log Link Part adjusted for Marginalized Model ####
  
  Xgp  = X
  mugp = exp(Xgp %*% mupars - log(1 - logi0)) ## - log(1 - logi0) is the adjustment ##
  
  #### Log Likelihoods (LLK) for Logit and Log Part ####
  
  logliklogit = log(logi0 + exp(log(1 - logi0) - mugp))
  loglikgp  = log(1 - logi0) + dgenpois2(y, meanpar = mugp, disppar = theta, log = TRUE)
  
  #### Estimating Joint LLK for Zero and Count Parts ####
  
  y.zero = y == 0  
  y.count = y > 0
  loglik = sum(logliklogit[y.zero]) + sum(loglikgp[y.count])
  return(-loglik + lambda * (par %*% par))
}

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

######################################################################
############################## Main Function #########################
######################################################################

### MLE ###
MAZIC.MLE = function(y, x, startVals, opt.method = c('Nelder-Mead', 'BFGS'),
                     mod = c("marzip","marzinb", "marzigp"), verbose = TRUE)
{
  
  #### Estimating Model Parameters ####
  
  fn = eval(parse(text = paste('reg', mod, sep = '.')))
  fit = optim(par = startVals,
              fn = fn,
              X = x,
              y = y,
              method = opt.method,
              control = list(maxit = 5000, reltol = 1e-12),
              hessian = TRUE)
  
  #### Gathering Estimates ####
  
  npar = length(fit$par)
  n = length(y)
  B  = fit$par
  se = tryCatch({
    se = suppressWarnings(sqrt(diag(solve((-fit$hessian)))))
  },
  error = function(e){
    return(NULL)
  })
  if(is.null(se)){
    se = rep(NA, length(B))
    Z  = rep(NA, length(B))
    p  = rep(NA, length(B))
  }else{
    Z  = B/se
    p  = pnorm(abs(Z), lower = FALSE)*2
  }
  # logLK = fit$value
  # aic = -(2 * logLK) + (2 * npar)
  # bic = (-2 * logLK) + (log(n) * npar)
  
  #### Summarizing Outputs ####
  
  resTable = data.frame(Estimate = B, Std.Error = se, z.value = Z, Pval = p)
  resTable.coeff.zero = resTable[grep('zero', rownames(resTable)),]
  resTable.coeff.mean = resTable[grep('mean', rownames(resTable)),]
  resTable.coeff = data.frame(rbind(resTable.coeff.zero, resTable.coeff.mean))
  if(nrow(resTable.coeff) == nrow(resTable)){
    add.parsEstimates = NULL
  }else{
    ex.pars = resTable[!rownames(resTable) %in% rownames(resTable.coeff),]
    add.parsEstimates = ex.pars[, c('Estimate', 'Std.Error')]
  }
  
  #### Extracting Variance Covariance Matrix and beta MLE estimates ####
  
  vcov = tryCatch(solve(as.matrix(-fit$hessian)),
                  error=function(e) {
                    warning(e$message, call=FALSE)
                    k <- nrow(as.matrix(fit$hessian))
                    return(matrix(NA, k, k))
                  })
  beta = resTable.coeff$Estimate
  names(beta) = rownames(resTable.coeff)
  mles = resTable$Estimate
  names(mles) = rownames(resTable)
  
  if(verbose){
    message('Snapshot of the MLE estimates')
    print(mles)
    message('Snapshot of the variance covariance matrix of MLE estimates')
    print(vcov[1:5, 1:5])
  }
  
  #### Giving Back Likelihood and Variance Covariance matrix ####
  
  fintab = list(Regression.tab = resTable, vcov = vcov, beta.estimate = beta)
  return(fintab)
}

### Ridge ###
MAZIC.ridge = function(y, x, startVals, opt.method = c('Nelder-Mead', 'BFGS'),
                     mod = c("marzip","marzinb", "marzigp"), lambda, verbose = TRUE)
{
  
  #### Estimating Model Parameters ####
  
  fn = eval(parse(text = paste('reg', mod, sep = '.')))
  fit = optim(par = startVals,
              fn = fn,
              X = x,
              y = y,
              method = opt.method,
              control = list(maxit = 5000, reltol = 1e-12),
              lambda = 0.001,
              hessian = TRUE)
  
  #### Gathering Estimates ####
  
  npar = length(fit$par)
  n = length(y)
  B  = fit$par
  se = tryCatch({
    se = suppressWarnings(sqrt(diag(solve((-fit$hessian)))))
  },
  error = function(e){
    return(NULL)
  })
  if(is.null(se)){
    se = rep(NA, length(B))
    Z  = rep(NA, length(B))
    p  = rep(NA, length(B))
  }else{
    Z  = B/se
    p  = pnorm(abs(Z), lower = FALSE)*2
  }
  # logLK = fit$value
  # aic = -(2 * logLK) + (2 * npar)
  # bic = (-2 * logLK) + (log(n) * npar)
  
  #### Summarizing Outputs ####
  
  resTable = data.frame(Estimate = B, Std.Error = se, z.value = Z, Pval = p)
  resTable.coeff.zero = resTable[grep('zero', rownames(resTable)),]
  resTable.coeff.mean = resTable[grep('mean', rownames(resTable)),]
  resTable.coeff = data.frame(rbind(resTable.coeff.zero, resTable.coeff.mean))
  if(nrow(resTable.coeff) == nrow(resTable)){
    add.parsEstimates = NULL
  }else{
    ex.pars = resTable[!rownames(resTable) %in% rownames(resTable.coeff),]
    add.parsEstimates = ex.pars[, c('Estimate', 'Std.Error')]
  }
  
  #### Extracting Variance Covariance Matrix and beta MLE estimates ####
  
  vcov = tryCatch(solve(as.matrix(-fit$hessian)),
                  error=function(e) {
                    warning(e$message, call=FALSE)
                    k <- nrow(as.matrix(fit$hessian))
                    return(matrix(NA, k, k))
                  })
  beta = resTable.coeff$Estimate
  names(beta) = rownames(resTable.coeff)
  mles = resTable$Estimate
  names(mles) = rownames(resTable)
  
  if(verbose){
    message('Snapshot of the MLE estimates')
    print(mles)
    message('Snapshot of the variance covariance matrix of MLE estimates')
    print(vcov[1:5, 1:5])
  }
  
  #### Giving Back Likelihood and Variance Covariance matrix ####
  
  fintab = list(Regression.tab = resTable, vcov = vcov, beta.estimate = beta)
  return(fintab)
}

####

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
  
  fitlass = glmnet::cv.glmnet(cov.star, y.star)
  lam = fitlass$lambda.min
  fitlass = glmnet::glmnet(cov.star, y.star, lambda = lam)
  
  #### Collecting Penalized Results ####
  
  coef = as.matrix(coef(fitlass, s = lam))
  colnames(coef) = 'Penalized.Estimates'
  tLL = fitlass$nulldev - deviance(fitlass)
  k = fitlass$df
  n = fitlass$nobs
  aic = -2*tLL+2*k
  bic = log(n)*k - tLL
  
  #### Giving Back Intercept Terms ####
  
  coef = as.data.frame(coef)
  coef = coef[!rownames(coef) %in% '(Intercept)',, drop = FALSE]
  tmp = sigma.12 %*% solve(sigma.22) %*% (beta - as.vector(coef$Penalized.Estimates))
  coef.zero = coef[grep('zero', rownames(coef)), , drop = FALSE]
  coef.mean = coef[grep('mean', rownames(coef)), , drop = FALSE]
  Penalized.Coeff.zero = c(tmp[1, 1], coef.zero$Penalized.Estimates)
  names(Penalized.Coeff.zero)[2:length(Penalized.Coeff.zero)] = rownames(coef.zero)
  Penalized.Coeff.mean = c(tmp[2, 1], coef.mean$Penalized.Estimates)
  names(Penalized.Coeff.mean)[2:length(Penalized.Coeff.mean)] = rownames(coef.mean)
  All.Penalized.Coeffs = list(zero = Penalized.Coeff.zero, mean = Penalized.Coeff.mean)
  
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