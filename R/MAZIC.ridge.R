#' Ridge Regression estimates
#'
#' @description The function returns the ridge estimates of each covariate
#'
#' @param A vector of response, a matrix of covariates, starting values, optimization method,
#' names of the distribution, and verbose option
#'
#' @return The ridge estimates
#' @export

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