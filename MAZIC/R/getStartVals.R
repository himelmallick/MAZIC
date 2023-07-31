#' Obtaining starting values for Optim
#'
#' @description The function creates starting values corresponding to the Marginalized Distributions for feed into Optim function
#'
#' @param A vector of response, a matrix of covariates, start values if exist, type of Marginalized Distribution
#'
#' @return A vector of start values
#' @export


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