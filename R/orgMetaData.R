#' Organize Design Matrix
#'
#' @description The function properly organize design matrix for input to MAZIC L1 function
#'
#' @param mat A dataframe and the formula of the model
#'
#' @return An organized design matrix
#' @export

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