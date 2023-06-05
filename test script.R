#### Load Libraries ####

library(haven)
library(pscl)
library(mcount)
library(zic)
library(MASS)
library(VGAM)

#### Set Path ####

dir = 'C:\\Users\\prith\\Desktop\\MAZIC'
setwd(dir)

source('MAZIC2.0.R', verbose = T)
#### Load Data ####

data(docvisits)
data = docvisits
rm(docvisits)
message('Loading German Health Care Data with dimensions')

## Loading German Health Care Data with dimensions

dim(data)

## [1] 1812   23

#### Setting up Regression Scheme ####

resp = 'docvisits'
varnames = 'all'
if(varnames == 'all'){
  covnames = names(data[, !names(data) %in% resp])
}else{
  covnames = varnames
}
formula = as.formula(paste(resp, '~', paste(covnames, collapse = "+"), sep = ''))
curr_mod = 'marzinb'
opt.method = 'Nelder-Mead'
#startVals = TRUE

#### Applying MAZIC Function ####

fit = MAZIC.L1(data = data, formula = formula, startVals = TRUE, mod = curr_mod, opt.method = opt.method, verbose = TRUE)

