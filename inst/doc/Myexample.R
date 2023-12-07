## -----------------------------------------------------------------------------
library(SA23204174)
library(MASS)
library(corrplot)
library(gridExtra)
library(stargazer)
library(MCMCpack)
library(matrixcalc)
library(MASS)
n_obs = 200
n_coeff = 5
rho1 = 0.2
rho2 = 0.3
rho3 = 0.5

true_proportions = c(rho1, rho2, rho3)*n_obs

intercept_beta = 0.5
beta1 = c(1,-2,3,2,4)
beta2 = c(-1,1,1,-2,-2)
beta3 = c(-1,-2,-3,2,-1)

# build covariance matrix
covariance_matrix = matrix(0.5, nrow = length(beta1), ncol = (length(beta1)))
for(i in 1:length(beta1))
{
  for(j in 1:length(beta1))
  {
    covariance_matrix[i,j] = covariance_matrix[i,j]^abs(i - j)
  }
}

dataset = data.matrix(data.frame(mvrnorm(n = n_obs,mu = rep(0,length(beta1)),Sigma = covariance_matrix)))
dataset = cbind(1, dataset)
colnames(dataset)[1] = 'intercept'

beta1 = c(intercept_beta,beta1)
beta2 = c(intercept_beta,beta2)
beta3 = c(intercept_beta,beta3)

y = c()
# beta1
for(i in 1:true_proportions[1])
{
  new_y = dataset[i,]%*%beta1
  y = c(y, new_y)
}

# beta2
for(i in (true_proportions[1]+1): (true_proportions[1]+true_proportions[2]))
{
  new_y = dataset[i,]%*%beta2
  y = c(y, new_y)
}

# beta3
for(i in (1 + true_proportions[1]+true_proportions[2]):n_obs)
{
  new_y = dataset[i,]%*%beta3
  y = c(y, new_y)
}

# add error
for(i in 1:length(y))
{
  y[i] = rnorm(1,0,1) + y[i]
}

TIMES = 1000

M = 3
X = as.matrix(dataset)
X = as.matrix(X[,order(colnames(X))])

## -----------------------------------------------------------------------------
initialized = init(X, y, M)

theta_array= initialized$theta
#beta_test1 = as.matrix(beta1)
#row.names(beta_test1) = row.names(theta_array$beta[[1]])

#beta_test2 = as.matrix(beta2)
#row.names(beta_test2) = row.names(theta_array$beta[[2]])

#theta_array$beta[[1]] =  beta_test1
#theta_array$beta[[2]] = beta_test2

z_array = initialized$matrix_z
history_theta = list()
history_z  = list()


## -----------------------------------------------------------------------------
for(i in 1:TIMES)
{
  new_z = UpdateArrayZ(theta_array,y,X,M)
  history_z[[i]] = new_z
  z_array=new_z
  #print('rho')
  theta_array = UpdateRho(theta_array,z_array,M)
  #print('subpop')
  covariates = UpdateSubPop(theta_array,z_array, X,M)
  #print('W_m inv')
  theta_array = UpdateW_inv(theta_array, covariates, M)
  #print('sigma')
  #print(theta_array$sigma)
  theta_array = UpdateSigma(theta_array, z_array,covariates,M,y)
  #print('beta')
  theta_array = UpdateBeta(z_array,covariates,theta_array,M,y)
  #print(theta_array$sigma)
  #print('r')
  #theta_array = UpdateR(theta_array, z_array, X,y)
  #print(theta_array$r)
  history_theta[[i]] = theta_array
}

## -----------------------------------------------------------------------------
n_covariates = 5
summary_parameters = data.frame(matrix(ncol = n_covariates+1, nrow = M))
colnames(summary_parameters) = row.names(theta_array$old_beta[[1]])

beta_estimates = list()
for(idx_m in 1:M){
  for(theta in history_theta){
    beta_estimates[[toString(idx_m)]] = rbind(beta_estimates[[toString(idx_m)]], theta$beta[[idx_m]])
  }
}

## -----------------------------------------------------------------------------
for(idx_m in 1:length(beta_estimates))
{
  print('Subpopulation n:')
  print(idx_m)
  par(mfrow=c(2,3))
  list_covariates = sort(unique(row.names(beta_estimates[[idx_m]])))
  for(one_cov in list_covariates)
  {
    values = beta_estimates[[idx_m]][which(row.names(beta_estimates[[idx_m]])==one_cov),]
    one_median = median(values)
    summary_parameters[idx_m,one_cov] = one_median
    plot(values,type = 'l', main = one_cov)
  }
}

## -----------------------------------------------------------------------------
library(SA23204174)
p=MyEM(156,185,95)

