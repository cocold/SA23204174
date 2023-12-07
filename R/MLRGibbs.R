#' @title A dataset used for illustration.
#' @name data
#' @description This dataset is used to fit models.
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' summary(data)
#' }
NULL

#' @title Initialize rho
#' @description Initialize rho
#' @param M the number of subpopulations
#' @return a vector of size
#' @examples
#' \dontrun{
#' vector_rho_0 <- init_rho(3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
init_rho = function(M){
  arr<-rep(M,M)
  vector_rho_0 = rDirichlet(1,arr)
  return(vector_rho_0)
}

#' @title Initialize r
#' @description Initialize r
#' @param X the independent variables
#' @param M the number of subpopulation
#' @return a matrix 
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' vector_r_0 <- init_r(X,3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
init_r = function(X, M){
  # r - active variables
  vec_bern = c()
  for(x in 1:M)
  {
    vec_bern = rbind(vec_bern, rbinom(dim(X)[2], size  = 1,prob = 1))
  }
  vector_r_0 = vec_bern
  colnames(vector_r_0) = colnames(X)
  return(vector_r_0)
}

#' @title Initialize z
#' @description Initialize z
#' @param X the independent variables
#' @param vector_rho_0 the rho vector
#' @return a matrix 
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' rho<-init_rho(3)
#' z<-ini_z(X,rho)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
init_z = function(X, vector_rho_0)
{
  sample_z = rmultinom(dim(X)[1], 1, prob = vector_rho_0)
  matrix_z_0 = t(sample_z)
  return(matrix_z_0)
}

#' @title Initialize theta
#' @description Initialize theta
#' @param X the independent variables
#' @param M the number of subpopulation
#' @param matrix_z_0 the z matrix
#' @param matrix_r_0 the r matrix
#' @param vector_rho_0 the rho vector
#' @param y the dependent variable
#' @return return an array of beta betahat rho sigma oldbeta r 
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' rho <- init_rho(3)
#' r<-init_r(X,3)
#' z<-ini_z(X,rho)
#' theta<-init_theta(3,z,r,rho,X)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
##Initialize Beta, Beta-hat, sigma and w_m_inv
init_theta = function(M,matrix_z_0,matrix_r_0,vector_rho_0,X,y)
{
  # initialize w_m,it gives the inverse of the matrix w_m
  w_m_inv_0 = list()
  for(idx_m in 1:M)
  {
    # define X_m(r_m)
    x_m_rm_0 = as.matrix(X[which(matrix_z_0[,idx_m] == 1),which(matrix_r_0[idx_m,] == 1),drop=FALSE])
    # check if the matrx X_m(r_m)'*X_m(r_m) is singular or not
    if(is.singular.matrix(t(x_m_rm_0)%*%x_m_rm_0))
    {
      # w_m(r_m)^(-1) = (X_m(r_m)'*X_m(r_m) + lambda_m*I)
      # lambda_m = 1/p (p is the number of covariates)
      one_w_m_inv = t(x_m_rm_0)%*%x_m_rm_0 + (1/dim(x_m_rm_0)[2])*diag(dim(x_m_rm_0)[2])
    }
    else
    {
      # w_m(r_m)^(-1) = (X_m(r_m)'*X_m(r_m)) 
      one_w_m_inv = t(x_m_rm_0)%*%x_m_rm_0
    }
    # add inverse of w_m(r_m) to vector
    w_m_inv_0[[idx_m]] = one_w_m_inv
  }
  # initialize beta_hat, beta and sigma
  vector_beta_hat_0 = list()
  vector_beta_0 = list()
  vector_beta_o_0=list()
  vector_sigma_0 = list()
  for(idx_m in 1:M)
  {
    # beta_hat = [X_m(r_m)'*X_m(r_m)]^(-1) * X_m(r_m)'*Y_m
    x_m_rm_0 = as.matrix(X[which(matrix_z_0[,idx_m] == 1),which(matrix_r_0[idx_m,] == 1),drop=FALSE])
    w_m = solve(w_m_inv_0[[idx_m]])
    selected_y = y[which(matrix_z_0[,idx_m] == 1), drop=F]
    beta_hat_0 = w_m%*%t(x_m_rm_0)%*%selected_y
    # now build the variance of beta, beta_variance = g_m * sigma_m^2 * w_m(r_m)
    # g_m
    g_m = sum(matrix_z_0[,idx_m])
    # sigma_0 - selected value equals to 1 in order to avoid 
    sigma_0 = 1 #rinvgamma(1, 0.1, 0.1)
    # beta_variance 
    var_beta_m_0 = g_m*sigma_0*solve(w_m_inv_0[[idx_m]])
    # sample of beta_m(r_m)
    beta_m_0 = mvrnorm(n = 1, mu = beta_hat_0, Sigma = var_beta_m_0)
    beta_m_o_0= array(rep(0,dim(X)[2]),dim=c(dim(X)[2],1),dimnames = list(colnames(X)))         
    # update arrays of beta, beta_hat and sigma
    vector_beta_hat_0[[idx_m]] = as.matrix(beta_hat_0)
    vector_beta_0[[idx_m]] = as.matrix(beta_m_0)
    for (cov_name in rownames(vector_beta_0[[idx_m]]))
    {
      beta_m_o_0[cov_name,1]=vector_beta_0[[idx_m]][cov_name,1]
    }
    vector_beta_o_0[[idx_m]]=as.matrix(beta_m_o_0)
    vector_sigma_0[[idx_m]] = sigma_0
  }
  # theta initialization
  theta = list('beta' = vector_beta_0,'beta_hat' = vector_beta_hat_0,'sigma' = vector_sigma_0,'rho' = vector_rho_0,'r' = matrix_r_0,'old_beta' = vector_beta_o_0,'w_m_inv' = w_m_inv_0)
}

#' @title Initialize all the parameters
#' @description Initialize all the parameters
#' @param X the independent variables
#' @param y the dependent variables
#' @param M the number of subpopulation
#' @return an array of theta and z
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
init= function(X, y, M)
{
  vector_rho_0 = init_rho(M)
  matrix_r_0 = init_r(X,M)
  matrix_z_0 = init_z(X,vector_rho_0)
  theta = init_theta(M,matrix_z_0,matrix_r_0,vector_rho_0,X,y)
  return(list('theta' = theta, 'matrix_z' = matrix_z_0))
}

#' @title Update Z:only one scalar at a time
#' @description Update Z:only one scalar at a time
#' @param Xm the independent variables
#' @param ym the dependent variables
#' @param M the number of subpopulation
#' @param theta_array the defined theta above
#' @return an scalar of z
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' GZ<-GetPosteriorZ(y[1],X[1,],theta_array,3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
GetPosteriorZ = function(ym,Xm, theta_array, M)
{
  num = c()
  for(idx_m in 1:M)
  {
    selected_beta = theta_array[['beta']][[idx_m]]
    label_beta = rownames(selected_beta)
    one_x = Xm[label_beta]
    mu_normal = one_x%*%selected_beta
    sigma_m = theta_array[['sigma']][[idx_m]]
    f_norm = dnorm(ym, mean =  mu_normal, sd = sqrt(sigma_m))
    rho_m = theta_array[['rho']][[idx_m]]
    prod_f_rho = rho_m*f_norm
    num = c(num, prod_f_rho)
  }
  if (sum(num)==0)
  {
    prob_z=rep(1/M,M)
  }
  else
    prob_z = num/sum(num)
  return(sample(x = 1:M, size=1, prob = prob_z))
}

#' @title Update the whole vector
#' @description Update the whole vector
#' @param X the independent variables
#' @param y the dependent variables
#' @param M the number of subpopulation
#' @param theta_array the defined theta above
#' @return an matrix of z
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
UpdateArrayZ = function(theta_array,y,X,M)
{
  vector_z = c()
  for(idx_obs in 1:dim(X)[1])
  {
    ym = y[idx_obs]
    Xm = X[idx_obs,]
    z_idx_m<-rep(0,M)
    z_idx_m[GetPosteriorZ(ym,Xm, theta_array, M)]<-1
    vector_z = c(vector_z, z_idx_m)
  }
  matrix_z<-matrix(vector_z,nrow=dim(X)[1],ncol=M,byrow=TRUE)
  return(matrix_z)
}

#' @title Update the rho vector
#' @description Update the rho vector
#' @param matrix_z the z matrix
#' @param M the number of subpopulation
#' @param theta_array the defined theta above
#' @return updated theta_array
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' theta_array<-UpdateRho(theta_array,new_Z,3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
UpdateRho = function(theta_array,matrix_z,M)
{
  vector_n_proportion = colSums(matrix_z) + M
  vector_rho = rDirichlet(1,vector_n_proportion)
  theta_array$rho=vector_rho
  return(theta_array)
}

#' @title Update the subpopulation
#' @description Update the subpopulation
#' @param matrix_z the z matrix
#' @param M the number of subpopulation
#' @param theta_array the defined theta above
#' @param X the independent variables
#' @return updated theta_array
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' covariates<-UpdateSubPop(theta_array,new_Z,X,3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
UpdateSubPop = function(theta_array,matrix_z,X,M)
{
  new_x = list()
  for(idx_m in 1:M)
  {
    one_for_subpop = X[which(matrix_z[,idx_m] == 1),which(theta_array[['r']][idx_m,] == 1),drop=FALSE]
    new_x[[idx_m]] = one_for_subpop
  }
  return(new_x)
}

#' @title Update the w_inv
#' @description Update the w_inv
#' @param covariates the new x subpopulation
#' @param M the number of subpopulation
#' @param theta_array the defined theta above
#' @return updated theta_array
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' covariates<-UpdateSubPop(theta_array,new_Z,X,3)
#' theta_array<-UpdateW_inv(theta_array, covariates, 3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
UpdateW_inv = function(theta_array, covariates, M){
  for(idx_m in 1:M){
    x_product = t(covariates[[idx_m]])%*%covariates[[idx_m]]
    if(is.singular.matrix(x_product)){
      lambda_m = (1/dim(covariates[[idx_m]])[2])
      identity_matrix = diag(dim(covariates[[idx_m]])[2])
      one_w_m_inv = x_product + lambda_m*identity_matrix
    }
    else{
      one_w_m_inv = x_product
    }
    theta_array$w_m_inv[[idx_m]] = one_w_m_inv
  }
  return(theta_array)
}

#' @title Update the betahat
#' @description Update the betahat
#' @param covariates the new x subpopulation
#' @param matrix_z the z matrix
#' @param theta_array the defined theta above
#' @param y the dependent variable
#' @return updated theta_array
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' covariates<-UpdateSubPop(theta_array,new_Z,X,3)
#' vector_beta_hat<-UpdateBetaHat(theta_array, new_Z, covariates)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
UpdateBetaHat = function(theta_array, matrix_z, covariates,y)
{
  vector_beta_hat = list()
  for(idx_m in 1:length(theta_array$beta))
  {
    w_m = solve(theta_array$w_m_inv[[idx_m]])
    y_m = y[which(matrix_z[,idx_m] == 1)]
    beta_hat = w_m%*%t(covariates[[idx_m]])%*%y_m
    vector_beta_hat[[idx_m]] = beta_hat
  }
  return(vector_beta_hat)
}

#' @title Update sigma
#' @description Update sigma
#' @param covariates the new x subpopulation
#' @param matrix_z the z matrix
#' @param theta_array the defined theta above
#' @param M the number of subpopulations
#' @param y the dependent variable
#' @return updated theta_array
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' covariates<-UpdateSubPop(theta_array,new_Z,X,3)
#' theta_array<-UpdateSigma(theta_array, new_Z, covariates,3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
UpdateSigma = function(theta_array,matrix_z,covariates,M,y)
{
  # update a_m, b_m and sigma_m
  vector_q = rowSums(theta_array[['r']])
  vector_n_proportion = colSums(matrix_z)
  #update a 
  vector_a = vector_n_proportion + vector_q + 0.001
  # update beta-hat
  theta_array[['beta_hat']] = UpdateBetaHat(theta_array, matrix_z, covariates,y)
  # update_b
  vector_b = c()
  for(idx_m in 1:length(theta_array$beta))
  {
    selected_beta = theta_array[['beta']][[idx_m]]
    x_m_beta_m = covariates[[idx_m]]%*%selected_beta
    num_b = y[which(matrix_z[,idx_m] == 1)] - x_m_beta_m
    num_b = t(num_b)%*%num_b
    selected_beta_hat = theta_array[['beta_hat']][[idx_m]]
    diff_beta_beta_hat = selected_beta - selected_beta_hat
    w_m_inv = theta_array$w_m_inv[[idx_m]]
    frac_num = t(diff_beta_beta_hat)%*%w_m_inv%*%diff_beta_beta_hat
    second_term = frac_num/vector_n_proportion[[idx_m]]
    one_b = num_b + second_term + 0.001
    vector_b = c(vector_b, one_b)
  }
  # sample sigma_m
  vector_sigma = c()
  for(idx_m in 1:length(theta_array$beta))
  {
    sigma_sample = 1/rgamma(n = 1, shape = vector_a[idx_m]/2, scale = vector_b[idx_m]/2)
    vector_sigma = c(vector_sigma, sigma_sample)
  }
  theta_array$sigma = vector_sigma
  return(theta_array)
}

#' @title Update Beta
#' @description Update Beta
#' @param covariates the new x subpopulation
#' @param matrix_z the z matrix
#' @param theta_array the defined theta above
#' @param M the number of subpopulations
#' @param y the dependent variable
#' @return updated theta_array
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' covariates<-UpdateSubPop(theta_array,new_Z,X,3)
#' theta_array<-UpdateBeta(theta_array,covariates,new_Z,3)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
UpdateBeta = function(matrix_z,covariates,theta_array,M,y)
{
  # initialization: Omega, n_m, w_m
  vector_omega = list()
  # update beta
  for(idx_m in 1:length(theta_array$beta))
  {
    proportion = sum(matrix_z[,idx_m])
    cov_multiplied = t(covariates[[idx_m]])%*%covariates[[idx_m]]
    w_m_inv = theta_array$w_m_inv[[idx_m]]
    num_omega = proportion*cov_multiplied + w_m_inv
    sigma_m = theta_array$sigma[[idx_m]]
    one_inverse_omega = num_omega/(proportion*sigma_m)
    omega = solve(one_inverse_omega)
    vector_omega[[idx_m]] = omega
  }
  theta_array$omega = vector_omega
  # mu
  vector_mu = list()
  for(idx_m in 1:length(theta_array$beta))
  {
    proportion = sum(matrix_z[,idx_m])
    x_y_multiplied = t(covariates[[idx_m]])%*%y[which(matrix_z[,idx_m] == 1)]
    first_term = proportion*x_y_multiplied
    w_m_inv = theta_array$w_m_inv[[idx_m]]
    beta_hat = theta_array$beta_hat[[idx_m]]
    second_term = w_m_inv%*%beta_hat
    one_omega = vector_omega[[idx_m]]
    sigma_m = theta_array$sigma[[idx_m]]
    mu = one_omega%*%(first_term + second_term)
    mu = mu/(proportion*sigma_m)
    vector_mu[[idx_m]] = mu
  }
  # sample of beta
  vector_beta = list()
  for(idx_m in 1:length(theta_array$beta))
  {
    estimate = as.matrix(mvrnorm(n=1, mu =vector_mu[[idx_m]],Sigma = vector_omega[[idx_m]])) 
    vector_beta[[idx_m]] = estimate
    for(cov_name in rownames(estimate))
    {
      theta_array$old_beta[[idx_m]][cov_name,1] = estimate[cov_name,1]
    }
  }
  #print('Vector beta')
  #print(vector_beta)
  theta_array$beta = vector_beta
  return(theta_array)
}

#' @title Generate LogLike of r
#' @description Generate LogLike of r
#' @param idx_m the mth model of the subpopulation
#' @param matrix_z the z matrix
#' @param y the dependent variable
#' @param theta_array the defined theta above
#' @param X the independent variables
#' @param single_label the chosen variable to be decided whether to be included in the mth model
#' @return updated theta_array
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' prob<-ComputeProb(theta_array,new_Z, 1,1,X)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
ComputeProb = function(theta_array,matrix_z, idx_m,single_label,X,y)
{
  # compute the prob of having active a specific variable in a sub-population
  active_covariates = labels(which(theta_array$r[idx_m,] == 1))
  if(single_label %in% active_covariates)
  {
    ### 1st case: we can compute easily l_n(0) and l_n(1) 
    # create loglike with 0 and 1
    mu_vector = c()
    for(idx_obs in which(matrix_z[,idx_m] == 1))
    {
      specific_covariates = X[idx_obs,active_covariates ,drop = F]
      beta_m = theta_array$beta[[idx_m]]
      mu = specific_covariates%*%beta_m
      mu_vector = c(mu_vector, mu)
    }
    one_y = y[which(matrix_z[,idx_m] == 1)]
    sigma_m = theta_array$sigma[[idx_m]]
    prob_1 = dnorm(x = one_y,mean = mu_vector,sd=sqrt(sigma_m))
    prob_1 = sum(log(prob_1))
    # prob0
    without_label = active_covariates[active_covariates!= single_label]
    mu_vector = c()
    for(idx_obs in which(matrix_z[,idx_m] == 1))
    {
      specific_covariates = X[idx_obs,without_label ,drop = F]
      beta_m_without_j = theta_array$beta[[idx_m]][without_label,]
      mu = specific_covariates%*%beta_m_without_j
      mu_vector = c(mu_vector, mu)
    }
    one_y = y[which(matrix_z[,idx_m] == 1)]
    sigma_m = theta_array$sigma[[idx_m]]
    prob_0 = dnorm(x = one_y,mean = mu_vector,sd=sqrt(sigma_m))
    prob_0 = sum(log(prob_0))
    probab = c(prob_0, prob_1)
    res = probab
  }
  else
  {
    ### 2nd case: we need to re-activate the variable, 
    ### using the last non-zero value computed before 
    #prob0
    mu_vector = c()
    for(idx_obs in which(matrix_z[,idx_m] == 1))
    {
      specific_covariates = X[idx_obs,active_covariates ,drop = F]
      beta_m = theta_array$beta[[idx_m]]
      mu = specific_covariates%*%beta_m
      mu_vector = c(mu_vector, mu)
    }
    one_y = y[which(matrix_z[,idx_m] == 1)]
    sigma_m = theta_array$sigma[[idx_m]]
    prob_0 = dnorm(x = one_y,  mean = mu_vector,sd=sqrt(sigma_m))
    prob_0 = sum(log(prob_0))
    # prob1
    mu_vector = c()
    # add the deactivated covariate j, taking its last computed value (beta*) 
    beta_m = theta_array$beta[[idx_m]]
    last_value_computed = as.matrix(theta_array$old_beta[[idx_m]][single_label,])
    reactivated = rbind(beta_m, last_value_computed)
    # order by names
    reactivated = as.matrix(reactivated[order(rownames(reactivated)),])
    labels_reactivated = rownames(reactivated)
    for(idx_obs in which(matrix_z[,idx_m] == 1))
    {
      mu = X[idx_obs,labels_reactivated,drop = F]%*%reactivated
      mu_vector = c(mu_vector, mu)
    }
    one_y = y[which(matrix_z[,idx_m] == 1)]
    sigma_m = theta_array$sigma[[idx_m]]
    prob_1 = dnorm(x = one_y,mean = mu_vector,sd=sqrt(sigma_m))
    prob_1 = sum(log(prob_1))
    # define the probabilities
    probab = c(prob_0, prob_1)
    res = probab
  }
  return(res)
}

#' @title update r
#' @description update r
#' @param matrix_z the z matrix
#' @param theta_array the defined theta above
#' @param X the independent variable
#' @param y the dependent variable
#' @return updated theta_array
#' @examples
#' \dontrun{
#' X<-matrix(rep(0,100),nrow=10,ncol=10)
#' y<-rep(1,100)
#' initialized<-init(3,X,y)
#' theta_array<-initialized$theta
#' new_Z<-UpdateArrayZ(theta_array,y,X,3)
#' prob<-UpdateR(theta_array,new_Z)
#' }
#' @import DIRECT
#' @import mcmc
#' @export
UpdateR = function(theta_array, matrix_z,X,y){
  # initialize r
  old_r_array = theta_array$r
  new_r_array = theta_array$r
  n_subpop = dim(theta_array$r)[1]
  for(idx_m in 1:n_subpop)
  {
    # append the likelihood in a matrix
    # compute probabilities and sample r
    for(single_label in colnames(theta_array$r))
    {
      if(single_label!='intercept')
      {
        likelihood = ComputeProb(theta_array,matrix_z, idx_m,single_label,X,y)
        one_prob = 1/(1 + exp(likelihood[1] - likelihood[2]))
        probability = c(1- one_prob, one_prob)
        one_r = sample(x = 0:1, size = 1, prob = probability)
        new_r_array[idx_m,single_label] = one_r
        # update label (keep (1) or discard (0))
      }
    }
    ### at this point we need to: activate/deactive the covariates (changing beta values)
    # reactivate beta by r
    new_active = c()
    for(single_label in colnames(theta_array$r))
    {
      if((old_r_array[idx_m,single_label] == 0) & (new_r_array[idx_m,single_label] == 1))
      {
        new_active = c(new_active, single_label)
        beta_m = theta_array$beta[[idx_m]]
        last_value = as.matrix(theta_array$old_beta[[idx_m]][single_label,] )
        reactivated = rbind(beta_m, last_value)
        reactivated = as.matrix(reactivated[order(rownames(reactivated)),])
        theta_array$beta[[idx_m]] = reactivated
      }
    }
  }
  # deactivate the beta by r
  theta_array$r = new_r_array
  for(idx_m in 1:length(theta_array$beta))
  {
    active_labels = names(which(theta_array$r[idx_m,]==1))
    theta_array$beta[[idx_m]] = as.matrix(theta_array$beta[[idx_m]][active_labels,])
  }
  return(theta_array)
}
#' @import microbenchmark
#' @import TSA
#' @import bootstrap
#' @import boot
#' @import DAAG
#' @import coda
#' @import gridExtra
#' @import stargazer
#' @import matrixcalc
#' @import MASS
#' @import mcmc
#' @import DIRECT
#' @import readxl
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma dnorm rbinom rmultinom 
#' @useDynLib SA23204174
NULL
