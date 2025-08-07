wasp_bayes<-function(n_subset, mesh_size=100, n_samp=100, Z, X, kmbayes_res,iter,... ) {


  for (i in 1:n_subset) {                          #####sample from each subset
    varname<-as.character(paste0('post_samp_',i))
    assign(varname, kmbayes_res[[i]])
    max_iter=length(get(varname)$lambda)
    assign(as.character(paste0('samp_',i)),cbind(get(varname)$lambda[(max_iter-n_samp+1):max_iter], get(varname)$beta[(max_iter-n_samp+1):max_iter,], get(varname)$sigsq.eps[(max_iter-n_samp+1):max_iter],
                                                 get(varname)$r[(max_iter-n_samp+1):max_iter,]))
  }

  samp_total<-rbind(samp_1, samp_2)               ######sample combined
  if (n_subset>=3) {
    for (i in 3:n_subset) {
      varname<-as.character(paste0('samp_',i))
      samp_total<-rbind(samp_total,get(varname))
    }
  }

  overallmin<-apply(samp_total, 2, min)
  overallmax<-apply(samp_total, 2, max)

  ######number of atoms in WASP
  D=list(rep(NA,n_subset))
  ####kappa####
  overall_lambda<-overallmin[1]+seq(1,mesh_size,1)*(overallmax[1]-overallmin[1])/mesh_size
  for (i in 1:n_subset) {

    D[[i]]<-matrix(0,length(overall_lambda),n_samp)
    samp<-as.character(paste0('samp_',i))
    D[[i]]<-(fields::rdist(overall_lambda,get(samp)[,1]))^2
  }

  Res_lambda<-overall_dens(n_subset, n_samp, overall_lambda, D,iter, ...)
  dens_lambda<-Res_lambda$dens
  wasp_samp_lambda<-Res_lambda$samp
  message('lambda done')
  idx=1


  ####beta####
  dens_beta<-list(rep(NA,ncol(post_samp_1$beta)))
  wasp_samp_beta<-list(rep(NA,ncol(post_samp_1$beta)))
  for (j in 1:ncol(post_samp_1$beta)) {
    overall_beta<-overallmin[idx+j]+seq(1,mesh_size,1)*
      (overallmax[idx+j]-overallmin[idx+j])/mesh_size
    for (i in 1:n_subset) {
      D[[i]]<-matrix(0,length(overall_beta),n_samp)
      samp<-as.character(paste0('samp_',i))
      D[[i]]<-(fields::rdist(overall_beta,get(samp)[,(idx+j)]))^2
    }

    Res_beta<-overall_dens(n_subset, n_samp, overall_beta, D,iter,...)
    dens_beta[[j]]<-Res_beta$dens
    wasp_samp_beta[[j]]<-Res_beta$samp
  }

  idx=idx+ncol(post_samp_1$beta)

  message('beta done')


  ####sigmasq####
  n_sigma<-ifelse(is.null(ncol(post_samp_1$sigsq.eps)),1,ncol(post_samp_1$sigsq.eps))
  dens_sigmasq<-list(rep(NA,n_sigma))
  wasp_samp_sigmasq<-list(rep(NA,n_sigma))
  for (j in 1:n_sigma) {
    overall_sigmasq<-overallmin[idx+j]+seq(1,mesh_size,1)*
      (overallmax[idx+j]-overallmin[idx+j])/mesh_size
    for (i in 1:n_subset) {
      D[[i]]<-matrix(0,length(overall_sigmasq),n_samp)
      samp<-as.character(paste0('samp_',i))
      D[[i]]<-(fields::rdist(overall_sigmasq,get(samp)[,(idx+j)]))^2
    }

    Res_sigmasq<-overall_dens(n_subset, n_samp, overall_sigmasq, D,iter,...)
    dens_sigmasq[[j]]<-Res_sigmasq$dens
    wasp_samp_sigmasq[[j]]<-Res_sigmasq$samp
  }

  idx=idx+n_sigma
  message('sigmasq done')


  ####r####
  dens_r<-list(rep(NA,ncol(post_samp_1$r)))
  wasp_samp_r<-list(rep(NA,ncol(post_samp_1$r)))
  for (j in 1:ncol(post_samp_1$r)) {
    overall_r<-overallmin[idx+j]+seq(1,mesh_size,1)*
      (overallmax[idx+j]-overallmin[idx+j])/mesh_size
    for (i in 1:n_subset) {

      D[[i]]<-matrix(0,length(overall_r),n_samp)
      samp<-as.character(paste0('samp_',i))
      D[[i]]<-(fields::rdist(overall_r,get(samp)[,(idx+j)]))^2
    }

    Res_r<-overall_dens(n_subset, n_samp, overall_r, D,iter,...)
    dens_r[[j]]<-Res_r$dens
    wasp_samp_r[[j]]<-Res_r$samp

  }

  idx=idx+ncol(post_samp_1$r)
  message('r done')


  wasp_samp<-list()
  wasp_samp$lambda<-wasp_samp_lambda
  wasp_samp$beta<-matrix(unlist(wasp_samp_beta),nrow=iter,ncol=ncol(kmbayes_res[[1]]$beta))
  wasp_samp$sigsq.eps<-unlist(wasp_samp_sigmasq)
  wasp_samp$r<-matrix(unlist(wasp_samp_r),nrow = iter,ncol = ncol(kmbayes_res[[1]]$r))

  return(list(dens_lambda=dens_lambda,dens_beta=dens_beta, dens_sigmasq=dens_sigmasq,
              dens_r=dens_r,
              wasp_samp=wasp_samp))
}

overall_dens<-function(n_subset, n_samp, overall_samp, D,iter,...) {
  # constants
  K  = n_subset
  Ni = n_samp
  N  = length(overall_samp)
  nx = N * (N+1)
  mx = K * N + N + 1
  In = diag(x=1, nrow=N, ncol = N)
  En =matrix(rep(1,N),nrow=1, ncol = N)

  # Generate matrix A0.
  A0  = matrix()
  for (p in 1:K) {
    Rp_i=In/Ni
    Rp=t(rep(1,Ni))%x%Rp_i
    if (p ==1) {A0  = Rp }
    else {A0  = bdiag(A0, Rp)}
  }


  A00=-1*In
  A00=(rep(1,K))%x%A00

  A0 = cbind(A00,A0)
  A0 = drop0(A0)
  b0 = matrix(0, nrow(A0), 1)

  # Generate matrix B from simplex constraints.
  B = matrix()
  for (p in 1:(Ni*K+1)) {
    if (p ==1) {B  = En }
    else {B = bdiag(B, En)}
  }
  # The hold matrix C.

  A = rbind(A0,B)
  A = drop0(A)
  # Generate the right hand size vector b.
  b = rbind(matrix(0, nrow = K*N, ncol=1), matrix(1, nrow=Ni*K + 1, ncol=1))
  b = drop0(b)

  # Generate the cost vector
  costVec=matrix(0, nrow = N, ncol=1)
  for (i in 1:n_subset) {
    costCell=D[[i]]/dim(D[[i]])[2]
    costVec=cbind(costVec,costCell)
  }

  c = drop0(costVec)

  #####using the LP solver
  model=list()
  model$A=A
  model$obj=as.matrix(c)
  model$rhs=as.matrix(b)
  model$modelsense='min'
  model$sense='='

  param=list()
  param$method = 2
  param$Presolve = 2
  param$Crossover = 0
  param$outputflag = 0

  result<- gurobi(model, param)
  #Define this soft-thresholding operator to remove small elements.
  softThresOper = function(x,t){
    a=sign(x)*max(abs(x) - t, 0)
    return(a)}

  # Recover the solution.
  aOptSol = as.matrix(result$x)[1:N, 1];
  xRest   = as.matrix(result$x)[(N+1):length(result$x)];
  for (p in 1:n_subset) {
    tOptSol = matrix(xRest[1:(N*Ni)], N, Ni)/Ni;
    xRest      = xRest[(N*Ni+1):length(xRest)];
    tOptSol = softThresOper(tOptSol, 1e-10);
  }

  ###calculate the overall distribution density####
  pr <- as.numeric(aOptSol)
  wasp_samp<- overall_samp[sample(1:nrow(as.matrix(overall_samp)), iter, replace = TRUE, prob = pr)]
  rr<-range(wasp_samp)
  bw<-max(diff(as.numeric(overall_samp)))
  if (bw<=0) {bw=1}
  dens<-bkde(wasp_samp,bandwidth=bw, range.x=rr)

  return(list(dens=dens, samp=wasp_samp))
}


makeKpart_wasp <- function(r, Z1, Z2 = NULL) {
  Z1r <- sweep(Z1, 2, sqrt(r), "*")
  if (is.null(Z2)) {
    Z2r <- Z1r
  } else {
    Z2r <- sweep(Z2, 2, sqrt(r), "*")
  }
  Kpart <- fields::rdist(Z1r, Z2r)^2
  Kpart <- exp(-Kpart)           ##K_Z,r=exp(-sum_M(r_m*(Z_im-Z_jm)^2))
  Kpart
}

####calculate V-matrix####
makeVmatrix_wasp <-function(r,lambda,Z){
  Kpart<-makeKpart_wasp(r,Z)
  V<-diag(1,nrow(Z),nrow(Z))+lambda*Kpart    ##V_lambda,Z,r=I_n+lambda*K_Z,r
  V
}

####calculate v-matrix inverse####
makeVinv_wasp<-function(r,lambda,Z){
  V<-makeVmatrix_wasp(r,lambda,Z)
  cholV <- chol(V)            ##V is PSD thus use Choleski decomp to calculate inverse
  Vinv <- chol2inv(cholV)
  logdetVinv<--2*sum(log(diag(cholV)))     ##need log transform since the detVinv is too large to calculate in MH step.
  detVinv<-exp(logdetVinv)
  Vcomp<-list(Vinv=Vinv,detVinv=detVinv,logdetVinv=logdetVinv)
  Vcomp
}


H.update_wasp<-function(samp_post,y,X,Z) {  ##h|beta,sigmasq,lambda,r,theta,y~N(lambda*K*Vinv(y-xbeta),sigmasq*lambda*K*Vinv)
  r<-ifelse(is.na(apply(samp_post$r,2,mean)),0,apply(samp_post$r,2,mean))
  Vcomp<-makeVinv_wasp(r,mean(samp_post$lambda),Z)
  Kpart<-makeKpart_wasp(r,Z)
  lambdaK<-mean(samp_post$lambda)*Kpart
  lambdaKV<-lambdaK %*% Vcomp$Vinv
  hmean<-lambdaKV %*% (y-X%*%colMeans(as.matrix(samp_post$beta)))
  hvar<-mean(samp_post$sigsq.eps)*lambdaKV
  sq<-try(chol(crossprod(Kpart,Vcomp$Vinv)),silent = TRUE)
  if (class(sq)[1]=='try-error'){
    sigsvd <- svd(crossprod(Kpart,Vcomp$Vinv))
    sq <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  }
  hstd<-sqrt(mean(samp_post$sigsq.eps)*mean(samp_post$lambda))*sq

  h.star<-hmean+crossprod(hstd,rnorm(length(hmean)))
  Hcomp<-h.star
  Hcomp
}


####generating hnew from Znew####
hnew.update_wasp<-function(samp_post,y,X,Z,Znew){
  ntot<-nrow(Z)+nrow(Znew)

  r<-ifelse(is.na(apply(samp_post$r,2,mean)),0,apply(samp_post$r,2,mean))
  Kpartmar<-makeKpart_wasp(r,Znew)
  Kpartcomb<-makeKpart_wasp(r,Z,Znew)
  Vcomp<-makeVinv_wasp(r,mean(samp_post$lambda),Z)
  laKpartcombVinv<-mean(samp_post$lambda)*t(Kpartcomb) %*% Vcomp$Vinv          ##hnew|beta,sigmasq,lambda,r,theta,y~N(lambda*t(Kznew,z)*Vinv*(y-xbeta),sigmasq*lambda(Kznew-lambda*t(Kznew,z)*Vinv*Kznew,z))
  hnewmean<-laKpartcombVinv %*% (y-X%*%as.matrix(colMeans(samp_post$beta)))
  hnewvar<-mean(samp_post$sigsq.eps)*mean(samp_post$lambda)*(Kpartmar-laKpartcombVinv%*%Kpartcomb)
  hnewstd<-chol(hnewvar)

  hnew<-hnewmean+crossprod(hnewstd,rnorm(nrow(Znew)))
  hnew
}


##################################################################################################################
######################################run wasp kmbayes############################################################
##################################################################################################################

#' Fit scalable Bayesian kernel machine regression
#'
#' Fits the Bayesian kernel machine regression (BKMR) model using Markov chain Monte Carlo (MCMC) methods with divide-and-conquer technique.
#'
#' @export
#'
#' @param y a vector of outcome data of length \code{n}.
#' @param Z an \code{n}-by-\code{M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n}-by-\code{K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
#' @param n_subset number of subsets to split the data
#' @param iter number of iterations to run the sampler
#' @param family a description of the error distribution and link function to be used in the model. Currently implemented for \code{gaussian} and \code{binomial} families.
#' @param id optional vector (of length \code{n}) of grouping factors for fitting a model with a random intercept. If NULL then no random intercept will be included.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate diagnostic information during the model fitting.
#' @param Znew optional matrix of new predictor values at which to predict \code{h}, where each row represents a new observation. This will slow down the model fitting, and can be done as a post-processing step using \code{\link{SamplePred}}
#' @param file_path working directory for saving the resulting rds files.
#' @param starting.values list of starting values for each parameter. If not specified default values will be chosen.
#' @param control.params list of parameters specifying the prior distributions and tuning parameters for the MCMC algorithm. If not specified default values will be chosen.
#' @param varsel TRUE or FALSE: indicator for whether to conduct variable selection on the Z variables in \code{h}
#' @param groups optional vector (of length \code{M}) of group indictors for fitting hierarchical variable selection if varsel=TRUE. If varsel=TRUE without group specification, component-wise variable selections will be performed.
#' @param knots optional matrix of knot locations for implementing the Gaussian predictive process of Banerjee et al (2008). Currently only implemented for models without a random intercept.
#' @param ztest optional vector indicating on which variables in Z to conduct variable selection (the remaining variables will be forced into the model).
#' @param rmethod for those predictors being forced into the \code{h} function, the method for sampling the \code{r[m]} values. Takes the value of 'varying' to allow separate \code{r[m]} for each predictor; 'equal' to force the same \code{r[m]} for each predictor; or 'fixed' to fix the \code{r[m]} to their starting values
#' @param est.h TRUE or FALSE: indicator for whether to sample from the posterior distribution of the subject-specific effects h_i within the main sampler. This will slow down the model fitting.
#'
#' @return an object of class "bkmrfit", which has the associated methods:
#' \itemize{
#'   \item \code{\link{print}} (i.e., \code{\link{print.bkmrfit}})
#'   \item \code{\link{summary}} (i.e., \code{\link{summary.bkmrfit}})
#' }
#'
#' @import utils foreach doParallel parallel

kmbayes_Wasp<-function(Z,X,y,n_subset,n_samp=200, iter=1000, parallel=FALSE, n_cores=NULL, varsel=FALSE, est.h=TRUE, Znew=NULL,knots=NULL, file_path=NULL,save_loc=FALSE,...){
  X=as.matrix(X)

  kmbayes_res<-list()
  m=floor(nrow(Z)/ceiling(n_subset))
  mod<-nrow(Z)-m*ceiling(n_subset)

  time1=Sys.time()
  kmbayes_res_list=list()


  # parallel::detectCores()
	# n.cores <- ceiling(parallel::detectCores()/2)
	# my.cluster <- parallel::makeCluster(n.cores)
	# doParallel::registerDoParallel(cl = my.cluster)
	# ---- PARALLELIZATION SETUP ----
	if (parallel) {
		if (is.null(n_cores)) {
			n_cores <- ceiling(parallel::detectCores()/2)
		}
		nphys <- parallel::detectCores()
		if (n_cores > nphys) {
			warning(sprintf("n_cores(%d) > detectCores(%d)", n_cores, nphys))
		}
		
		cl <- parallel::makeCluster(n_cores)   # default: PSOCK
		on.exit({
		  parallel::stopCluster(cl)   # 1. Shut down the cluster and free its system resources
		  foreach::registerDoSEQ()    # 2. Reset the foreach backend to sequential execution
		}, add = TRUE)

		
		doParallel::registerDoParallel(cl)
		
		message("Parallelization started: ", foreach::getDoParName(),
										", workers = ", foreach::getDoParWorkers())
	}

  kmbayes_res_list<-foreach(i =1:ceiling(n_subset))%dopar%{
    if (i <=mod){
      n_f=i*m+i
      Z_i=Z[((i-1)*m+i):n_f,]
      X_i=as.matrix(X[((i-1)*m+i):n_f,])
      y_i=y[((i-1)*m+i):n_f]
      indx=seq(((i-1)*m+i),n_f,1)-((i-1)*m+i)+1
    }
    else{n_f=min(i*m+mod, nrow(Z))
    Z_i=Z[((i-1)*m+mod+1):n_f,]
    X_i=as.matrix(X[((i-1)*m+mod+1):n_f,])
    y_i=y[((i-1)*m+mod+1):n_f]
    indx=seq(((i-1)*m+mod+1),n_f,1)-((i-1)*m+mod)
    }

    knots100<-try(fields::cover.design(Z_i,nd=100)$design,silent = TRUE)

    if(is.null(knots)){knots100=NULL}
    if (class(knots100)[1]=='try-error'){knots100=NULL}
    kmbayes_res<-kmbayes(y = y_i, Z = Z_i, X = X_i, Znew=Znew, iter = iter, knots=knots100, control.param=list(lambda.jump=rep(2.8,1),mu.lambda=rep(15,1),sigma.lambda=rep(8,1)),verbose = F, varsel = varsel,sketch.type = 'sub_sampling',p=ceiling(nrow(Z)/n_subset),n_subset=n_subset, indx=indx,Snum=1,est.h = FALSE,h.post.type = h.post.type,est.last.h = 1000)

    if (save_loc==T){
      saveRDS(kmbayes_res, file = paste0(file_path,'/wasp_res_',i,'.RDS'))
      cat(paste('NUMBER ',i,' subset is completed, system time is: ', difftime(Sys.time(),time1)),file=paste0(file_path,'/wasp_res.txt'),sep='\n',append = TRUE)
    }


    return(kmbayes_res)
  }
  # stopCluster(cl = my.cluster)
  message(paste('Subset MCMC is completed, system time is: ', round(difftime(Sys.time(),time1, units = "mins"),2),'mins'))

  # for (i in 1:ceiling(n_subset)){
  #   if (i <=mod){
  #     n_f=i*m+i
  #     Z_i=Z[((i-1)*m+i):n_f,]
  #     X_i=as.matrix(X[((i-1)*m+i):n_f,])
  #     y_i=y[((i-1)*m+i):n_f]
  #     indx=seq(((i-1)*m+i),n_f,1)-((i-1)*m+i)+1
  #   }
  #   else{n_f=min(i*m+mod, nrow(Z))
  #   Z_i=Z[((i-1)*m+mod+1):n_f,]
  #   X_i=as.matrix(X[((i-1)*m+mod+1):n_f,])
  #   y_i=y[((i-1)*m+mod+1):n_f]
  #   indx=seq(((i-1)*m+mod+1),n_f,1)-((i-1)*m+mod)
  #   }
  #
  #   knots100<-try(fields::cover.design(Z_i,nd=100)$design,silent = TRUE)
  #
  #   if(is.null(knots)){knots100=NULL}
  #   if (class(knots100)[1]=='try-error'){knots100=NULL}
  #   kmbayes_res<-kmbayes(y = y_i, Z = Z_i, X = X_i, Znew=Znew, iter = iter, knots=knots100, control.param=list(lambda.jump=rep(2.8,1),mu.lambda=rep(15,1),sigma.lambda=rep(8,1)),verbose = F, varsel = varsel,sketch.type = 'sub_sampling',p=ceiling(nrow(Z)/n_subset),n_subset=n_subset, indx=indx,Snum=1,est.h = FALSE,h.post.type = h.post.type,est.last.h = 1000)
  #
  #   if (save_loc==T){
  #     saveRDS(kmbayes_res, file = paste0(file_path,'/wasp_res_',i,'.RDS'))
  #     cat(paste('NUMBER ',i,' subset is completed, system time is: ', difftime(Sys.time(),time1)),file=paste0(file_path,'/wasp_res.txt'),sep='\n',append = TRUE)
  #   }
    return(kmbayes_res_list)

  }









linear_prog<-function(kmbayes_res=kmbayes_res,n_subset=n_subset,iter=iter,n_samp=n_samp,est.h=est.h,y=y,X=X,Z=Z)  {

  if (n_subset==1) {
    samp<-list()
    samp$lambda<-as.matrix(kmbayes_res[[1]]$lambda)
    samp$beta<-as.matrix(kmbayes_res[[1]]$beta)
    samp$sigsq.eps<-kmbayes_res[[1]]$sigsq.eps
    samp$r<-as.matrix(kmbayes_res[[1]]$r)



    Res_final<-list()
    Res_final<-samp
    #Res_final$lambda_mean<-mean(samp$lambda)
    #Res_final$lambda_std<-sd(samp$lambda)
    #Res_final$beta_mean<-colMeans(as.matrix(samp$beta))
    #Res_final$beta_std<-apply(as.matrix(samp$beta),2,sd)
    #Res_final$sigmasq_mean<-mean(samp$sigmasq)
    #Res_final$sigmasq_std<-sd(samp$sigmasq)
    #Res_final$r_mean<-colMeans(samp$r)
    #Res_final$r_std<-apply(samp$r,2,sd)
    #Res_final$h_hat<-h_hat_wasp
  }
  else{



    den<-wasp_bayes(ceiling(n_subset),60,n_samp=n_samp,Z,X,kmbayes_res=kmbayes_res,iter=iter)
    #h_hat_wasp<-H.update_wasp(den$wasp_samp,y,X,Z)

    Res_final<-den$wasp_samp
    #Res_final<-list()
    #Res_final$lambda_mean<-mean(den$wasp_samp$lambda)
    #Res_final$lambda_std<-sd(den$wasp_samp$lambda)
    #Res_final$beta_mean<-colMeans(den$wasp_samp$beta)
    #Res_final$beta_std<-apply(den$wasp_samp$beta,2,sd)
    #Res_final$sigmasq_mean<-mean(den$wasp_samp$sigmasq)
    #Res_final$sigmasq_std<-sd(den$wasp_samp$sigmasq)
    #Res_final$r_mean<-colMeans(den$wasp_samp$r)
    #Res_final$r_std<-apply(den$wasp_samp$r,2,sd)
    #Res_final$h_hat<-h_hat_wasp

  }

  if (est.h){
    h_hat_wasp<-H.update_wasp(Res_final,y=y,X=X,Z=Z)
    h_hat_mean<-h_hat_wasp
    Res_final$h_hat <- h_hat_wasp
  }


  message(paste0('WASP done'))
 return(Res_final)
}



#' Fit scalable Bayesian kernel machine regression
#'
#' Fits the Bayesian kernel machine regression (BKMR) model using Markov chain Monte Carlo (MCMC) methods with divide-and-conquer technique.
#'
#'
#' @param y a vector of outcome data of length \code{n}.
#' @param Z an \code{n}-by-\code{M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n}-by-\code{K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
#' @param n_subset number of subsets to split the data
#' @param n_samp number of samples drew from final MCMC
#' @param iter number of iterations to run the sampler
#' @param family a description of the error distribution and link function to be used in the model. Currently implemented for \code{gaussian} and \code{binomial} families.
#' @param id optional vector (of length \code{n}) of grouping factors for fitting a model with a random intercept. If NULL then no random intercept will be included.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate diagnostic information during the model fitting.
#' @param Znew optional matrix of new predictor values at which to predict \code{h}, where each row represents a new observation. This will slow down the model fitting, and can be done as a post-processing step using \code{\link{SamplePred}}
#' @param file_path working directory for saving the resulting rds files.
#' @param linprog TRUE or FALSE: flag indicating whether to preform linear programming for combining samples from subset MCMC results.
#' @param starting.values list of starting values for each parameter. If not specified default values will be chosen.
#' @param control.params list of parameters specifying the prior distributions and tuning parameters for the MCMC algorithm. If not specified default values will be chosen.
#' @param varsel TRUE or FALSE: indicator for whether to conduct variable selection on the Z variables in \code{h}
#' @param groups optional vector (of length \code{M}) of group indictors for fitting hierarchical variable selection if varsel=TRUE. If varsel=TRUE without group specification, component-wise variable selections will be performed.
#' @param knots optional matrix of knot locations for implementing the Gaussian predictive process of Banerjee et al (2008). Currently only implemented for models without a random intercept.
#' @param ztest optional vector indicating on which variables in Z to conduct variable selection (the remaining variables will be forced into the model).
#' @param rmethod for those predictors being forced into the \code{h} function, the method for sampling the \code{r[m]} values. Takes the value of 'varying' to allow separate \code{r[m]} for each predictor; 'equal' to force the same \code{r[m]} for each predictor; or 'fixed' to fix the \code{r[m]} to their starting values
#' @param est.h TRUE or FALSE: indicator for whether to sample from the posterior distribution of the subject-specific effects h_i within the main sampler. This will slow down the model fitting.
#' @param save_loc TRUE or FALSE indicator for whether to save intemediate or final result on local machine.
#' @param parallel Bool, whether to do parallel tasks. 
#' @param n_cores Number of cores for parallel tasks
#' @return an object of class "bkmrfit", which has the associated methods:
#' \itemize{
#'   \item \code{\link{print}} (i.e., \code{\link{print.bkmrfit}})
#'   \item \code{\link{summary}} (i.e., \code{\link{summary.bkmrfit}})
#' }
#'
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 1000, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#'
#' ## Fit model with 5 subsests
#' set.seed(111)
#' fitkm <- skmbayes(y = y, Z = Z, X = X, iter = 1000, n_subset=5, verbose = FALSE, varsel = FALSE)
#'
#' @import  mvtnorm foreach doParallel lpSolve Matrix KernSmooth nlme magrittr fields slam parallel
#' @export
#'



skmbayes <- function(Z,X,y,n_subset=1,n_samp=200, iter=1000, parallel=FALSE, n_cores=NULL, varsel=FALSE, est.h=TRUE, Znew=NULL,file_path=NULL,save_loc=FALSE, ...){

  time1 <- Sys.time()
  if (is.null(file_path)){file_path=getwd()}

	kmbayes_res<-kmbayes_Wasp(Z=Z,X=X,y=y,n_subset=n_subset,n_samp=n_samp, iter=iter, parallel=parallel, n_cores=n_cores, varsel=varsel, est.h=est.h, Znew=Znew,file_path=file_path,save_loc=save_loc,...)

  # if (linprog==TRUE){
  #   if(save_loc==T){
  #     for (j in 1:n_subset){
  #       kmbayes_res[[j]]<-readRDS(file = paste0(file_path,'/wasp_res_',j,'.RDS'))
  #     }
  #   }

  #   N_part <- floor(sqrt(n_subset))
  #   res_split <- split(seq(1,n_subset,1), ceiling(seq_along(seq(1,n_subset,1))/N_part))

  #   parallel::detectCores()
  #   n.cores <- ceiling(parallel::detectCores()/2)
  #   my.cluster <- parallel::makeCluster(n.cores)
  #   doParallel::registerDoParallel(cl = my.cluster)

  #   Res_final_split<-foreach(i=1:ceiling(n_subset/N_part))%dopar%{
  #     sel <- res_split[[i]]
  #     kmbayes_res_split<-kmbayes_res[sel]

  #     Res_final_s<-linear_prog(kmbayes_res=kmbayes_res_split,n_subset=length(sel),iter=iter,n_samp=n_samp,Z=Z,X=X,y=y,est.h = est.h)

  #     if (save_loc==T){
  #       saveRDS(Res_final_s, file = paste0(file_path,'/comb_',i,'.RDS'))
  #     }
  #     return(Res_final_s)
  #   }
  #   stopCluster(cl = my.cluster)

  #   message(paste('Split linear programming is completed, system time is: ', round(difftime(Sys.time(),time1, units = "mins"),2),'mins'))

  #   Res_final<-linear_prog(kmbayes_res=Res_final_split,n_subset=ceiling(n_subset/N_part),iter=iter,n_samp=n_samp,Z=Z,X=X,y=y,est.h = est.h)

  #   message(paste('FBKMR is completed, system time is: ', round(difftime(Sys.time(),time1, units = "mins"),2),'mins'))

  #       time2 <- Sys.time()

  #       Res_final$X <- X
  #       Res_final$Z <- Z
  #       Res_final$y <- y
  #       Res_final$n_subset <- n_subset
  #       Res_final$iter <- iter
  #       Res_final$est.h <-est.h
  #       Res_final$varsel <- varsel
  #       Res_final$lambda<-as.matrix(Res_final$lambda)
  #       Res_final$time1 <-time1
  #       Res_final$time2 <-time2

  #       if(save_loc==T){saveRDS(Res_final, file = paste0(file_path,'/wasp_res_final','.RDS'))}
  #       class(Res_final) <- c("bkmrfit", class(Res_final))
  #       return(Res_final)
  # }
  #if(linprog==F){
    message(paste('FBKMR is completed without overall result, system time is: ', round(difftime(Sys.time(),time1, units = "mins"),2),'mins'))
    return(kmbayes_res)
  #}


  }


