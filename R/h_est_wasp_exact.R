linear_prog_h<-function(h_est, n_subset, mesh_size=100, n_samp=100){
  library(lpSolve)
  library(Matrix)

  for (i in 1:n_subset) {                          #####sample from each subset
    varname<-as.character(paste0('post_samp_',i))
    samp<-matrix(NA,nrow=n_samp,ncol=length(h_est[[1]]$est))
    for (j in 1:length(h_est[[i]]$est)){
      samp[,j]<-rnorm(n_samp,mean=h_est[[i]]$est[j],sd=h_est[[i]]$se[j])
    }

    assign(as.character(paste0('samp_',i)),samp)
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


  ####h####
  dens_h<-list(rep(NA,length(h_est[[1]]$est)))
  wasp_samp_r<-list(rep(NA,length(h_est[[1]]$est)))
  for (j in 1:length(h_est[[1]]$est)) {
    overall_r<-overallmin[j]+seq(1,mesh_size,1)*
      (overallmax[j]-overallmin[j])/mesh_size
    for (i in 1:n_subset) {

      D[[i]]<-matrix(0,length(overall_r),n_samp)
      samp<-as.character(paste0('samp_',i))
      D[[i]]<-(fields::rdist(overall_r,get(samp)[,j]))^2
    }

    Res_r<-overall_dens_h(n_subset, n_samp, overall_r, D)
    wasp_samp_r[[j]]<-Res_r$samp

    if (j%%(length(h_est[[1]]$est)/10)==0) {
      message(paste((j/(length(h_est[[1]]$est)/10)*10),'% of h has been estimated.'))}
  }

  print('h done')


  wasp_samp<-list()
  samp_wasp<-matrix(unlist(wasp_samp_r),nrow = 1000,ncol = length(h_est[[1]]$est))
  wasp_samp$est<-colMeans(samp_wasp)
  wasp_samp$se<-apply(samp_wasp,2,sd)
  rm(samp_wasp)

  return(list(wasp_samp=wasp_samp))
}

overall_dens_h<-function(n_subset, n_samp, overall_samp, D) {
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

  # result <- gurobi(model, param)
  # aOptSol = as.matrix(result$x)[1:N, 1]

  #Define this soft-thresholding operator to remove small elements.
  softThresOper = function(x,t){
    a=sign(x)*max(abs(x) - t, 0)
    return(a)}

  # Recover the solution.
  # aOptSol = as.matrix(result$x)[1:N, 1];
  xRest   = as.matrix(result$x)[(N+1):length(result$x)];
  for (p in 1:n_subset) {
    tOptSol = matrix(xRest[1:(N*Ni)], N, Ni)/Ni;
    xRest      = xRest[(N*Ni+1):length(xRest)];
    tOptSol = softThresOper(tOptSol, 1e-10);
  }

  ###calculate the overall distribution density####
  pr <- as.numeric(aOptSol)
  wasp_samp<- overall_samp[sample(1:nrow(as.matrix(overall_samp)), 1000, replace = TRUE, prob = pr)]


  return(list(samp=wasp_samp))
}


