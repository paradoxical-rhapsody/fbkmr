#' Compute the posterior mean and variance of \code{h} at a new predictor values
#'
#' @inheritParams skmbayes
#' @param fit An object containing the results returned by a the \code{kmbayes} function
#' @param Znew matrix of new predictor values at which to predict new \code{h}, where each row represents a new observation. If set to NULL then will default to using the observed exposures Z.
#' @param sel selects which iterations of the MCMC sampler to use for inference; see details
#' @export
postmeanhnew_wasp<-function(fit,X,y,Z, Znew = NULL, sel = NULL) {
  if (!is.null(Znew)) {
    if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
    if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)
  }

  len <- length(fit$lambda)
  sigsq.eps <- mean(fit$sigsq.eps[(len/2+1):len])
  r <-colMeans(fit$r[(len/2+1):len,])
  beta <- colMeans(fit$beta[(len/2+1):len,])
  lambda <-mean(fit$lambda[(len/2+1):len])

  K<-makeKpart_wasp(r, Z)
  V <- diag(1, nrow(Z), nrow(Z)) + lambda*K
  cholV <- chol(V)
  Vinv <- chol2inv(cholV)

  if (!is.null(Znew)) {
    n0 <- nrow(Z)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    Kmat <- makeKpart_wasp(r, rbind(Z, Znew))
    Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
    Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
    Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]

    lamK10Vinv <- lambda*Kmat10 %*% Vinv
    postvar <- lambda*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
    postmean <- lamK10Vinv %*% (y - X%*%beta)
  } else {
    lamKVinv <- lambda*K%*%Vinv
    postvar <- lambda*sigsq.eps*(K - lamKVinv%*%K)
    postmean <- lamKVinv %*% (y - X%*%beta)
  }
  ret <- list(postmean = drop(postmean), postvar = postvar)
  ret
}

