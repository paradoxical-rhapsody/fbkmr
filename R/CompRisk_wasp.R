##############Overall risks###############################
postmeanhnew_wasp<-function(fit,X,y,Z, Znew = NULL, sel = NULL) {


  if (!is.null(Znew)) {
    if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
    if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)
  }

  len <- length(fit$lambda)
  sigsq.eps <- mean(fit$sigsq.eps[(len/2+1):len])
  r <-colMeans(fit$r[(len/2+1):len,])
  beta <- colMeans(as.matrix(fit$beta[(len/2+1):len,]))
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

########function##############
riskSummary.approx_wasp <- function(point1, point2, preds.fun, n_subset=1,...) {
  cc <- c(-1, 1)
  newz <- rbind(point1, point2)
  preds <- preds.fun(newz, ...)
  diff <- drop(cc %*% preds$postmean)
  pred.var<-preds$postvar/n_subset
  diff.sd <- drop(sqrt(cc %*% pred.var %*% cc))
  c(est = diff, sd = diff.sd)
}




#' Calculate overall risk summaries
#'
#' Compare estimated \code{h} function when all predictors are at a particular quantile to when all are at a second fixed quantile
#' @inheritParams skmbayes
#' @inheritParams postmeanhnew_wasp
#' @inherit postmeanhnew_wasp details
#' @param qs vector of quantiles at which to calculate the overall risk summary
#' @param q.fixed a second quantile at which to compare the estimated \code{h} function
#' @export
#'
OverallRiskSummaries_wasp <- function(fit, y = NULL, Z = NULL, X = NULL, Z.q=NULL, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "approx", sel = NULL, n_subset=1) {


    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X


  point1 <- apply(Z.q, 2, quantile, q.fixed)
  if (method %in% c("approx", "exact")) {
    preds.fun <- function(znew) postmeanhnew_wasp(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
    riskSummary <- riskSummary.approx_wasp
  } else {
    stop("method must be one of c('approx', 'exact')")
  }
  risks.overall <- t(sapply(qs, function(quant) riskSummary(point1 = point1, point2 = apply(Z.q, 2, quantile, quant), preds.fun = preds.fun, n_subset=n_subset)))
  risks.overall <- data.frame(quantile = qs, risks.overall)
}



########################Singular Risk################
##########function####################
VarRiskSummary_wasp <- function(whichz = 1, fit, y = NULL, Z = NULL, X = NULL, Z.q=NULL, qs.diff = c(0.25, 0.75), q.fixed = 0.5, method = "approx", sel = NULL, n_subset=1,...) {


    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X


  point2 <- point1 <- apply(Z.q, 2, quantile, q.fixed)
  point2[whichz] <- apply(Z.q[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
  point1[whichz] <- apply(Z.q[, whichz, drop = FALSE], 2, quantile, qs.diff[1])
  # point1 <- makePoint(whichz, Z, qs.diff[1], q.fixed)
  # point2 <- makePoint(whichz, Z, qs.diff[2], q.fixed)
  if (method %in% c("approx", "exact")) {
    preds.fun <- function(znew) postmeanhnew_wasp(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
    riskSummary <- riskSummary.approx_wasp
  }  else {
    stop("method must be one of c('approx', 'exact')")
  }
  riskSummary(point1 = point1, point2 = point2, preds.fun = preds.fun, ...)
}

#' Single Variable Risk Summaries
#'
#' Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile)
#'
#' @inheritParams skmbayes
#' @inheritParams ExtractEsts
#' @inheritParams OverallRiskSummaries_wasp
#' @inherit postmeanhnew_wasp details
#' @param qs.diff vector indicating the two quantiles \code{q_1} and \code{q_2} at which to compute \code{h(z_{q2}) - h(z_{q1})}
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in \code{Z}
#' @param z.names optional vector of names for the columns of \code{z}
#' @param ... other argumentd to pass on to the prediction function
#' @param which.z vector indicating which variables (columns of \code{Z}) for which the summary should be computed
#' @export
SingVarRiskSummaries_wasp <- function(fit, y = NULL, Z = NULL, X = NULL, Z.q=NULL, which.z = 1:ncol(Z), qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), method = "approx", sel = NULL, z.names = colnames(Z), n_subset=1, ...) {


    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X


  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))

  df <- dplyr::tibble()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk <- VarRiskSummary_wasp(whichz = which.z[j], fit = fit, y = y, Z = Z, X = X, Z.q=Z.q, qs.diff = qs.diff, q.fixed = q.fixed[i], method = method, sel = sel, n_subset=n_subset,...)
      df0 <- dplyr::tibble(q.fixed = q.fixed[i], variable = z.names[j], est = risk["est"], sd = risk["sd"])
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df <- dplyr::mutate_(df, variable = ~factor(variable, levels = z.names[which.z]), q.fixed = ~as.factor(q.fixed))
  attr(df, "qs.diff") <- qs.diff
  df
}


