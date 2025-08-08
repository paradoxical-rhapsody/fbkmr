##############Overall risks###############################

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
#' @param n.cores number of cores for parallel computing
#' @param qs vector of quantiles at which to calculate the overall risk summary
#' @param q.fixed a second quantile at which to compare the estimated \code{h} function
#' @import utils foreach doParallel parallel
#' @export
#'
OverallRiskSummaries_wasp <- function(fit, y = NULL, Z = NULL, X = NULL, parallel = FALSE, n_cores = NULL, Z.q=NULL, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "approx", sel = NULL, n_subset=1, ...) {
  if (class(fit)[1]=='bkmrfit'){fit<-list(fit)}
  n_subset=length(fit)
  risks.overall<-list()

  # parallel::detectCores()
	# 
	# if(is.null(n.cores)){n.cores <- ceiling(parallel::detectCores()/2)}
	# my.cluster <- parallel::makeCluster(n.cores)
	# doParallel::registerDoParallel(cl = my.cluster)
	# ---- PARALLELIZATION SEPUP ----
	if (parallel) {
		if (is.null(n_cores)) {
			n_cores <- ceiling(parallel::detectCores()/2)
		}
		nphys <- parallel::detectCores()
		if (n_cores > nphys) {
			warning(sprintf("n_cores(%d) > detectCores(%d)", n_cores, nphys))
		}
		
		cl <- parallel::makeCluster(n_cores)
		on.exit({
		  parallel::stopCluster(cl)   # 1. Shut down the cluster and free its system resources
		  foreach::registerDoSEQ()    # 2. Reset the foreach backend to sequential execution
		}, add = TRUE)
		
		doParallel::registerDoParallel(cl)
		
		message("Parallelization started: ", foreach::getDoParName(),
										", workers = ", foreach::getDoParWorkers())
	}

  risks.overall<-foreach(j=1: n_subset)%dopar%{
    if (inherits(fit[[j]], "bkmrfit")) {
      y <- fit[[j]]$y
      Z <- fit[[j]]$Z
      X <- fit[[j]]$X
    }

    point1 <- apply(Z, 2, quantile, q.fixed)
    if (method %in% c("approx", "exact")) {
      preds.fun <- function(znew) postmeanhnew_wasp(fit = fit[[j]], y = y, Z = Z, X = X, Znew = znew, sel = sel)
      riskSummary <- riskSummary.approx_wasp
    } else {
      stop("method must be one of c('approx', 'exact')")
    }
    df <- t(sapply(qs, function(quant) riskSummary.approx_wasp(point1 = point1, point2 = apply(Z, 2, quantile, quant), preds.fun = preds.fun, n_subset=n_subset)))
    df <- data.frame(quantile = qs, df)
    return( df)
  }

  risks.overall_overall <- risks.overall[[1]]
				   
  # stopCluster(cl = my.cluster)
  # if (parallel) stopCluster(cl)

  if(n_subset>1){##median combination of normal distributions is taking average of their mean and std
    for (i in 2:n_subset){
      risks.overall_overall$est<-risks.overall_overall$est+risks.overall[[i]]$est
      risks.overall_overall$sd<-risks.overall_overall$sd+risks.overall[[i]]$sd
    }
    risks.overall_overall$est<-risks.overall_overall$est/n_subset
    risks.overall_overall$sd<-risks.overall_overall$sd/n_subset
  }
  risks.overall_overall

}






########################Singular Risk################
##########function####################
VarRiskSummary_wasp <- function(whichz = 1, fit, y = NULL, Z = NULL, X = NULL, Z.q=NULL, qs.diff = c(0.25, 0.75), q.fixed = 0.5, method = "approx", sel = NULL, n_subset=1,...) {


    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X


  point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
  point2[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
  point1[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile, qs.diff[1])
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
#' @param n.cores number of cores for parallel computing
#' @import utils foreach doParallel parallel
#' @export
SingVarRiskSummaries_wasp <- function(fit, y = NULL, Z = NULL, X = NULL, parallel = FALSE, n_cores = NULL, Z.q=NULL, which.z = 1:ncol(Z), qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), method = "approx", sel = NULL, z.names = colnames(Z), n_subset=1, n.cores=NULL,...) {

  if (class(fit)[1]=='bkmrfit'){fit<-list(fit)}
  n_subset=length(fit)
  risks.singvar<-list()

  # parallel::detectCores()
	# if(is.null(n.cores)){n.cores <- ceiling(parallel::detectCores()/2)}
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

  risks.singvar<-foreach(j=1: n_subset)%dopar%{
    if (inherits(fit[[j]], "bkmrfit")) {
      y <- fit[[j]]$y
      Z <- fit[[j]]$Z
      X <- fit[[j]]$X
    }

    if (is.null(z.names)) {
      z.names <- paste0("z", 1:ncol(Z))
    }
    df <- dplyr::tibble()
    for(i in seq_along(q.fixed)) {
      for(k in seq_along(which.z)) {
        risk <- VarRiskSummary_wasp(whichz = which.z[k], fit = fit[[j]], y = y, Z = Z, X = X, qs.diff = qs.diff, q.fixed = q.fixed[i], method = method, sel = sel, n_subset=n_subset)
        df0 <- dplyr::tibble(q.fixed = q.fixed[i], variable = z.names[k], est = risk["est"], sd = risk["sd"])
        df <- dplyr::bind_rows(df, df0)
      }
    }
    df <- dplyr::mutate_(df, variable = ~factor(variable, levels = z.names[which.z]), q.fixed = ~as.factor(q.fixed))
    attr(df, "qs.diff") <- qs.diff

    return(df)
  }
  # stopCluster(cl = my.cluster)


  risks.singvar_overall<-risks.singvar[[1]]
  if(n_subset>1){##median combination of normal distributions is taking average of their mean and std
    for (i in 2:n_subset){
      risks.singvar_overall$est<-risks.singvar_overall$est+risks.singvar[[i]]$est
      risks.singvar_overall$sd<-risks.singvar_overall$sd+risks.singvar[[i]]$sd
    }
    risks.singvar_overall$est<-risks.singvar_overall$est/n_subset
    risks.singvar_overall$sd<-risks.singvar_overall$sd/n_subset
  }
  risks.singvar_overall


}






