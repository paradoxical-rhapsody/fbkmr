make_r_params_comp <- function(r.params, whichcomp) {
  for(i in seq_along(r.params)) {
    if(length(r.params[[i]]) > 1) {
      r.params[[i]] <- r.params[[i]][whichcomp]
    }
  }
  r.params
}

set.r.params <- function(r.prior, comp, r.params) {
	if(r.prior == "gamma") {
		if(length(r.params$mu.r) > 1) r.params$mu.r <- r.params$mu.r[comp]
		if(length(r.params$sigma.r) > 1) r.params$sigma.r <- r.params$sigma.r[comp]
		if(length(r.params$r.jump1) > 1) r.params$r.jump1 <- r.params$r.jump1[comp]
		if(length(r.params$r.jump2) > 1) r.params$r.jump2 <- r.params$r.jump2[comp]
	}
	if(r.prior %in% c("unif", "invunif")) {
		if(length(r.params$r.a) > 1) r.params$r.a <- r.params$r.a[comp]
		if(length(r.params$r.b) > 1) r.params$r.b <- r.params$r.b[comp]
		if(length(r.params$r.jump2) > 1) r.params$r.jump2 <- r.params$r.jump2[comp]
	}
	r.params
}

set.r.MH.functions <- function(r.prior) {
	if(r.prior == "gamma") {
		# r.params <- list(mu.r, sigma.r, r.muprop, r.jump1, r.jump2)
		rprior.logdens <- function(x, r.params) {
			mu.r <- r.params$mu.r
			sigma.r <- r.params$sigma.r
			dgamma(x, shape=mu.r^2/sigma.r^2, rate=mu.r/sigma.r^2, log=TRUE)
		}
		rprop.gen1 <- function(r.params) {
			r.muprop <- r.params$r.muprop
			r.jump <- r.params$r.jump1
			rgamma(1, shape=r.muprop^2/r.jump^2, rate=r.muprop/r.jump^2)
		}
		rprop.logdens1 <- function(x, r.params) {
			r.muprop <- r.params$r.muprop
			r.jump <- r.params$r.jump1
			dgamma(x, shape=r.muprop^2/r.jump^2, rate=r.muprop/r.jump^2, log=TRUE)
		}
		rprop.gen2 <- function(current, r.params) {
			r.jump <- r.params$r.jump2
			rgamma(1, shape=current^2/r.jump^2, rate=current/r.jump^2)
		}
		rprop.logdens2 <- function(prop, current, r.params) {
			r.jump <- r.params$r.jump2
			dgamma(prop, shape=current^2/r.jump^2, rate=current/r.jump^2, log=TRUE)
		}
		rprop.gen <- function(current, r.params) {
			r.jump <- r.params$r.jump
			rgamma(1, shape=current^2/r.jump^2, rate=current/r.jump^2)
		}
		rprop.logdens <- function(prop, current, r.params) {
			r.jump <- r.params$r.jump
			dgamma(prop, shape=current^2/r.jump^2, rate=current/r.jump^2, log=TRUE)
		}
	}

	if(r.prior == "invunif") {
		# r.params <- list(r.a, r.b, r.jump2)
		rprior.logdens <- function(x, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			ifelse(1/r.b <= x & x <= 1/r.a, -2*log(x) - log(r.b - r.a), log(0))
		}
		rprop.gen1 <- function(r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			1/runif(1, r.a, r.b)
		}
		rprop.logdens1 <- function(x, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			ifelse(1/r.b <= x & x <= 1/r.a, -2*log(x) - log(r.b - r.a), log(0))
		}
		rprop.gen2 <- function(current, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			r.jump <- r.params$r.jump2
			truncnorm::rtruncnorm(1, a = 1/r.b, b = 1/r.a, mean = current, sd = r.jump)
		}
		rprop.logdens2 <- function(prop, current, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			r.jump <- r.params$r.jump2
			log(truncnorm::dtruncnorm(prop, a = 1/r.b, b = 1/r.a, mean = current, sd = r.jump))
		}
		rprop.gen <- function(current, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			r.jump <- r.params$r.jump
			truncnorm::rtruncnorm(1, a = 1/r.b, b = 1/r.a, mean = current, sd = r.jump)
		}
		rprop.logdens <- function(prop, current, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			r.jump <- r.params$r.jump
			log(truncnorm::dtruncnorm(prop, a = 1/r.b, b = 1/r.a, mean = current, sd = r.jump))
		}
	}

	if(r.prior == "unif") {
		# r.params <- list(r.a, r.b, r.jump2)
		rprior.logdens <- function(x, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			dunif(x, r.a, r.b, log=TRUE)
		}
		rprop.gen1 <- function(r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			runif(1, r.a, r.b)
		}
		rprop.logdens1 <- function(x, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			dunif(x, r.a, r.b, log=TRUE)
		}
		rprop.gen2 <- function(current, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			r.jump <- r.params$r.jump2
			truncnorm::rtruncnorm(1, a = r.a, b = r.b, mean = current, sd = r.jump)
		}
		rprop.logdens2 <- function(prop, current, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			r.jump <- r.params$r.jump2
			log(truncnorm::dtruncnorm(prop, a = r.a, b = r.b, mean = current, sd = r.jump))
		}
		rprop.gen <- function(current, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			r.jump <- r.params$r.jump
			truncnorm::rtruncnorm(1, a = r.a, b = r.b, mean = current, sd = r.jump)
		}
		rprop.logdens <- function(prop, current, r.params) {
			r.a <- r.params$r.a
			r.b <- r.params$r.b
			r.jump <- r.params$r.jump
			log(truncnorm::dtruncnorm(prop, a = r.a, b = r.b, mean = current, sd = r.jump))
		}
	}

	list(rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
}
############################################################################################################################
sketching <- function(sketch.type,n,p,indx,Snum){
  if (sketch.type == 'gaussian') { #p <- floor(n^(1/5));;
    I <- diag(p); indx <- NULL
    S <- matrix(0,n,p)
    for (i in c(1:Snum)) S <- S + t(mvrnorm(p,rep(0,n),diag(n)/p))
    S <- S/Snum
  }else if(sketch.type == 'sub_sampling'){ # sketch.type == 'sub sampling'
    
    if(is.null(indx)) {indx <- sample.int(n,size=p,replace=F)}
    I_n <- diag(n) # This code could be made more efficient
    I <- diag(p); S <- sqrt(n/p)*I_n[,indx]
  } else {
    stop("sketch.type must be either None, gaussian or sub sampling")
  }
  return(list(S=S,I=I,indx=indx))
}

############################################################################################################################