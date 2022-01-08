.ridgeLMold <- function(Y, X, lambda, target){

	## ---------------------------------------------------------------------
	## Targeted ridge estimation of the linear regression parameter
	## ---------------------------------------------------------------------
	## Arguments
	## Y      : response vector
	## X      : design matrix
	## lambda : penalty parameter
	## target : target for regression parameter
	## ---------------------------------------------------------------------
	## Value:
	## The ridge estimate of the regression parameter, a numeric.
	## ---------------------------------------------------------------------
	## Authors     : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	if (ncol(X) == 1){
		# if there is only a single covariate
		return(solve(crossprod(X, X) + rep(lambda, ncol(X))) %*% 
		       (crossprod(X, Y) + lambda * target))
	}
	if (nrow(X) >= ncol(X)){
		# low-dimensional case
		return(crossprod(solve(crossprod(X, X) + diag(rep(lambda, ncol(X)))),
 		                 (crossprod(X, Y) + lambda * target)))
	}
	if (nrow(X) < ncol(X)){
		# high-dimensional case
		XXT       <- tcrossprod(X) / lambda
		diag(XXT) <- diag(XXT) + 1
		XYpT      <- (crossprod(X, Y) + lambda * target)
		return((XYpT / lambda - lambda^(-2) * crossprod(X, t(crossprod(tcrossprod(X, t(XYpT)), solve(XXT)))))[,1])
	}
}




.ridgeBLMold <- function(Y, X, lambda, target, minSuccDiff, maxIter){

	## ---------------------------------------------------------------------
	## Targeted ridge estimation of the logistic regression parameter
	## ---------------------------------------------------------------------
	## Arguments 
	## Y           : A numeric being the response vector comprising zero's 
	##               and ones only.
	## X           : A matrix, the design matrix. The number of rows should
	##               match the number of elements of Y.
	## lambda      : A numeric, the ridge penalty parameter.
	## target      : A numeric towards which the estimate is shrunken.
	## minSuccDiff : A numeric, the minimum distance between two successive 
	##               estimates to be achieved.
	## maxIter     : A numeric specifying the maximum number of iterations.
	## ---------------------------------------------------------------------
	## Value:
	## The ridge estimate of the regression parameter, a numeric.
	## ---------------------------------------------------------------------
	## Authors     : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# initiate
	loglikPrev <- -10^(10)
	if (ncol(X) >= nrow(X)){ 
		XXT     <- tcrossprod(X) / lambda
		Xtarget <- t(tcrossprod(target, X))
	}
	lp <- rep(0, length(Y))
	for (iter in 1:maxIter){
		# calculate the weights
		Ypred <- 1 / (1 + exp(-lp))
		W0    <- Ypred * (1 - Ypred)
		if (min(W0) <= .Machine$double.eps){
			W0[which(W0 < .Machine$double.eps)] <- .Machine$double.eps
		}

		# distinguish between the high- and low-dimensional cases
		if (ncol(X) >= nrow(X)){ 
			# obtain part one of the estimator
			Z <- lp + (Y - Ypred)/W0

			# now obtain the IRWLS update efficiently
			diag(XXT) <- diag(XXT) + 1/W0
			slh       <- crossprod(solve(XXT), Z-Xtarget)
			diag(XXT) <- diag(XXT) - 1/W0
			lp        <- crossprod(XXT, slh)
			penalty   <- sum(lp * slh) / 2
			lp        <- Xtarget + lp
			loglik    <- sum(Y * lp - log(1 + exp(lp))) - penalty

		} else {
			# obtain part one of the estimator
			XWZpT <- lambda * target + crossprod(X, W0 * lp + Y - Ypred)[,1]

			# evaluate X %*% X^t and add the reciprocal weight to 
			# its diagonal
			XTX       <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
			diag(XTX) <- diag(XTX) + lambda 
	
			# now obtain the IRWLS update efficiently
			betas <- crossprod(solve(XTX), XWZpT)
			lp    <- tcrossprod(X, t(betas))

			# evaluate the loglikelihood
			loglik <- sum(Y * lp - log(1 + exp(lp))) - 0.5 * lambda * sum((betas - target)^2)
		}

		# assess convergence
		if (abs(loglik - loglikPrev) < minSuccDiff){ 
			break 
		} else {
			loglikPrev <- loglik
		}
	}
	if (ncol(X) >= nrow(X)){ 
		betas <- target + crossprod(X, slh) / lambda
	}
	return(betas=betas)
}


