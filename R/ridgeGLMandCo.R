##------------------------------------------------------------------------------
## Functions for the ridge signal and error precision estimators from studies with
## replications
##------------------------------------------------------------------------------

.ridgeLM <- function(Y, X, lambda, target){

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
		return(solve(crossprod(X, X) + diag(rep(lambda, ncol(X)))) %*% 
		       (crossprod(X, Y) + lambda * target))
	}
	if (nrow(X) > ncol(X)){
		# high-dimensional case
		XXT       <- tcrossprod(X) / lambda
		diag(XXT) <- diag(XXT) + 1
		XYpT      <- (crossprod(X, Y) + lambda * target)
		return((XYpT / lambda - lambda^(-2) * crossprod(crossprod(X, (crossprod(solve(XXT), X))), XYpT))[,1])
	}
}



.ridgeBLMtemp <- function(Y, X, lambda, target, minSuccDiff, maxIter, output){

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
	## output      : A character, "light" or "heavy" (or anything else), 
	##               that specifies whether only the regression estimate
	##               should be returned.
	## ---------------------------------------------------------------------
	## Value, if output="light":
	## The ridge estimate of the regression parameter, a numeric.
	## ---------------------------------------------------------------------
	## Authors     : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# initiate
	betas      <- rep(0, ncol(X))
	loglikPrev <- -10^(10)
	for (iter in 1:maxIter){
		# evaluate the linear predictor
		lp0 <- as.numeric(tcrossprod(betas, X))

		# assess convergence
		loglik <- sum(Y * lp0 - log(1 + exp(lp0))) + 0.5 * lambda * sum((betas - target)^2)
		if (abs(loglik - loglikPrev) < minSuccDiff){ 
			break 
		} else {
			loglikPrev <- loglik
		}
	
		# calculate the weights
		W0 <- exp(lp0) / (1 + exp(lp0))^2

		# obtain part one of the estimator
		XWZpT <- lambda * target + crossprod(X, W0 * lp0 + 
		            	           (Y - exp(lp0) / (1 + exp(lp0))))  

		# distinguish between the high- and low-dimensional cases
		if (ncol(X) >= nrow(X)){ 
			if (min(W0) > .Machine$double.eps){
				# evaluate X %*% X^t and add the reciprocal 
				# weight to its diagonal
				XXT       <- tcrossprod(X) / lambda
				diag(XXT) <- diag(XXT) + 1/W0 
		
				# now obtain the IRWLS update efficiently
				betas  <- (XWZpT / lambda - lambda^(-2) * 
				          crossprod(crossprod(X, (crossprod(solve(XXT), X))), XWZpT))[,1]
			} else {		
				# evaluate X %*% X^t and add the reciprocal 
				# weight to its diagonal
				X         <- sweep(X, 1, sqrt(W0), FUN="*")
				XXT       <- tcrossprod(X) / lambda
				diag(XXT) <- diag(XXT) + 1 
	
				# now obtain the IRWLS update efficiently
				betas  <- (XWZpT / lambda - lambda^(-2) * 
				         crossprod(crossprod(X, (crossprod(solve(XXT), X))), XWZpT))[,1]
				X      <- sweep(X, 1, sqrt(W0), FUN="/")
			}
		} else {
			# evaluate X %*% X^t and add the reciprocal weight to 
			# its diagonal
			XTX       <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
			diag(XTX) <- diag(XTX) + lambda 
	
			# now obtain the IRWLS update efficiently
			betas <- crossprod(solve(XTX), XWZpT)[,1]
		}
	}
	if (output == "light"){
		return(betas=betas)
	} else {
		return(list(betas=betas, loglik=loglik, nIter=iter))
	}
}


.ridgeBLM <- function(Y, X, lambda, target, minSuccDiff, maxIter){

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
	betas      <- rep(0, ncol(X))
	loglikPrev <- -10^(10)
	for (iter in 1:maxIter){
		# evaluate the linear predictor
		lp0 <- as.numeric(tcrossprod(betas, X))

		# assess convergence
		loglik <- sum(Y * lp0 - log(1 + exp(lp0))) + 0.5 * lambda * sum((betas - target)^2)
		if (abs(loglik - loglikPrev) < minSuccDiff){ 
			break 
		} else {
			loglikPrev <- loglik
		}
	
		# calculate the weights
		W0 <- exp(lp0) / (1 + exp(lp0))^2

		# obtain part one of the estimator
		XWZpT <- lambda * target + crossprod(X, W0 * lp0 + 
		            	           (Y - exp(lp0) / (1 + exp(lp0))))  

		# distinguish between the high- and low-dimensional cases
		if (ncol(X) >= nrow(X)){ 
			if (min(W0) > .Machine$double.eps){
				# evaluate X %*% X^t and add the reciprocal 
				# weight to its diagonal
				XXT       <- tcrossprod(X) / lambda
				diag(XXT) <- diag(XXT) + 1/W0 
		
				# now obtain the IRWLS update efficiently
				betas  <- (XWZpT / lambda - lambda^(-2) * 
				          crossprod(crossprod(X, (crossprod(solve(XXT), X))), XWZpT))[,1]
			} else {		
				# evaluate X %*% X^t and add the reciprocal 
				# weight to its diagonal
				X         <- sweep(X, 1, sqrt(W0), FUN="*")
				XXT       <- tcrossprod(X) / lambda
				diag(XXT) <- diag(XXT) + 1 
	
				# now obtain the IRWLS update efficiently
				betas  <- (XWZpT / lambda - lambda^(-2) * 
				         crossprod(crossprod(X, (crossprod(solve(XXT), X))), XWZpT))[,1]
				X      <- sweep(X, 1, sqrt(W0), FUN="/")
			}
		} else {
			# evaluate X %*% X^t and add the reciprocal weight to 
			# its diagonal
			XTX       <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
			diag(XTX) <- diag(XTX) + lambda 
	
			# now obtain the IRWLS update efficiently
			betas <- crossprod(solve(XTX), XWZpT)[,1]
		}
	}
	return(betas=betas)
}



ridgeGLM <- function(Y, X, lambda, target=rep(0, ncol(X)), model, minSuccDiff=10^(-10), maxIter=100){

	## ---------------------------------------------------------------------
	## Targeted ridge estimation of the regression parameter of the 
	## generalized linear model (GLM)
	## ---------------------------------------------------------------------
	## Arguments 
	## Y           : A numeric being the response vector.
	## X           : A matrix, the design matrix. The number of rows should
	##               match the number of elements of Y.
	## lambda      : A numeric, the ridge penalty parameter.
	## target      : A numeric towards which the estimate is shrunken.
	## model       : A character indicating which generalized linear model
	##               model instance is to be fitted.
	## minSuccDiff : A numeric, the minimum distance between two successive 
	##               estimates to be achieved.
	## maxIter     : A numeric specifying the maximum number of iterations.
	## ---------------------------------------------------------------------
	## Value:
	## The ridge estimate of the regression parameter, a numeric.
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	if (model == "linear"){
		return(.ridgeLM(Y, X, lambda, target))
	}
	if (model == "logistic"){
		return(.ridgeBLM(Y, X, lambda, target, minSuccDiff, maxIter))
	}
}




ridgeGLMmultiT <- function(Y, X, lambdas, targetMat, model, minSuccDiff=10^(-10), maxIter=100){

	## ---------------------------------------------------------------------
	## Multi-targeted ridge estimation of the regression parameter of the 
	## generalized linear model (GLM).
	## ---------------------------------------------------------------------
	## Arguments 
	## Y           : A numeric being the response vector.
	## X           : A matrix, the design matrix. The number of rows should
	##               match the number of elements of Y.
	## lambdas     : A numeric, vector of penalty parameters, one per target
	## targetMat   : A matrix with targets for the regression parameter 
	##               as columns
	## model       : A character indicating which generalized linear model
	##               model instance is to be fitted.
	## minSuccDiff : A numeric, the minimum distance between two successive 
	##               estimates to be achieved.
	## maxIter     : A numeric specifying the maximum number of iterations.
	## ---------------------------------------------------------------------
	## Value:
	## The ridge estimate of the regression parameter, a numeric.
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# create weighted average of targets
	mixedTarget <- rep(0, ncol(X))
	lambdaTotal <- sum(lambdas)
	for (k in 1:length(lambdas)){
		mixedTarget <- mixedTarget + (lambdas[k] / lambdaTotal) * targetMat[,k]
	}

	# evaluate the regression estimator
	if (model == "linear"){
		return(.ridgeLM(Y, X, lambdaTotal, mixedTarget))
	}
	if (model == "logistic"){
		return(.ridgeBLM(Y, X, lambdaTotal, mixedTarget, minSuccDiff, maxIter))
	}
}

.loglikBLM <- function(Y, X, betas){

	## ---------------------------------------------------------------------
	## Function calculates the loglikelihood of the logistic regression model
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## X       : design matrix
	## betas   : regression parameter
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the loglikelihood of the logistic regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# evaluate the linear predictor
	lp0 <- as.numeric(tcrossprod(betas, X))
	return(-sum(Y * lp0 - log(1 + exp(lp0))))
}

.loglikLM <- function(Y, X, betas){

	## ---------------------------------------------------------------------
	## Function calculates the loglikelihood of the linear regression model
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## X       : design matrix
	## betas   : regression parameter
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the loglikelihood of the linear regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	return(length(Y) * (log(2 * pi * sum((Y - X %*% betas)^2) / length(Y)) + 1) / 2)
}

.sosLM <- function(Y, X, betas){

	## ---------------------------------------------------------------------
	## Function calculates the sum-of-squares of the linear regression model
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## X       : design matrix
	## betas   : regression parameter
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the sum-of-squares of the linear regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	return(sum((Y - X %*% betas)^2))
}

.kcvLMlossPleqN <- function(lambda, Y, X, folds, target, loss="loglik"){

	## ---------------------------------------------------------------------
	## Function yields the cross-validated loss of the ridge
	## regression estimator with a nonnull target
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## X       : design matrix
	## folds   : list with folds
	## target  : shrinkage target for regression parameter
	## loss    : loss to be used in CV, either "sos" (sum-of-squares) or 
	##           "loglik" (loglikelihood)
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the linear regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# evaluate loss per fold
	cvLoss <- 0
	for (k in 1:length(folds)){
		if (loss=="sos"){
			cvLoss <- cvLoss + .sosLM(Y[folds[[k]]], 
	                                          X[folds[[k]],,drop=FALSE], 
	                                          ridgeGLM(Y[-folds[[k]]], 
		                                           X[-folds[[k]],,drop=FALSE], 
	                                                   lambda, target, "linear"))
		}
		if (loss=="loglik"){
			cvLoss <- cvLoss + .loglikLM(Y[folds[[k]]], 
	                                             X[folds[[k]],,drop=FALSE], 
	                                             ridgeGLM(Y[-folds[[k]]], 
		                                              X[-folds[[k]],,drop=FALSE], 
	                                                      lambda, target, "linear"))
		}

	}
	
	# average over the folds
	return(cvLoss / length(folds))
}

.kcvLMlossPgeqN <- function(lambda, Y, X, folds, target, loss){

	## ---------------------------------------------------------------------
	## Function yields the cross-validated loglikelihood of the ridge
	## regression estimator with a nonnull target of the logistic regression
	## model
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## X       : design matrix
	## folds   : list with folds
	## target  : shrinkage target for regression parameter
	## ...     : other arguments passed on to the 'ridgeGLM'-function
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the logistic regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# evaluate loss per fold
	cvLoss <- 0
	for (k in 1:length(folds)){
		Yhat   <- t(tcrossprod(target, X[folds[[k]],])) + 
		          tcrossprod(X[folds[[k]],], X[-folds[[k]],]) %*% 
		          solve(diag(rep(lambda, nrow(X)-length(folds[[k]]))) + 
		                X[-folds[[k]],] %*% t(X[-folds[[k]],])) %*% 
		          (matrix(Y[-folds[[k]]], ncol=1) - 
		           t(tcrossprod(target, X[-folds[[k]],])))
		if (loss == "loglik"){
			cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
		                           (log(2 * pi * sum((Y[folds[[k]]] - Yhat)^2) / 
			                   length(folds[[k]])) + 1) / 2
		}
		if (loss == "sos"){
			cvLoss <- cvLoss + sum((Y[folds[[k]]] - Yhat)^2)
		}
	}
	
	# average over the folds
	return(cvLoss / length(folds))
}



.kcvBLMlossPleqN_old <- function(lambda, Y, X, folds, target, ...){

	## ---------------------------------------------------------------------
	## Function yields the cross-validated loglikelihood of the ridge
	## regression estimator with a nonnull target of the logistic regression
	## model
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## X       : design matrix
	## folds   : list with folds
	## target  : shrinkage target for regression parameter
	## ...     : other arguments passed on to the 'ridgeGLM'-function
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the logistic regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# evaluate loss per fold
	cvLoss <- 0
	for (k in 1:length(folds)){
		cvLoss <- cvLoss + .loglikLM(Y[folds[[k]]], 
                                             X[folds[[k]],,drop=FALSE], 
                                             ridgeGLM(Y[-folds[[k]]], 
	                                             X[-folds[[k]],,drop=FALSE], 
                                                     lambda, target, "logistic", ...))
	}
	
	# average over the folds
	return(cvLoss / length(folds))
}



.kcvBLMloss <- function(lambda, Y, X, folds, target, minSuccDiff, maxIter){

	## ---------------------------------------------------------------------
	## Function yields the cross-validated loglikelihood of the ridge
	## regression estimator with a nonnull target of the logistic regression
	## model
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## X       : design matrix
	## folds   : list with folds
	## target  : shrinkage target for regression parameter
	## ...     : other arguments passed on to the 'ridgeGLM'-function
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the logistic regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# initiate
	cvLoss <- 0
	for (k in 1:length(folds)){
		# convergence criterion
		loglikPrev <- -10^(10)

		# evaluate the initial linear predictor
		lpOld <- as.numeric(tcrossprod(target, X[-folds[[k]],,drop=FALSE]))

		for (iter in 1:maxIter){
			# assess convergence
			loglik <- sum(Y[-folds[[k]]] * lpOld - log(1 + exp(lpOld))) 
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
	
			# calculate the weights
			W0 <- exp(lpOld) / (1 + exp(lpOld))^2

			# obtain the adjusted response
			WZ <- W0 * lpOld + (Y[-folds[[k]]] - exp(lpOld) / (1 + exp(lpOld)))  

			# distinguish between the high- and low-dimensional cases
			if (ncol(X) >= nrow(X[-folds[[k]],])){ 
				if (min(W0) > .Machine$double.eps){
					# evaluate X %*% X^t and add the reciprocal 
					# weight to its diagonal
					XXT       <- tcrossprod(X[-folds[[k]],,drop=FALSE]) / lambda
					diag(XXT) <- diag(XXT) + 1/W0 

					# obtain the IRWLS update of the linear predictors efficiently
					XXTZpT <- tcrossprod(solve(XXT), (1/W0) * (WZ - W0 * tcrossprod(target, X[-folds[[k]],])))
					lpAll  <- t(tcrossprod(target, X)) +
					        (1/lambda) * tcrossprod(X, X[-folds[[k]],]) %*% XXTZpT
				} else {		
					# evaluate X %*% X^t and add the reciprocal 
					# weight to its diagonal
					XXT       <- tcrossprod(sweep(X[-folds[[k]], ], 1, sqrt(W0), FUN="*")) / lambda
					diag(XXT) <- diag(XXT) + 1 

					# obtain the IRWLS update of the linear predictors efficiently
					XXTZpT <- tcrossprod(solve(XXT), sqrt(1/W0) * (WZ - W0 * tcrossprod(target, sweep(X[-folds[[k]], ], 1, sqrt(W0), FUN="*"))))
					lpAll  <- t(tcrossprod(target, X)) +
					         (1/lambda) * tcrossprod(X, sweep(X[-folds[[k]], ], 1, sqrt(W0), FUN="*")) %*% XXTZpT
				}
			} else {
				# evaluate X %*% X^t and add the reciprocal weight to 
				# its diagonal
				XWZpT     <- lambda * target + crossprod(X[-folds[[k]], ], WZ)  
				XTX       <- crossprod(sweep(X[-folds[[k]], ], 1, sqrt(W0), FUN="*"))
				diag(XTX) <- diag(XTX) + lambda 
		
				# now obtain the IRWLS update of the linear predictor efficiently
				lpAll  <- as.numeric(tcrossprod(crossprod(solve(XTX), XWZpT)[,1], X))
			}
			lpOld  <- lpAll[-folds[[k]]]
			lpNew  <- lpAll[ folds[[k]]]
		}
		cvLoss <- cvLoss + sum(Y[folds[[k]]] * lpNew - log(1 + exp(lpNew))) 
	}

	# average over the folds
	return(-cvLoss / length(folds))
}




.optPenaltyLM.kCVauto <- function(Y, X, lambdaInit, folds, target, loss){

	## ---------------------------------------------------------------------
	## Function finds the optimal penalty parameter of the targeted ridge 
	## regression estimator. The optimum is defined as the minimizer of the 
	## cross-validated loss associated with the estimator. 
	## --------------------------------------------------------------------- 
	## Arguments
	## Y          : response vector
	## X          : design matrix
	## lambdaInit : initial guess of the optimal penalty parameter
	## folds      : list with folds
	## target     : shrinkage target for regression parameter
	## loss       : loss to be used in CV, either "sos" (sum-of-squares) or 
	##              "loglik" (loglikelihood)
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the optimal cross-validated penalty parameter of the
	## ridge estimator of the linear regression model parameter
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# search by the Brent algorithm
	if (nrow(X) >= ncol(X)){
		optLambda <- optim(lambdaInit, .kcvLMlossPleqN, Y=Y, X=X, folds=folds, 
		                   target=target, loss=loss, method="Brent", 
		                   lower=10^(-10), upper=10^(10))$par
	}
	if (nrow(X) <ncol(X)){
		optLambda <- optim(lambdaInit, .kcvLMlossPgeqN, Y=Y, X=X, folds=folds, 
		                   target=target, loss=loss, method="Brent", 
		                   lower=10^(-10), upper=10^(10))$par
	}
	return(optLambda)
}



.optPenaltyBLM.kCVauto <- function(Y, X, lambdaInit, folds, target, minSuccDiff, 
                                   maxIter){

	## ---------------------------------------------------------------------
	## Function finds the optimal penalty parameter of the targeted ridge 
	## regression estimator of the logistic regression model parameter. The 
	## optimum is defined as the minimizer of the cross-validated loss 
	## associated with the estimator. 
	## --------------------------------------------------------------------- 
	## Arguments
	## Y           : response vector
	## X           : design matrix
	## lambdaInit  : initial guess of the optimal penalty parameter
	## folds       : list with folds
	## target      : shrinkage target for regression parameter
	## minSuccDiff : A numeric, the minimum distance between two successive 
	##               estimates to be achieved.
	## maxIter     : A numeric specifying the maximum number of iterations.
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the optimal cross-validated penalty parameter of the
	## ridge estimator of the logistic regression model parameter
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# search by the Brent algorithm
	optLambda <- optim(lambdaInit, .kcvBLMloss, Y=Y, X=X, folds=folds, 
	        	           target=target, method="Brent", 
	        	           minSuccDiff=minSuccDiff, maxIter=maxIter, 
	        	           lower=10^(-10), upper=10^(10))$par
	return(optLambda)
}


optPenaltyGLM.kCVauto <- function(Y, X, lambdaInit, fold=nrow(X), stratified=TRUE,
                                  target=rep(0, ncol(X)), model, loss="loglik",
                                  minSuccDiff=10^(-5), maxIter=100){

	## ---------------------------------------------------------------------
	## Function finds the optimal penalty parameter of the targeted ridge 
	## regression estimator of the generalized linear model parameter. The 
	## optimum is defined as the minimizer of the cross-validated loss 
	## associated with the estimator. 
	## --------------------------------------------------------------------- 
	## Arguments
	## Y           : response vector
	## X           : design matrix
	## lambdaInit  : initial guess of the optimal penalty parameter
	## fold        : list with folds
	## stratified  : boolean, should the data by stratified such that 
	##               range of the response is comparable among folds
	## target      : shrinkage target for regression parameter
	## model       : A character indicating which generalized linear model
	##               model instance is to be fitted.
	## loss        : loss to be used in CV, either "sos" (sum-of-squares) or 
	##              "loglik" (loglikelihood)
	## minSuccDiff : A numeric, the minimum distance between two successive 
	##               estimates to be achieved.
	## maxIter     : A numeric specifying the maximum number of iterations.
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the optimal cross-validated penalty parameter of the
	## ridge estimator of the generalized linear regression model parameter
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# create folds
	if (!stratified){
	        fold    <- max(min(ceiling(fold), nrow(X)), 2)
        	fold    <- rep(1:fold, ceiling(nrow(X)/fold))[1:nrow(X)]
        	shuffle <- sample(1:nrow(X), nrow(X))
        	folds   <- split(shuffle, as.factor(fold))
	}
	if (stratified){
		if (model == "linear"){	
			idSorted <- c(order(Y), rep(NA, ceiling(length(Y) / fold) * fold - length(Y)))
			folds    <- matrix(idSorted, nrow=fold)
			folds    <- apply(folds, 2, function(Z){ Z[sample(1:length(Z), length(Z))] })
			folds    <- split(as.numeric(folds), 
			                  as.factor(matrix(rep(1:fold, ncol(folds)), nrow=fold)))
			folds    <- lapply(folds, function(Z){ Z[!is.na(Z)] })
		}
		if (model == "logistic"){		
			idSorted <- c(which(Y == 0)[sample(1:sum(Y==0), sum(Y==0))], 
			              which(Y == 1)[sample(1:sum(Y==1), sum(Y==1))])
			idSorted <- c(idSorted, rep(NA, ceiling(length(Y) / fold) * fold - length(Y)))
			folds    <- matrix(idSorted, nrow=fold)
			folds   <- split(as.numeric(folds), 
			                 as.factor(matrix(rep(1:fold, ncol(folds)), nrow=fold)))
			folds   <- lapply(folds, function(Z){ Z[!is.na(Z)] })
		}		
	}

	# search by the Brent algorithm
	if (model == "linear"){
		return(.optPenaltyLM.kCVauto(Y, X, lambdaInit, folds, target, loss))
	}
	if (model == "logistic"){
		return(.optPenaltyBLM.kCVauto(Y, X, lambdaInit, folds, target, 
		                              minSuccDiff, maxIter))
	}
}


###############################################################################

.kcvLMlossPleqN_multiT <- function(lambdas, Y, X, folds, targetMat, loss="loglik"){

	## ---------------------------------------------------------------------
	## Function yields the cross-validated loss of the ridge
	## regression estimator with multiple nonnull targets
	## --------------------------------------------------------------------- 
	## Arguments
	## lambdas   : A numeric, vector of penalty parameters, one per target
	## Y         : response vector
	## X         : design matrix
	## folds     : list with folds
	## targetMat : A matrix with targets for the regression parameter 
	##               as columns
	## loss      : loss to be used in CV, either "sos" (sum-of-squares) or 
	##             "loglik" (loglikelihood)
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the linear regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# create weighted average of targets
	mixedTarget <- rep(0, ncol(X))
	lambdaTotal <- sum(lambdas)
	for (k in 1:length(lambdas)){
		mixedTarget <- mixedTarget + (lambdas[k] / lambdaTotal) * targetMat[,k]
	}

	# evaluate loss per fold
	cvLoss <- 0
	for (k in 1:length(folds)){
		if (loss=="sos"){
			cvLoss <- cvLoss + .sosLM(Y[folds[[k]]], 
	                                          X[folds[[k]],,drop=FALSE], 
	                                          ridgeGLM(Y[-folds[[k]]], 
		                                          X[-folds[[k]],,drop=FALSE], 
	                                                  lambdaTotal, mixedTarget, "linear"))
		}
		if (loss=="loglik"){
			cvLoss <- cvLoss + .loglikLM(Y[folds[[k]]], 
	                                             X[folds[[k]],,drop=FALSE], 
	                                             ridgeGLM(Y[-folds[[k]]], 
		                                             X[-folds[[k]],,drop=FALSE], 
	                                                     lambdaTotal, mixedTarget, "linear"))
		}

	}
	
	# average over the folds
	return(cvLoss / length(folds))
}



.kcvLMlossPgeqN_multiT <- function(lambdas, Y, X, folds, targetMat, loss="loglik"){

	## ---------------------------------------------------------------------
	## Function yields the cross-validated loss of the ridge
	## regression estimator with multiple nonnull targets
	## --------------------------------------------------------------------- 
	## Arguments
	## lambdas   : A numeric, vector of penalty parameters, one per target
	## Y         : response vector
	## X         : design matrix
	## folds     : list with folds
	## targetMat : A matrix with targets for the regression parameter 
	##               as columns
	## loss      : loss to be used in CV, either "sos" (sum-of-squares) or 
	##             "loglik" (loglikelihood)
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the linear regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# create weighted average of targets
	mixedTarget <- rep(0, ncol(X))
	lambdaTotal <- sum(lambdas)
	for (k in 1:length(lambdas)){
		mixedTarget <- mixedTarget + (lambdas[k] / lambdaTotal) * targetMat[,k]
	}

	# evaluate loss per fold
	cvLoss <- 0
	for (k in 1:length(folds)){
		Yhat   <- t(tcrossprod(mixedTarget, X[folds[[k]],])) + 
		          tcrossprod(X[folds[[k]],], X[-folds[[k]],]) %*% 
		          solve(diag(rep(lambdaTotal, nrow(X)-length(folds[[k]]))) + 
		                X[-folds[[k]],] %*% t(X[-folds[[k]],])) %*% 
		          (matrix(Y[-folds[[k]]], ncol=1) - 
		           t(tcrossprod(mixedTarget, X[-folds[[k]],])))
		if (loss == "loglik"){
			cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
		                           (log(2 * pi * sum((Y[folds[[k]]] - Yhat)^2) / 
			                   length(folds[[k]])) + 1) / 2
		}
		if (loss == "sos"){
			cvLoss <- cvLoss + sum((Y[folds[[k]]] - Yhat)^2)
		}
	}
	
	# average over the folds
	return(cvLoss / length(folds))
}



.kcvBLMloss_multiT <- function(lambdas, Y, X, folds, targetMat, minSuccDiff, maxIter){

	## ---------------------------------------------------------------------
	## Function yields the cross-validated loss of the ridge logistic
	## regression estimator with multiple nonnull targets
	## --------------------------------------------------------------------- 
	## Arguments
	## lambdas   : A numeric, vector of penalty parameters, one per target
	## Y         : response vector
	## X         : design matrix
	## folds     : list with folds
	## targetMat : A matrix with targets for the regression parameter 
	##             as columns
	## ...       : other arguments passed on to the 'ridgeGLM'-function
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the logistic regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# create weighted average of targets
	mixedTarget <- rep(0, ncol(X))
	lambdaTotal <- sum(lambdas)
	for (k in 1:length(lambdas)){
		mixedTarget <- mixedTarget + (lambdas[k] / lambdaTotal) * targetMat[,k]
	}

	# calculate cross-validated loss
	return(.kcvBLMloss(lambdaTotal, Y, X, folds, mixedTarget, minSuccDiff, maxIter))
}





.optPenaltyLMmultT.kCVauto <- function(Y, X, lambdas, folds, 
                                      targetMat, loss=loss){

	## --------------------------------------------------------------------
	## Function finds the optimal penalty parameter of the targeted ridge 
	## regression estimator of the logistic regression model parameter. The 
	## optimum is defined as the minimizer of the cross-validated loss 
	## associated with the estimator. 
	## --------------------------------------------------------------------- 
	## Arguments
	## Y           : response vector
	## X           : design matrix
	## lambdaInit  : vector of penalty parameters
	## folds       : list with folds
	## targetMat   : A matrix with targets for the regression parameter 
	##               as columns
	## loss        : loss to be used in CV, either "sos" (sum-of-squares) or 
	##              "loglik" (loglikelihood). For model="linear" only.
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the optimal cross-validated penalty parameter of the
	## ridge estimator of the logistic regression model parameter
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# search for optimal parameters
	if (nrow(X) >= ncol(X)){
		lambdasOpt <- constrOptim(lambdas, 
					  .kcvLMlossPleqN_multiT, 
					  grad=NULL, 
				          Y=Y, 
				          X=X, 
					  targetMat=targetMat, 
					  folds=folds, 
	        	                  loss=loss,
					  ui=diag(length(lambdas)), 
					  ci=rep(0, length(lambdas)))$par
	}
	if (nrow(X) < ncol(X)){
		lambdasOpt <- constrOptim(lambdas, 
					  .kcvLMlossPgeqN_multiT, 
					  grad=NULL, 
				          Y=Y, 
				          X=X, 
					  targetMat=targetMat, 
					  folds=folds, 
	        	                  loss=loss,
					  ui=diag(length(lambdas)), 
					  ci=rep(0, length(lambdas)))$par
	}

	return(lambdasOpt)
}

.optPenaltyBLMmultT.kCVauto <- function(Y, X, lambdas, folds, targetMat,
                                        minSuccDiff, maxIter){

	## --------------------------------------------------------------------
	## Function finds the optimal penalty parameter of the targeted ridge 
	## regression estimator of the logistic regression model parameter. The 
	## optimum is defined as the minimizer of the cross-validated loss 
	## associated with the estimator. 
	## --------------------------------------------------------------------- 
	## Arguments
	## Y           : response vector
	## X           : design matrix
	## lambdaInit  : vector of penalty parameters
	## folds       : list with folds
	## targetMat   : A matrix with targets for the regression parameter 
	##               as columns
	## minSuccDiff : A numeric, the minimum distance between two successive 
	##               estimates to be achieved.
	## maxIter     : A numeric specifying the maximum number of iterations.
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the optimal cross-validated penalty parameter of the
	## ridge estimator of the logistic regression model parameter
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# search for optimal parameters
	lambdasOpt <- constrOptim(lambdas, 
				  .kcvBLMloss_multiT, 
				  grad=NULL, 
			          Y=Y, 
			          X=X, 
				  targetMat=targetMat, 
				  folds=folds, 
	                          minSuccDiff=minSuccDiff,
	                          maxIter=maxIter,                          
				  ui=diag(length(lambdas)), 
				  ci=rep(0, length(lambdas)))$par
	return(lambdasOpt)
}



optPenaltyGLMmultiT.kCVauto <- function(Y, X, lambdasInit, fold=nrow(X), 
                                        stratified=TRUE,
                                        targetMat, model, loss="loglik",
                                        minSuccDiff=10^(-5), maxIter=100){

	## --------------------------------------------------------------------
	## Function finds the optimal penalty parameter of the targeted ridge 
	## regression estimator of the generalized linear model parameter. The 
	## optimum is defined as the minimizer of the cross-validated loss 
	## associated with the estimator. 
	## --------------------------------------------------------------------- 
	## Arguments
	## Y           : response vector
	## X           : design matrix
	## lambdaInit  : initial guess of the optimal penalty parameter
	## fold        : number of cross-validation folds
	## stratified  : boolean, should the data by stratified such that 
	##               range of the response is comparable among folds
	## target      : shrinkage target for regression parameter
	## model       : A character indicating which generalized linear model
	##               model instance is to be fitted.
	## loss        : loss to be used in CV, either "sos" (sum-of-squares) or 
	##              "loglik" (loglikelihood). For model="linear" only.
	## minSuccDiff : A numeric, the minimum distance between two successive 
	##               estimates to be achieved.
	## maxIter     : A numeric specifying the maximum number of iterations.
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the optimal cross-validated penalty parameter of the
	## ridge estimator of the generalized linear regression model parameter
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------


	# create folds
	if (!stratified){
	        fold    <- max(min(ceiling(fold), nrow(X)), 2)
        	fold    <- rep(1:fold, ceiling(nrow(X)/fold))[1:nrow(X)]
        	shuffle <- sample(1:nrow(X), nrow(X))
        	folds   <- split(shuffle, as.factor(fold))
	}
	if (stratified){
		if (model == "linear"){	
			idSorted <- c(order(Y), rep(NA, ceiling(length(Y) / fold) * fold - length(Y)))
			folds    <- matrix(idSorted, nrow=fold)
			folds    <- apply(folds, 2, function(Z){ Z[sample(1:length(Z), length(Z))] })
			folds    <- split(as.numeric(folds), 
			                  as.factor(matrix(rep(1:fold, ncol(folds)), nrow=fold)))
			folds    <- lapply(folds, function(Z){ Z[!is.na(Z)] })
		}
		if (model == "logistic"){		
			idSorted <- c(which(Y == 0)[sample(1:sum(Y==0), sum(Y==0))], 
			              which(Y == 1)[sample(1:sum(Y==1), sum(Y==1))])
			idSorted <- c(idSorted, rep(NA, ceiling(length(Y) / fold) * fold - length(Y)))
			folds    <- matrix(idSorted, nrow=fold)
			folds   <- split(as.numeric(folds), 
			                 as.factor(matrix(rep(1:fold, ncol(folds)), nrow=fold)))
			folds   <- lapply(folds, function(Z){ Z[!is.na(Z)] })
		}		
	}

	# search by the Brent algorithm
	if (model == "linear"){
		return(.optPenaltyLMmultT.kCVauto(Y, X, lambdasInit, folds, targetMat, loss))
	}
	if (model == "logistic"){
		return(.optPenaltyBLMmultT.kCVauto(Y, X, lambdasInit, folds, targetMat, 
		                                   minSuccDiff, maxIter))
	}
}



