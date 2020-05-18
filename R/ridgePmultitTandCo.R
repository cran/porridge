.kcvlTmult <- function(lambda, Y, targetList, folds){
	## ---------------------------------------------------------------------
	## Function that calculates a cross-validated negative log-likelihood
	## score for the specified penalty parameter vector.
	## ---------------------------------------------------------------------
	## Arguments 
	## - lambda     : vector of penalty parameter values, one per target
	## - Y          : (raw) Data matrix, variables in columns
	## - targetList : target (precision terms) for Type I estimators,
	##                default = default.target(covML(Y))
	## - folds      : list with the cross-validation sample splits
	## ---------------------------------------------------------------------
	## Value
	## A numeric, the cross-validated negative loglikelihood
	## ---------------------------------------------------------------------
	## Authors  : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# create cv-storage variable
	cvLL <- 0
	for (f in 1:length(folds)){
		S    <- crossprod(Y[-folds[[f]], , drop=FALSE]) / nrow(Y[-folds[[f]], , drop=FALSE])
		Phat <- .armaRidgePmultiT(S, lambda, targetList)
		cvLL <- cvLL -determinant(Phat)$modulus[1] + sum(crossprod(Y[folds[[f]], , drop=FALSE])*Phat) / length(folds[[f]])
	}
	return(cvLL / length(folds))
}



ridgePmultiT <- function(S, lambda, targetList){
	## ---------------------------------------------------------------------
	## Ridge precision matrix estimation with multiple targets
	## ---------------------------------------------------------------------
	## Arguments 
	## S		: a sample covariance matrix
        ## lambda       : vector of penalty parameters
        ## targetList   : target matrices towards the precision matrices are potentially shrunken
	## ---------------------------------------------------------------------
	## Value
	## A matrix, a ridge precision matrix estimate
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(S, "matrix")){ stop("Input (S) should be a matrix") }
	if (as.character(class(lambda)) != "numeric"){ stop("Input (lambda) is of wrong class.")  }
	if (length(lambda) <= 1){ stop("Input (lambda) is of wrong length.") }
	if (any(is.na(lambda))){ stop("Input (lambda) contains a nonpositive number.") }
	if (any(lambda <= 0)){ stop("Input (lambda) contains a nonpositive number.") }
	if (!is(targetList, "list")){ stop("Input (targetList) is of wrong class.") }
	if (length(targetList) != length(lambda)){ stop("Arguments (lambda & targetList) do not match") }
	for (k in 1:length(targetList)){
		if (!is(targetList[[k]], "matrix")){ 
			stop("Elements of input (targetList) of wrong class") 
		}		
		if (dim(S)[2] != nrow(targetList[[k]])){ 
			stop("Dimensions of soms elements of input (targetList) do not match that of other input (S).") 
		}
		if (dim(S)[2] != ncol(targetList[[k]])){ 
			stop("Dimensions of some elements of input (targetList) do not match that of other input (S).") 
		}
	}

	# estimate and return ridge precision matrix estimate
	return(.armaRidgePmultiT(S, lambda, targetList))
}


optPenaltyPmultiT.kCVauto <- function(Y, 
                                      lambdaMin, 
                                      lambdaMax,
                                      lambdaInit=(lambdaMin+lambdaMax)/2,
                                      fold=nrow(Y), 
                                      targetList){

	## ---------------------------------------------------------------------
	## Cross-validation for choosing the penalty parameter of 
	## the ridge precision matrix estimator with multiple targets
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
        ## lambdaMin    : vector of minimum values of the penalty parameters
        ## lambdaMax    : vector of maximum value of the penalty parameters
        ## lambdaInit   : initial value of penalty parameters for starting optimization
        ## fold         : cross-validation fold, default gives LOOCV
        ## targetList   : target matrix towards the precision matrices are potentially shrunken
	## ---------------------------------------------------------------------
	## Value
	## cvLL		: negative cross-validated loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){ stop("Input (Y) should be a matrix") }
	if (as.character(class(lambdaInit)) != "numeric"){ stop("Input (lambdaInit) is of wrong class.")  }
	if (length(lambdaInit) <= 1){ stop("Input (lambdaInit) is of wrong length.") }
	if (any(is.na(lambdaInit))){ stop("Input (lambdaInit) contains a nonpositive number.") }
	if (any(lambdaInit <= 0)){ stop("Input (lambdaInit) contains a nonpositive number.") }
	if (as.character(class(lambdaMin)) != "numeric"){ stop("Input (lambdaMin) is of wrong class.")  }
	if (length(lambdaMin) <= 1){ stop("Input (lambdaMin) is of wrong length.") }
	if (any(is.na(lambdaMin))){ stop("Input (lambdaMin) contains a nonpositive number.") }
	if (any(lambdaMin <= 0)){ stop("Input (lambdaMin) contains a nonpositive number.") }
	if (as.character(class(lambdaMax)) != "numeric"){ stop("Input (lambdaMax) is of wrong class.")  }
	if (length(lambdaMax) <= 1){ stop("Input (lambdaMax) is of wrong length.") }
	if (any(is.na(lambdaMax))){ stop("Input (lambdaMax) contains a nonpositive number.") }
	if (any(lambdaMax <= 0)){ stop("Input (lambdaMax) contains a nonpositive number.") }
	if (class(fold) != "numeric" & class(fold) != "integer"){ stop("Input (fold) is of wrong class") }
	if ((fold <=  1) | (fold > nrow(Y))){ stop("Input (fold) out of range") }
	if (!is(targetList, "list")){ stop("Input (targetList) is of wrong class.") }
	if (length(targetList) != length(lambdaInit)){ stop("Arguments (lambdaInit & targetList) do not match") }
	for (k in 1:length(targetList)){
		if (!is(targetList[[k]], "matrix")){ 
			stop("Elements of input (targetList) of wrong class") 
		}		
		if (dim(Y)[2] != nrow(targetList[[k]])){ 
			stop("Dimensions of soms elements of input (targetList) do not match that of other input (Y).") 
		}
		if (dim(Y)[2] != ncol(targetList[[k]])){ 
			stop("Dimensions of some elements of input (targetList) do not match that of other input (Y).") 
		}
	}


	# make k-folds as list
	fold    <- max(min(ceiling(fold), nrow(Y)), 2)
	fold    <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
	shuffle <- sample(1:nrow(Y), nrow(Y))
	folds   <- split(shuffle, as.factor(fold))

	# determine optimal value of ridge penalty parameter
	optLambda <- optim(lambdaInit, 
	                   .kcvlTmult, 
	                   method="L-BFGS-B", 
	                   lower=lambdaMin,
	                   upper=lambdaMax, 
	                   Y=Y, 
	                   folds=folds,
	                   targetList=targetList)$par

	# return
	return(optLambda)
}



