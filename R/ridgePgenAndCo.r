ridgePgen <- function(S, lambda, target, nInit=100, minSuccDiff=10^(-10)){
	if (nrow(S) != ncol(S)){ stop("S should be a square matrix") }
	if (nrow(lambda) != ncol(lambda)){ stop("lambda should be a square matrix") }
	if (nrow(target) != ncol(target)){ stop("target should be a square matrix") }
	if (!isSymmetric(S)){ stop("S should be a symmetric matrix") }
	if (!isSymmetric(lambda)){ stop("lambda should be a symmetric matrix") }
	if (!isSymmetric(target)){ stop("target should be a symmetric matrix") }
	if (as.character(class(minSuccDiff)) != "numeric"){ stop("Input (llDiff) is of wrong class.") }
	if (length(minSuccDiff) != 1){ stop("Input (minSuccDiff) is of wrong length.") }
	if (is.na(minSuccDiff)){ stop("Input (minSuccDiff) is not a positive number.") }
	if (minSuccDiff <= 0){ stop("Input (minSuccDiff) is not a positive number.") }
	if (as.character(class(nInit)) != "numeric" & as.character(class(nInit)) !=  "integer"){ stop("Input (nInit) is of wrong class.") }
	if (length(nInit) != 1){ stop("Input (nInit) is of wrong length.") }
	if (is.na(nInit)){ stop("Input (nInit) is not a positive integer.") }
	if (nInit < 0){ stop("Input (nInit) is not a positive integer.") }
	if (nInit%%1 != 0){ stop("Input (nInit) is not a positive integer.") }
 	
	# return results
	return(.armaRidgePgen(S, lambda, target, nInit, minSuccDiff))
}



ridgePgen.kCV <- function(lambda,
                          Y, 
                          fold=nrow(Y),
                          target, 
                          nInit=100, 
                          minSuccDiff=10^(-5)){ 
	## ---------------------------------------------------------------------
	## Function that calculates the k-fold cross-validated loglikelihood 
	## ---------------------------------------------------------------------
	## Arguments:    
	## -> lambda         : penalty matrix.  
	## -> Y              : data matrix with rows and columns containing samples
	##                     and variables, respectively
	## -> target         : target matrix for the precision 
	## -> folds          : cross-validation sample splits
	## -> nInit	         : maximum number of iterations. Passed on to the 
	##                     ridgePgen-function.
	## -> minSuccDiff    : minimum successive difference (in the penalized 
	##                     loglikehood) to be achieved. Passed on to the 
	##                     ridgePgen-function.
	## ---------------------------------------------------------------------
	## Value:
	## -> cvLoglik       : cross-validated loglikelihood.
	## ---------------------------------------------------------------------
	## Authors           : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (!is(lambda, "matrix")){  
		stop("Input (lambdas) is of wrong class") 
	}
	if (nrow(lambda) != ncol(lambda)){ 
	    stop("lambda should be a square matrix") 
	}
	if (!isSymmetric(lambda)){ 
	    stop("lambda should be a symmetric matrix") 
	}
	if (any(lambda <= 0)){ 
		stop("Input (lambda) must be positive") 
	}
	if (class(fold) != "numeric" & class(fold) != "integer"){
		stop("Input (fold) is of wrong class") 
	}
	if ((fold <=  1) | (fold > nrow(Y))){ 
		stop("Input (fold) out of range") 
	}
	if (!isSymmetric(target)){ 
	    stop("target should be a symmetric matrix") 
	}
	if (as.character(class(minSuccDiff)) != "numeric"){ 
	        stop("Input (minSuccDiff) is of wrong class.") 
	}
	if (length(minSuccDiff) != 1){ 
	    stop("Input (minSuccDiff) is of wrong length.") 
    }
	if (is.na(minSuccDiff)){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (minSuccDiff <= 0){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (as.character(class(nInit)) != "numeric" & as.character(class(nInit)) !=  "integer"){ 
	    stop("Input (nInit) is of wrong class.") 
	}
	if (length(nInit) != 1){ 
	    stop("Input (nInit) is of wrong length.") 
	}
	if (is.na(nInit)){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit < 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit%%1 != 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}

	# make k-folds as list
	fold         <- max(min(ceiling(fold), nrow(Y)), 2)
	fold         <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
	shuffle      <- sample(1:nrow(Y), nrow(Y))
	folds        <- split(shuffle, as.factor(fold))
	names(folds) <- NULL

	# optimize CV-loss
	cvLoglik <- .armaKcvlossR(lambda, 
                              Y=Y, 
                              target=target, 
                              folds=folds,                               
                              nInit=nInit, 
                              minSuccDiff=minSuccDiff)

	# return
	return(cvLoglik)
}


ridgePgen.kCV.banded <- function(lambda,
                                 Y, 
                                 fold=nrow(Y),
                                 target, 
                                 zeros=matrix(nrow=0, ncol=2), 
                                 penalize.diag=TRUE, 
                                 nInit=100, 
                                 minSuccDiff=10^(-5)){ 
	## ---------------------------------------------------------------------
	## Calculation of the k-fold cross-validated negative (!) loglikelihood, 
	## with a penalty structure that encourages a banded precision.
	## ---------------------------------------------------------------------
	## Arguments:    
	## -> lambda    : numeric with penalty parameter value
	## -> Y              : data matrix with rows and columns containing samples
	##                     and variables, respectively
	## -> target         : target matrix for the precision 
	## -> folds          : cross-validation sample splits
	## -> zerosR         : row-index of zero precision elements.
	## -> zerosC         : column-index of zero precision elements.
	## -> penalize_diag  : logical indicating whether the diagonal should 
	##                     be penalized    
	## -> nInit	         : maximum number of iterations. Passed on to the 
	##                     ridgePgen-function.
	## -> minSuccDiff    : minimum successive difference (in the penalized 
	##                     loglikehood) to be achieved. Passed on to the 
	##                     ridgePgen-function.
	## ---------------------------------------------------------------------
	## Value:
	## -> cvLoglik       : cross-validated negative loglikelihood.
	## ---------------------------------------------------------------------
	## Authors           : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(lambda) != "numeric"){ 
		stop("Input (lambda) is of wrong class") 
	}
	if (length(lambda) < 1){ 
		stop("Input (lambda) must be a numeric") 
	}
	if (any(lambda <= 0)){ 
		stop("Input (lambda) must be positive") 
	}
	if (class(fold) != "numeric" & class(fold) != "integer"){
		stop("Input (fold) is of wrong class") 
	}
	if ((fold <=  1) | (fold > nrow(Y))){ 
		stop("Input (fold) out of range") 
	}
	if (!isSymmetric(target)){ 
	    stop("target should be a symmetric matrix") 
	}
	if (as.character(class(minSuccDiff)) != "numeric"){ 
	        stop("Input (minSuccDiff) is of wrong class.") 
	}
	if (length(minSuccDiff) != 1){ 
	    stop("Input (minSuccDiff) is of wrong length.") 
    }
	if (is.na(minSuccDiff)){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (minSuccDiff <= 0){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (as.character(class(nInit)) != "numeric" & as.character(class(nInit)) !=  "integer"){ 
	    stop("Input (nInit) is of wrong class.") 
	}
	if (length(nInit) != 1){ 
	    stop("Input (nInit) is of wrong length.") 
	}
	if (is.na(nInit)){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit < 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit%%1 != 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (!is(zeros, "matrix")){  
	    stop("Input (zeros) is of wrong class.") 
	}    
	if(ncol(zeros) != 2){ 
	        stop("Wrong dimensions of the (zeros).") 
	} 
	if (as.character(class(penalize.diag)) != "logical"){ 
	    stop("Input (penalize.diag) is of wrong class.") 
	}    

	# make k-folds as list
	fold         <- max(min(ceiling(fold), nrow(Y)), 2)
	fold         <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
	shuffle      <- sample(1:nrow(Y), nrow(Y))
	folds        <- split(shuffle, as.factor(fold))
	names(folds) <- NULL

	# optimize CV-loss
	cvLoglik <- .armaKCVlossR_banded(lambda, 
                                     Y=Y, 
                                     folds=folds, 
                                     target=target, 
                                     penalize_diag=penalize.diag,
                                     zerosR=zeros[,1]-1, 
                                     zerosC=zeros[,2]-1, 
                                     nInit=nInit, 
                                     minSuccDiff=minSuccDiff)

	# return
	return(cvLoglik)
}



ridgePgen.kCV.groups <- function(lambdaGrps,
                                 Y, 
                                 fold=nrow(Y),
                                 groups,
                                 target, 
                                 zeros=matrix(nrow=0, ncol=2), 
                                 penalize.diag=TRUE, 
                                 nInit=100, 
                                 minSuccDiff=10^(-5)){ 
	## ---------------------------------------------------------------------
	## Calculation of the k-fold cross-validated negative (!) loglikelihood, 
	## assuming that variates are grouped and penalized group-wise.
	## ---------------------------------------------------------------------
	## Arguments:    
	## -> lambdaGrps    : vector with penalty parameter values, one per group.
	##                     values should be specified in the same order as the 
	##                     first appearance of a group representative.  
	## -> Y              : data matrix with rows and columns containing samples
	##                     and variables, respectively
	## -> target         : target matrix for the precision 
	## -> folds          : cross-validation sample splits
	## -> groups         : vector indicating to which group a variate belongs.
	## -> zerosR         : row-index of zero precision elements.
	## -> zerosC         : column-index of zero precision elements.
	## -> penalize_diag  : logical indicating whether the diagonal should 
	##                     be penalized    
	## -> nInit	     : maximum number of iterations. Passed on to the 
	##                     ridgePgen-function.
	## -> minSuccDiff    : minimum successive difference (in the penalized 
	##                     loglikehood) to be achieved. Passed on to the 
	##                     ridgePgen-function.
	## ---------------------------------------------------------------------
	## Value:
	## -> cvLoglik       : cross-validated negative loglikelihood.
	## ---------------------------------------------------------------------
	## Authors           : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(lambdaGrps) != "numeric"){ 
		stop("Input (lambdaGrps) is of wrong class") 
	}
	if (length(lambdaGrps) < 1){ 
		stop("Input (lambdaGrps) must be a numeric") 
	}
	if (any(lambdaGrps <= 0)){ 
		stop("Input (lambdaGrps) must be positive") 
	}
	if (class(fold) != "numeric" & class(fold) != "integer"){
		stop("Input (fold) is of wrong class") 
	}
	if ((fold <=  1) | (fold > nrow(Y))){ 
		stop("Input (fold) out of range") 
	}
	if (class(groups) != "numeric" & class(groups) != "integer" & class(groups) != "character"){
		stop("Input (groups) is of wrong class") 
	}
	if ((length(groups) <=  1) | (length(groups) > ncol(Y))){ 
		stop("Input (groups) out of range") 
	}
	if (!isSymmetric(target)){ 
	    stop("target should be a symmetric matrix") 
	}
	if (as.character(class(minSuccDiff)) != "numeric"){ 
	        stop("Input (minSuccDiff) is of wrong class.") 
	}
	if (length(minSuccDiff) != 1){ 
	    stop("Input (minSuccDiff) is of wrong length.") 
    }
	if (is.na(minSuccDiff)){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (minSuccDiff <= 0){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (as.character(class(nInit)) != "numeric" & as.character(class(nInit)) !=  "integer"){ 
	    stop("Input (nInit) is of wrong class.") 
	}
	if (length(nInit) != 1){ 
	    stop("Input (nInit) is of wrong length.") 
	}
	if (is.na(nInit)){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit < 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit%%1 != 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (!is(zeros, "matrix")){  
	    stop("Input (zeros) is of wrong class.") 
	}    
	if(ncol(zeros) != 2){ 
	        stop("Wrong dimensions of the (zeros).") 
	} 
	if (as.character(class(penalize.diag)) != "logical"){ 
	    stop("Input (penalize.diag) is of wrong class.") 
	}    

	# make k-folds as list
	fold         <- max(min(ceiling(fold), nrow(Y)), 2)
	fold         <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
	shuffle      <- sample(1:nrow(Y), nrow(Y))
	folds        <- split(shuffle, as.factor(fold))
	names(folds) <- NULL

	# optimize CV-loss
	cvLoglik <- .armaKCVlossR_groups(lambdaGrps, 
                                     Y=Y, 
                                     folds=folds, 
                                     groups=groups, 
                                     target=target, 
                                     penalize_diag=penalize.diag,
                                     zerosR=zeros[,1]-1, 
                                     zerosC=zeros[,2]-1, 
                                     nInit=nInit, 
                                     minSuccDiff=minSuccDiff)

	# return
	return(cvLoglik)
}



optPenaltyPgen.kCVauto.banded <- function(Y, 
                                      lambdaMin, 
                                      lambdaMax, 
                                      lambdaInit=(lambdaMin + lambdaMax)/2,
                                      fold=nrow(Y),
                                      target, 
                                      zeros=matrix(nrow=0, ncol=2), 
                                      penalize.diag=TRUE, 
                                      nInit=100, 
                                      minSuccDiff=10^(-5)){ 
	## ---------------------------------------------------------------------
	## Function that determines the optimal penalty parameters through 
	## maximization of the k-fold cross-validated log-likelihood score, 
	## assuming that variates are grouped and penalized group-wise.
	## ---------------------------------------------------------------------
	## Arguments:    
	## -> lambdas        : vector with penalty parameter values, one per group.
	##                     values should be specified in the same order as the 
	##                     first appearance of a group representative.  
	## -> Y              : data matrix with rows and columns containing samples
	##                     and variables, respectively
	## -> target         : target matrix for the precision 
	## -> folds          : cross-validation sample splits
	## -> zerosR         : row-index of zero precision elements.
	## -> zerosC         : column-index of zero precision elements.
	## -> penalize_diag : logical indicating whether the diagonal should 
	##                     be penalized    
	## -> nInit	         : maximum number of iterations. Passed on to the 
	##                     ridgePgen-function.
	## -> minSuccDiff    : minimum successive difference (in the penalized 
	##                     loglikehood) to be achieved. Passed on to the 
	##                     ridgePgen-function.
	## ---------------------------------------------------------------------
	## Value:
	## -> optLambdas     : vector with optimal parameter values.
	## ---------------------------------------------------------------------
	## Authors           : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(lambdaMin) != "numeric"){ 
		stop("Input (lambdaMin) is of wrong class") 
	}
	if (length(lambdaMin) < 1){ 
		stop("Input (lambdaMin) must be a scalar") 
	}
	if (any(lambdaMin <= 0)){ 
		stop("Input (lambdaMin) must be positive") 
	}
	if (class(lambdaMax) != "numeric"){ 
		stop("Input (lambdaMax) is of wrong class") 
	}
	if (length(lambdaMax) < 1){ 
		stop("Input (lambdaMax) must be a scalar") 
	}
	if (any(lambdaMax <= lambdaMin)){ 
		stop("Input (lambdaMax) must be larger than lambdaMin") 
	}
	if (class(lambdaInit) != "numeric"){ 
		stop("Input (lambdaInit) is of wrong class") 
	}
	if (length(lambdaInit) < 1){ 
		stop("Input (lambdaInit) must be a scalar") 
	}
	if (any(lambdaInit <= lambdaMin)){ 
		stop("Input (lambdaInit) must be larger than lambdaMin") 
	}
	if (any(lambdaInit > lambdaMax)){ 
		stop("Input (lambdaInit) must be smaller than lambdaMax") 
	}
	if (class(fold) != "numeric" & class(fold) != "integer"){
		stop("Input (fold) is of wrong class") 
	}
	if ((fold <=  1) | (fold > nrow(Y))){ 
		stop("Input (fold) out of range") 
	}
	if (!isSymmetric(target)){ 
	    stop("target should be a symmetric matrix") 
	}
	if (as.character(class(minSuccDiff)) != "numeric"){ 
	        stop("Input (minSuccDiff) is of wrong class.") 
	}
	if (length(minSuccDiff) != 1){ 
	    stop("Input (minSuccDiff) is of wrong length.") 
    }
	if (is.na(minSuccDiff)){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (minSuccDiff <= 0){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (as.character(class(nInit)) != "numeric" & as.character(class(nInit)) !=  "integer"){ 
	    stop("Input (nInit) is of wrong class.") 
	}
	if (length(nInit) != 1){ 
	    stop("Input (nInit) is of wrong length.") 
	}
	if (is.na(nInit)){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit < 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit%%1 != 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (!is(zeros, "matrix")){  
	    stop("Input (zeros) is of wrong class.") 
	}    
	if(ncol(zeros) != 2){ 
	        stop("Wrong dimensions of the (zeros).") 
	} 
	if (as.character(class(penalize.diag)) != "logical"){ 
	    stop("Input (penalize.diag) is of wrong class.") 
	}    

	# make k-folds as list
	fold         <- max(min(ceiling(fold), nrow(Y)), 2)
	fold         <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
	shuffle      <- sample(1:nrow(Y), nrow(Y))
	folds        <- split(shuffle, as.factor(fold))
	names(folds) <- NULL

	# optimize CV-loss
	optLambda <- optim(lambdaInit, 
                       .armaKCVlossR_banded, 
                       gr=NULL, 
                       method="Brent",
                       Y=Y, 
                       lower=lambdaMin, 
                       upper=lambdaMax, 
                       folds=folds, 
                       target=target, 
                       penalize_diag=penalize.diag,
                       zerosR=zeros[,1]-1, 
                       zerosC=zeros[,2]-1, 
                       nInit=nInit, 
                       minSuccDiff=minSuccDiff)$par

	# return
	return(optLambda)
}



optPenaltyPgen.kCVauto.groups <- function(Y, 
                                      lambdaMin, 
                                      lambdaMax, 
                                      lambdaInit=(lambdaMin + lambdaMax)/2,
                                      fold=nrow(Y), 
                                      groups, 
                                      target, 
                                      zeros=matrix(nrow=0, ncol=2), 
                                      penalize.diag=TRUE, 
                                      nInit=100, 
                                      minSuccDiff=10^(-5)){ 
	## ---------------------------------------------------------------------
	## Function that determines the optimal penalty parameters through 
	## maximization of the k-fold cross-validated log-likelihood score, 
	## assuming that variates are grouped and penalized group-wise.
	## ---------------------------------------------------------------------
	## Arguments:    
	## -> lambdas        : vector with penalty parameter values, one per group.
	##                     values should be specified in the same order as the 
	##                     first appearance of a group representative.  
	## -> Y              : data matrix with rows and columns containing samples
	##                     and variables, respectively
	## -> target         : target matrix for the precision 
	## -> folds          : cross-validation sample splits
	## -> groups         : vector indicating to which group a variate belongs.
	## -> zerosR         : row-index of zero precision elements.
	## -> zerosC         : column-index of zero precision elements.
	## -> penalize_diag  : logical indicating whether the diagonal should 
	##                     be penalized    
	## -> nInit	         : maximum number of iterations. Passed on to the 
	##                     ridgePgen-function.
	## -> minSuccDiff    : minimum successive difference (in the penalized 
	##                     loglikehood) to be achieved. Passed on to the 
	##                     ridgePgen-function.
	## ---------------------------------------------------------------------
	## Value:
	## -> optLambdas     : vector with optimal parameter values.
	## ---------------------------------------------------------------------
	## Authors           : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(lambdaMin) != "numeric"){ 
		stop("Input (lambdaMin) is of wrong class") 
	}
	if (length(lambdaMin) < 1){ 
		stop("Input (lambdaMin) must be a scalar") 
	}
	if (any(lambdaMin <= 0)){ 
		stop("Input (lambdaMin) must be positive") 
	}
	if (class(lambdaMax) != "numeric"){ 
		stop("Input (lambdaMax) is of wrong class") 
	}
	if (length(lambdaMax) < 1){ 
		stop("Input (lambdaMax) must be a scalar") 
	}
	if (any(lambdaMax <= lambdaMin)){ 
		stop("Input (lambdaMax) must be larger than lambdaMin") 
	}
	if (class(lambdaInit) != "numeric"){ 
		stop("Input (lambdaInit) is of wrong class") 
	}
	if (length(lambdaInit) < 1){ 
		stop("Input (lambdaInit) must be a scalar") 
	}
	if (any(lambdaInit <= lambdaMin)){ 
		stop("Input (lambdaInit) must be larger than lambdaMin") 
	}
	if (any(lambdaInit > lambdaMax)){ 
		stop("Input (lambdaInit) must be smaller than lambdaMax") 
	}
	if (class(fold) != "numeric" & class(fold) != "integer"){
		stop("Input (fold) is of wrong class") 
	}
	if ((fold <=  1) | (fold > nrow(Y))){ 
		stop("Input (fold) out of range") 
	}
	if (class(groups) != "numeric" & class(groups) != "integer" & class(groups) != "character"){
		stop("Input (groups) is of wrong class") 
	}
	if ((length(groups) <=  1) | (length(groups) > ncol(Y))){ 
		stop("Input (groups) out of range") 
	}
	if (!isSymmetric(target)){ 
	    stop("target should be a symmetric matrix") 
	}
	if (as.character(class(minSuccDiff)) != "numeric"){ 
	        stop("Input (minSuccDiff) is of wrong class.") 
	}
	if (length(minSuccDiff) != 1){ 
	    stop("Input (minSuccDiff) is of wrong length.") 
	}
	if (is.na(minSuccDiff)){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (minSuccDiff <= 0){ 
	    stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (as.character(class(nInit)) != "numeric" & as.character(class(nInit)) !=  "integer"){ 
	    stop("Input (nInit) is of wrong class.") 
	}
	if (length(nInit) != 1){ 
	    stop("Input (nInit) is of wrong length.") 
	}
	if (is.na(nInit)){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit < 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit%%1 != 0){ 
	    stop("Input (nInit) is not a positive integer.") 
	}
	if (!is(zeros, "matrix")){  
	    stop("Input (zeros) is of wrong class.") 
	}    
	if(ncol(zeros) != 2){ 
	        stop("Wrong dimensions of the (zeros).") 
	} 
	if (as.character(class(penalize.diag)) != "logical"){ 
	    stop("Input (penalize.diag) is of wrong class.") 
	}    

	# make k-folds as list
	fold         <- max(min(ceiling(fold), nrow(Y)), 2)
	fold         <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
	shuffle      <- sample(1:nrow(Y), nrow(Y))
	folds        <- split(shuffle, as.factor(fold))
	names(folds) <- NULL

	# optimize CV-loss
	optLambda <- optim(lambdaInit, 
			   .armaKCVlossR_groups, 
			   gr=NULL, 
			   method="L-BFGS-B",
			   Y=Y, 
			   lower=lambdaMin, 
			   upper=lambdaMax, 
			   folds=folds, 
			   groups=groups, 
			   target=target, 
			   penalize_diag=penalize.diag,
			   zerosR=zeros[,1]-1, 
			   zerosC=zeros[,2]-1, 
			   nInit=nInit, 
			   minSuccDiff=minSuccDiff)$par

	# return
	return(optLambda)
}

