##------------------------------------------------------------------------------
## Functions for the ridge GGM mixture estimator
##------------------------------------------------------------------------------

ridgeGGMmixture <- function(Y, 
                            K, 
                            lambda, 
                            target,                                    
                            iWeights=matrix(sample(seq(0+1/nrow(Y), 
                                                       1-1/nrow(Y), 
                                                       by=1/(2*nrow(Y))), 
                                                   nrow(Y)*K, 
                                                   replace=TRUE), 
                                            nrow=nrow(Y), 
                                            ncol=K),
                            nInit=100,
                            minSuccDiff=10^(-10),
                            minMixProp=0.01){

	## ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## K		: number of mixture components.
        ## lambda       : ridge penalty parameter
        ## target       : target matrix towards the precision matrices are shrunken
	## iWeights	: sample-specific positive component weights 
	##                (may be larger than one as they are rescaled to sum to one)
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized loglikehood)
	##                to be achieved
	## minMixProp	: smallest mixing probability tolerated.
	## ---------------------------------------------------------------------
	## Value
	## A list-object with the following slots:
	## mu		: estimated mean vectors 
	## P		: estimated mixture precision matrices 
	## weights	: estimated component memberships.
	## penLL	: penalized loglikelihood
	## ---------------------------------------------------------------------
	## Authors      : Mehran Aflakparast & Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(K) != "numeric" & class(K) != "integer"){
		stop("Input (K) is of wrong class") 
	}
	if ((K <=  1) | (K > nrow(Y))){ 
		stop("Input (K) out of range") 
	}
	if (as.character(class(lambda)) != "numeric"){ 
		stop("Input (lambda) is of wrong class.") 
	}
	if (length(lambda) != 1){ 
		stop("Input (lambda) is of wrong length.") 
	}
	if (is.na(lambda)){ 
		stop("Input (lambda) is not a positive number.") 
	}
	if (lambda <= 0){ 
		stop("Input (lambda) is not a positive number.") 
	}
	if (!is(target, "matrix")){  
		stop("Input (target) is of wrong class.") 
	}
	if (dim(Y)[2] != nrow(target)){ 
		stop("Dimensions of input (target) do not match that of other input (Y).") 
	}
	if (dim(Y)[2] != ncol(target)){ 
		stop("Dimensions of input (target) do not match that of other input (Y).")  
	}
	if (!is(iWeights, "matrix")){  
		stop("Input (iWeights) is of wrong class.") 
	}
	if (dim(Y)[1] != nrow(iWeights)){ 
		stop("Dimensions of input (iWeights) do not match that of other input (Y).") 
	}
	if (K != ncol(iWeights)){ 
		stop("Dimensions of input (iWeights) do not match that of other input (K).")  
	}
	if (class(nInit) != "numeric" & class(nInit) != "integer"){
		stop("Input (nInit) is of wrong class") 
	}
	if (nInit <=  1){ 
		stop("Input (nInit) out of range") 
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
	if (length(minMixProp) != 1){ 
		stop("Input (minMixProp) is of wrong length.") 
	}
	if (is.na(minMixProp)){ 
		stop("Input (minMixProp) is not a positive number.") 
	}
	if (minMixProp <= 0){ 
		stop("Input (minMixProp) is not a positive number.") 
	}

	# estimate and return GGM mixture model parameters
	return(.armaRidgeGGMmixture(Y, 
	                            K, 
	                            lambda, 
	                            target, 
	                            iWeights, 
	                            nInit, 
	                            minSuccDiff, 
	                            minMixProp))
}


##------------------------------------------------------------------------------
##
## K-fold cross-validation for the GGM mixture setting
##
##------------------------------------------------------------------------------

optPenaltyGGMmixture.kCVauto <- function(Y, 
                                       K, 
                                       lambdaMin, 
                                       lambdaMax,
                                       lambdaInit=(lambdaMin+lambdaMax)/2,
                                       fold=nrow(Y), 
                                       target,                                    
                                       iWeights=matrix(sample(seq(0+1/nrow(Y), 
                                                                  1-1/nrow(Y), 
                                                                  by=1/(2*nrow(Y))), 
                                                              nrow(Y)*K, 
                                                              replace=TRUE), 
                                                       nrow=nrow(Y), 
                                                       ncol=K),
                                       nInit=100,
                                       minSuccDiff=10^(-10),
                                       minMixProp=0.01){

	## ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## K		: number of mixture components.
        ## lambdaMin    : minimum value penalty parameter
        ## lambdaMax    : maximum value penalty parameter
        ## lambdaInit   : initial value for lambda for starting optimization
        ## fold         : cross-validation fold, default gives LOOCV
        ## target       : target matrix towards the precision matrices are shrunken
	## iWeights	: sample-specific positive component weights 
	##                (may be larger than one as they are rescaled to sum to one)
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized loglikehood)
	##                to be achieved
	## minMixProp	: smallest mixing probability tolerated.
	## ---------------------------------------------------------------------
	## Value
	## cvLL		: negative cross-validated loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## ---------------------------------------------------------------------
	## Authors      : Mehran Aflakparast & Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(K) != "numeric" & class(K) != "integer"){
		stop("Input (K) is of wrong class") 
	}
	if ((K <=  1) | (K > nrow(Y))){ 
		stop("Input (K) out of range") 
	}
	if (class(lambdaMin) != "numeric"){ 
		stop("Input (lambdaMin) is of wrong class") 
	}
	if (length(lambdaMin) != 1){ 
		stop("Input (lambdaMin) must be a scalar") 
	}
	if (lambdaMin <= 0){ 
		stop("Input (lambdaMin) must be positive") 
	}
	if (class(lambdaMax) != "numeric"){ 
		stop("Input (lambdaMax) is of wrong class") 
	}
	if (length(lambdaMax) != 1){ 
		stop("Input (lambdaMax) must be a scalar") 
	}
	if (lambdaMax <= lambdaMin){ 
		stop("Input (lambdaMax) must be larger than lambdaMin") 
	}
	if (class(lambdaInit) != "numeric"){ 
		stop("Input (lambdaInit) is of wrong class") 
	}
	if (length(lambdaInit) != 1){ 
		stop("Input (lambdaInit) must be a scalar") 
	}
	if (lambdaInit <= lambdaMin){ 
		stop("Input (lambdaInit) must be larger than lambdaMin") 
	}
	if (lambdaInit > lambdaMax){ 
		stop("Input (lambdaInit) must be smaller than lambdaMax") 
	}
	if (class(fold) != "numeric" & class(fold) != "integer"){
		stop("Input (fold) is of wrong class") 
	}
	if ((fold <=  1) | (fold > nrow(Y))){ 
		stop("Input (fold) out of range") 
	}

	# make k-folds as list
	fold    <- max(min(ceiling(fold), nrow(Y)), 2)
	fold    <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
	shuffle <- sample(1:nrow(Y), nrow(Y))
	folds   <- split(shuffle, as.factor(fold))

	# determine optimal value of ridge penalty parameter
	optLambda <- optim(lambdaInit, 
	                   .armaKcvlGGMmixture, 
	                   method="Brent", 
	                   lower=lambdaMin,
	                   upper=lambdaMax, 
	                   Y=Y, 
	                   K=K, 
	                   folds=folds,
	                   target=target,
	                   iWeights=iWeights,
	                   nInit=nInit,
	                   minSuccDiff=minSuccDiff,
	                   minMixProp=minMixProp)$par

	# return
	return(optLambda)
}



