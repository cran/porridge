##------------------------------------------------------------------------------
## Functions for the ridge signal and error precision estimators from studies with
## replications
##------------------------------------------------------------------------------

ridgePrep <- function(Y, 
                      ids, 
              	      lambdaZ, 
              	      lambdaE, 		
                      targetZ=matrix(0, ncol(Y), ncol(Y)),
                      targetE=matrix(0, ncol(Y), ncol(Y)),
		      nInit=100,
                      minSuccDiff=10^(-10)){

	## ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a 'signal+noise' GGM
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size (including the 
	##                repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## lambdaZ      : ridge penalty parameter for the signal precision matrix
        ## lambdaE      : ridge penalty parameter for the error precision matrix
        ## targetZ      : target matrix towards the signal precision matrix is shrunken
        ## targetE      : target matrix towards the error precision matrix is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized loglikehood)
	##                to be achieved
	## ---------------------------------------------------------------------
	## Value
	## A list-object with the following slots:
	## Pz		: estimated signal precision matrix
	## Pe		: estimated error precision matrix  
	## penLL	: penalized loglikelihood
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(ids) != "numeric" & class(ids) != "integer"){
		stop("Input (ids) is of wrong class") 
	}
	if (as.character(class(lambdaZ)) != "numeric"){ 
		stop("Input (lambdaZ) is of wrong class.") 
	}
	if (length(lambdaZ) != 1){ 
		stop("Input (lambdaZ) is of wrong length.") 
	}
	if (is.na(lambdaZ)){ 
		stop("Input (lambdaZ) is not a positive number.") 
	}
	if (lambdaE <= 0){ 
		stop("Input (lambdaE) is not a positive number.") 
	}
	if (as.character(class(lambdaE)) != "numeric"){ 
		stop("Input (lambdaE) is of wrong class.") 
	}
	if (length(lambdaE) != 1){ 
		stop("Input (lambdaE) is of wrong length.") 
	}
	if (is.na(lambdaE)){ 
		stop("Input (lambdaE) is not a positive number.") 
	}
	if (lambdaE <= 0){ 
		stop("Input (lambdaE) is not a positive number.") 
	}
	if (!is(targetZ, "matrix")){  
		stop("Input (targetZ) is of wrong class.") 
	}
	if (dim(Y)[2] != nrow(targetZ)){ 
		stop("Dimensions of input (targetZ) do not match that of other input (Y).") 
	}
	if (dim(Y)[2] != ncol(targetZ)){ 
		stop("Dimensions of input (targetZ) do not match that of other input (Y).")  
	}
	if (!is(targetE, "matrix")){  
		stop("Input (targetE) is of wrong class.") 
	}
	if (dim(Y)[2] != nrow(targetE)){ 
		stop("Dimensions of input (targetE) do not match that of other input (Y).") 
	}
	if (dim(Y)[2] != ncol(targetE)){ 
		stop("Dimensions of input (targetE) do not match that of other input (Y).")  
	}
	if (class(nInit) != "numeric" & class(nInit) != "integer"){
		stop("Input (nInit) is of wrong class") 
	}
	if (nInit <  1){ 
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

	# estimate and return GGM mixture model parameters
	return(.armaRidgePrepEM(Y, 
	                        ids, 
	                        lambdaZ,
	                        lambdaE,  
	                        targetZ,
	                        targetE,  
	                        nInit, 
	                        minSuccDiff))
}


ridgePrepEdiag <- function(Y, 
                      ids, 
              	      lambdaZ, 
                      targetZ=matrix(0, ncol(Y), ncol(Y)),
		      nInit=100,
                      minSuccDiff=10^(-10)){

	## ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a 'signal+noise' GGM
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size (including the 
	##                repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## lambdaZ      : ridge penalty parameter for the signal precision matrix
        ## lambdaE      : ridge penalty parameter for the error precision matrix
        ## targetZ      : target matrix towards the signal precision matrix is shrunken
        ## targetE      : target matrix towards the error precision matrix is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized loglikehood)
	##                to be achieved
	## ---------------------------------------------------------------------
	## Value
	## A list-object with the following slots:
	## Pz		: estimated signal precision matrix
	## Pe		: estimated error precision matrix  
	## penLL	: penalized loglikelihood
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(ids) != "numeric" & class(ids) != "integer"){
		stop("Input (ids) is of wrong class") 
	}
	if (as.character(class(lambdaZ)) != "numeric"){ 
		stop("Input (lambdaZ) is of wrong class.") 
	}
	if (length(lambdaZ) != 1){ 
		stop("Input (lambdaZ) is of wrong length.") 
	}
	if (is.na(lambdaZ)){ 
		stop("Input (lambdaZ) is not a positive number.") 
	}
	if (!is(targetZ, "matrix")){  
		stop("Input (targetZ) is of wrong class.") 
	}
	if (dim(Y)[2] != nrow(targetZ)){ 
		stop("Dimensions of input (targetZ) do not match that of other input (Y).") 
	}
	if (dim(Y)[2] != ncol(targetZ)){ 
		stop("Dimensions of input (targetZ) do not match that of other input (Y).")  
	}
	if (class(nInit) != "numeric" & class(nInit) != "integer"){
		stop("Input (nInit) is of wrong class") 
	}
	if (nInit <  1){ 
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

	# estimate and return GGM mixture model parameters
	return(.armaRidgePrepEMdiag(Y, 
	                            ids, 
	                            lambdaZ,
	                            targetZ,
	                            nInit, 
	                            minSuccDiff))
}



##------------------------------------------------------------------------------
##
## K-fold cross-validation for the GGM mixture setting
##
##------------------------------------------------------------------------------

optPenaltyPrep.kCVauto <- function(Y, 
				   ids,
                                   lambdaInit,
                                   fold=nrow(Y), 
				   CVcrit,
				   splitting="stratified",
                                   targetZ=matrix(0, ncol(Y), ncol(Y)),
                                   targetE=matrix(0, ncol(Y), ncol(Y)),
                                   nInit=100,
                                   minSuccDiff=10^(-10)){

	## ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size (including the 
	##                repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## lambdaInit   : initial values for the penalty parameters when 
	##                starting optimization
        ## fold         : cross-validation fold, default gives LOOCV
	## CVcrit       : Cross-validation criterion to applied. Either "LL" (the 
	##                loglikelihood), "compLL" (complete loglikelihood), or
	##                "Qloss" (the quadratic loss).
	## splitting    : a character, either "replications", "samples" or "stratified",
	##                specifying either how the splits are to be formed: either 
	##                replications or samples are randomly divided over the "fold" 
	##                groups (first two options, respectively), or samples are
	##	          randomly divided over the groups but in stratified manner
	##                such that the total number of replications in each group is 
	##                roughly comparable.
        ## targetZ      : target matrix towards the signal precision matrix is shrunken
        ## targetE      : target matrix towards the error precision matrix is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized loglikehood)
	##                to be achieved
	## ---------------------------------------------------------------------
	## Value
	## optLambdas   : cross-validated optimal lambdas
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(ids) != "numeric" & class(ids) != "integer"){
		stop("Input (ids) is of wrong class") 
	}
	if (!is(targetZ, "matrix")){  
		stop("Input (targetZ) is of wrong class.") 
	}
	if (dim(Y)[2] != nrow(targetZ)){ 
		stop("Dimensions of input (targetZ) do not match that of other input (Y).") 
	}
	if (dim(Y)[2] != ncol(targetZ)){ 
		stop("Dimensions of input (targetZ) do not match that of other input (Y).")  
	}
	if (!is(targetE, "matrix")){  
		stop("Input (targetE) is of wrong class.") 
	}
	if (dim(Y)[2] != nrow(targetE)){ 
		stop("Dimensions of input (targetE) do not match that of other input (Y).") 
	}
	if (dim(Y)[2] != ncol(targetE)){ 
		stop("Dimensions of input (targetE) do not match that of other input (Y).")  
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
	if (length(lambdaInit) != 2){ 
		stop("Input (lambdaInit) must be a scalar") 
	}
	if (!all(lambdaInit > 0)){ 
		stop("Input (lambdaInit) must element-wise be positive") 
	}
	if (class(fold) != "numeric" & class(fold) != "integer"){
		stop("Input (fold) is of wrong class") 
	}
	if ((fold <=  1) | (fold > nrow(Y))){ 
		stop("Input (fold) out of range") 
	}
	if (class(CVcrit) != "character"){ 
		stop("Input (CVcrit) is of wrong class") 
	}
	if (length(intersect(CVcrit, c("LL", "compLL", "Qloss", "Qloss2" )))!=1){ 
		stop("Input (CVcrit) is wrongly specified") 
	}
	if (class(splitting) != "character"){ 
		stop("Input (splitting) is of wrong class") 
	}
	if (length(intersect(splitting, c("stratified", "replications", "samples")))!=1){ 
		stop("Input (splitting) is wrongly specified") 
	}

	# make k-folds as list
	if (splitting=="replications"){
		# folds over replications
		fold    <- max(min(ceiling(fold), nrow(Y)), 2)
		fold    <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
		shuffle <- sample(1:nrow(Y), nrow(Y))
		folds   <- split(shuffle, as.factor(fold))
	}
	if (splitting=="samples"){
		# folds over samples
		nSamples <- length(unique(ids))
		fold    <- max(min(ceiling(fold), nSamples), 2)
		fold    <- rep(1:fold, ceiling(nSamples/fold))[1:nSamples]
		shuffle <- sample(1:nSamples, nSamples)
		folds   <- split(shuffle, as.factor(fold))
		folds   <- lapply(folds, function(Z, ids){ which(ids %in% Z) }, ids=ids)
	}
	if (splitting=="stratified"){
		# folds over samples, but stratified for number of replications
		fold    <- max(min(ceiling(fold), nrow(Y)), 2)
		slh     <- order(table(ids))
		folds   <- rep(NA, ceiling(length(slh) / fold) * fold) 
		folds[1:length(slh)]   <- slh
		folds   <- matrix(folds, nrow=fold)
		folds   <- apply(folds, 2, function(Z){ Z[sample(1:length(Z), replace=FALSE)] })
		folds   <- split(as.numeric(folds), as.factor(matrix(rep(1:fold, ncol(folds)), nrow=fold)))
		folds   <- lapply(folds, function(Z, ids){ which(ids %in% Z[!is.na(Z)]) }, ids=ids)
	}

	# determine optimal value of ridge penalty parameter
	optLambdas <- constrOptim(lambdaInit, 
				  .armaRidgePrepKcvLL, 
				  grad=NULL, 
				  ui=diag(2), 
				  ci=rep(10^(-10), 2), 
				  Y=Y,
				  ids=ids-1, 
				  folds=folds,
				  CVcrit=CVcrit,
				  targetZ=targetZ,
				  targetE=targetE,
				  nInit=nInit,
				  minSuccDiff=minSuccDiff)$par

	# return
	return(optLambdas)
}


optPenaltyPrepEdiag.kCVauto <- function(Y, 
				   ids,
                                   lambdaInit,
                                   fold=nrow(Y), 
				   CVcrit,
				   splitting="stratified",
                                   targetZ=matrix(0, ncol(Y), ncol(Y)),
                                   nInit=100,
                                   minSuccDiff=10^(-10)){

	## ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size (including the 
	##                repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## lambdaInit   : initial values for the penalty parameter when 
	##                starting optimization
        ## fold         : cross-validation fold, default gives LOOCV
	## CVcrit       : Cross-validation criterion to applied. Either "LL" (the 
	##                loglikelihood), "compLL" (complete loglikelihood), or
	##                "Qloss" (the quadratic loss).
	## splitting    : a character, either "replications", "samples" or "stratified",
	##                specifying either how the splits are to be formed: either 
	##                replications or samples are randomly divided over the "fold" 
	##                groups (first two options, respectively), or samples are
	##	          randomly divided over the groups but in stratified manner
	##                such that the total number of replications in each group is 
	##                roughly comparable.
        ## targetZ      : target matrix towards the signal precision matrix is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized loglikehood)
	##                to be achieved
	## ---------------------------------------------------------------------
	## Value
	## optLambdas   : cross-validated optimal lambdas
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# input checks
	if (!is(Y, "matrix")){  
		stop("Input (Y) should be a matrix")
	}
	if (class(ids) != "numeric" & class(ids) != "integer"){
		stop("Input (ids) is of wrong class") 
	}
	if (!is(targetZ, "matrix")){  
		stop("Input (targetZ) is of wrong class.") 
	}
	if (dim(Y)[2] != nrow(targetZ)){ 
		stop("Dimensions of input (targetZ) do not match that of other input (Y).") 
	}
	if (dim(Y)[2] != ncol(targetZ)){ 
		stop("Dimensions of input (targetZ) do not match that of other input (Y).")  
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
	if (length(lambdaInit) != 1){ 
		stop("Input (lambdaInit) must be a scalar") 
	}
	if (!all(lambdaInit > 0)){ 
		stop("Input (lambdaInit) must element-wise be positive") 
	}
	if (class(fold) != "numeric" & class(fold) != "integer"){
		stop("Input (fold) is of wrong class") 
	}
	if ((fold <=  1) | (fold > nrow(Y))){ 
		stop("Input (fold) out of range") 
	}
	if (class(CVcrit) != "character"){ 
		stop("Input (CVcrit) is of wrong class") 
	}
	if (length(intersect(CVcrit, c("LL", "compLL", "Qloss", "Qloss2")))!=1){ 
		stop("Input (CVcrit) is wrongly specified") 
	}
	if (class(splitting) != "character"){ 
		stop("Input (splitting) is of wrong class") 
	}
	if (length(intersect(splitting, c("stratified", "replications", "samples")))!=1){ 
		stop("Input (splitting) is wrongly specified") 
	}

	# make k-folds as list
	if (splitting=="replications"){
		# folds over replications
		fold    <- max(min(ceiling(fold), nrow(Y)), 2)
		fold    <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
		shuffle <- sample(1:nrow(Y), nrow(Y))
		folds   <- split(shuffle, as.factor(fold))
	}
	if (splitting=="samples"){
		# folds over samples
		nSamples <- length(unique(ids))
		fold    <- max(min(ceiling(fold), nSamples), 2)
		fold    <- rep(1:fold, ceiling(nSamples/fold))[1:nSamples]
		shuffle <- sample(1:nSamples, nSamples)
		folds   <- split(shuffle, as.factor(fold))
		folds   <- lapply(folds, function(Z, ids){ which(ids %in% Z) }, ids=ids)
	}
	if (splitting=="stratified"){
		# folds over samples, but stratified for number of replications
		fold    <- max(min(ceiling(fold), nrow(Y)), 2)
		slh     <- order(table(ids))
		folds   <- rep(NA, ceiling(length(slh) / fold) * fold) 
		folds[1:length(slh)]   <- slh
		folds   <- matrix(folds, nrow=fold)
		folds   <- apply(folds, 2, function(Z){ Z[sample(1:length(Z), replace=FALSE)] })
		folds   <- split(as.numeric(folds), as.factor(matrix(rep(1:fold, ncol(folds)), nrow=fold)))
		folds   <- lapply(folds, function(Z, ids){ which(ids %in% Z[!is.na(Z)]) }, ids=ids)
	}

	# determine optimal value of ridge penalty parameter
	optLambdas <- optim(lambdaInit, 
			    .armaRidgePrepKcvLLdiag, 
			    method="Brent", 
			    lower=10^(-10), 
			    upper=10^(10),
			    Y=Y,
			    ids=ids-1, 
			    folds=folds,
			    CVcrit=CVcrit,
			    targetZ=targetZ,
			    nInit=nInit,
			    minSuccDiff=minSuccDiff)$par

	# return
	return(optLambdas)
}


