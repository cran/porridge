################################################################################
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Function to evaluate the degrees of freedom of the generalized ridge       ##
## estimator of the generalized linear model                                  ##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
################################################################################


ridgeGLMdof <- function(X, U=matrix(ncol=0, nrow=nrow(X)), lambda, lambdaG, 
                        Dg=matrix(0, ncol=ncol(X), nrow=ncol(X)), 
                        model="linear", linPred=rep(0,nrow(X))){

	## ---------------------------------------------------------------------
	## Function that evaluates ridge regression estimator with a regular and
	## generalized penalty. The fused penalty is specified by the matrix Dg.
	## --------------------------------------------------------------------- 
	## Arguments
	## X       : design matrix of the penalized covariates
	## U       : design matrix of the unpenalized covariates
	## lambda  : numeric, the regular ridge penalty parameter
	## lambdaG : numeric, the fused ridge penalty parameter
	## Dg      : nonnegative definite matrix that specifies the structure 
	##           of the generalized ridge penalty.
	## model   : ...
	## linPred : linear predictor
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the degrees of freedom of the generalized ridge 
	## regression estimate.
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

		if (model == "linear"){ 
			# evaluate intermediate matrices
			Dg <- rbind(cbind(crossprod(U), crossprod(U, X)), 
			            cbind(crossprod(X, U), 
			                  crossprod(X) + lambda * diag(ncol(X)) + 
			                                 lambdaG * Dg))
	
			# evaluate condition number (alternative implementation)
			# recCN <- 1/kappa(Hmiddle)
			# if (recCN < .Machine$double.eps){ 
			#	return(sum(diag(crossprod(ginv(Hmiddle, 
			#                                       tol=.Machine$double.eps),  
			#	                          crossprod(cbind(U, X))))))
			# } else {
			# 	return(sum(diag(solve(Hmiddle, crossprod(cbind(U, X))))))
			# }
	
			# evaluate condition number (alternative implementation)
			#          sum(diag(cbind(U, X) %*% 
			#              solve(rbind(cbind(t(U) %*% U, t(U) %*% X), 
			#                    cbind(t(X) %*% U, t(X) %*% X + 
			#                          lambda * diag(ncol(X)) + 
			#                          lambdaF * Df))) %*% 
			#                   rbind(t(U), t(X))))

			# return degrees of freedom
			return(sum(diag(qr.solve(Dg, crossprod(cbind(U, X)), 
			                         tol=.Machine$double.eps))))
		}
		
		if (model == "logistic"){
			# calculate weights
			W0 <- 1 / (1 + exp(-linPred))
			W0 <- as.numeric(W0 * (1 - W0))

			# evaluate intermediate matrices
			X  <- sweep(X, 1, sqrt(W0), FUN="*")
			if (ncol(U) > 0){
				U  <- sweep(U, 1, sqrt(W0), FUN="*")
				Pu <- U %*% solve(crossprod(U), t(U)) 
			}	
			A <- crossprod(t(X), solve(lambda * diag(ncol(X)) + 
			                           lambdaG * Dg, t(X)))
			
			# efficient trace calculation
			if (ncol(U) > 0){
				evs <- eigen((diag(nrow(X)) - Pu) %*% A %*% 
				             (diag(nrow(X)) - Pu), 
				             only.values=TRUE)$values
			} else {
				evs <- eigen(A, only.values=TRUE)$values
			}
				
			# return degrees of freedom
			return(ncol(U) + nrow(X) - sum(1/(1+evs)))
		}
}





################################################################################
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Function to generate the generalized ridge penalty matrix                  ##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
################################################################################

genRidgePenaltyMat <- function(pr, pc=pr, type="2dimA"){
	## ---------------------------------------------------------------------
	## Function produces a penalty parameter matrix to be used in the ridge
	## regression estimator.
	## --------------------------------------------------------------------- 
	## Arguments
	## pr : A positive integer, the row dimension of 2d covariate layout.
	## pc : A positive integer, the second dimension of 2d covariate layout.
	## type : a character specifiying the penalty matrix up to a factor.
	## with a two-dimensional covariate layout. The 
	## matrix is non-negative definite and specifies how neighboring 
	## parameters are to be fused.
	## ---------------------------------------------------------------------
	## Value:
	## A (non-negative definite) matrix.
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	if (type=="2dimA"){
		# prepare diagonal blocks
		diags            <- list(c(2, rep(3, pr-2), 2), rep(-1, pr-1))
		DdiagC           <- as.matrix(bandSparse(pr, k = -c(0:1), 
		                              diagonals=c(diags), symmetric=TRUE))
		colnames(DdiagC) <- rownames(DdiagC) <- NULL
		diags            <- list(c(3, rep(4, pr-2), 3), rep(-1, pr-1))
		DdiagM           <- as.matrix(bandSparse(pr, k = -c(0:1), 
		                              diagonals=c(diags), symmetric=TRUE))
		colnames(DdiagM) <- rownames(DdiagM) <- NULL

		# prepare diagonal blocks
		diags            <- list(rep(-1, pr*(pc-1)))
		D                <- as.matrix(bandSparse(pr*pc, k = -c(pr), 
		                              diagonals=c(diags), symmetric=TRUE))
		colnames(D)      <- rownames(D) <- NULL

		# fill the diagonal with the block
		for (k in 1:pc){
			if (k > 1 & k < pc){
				D[(k-1)*pr + 1:pr, (k-1)*pr + 1:pr] <- DdiagM
			} else {
				D[(k-1)*pr + 1:pr, (k-1)*pr + 1:pr] <- DdiagC
			}
		}
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (type=="1dim"){
		# prepare diagonal blocks
		diags <- list(c(1, rep(2, pr-2), 1), rep(-1, pr-1))
		D     <- as.matrix(bandSparse(pr, k = -c(0,1), 
		                              diagonals=c(diags), 
		                              symmetric=TRUE))
		colnames(D) <- rownames(D) <- NULL
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (type=="common"){
		D <- diag(pr) - matrix(1/pr, pr, pr)
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (type=="2dimD"){
		# prepare diagonal blocks
		diags            <- list(c(1, rep(2, pr-2), 1))
		DdiagC           <- as.matrix(bandSparse(pr, k = -c(0), 
		                              diagonals=c(diags), symmetric=TRUE))
		colnames(DdiagC) <- rownames(DdiagC) <- NULL
		diags            <- list(c(2, rep(4, pr-2), 2))
		DdiagM           <- as.matrix(bandSparse(pr, k = -c(0), 
		                              diagonals=c(diags), symmetric=TRUE))
		colnames(DdiagM) <- rownames(DdiagM) <- NULL

		# prepare diagonal blocks
		Bi <- c(0, rep(c(rep(-1, pr-1), 0), pc-1))
		Bo <- Bi[-c(1, length(Bi))]
		diags            <- list(Bi, Bo)
		D                <- as.matrix(bandSparse(pr*pc, k = -c(pr-1, pr+1), 
		                              diagonals=c(diags), symmetric=TRUE))
		colnames(D)      <- rownames(D) <- NULL

		# fill the diagonal with the block
		for (k in 1:pc){
			# if (k > 1 & k < min(pc, pr)){
			if (k > 1 & k < pc){
				D[(k-1)*pr + 1:pr, (k-1)*pr + 1:pr] <- DdiagM
			} else {
				D[(k-1)*pr + 1:pr, (k-1)*pr + 1:pr] <- DdiagC
			}
		}
	}
	return(D)
}


################################################################################
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Functions for generalized ridge GLM estimation                             ## 
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
################################################################################

.ridgeLM <- function(Y, X, U, lambda, lambdaG, Dg, target){
	## ---------------------------------------------------------------------
	## Function that evaluates ridge regression estimator with a regular and
	## generalized penalty. The fused penalty is specified by the matrix Dg.
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## X       : design matrix multiplied by the eigenvector matrix of the
	##           nonnegative definite matrix that specifies the structure 
	##           of spatial fused ridge penalty.
	## lambda  : numeric, the regular ridge penalty parameter
	## lambdaG : numeric, the fused ridge penalty parameter
	## Dg      : nonnegative definite matrix that specifies the structure 
	##           of the generalized ridge penalty.
	## target  : shrinkage target for regression parameter
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the generalized ridge regression estimate.
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	if (((max(abs(Dg)) > 0) && lambdaG > 0) && (ncol(U) > 0)){
		# generalized ridge penalization + unpenalized covariates

		if (nrow(X) >= ncol(X)){
			# efficient evaluation for n >= p

			# evaluate subexpressions of the estimator
			Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
			XTX  <- crossprod(X)
			tUTX <- crossprod(X, U)
			XTY  <- crossprod(X, Y) - crossprod(XTX, target)
			UTY  <- crossprod(U, Y) - crossprod(tUTX, target)
			UTU  <- crossprod(U)

			# evaluate unpenalized low-dim regression estimator
			gHat <- solve(UTU -  crossprod(solve(Dg + XTX, tUTX), tUTX), 
			              (UTY - crossprod(solve(Dg + XTX, XTY), tUTX)[1,]))

			# evaluate penalized high-dim regression estimator
			bHat <- target + solve(XTX + Dg, XTY - crossprod(t(tUTX), gHat))
		}

		if (nrow(X) < ncol(X)){
			# efficient evaluation for n < p

			# evaluate subexpressions of the estimator
			Dg     <- qr.solve(diag(rep(lambda, ncol(X))) + 
			                   lambdaG * Dg, t(X), tol=.Machine$double.eps)
			Y      <- Y - crossprod(t(X), target)
			XDXTpI <- crossprod(t(X), Dg) + diag(nrow(X))

			# evaluate unpenalized low-dim regression estimator
			gHat <- solve(crossprod(U, solve(XDXTpI, U)),
			              crossprod(U, solve(XDXTpI, Y)))

			# evaluate penalized high-dim regression estimator
			Y    <- Y - crossprod(t(U), gHat) 
 			bHat <- target + tcrossprod(Dg, t(solve(XDXTpI, Y)))[,1]
 		}
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (((max(abs(Dg)) > 0) && lambdaG > 0) && (ncol(U) == 0)){
		# generalized ridge penalization + no unpenalized covariates

		if (nrow(X) >= ncol(X)){
			# efficient evaluation for n >= p

			# evaluate subexpressions of the estimator
			Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
			XTX  <- crossprod(X)
			XTY  <- crossprod(X, Y) - crossprod(XTX, target)

			# evaluate penalized high-dim regression estimator
			bHat <- target + solve(XTX + Dg, XTY)
		}

		if (nrow(X) < ncol(X)){
			# efficient evaluation for n < p

			# evaluate subexpressions of the estimator
			Dg     <- qr.solve(diag(rep(lambda, ncol(X))) + 
			                   lambdaG * Dg, t(X), tol=.Machine$double.eps)
			Y      <- Y - crossprod(t(X), target)
			XDXTpI <- crossprod(t(X), Dg) + diag(nrow(X))

			# evaluate penalized high-dim regression estimator
 			bHat <- target + tcrossprod(Dg, t(solve(XDXTpI, Y)))[,1]
 		}
		gHat <- numeric()
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (((max(abs(Dg)) == 0) | lambdaG == 0) && (ncol(U) == 0)){
		# no generalized ridge penalization + no unpenalized covariates

		if (ncol(X) == 1){
			# if there is only a single covariate
			bHat <- solve(crossprod(X, X) + rep(lambda, ncol(X))) %*% 
			        (crossprod(X, Y) + lambda * target)
		}
		if (nrow(X) >= ncol(X)){
			# low-dimensional case
			bHat <- crossprod(solve(crossprod(X, X) + diag(rep(lambda, ncol(X)))),
	 		                 (crossprod(X, Y) + lambda * target))
		}
		if (nrow(X) < ncol(X)){
			# high-dimensional case
			XXT       <- tcrossprod(X) / lambda
			diag(XXT) <- diag(XXT) + 1
			XYpT      <- (crossprod(X, Y) + lambda * target)
			bHat      <- (XYpT / lambda - lambda^(-2) * 
			             crossprod(X, t(crossprod(tcrossprod(X, t(XYpT)), 
							      solve(XXT)))))[,1]
		}
		gHat <- numeric()
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (((max(abs(Dg)) == 0) | lambdaG == 0) && (ncol(U) > 0)){
		# no generalized ridge penalization + unpenalized covariates

		if (nrow(X) >= ncol(X)){
			# efficient evaluation for n >= p

			# evaluate subexpressions of the estimator
			tUTX        <- crossprod(X, U)
			XTXpI       <- crossprod(X) 
			XTY         <- crossprod(X, Y) - crossprod(XTXpI, target)
			UTY         <- crossprod(U, Y) - crossprod(tUTX, target)
			UTU         <- crossprod(U)
			diag(XTXpI) <- diag(XTXpI) + lambda

			# evaluate unpenalized low-dim regression estimator
			gHat <- solve(UTU - crossprod(solve(XTXpI, tUTX), tUTX), 
			             (UTY - as.numeric(crossprod(solve(XTXpI, XTY), tUTX))))

			# evaluate penalized high-dim regression estimator
			bHat <- target + solve(XTXpI, XTY - crossprod(t(tUTX), gHat))
		}

		if (nrow(X) < ncol(X)){
			# efficient evaluation for n < p

			# evaluate subexpressions of the estimator
			XXT       <- tcrossprod(X) / lambda
			diag(XXT) <- diag(XXT) + 1
			Y         <- Y - crossprod(t(X), target)

			# evaluate unpenalized low-dim regression estimator
			gHat <- solve(crossprod(U, solve(XXT, U)),
			              crossprod(U, solve(XXT, Y)))

			# evaluate penalized high-dim regression estimator
			Y    <- Y - crossprod(t(U), gHat) 
			bHat <- target + crossprod(X, solve(XXT, Y))/lambda
 		}
	}

	return(c(gHat, bHat))
}



.ridgeBLM <- function(Y, X, U, lambda, lambdaG, Dg, target, minSuccDiff, maxIter){

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

	if (((max(abs(Dg)) == 0) | lambdaG == 0) && (ncol(U) == 0)){
		# no generalized ridge penalization + no unpenalized covariates

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
				W0[which(W0 < .Machine$double.eps)] <- 
				              .Machine$double.eps
			}

			# distinguish between the high- and low-dimensional cases
			if (ncol(X) >= nrow(X)){ 
				# obtain part one of the estimator
				Z <- lp + (Y - Ypred)/W0

				# now obtain the IRWLS update efficiently
				diag(XXT) <- diag(XXT) + 1/W0
				slh       <- solve(XXT, Z-Xtarget)
				diag(XXT) <- diag(XXT) - 1/W0
				lp        <- crossprod(XXT, slh)
				penalty   <- sum(lp * slh) / 2
				lp        <- Xtarget + lp
				loglik    <- .loglikBLMlp(Y, lp) - penalty
			} else {
				# obtain part one of the estimator
				XWZpT <- lambda * target + 
				         as.numeric(crossprod(X, W0 * lp + Y - Ypred))

				# evaluate X %*% X^t and add the reciprocal weight to 
				# its diagonal
				XTX       <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
				diag(XTX) <- diag(XTX) + lambda 
	
				# now obtain the IRWLS update efficiently
				bHat <- solve(XTX, XWZpT)
				lp   <- tcrossprod(X, t(bHat))

				# evaluate the loglikelihood
				loglik <- .loglikBLMlp(Y, lp) - 
				          0.5 * lambda * sum((bHat - target)^2)
			}

			# assess convergence
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		if (ncol(X) >= nrow(X)){ 
			bHat <- target + crossprod(X, slh) / lambda
		}
		gHat <- numeric()
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (((max(abs(Dg)) != 0) && lambdaG != 0) && (ncol(U) == 0)){
		# generalized ridge penalization + no unpenalized covariates

		# initiate
		loglikPrev <- -10^(10)
		if (ncol(X) >= nrow(X)){ 
			# evaluate subexpressions of the estimator
			Dg          <- qr.solve(diag(rep(lambda, ncol(X))) + 
			                        lambdaG * Dg, t(X), 
			                        tol=.Machine$double.eps)
			XDXT        <- crossprod(t(X), Dg)
			diagXDXTorg <- diag(XDXT)
		} else{
			# evaluate subexpressions of the estimator
			Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
		}
		Xtarget  <- t(tcrossprod(target, X))
		lp       <- rep(0, length(Y))
		for (iter in 1:maxIter){
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lp))
			W0    <- Ypred * (1 - Ypred)
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- 
				              .Machine$double.eps
			}

			# distinguish between the high- and low-dimensional cases
			if (ncol(X) >= nrow(X)){ 
				# obtain part one of the estimator
				Z <- lp + (Y - Ypred)/W0

				# now obtain the IRWLS update efficiently
				diag(XDXT) <- diag(XDXT) + 1/W0
				slh        <- qr.solve(XDXT, Z-Xtarget, 
				                       tol=.Machine$double.eps)
				diag(XDXT) <- diagXDXTorg
				lp         <- crossprod(XDXT, slh)
				penalty    <- sum(lp * slh) / 2
				lp         <- Xtarget + lp
				loglik    <- .loglikBLMlp(Y, lp) - penalty
			} else {
				# obtain part one of the estimator
				XWZpT <- crossprod(Dg, target) + 
				         as.numeric(crossprod(X, W0 * lp + Y - Ypred))

				# evaluate X %*% X^t and add the reciprocal weight to 
				# its diagonal
				XTX <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
				XTX <- XTX + Dg
	
				# now obtain the IRWLS update efficiently
				bHat <- qr.solve(XTX, XWZpT, tol=.Machine$double.eps)
				lp   <- tcrossprod(X, t(bHat))

				# evaluate the loglikelihood
				loglik <- .loglikBLMlp(Y, lp) -
				          0.5 * sum(crossprod(Dg, bHat - target) *
				                    (bHat - target))
			}

			# assess convergence
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		if (ncol(X) >= nrow(X)){ 
			bHat <- target + as.numeric(tcrossprod(as.numeric(slh), Dg))
		}
		gHat <- numeric()
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (((max(abs(Dg)) == 0) | lambdaG == 0) && (ncol(U) > 0)){
		# no generalized ridge penalization + unpenalized covariates

		# initiate
		loglikPrev <- -10^(10)
		if (ncol(X) >= nrow(X)){ 
			XXT <- tcrossprod(X) / lambda
		}
		if (ncol(X) >= nrow(X)){ 
			tUTX <- crossprod(X, U)
		}
		Xtarget <- t(tcrossprod(target, X))
		lp      <- rep(0, length(Y))
		for (iter in 1:maxIter){
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lp))
			W0    <- as.numeric(Ypred * (1 - Ypred))
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- 
				              .Machine$double.eps
			}

			# distinguish between the high- and low-dimensional cases
			if (ncol(X) >= nrow(X)){ 
				# obtain part one of the estimator
				Z <- lp + (Y - Ypred)/W0

				# now obtain the IRWLS update efficiently
				diag(XXT) <- diag(XXT) + 1/W0
				slh       <- solve(XXT, Z-Xtarget)
				gHat      <- solve(crossprod(U, solve(XXT, U)),
					           crossprod(U, slh))
				slh       <- solve(XXT, Z - Xtarget - U %*% gHat)				
				diag(XXT) <- diag(XXT) - 1/W0
				lp        <- crossprod(XXT, slh)
				penalty   <- sum(lp * slh) / 2
				lp        <- Xtarget + lp + U %*% gHat
				loglik    <- .loglikBLMlp(Y, lp) - penalty
			}
			if (ncol(X) < nrow(X)){ 
				# adjusted response
				Z <- W0*lp + (Y - Ypred) - W0*Xtarget

				# evaluate subexpressions of the estimator
				XTXpI       <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
				tUTX        <- crossprod(X, sweep(U, 1, W0, FUN="*"))
				XTZ         <- crossprod(X, Z)
				UTZ         <- crossprod(U, Z)
				UTU         <- crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
				diag(XTXpI) <- diag(XTXpI) + lambda
	
				# evaluate unpenalized low-dim regression estimator
				gHat <- as.numeric(
				        solve(UTU - crossprod(solve(XTXpI, tUTX), tUTX), 
				             (UTZ - as.numeric(crossprod(solve(XTXpI, XTZ), 
				                    crossprod(X, sweep(U, 1, W0, FUN="*")))))))

				# obtain part one of the estimator
				XWZpT <- lambda * target + as.numeric(
				         crossprod(X, W0*lp + (Y - Ypred) - 
				         W0*as.numeric(tcrossprod(gHat, U))))

				# now obtain the IRWLS update efficiently
				bHat <- solve(XTXpI, XWZpT)
				lp   <- as.numeric(tcrossprod(X, t(bHat)) + 
				        as.numeric(tcrossprod(gHat, U)))

				# evaluate the loglikelihood
				loglik <- .loglikBLMlp(Y, lp) -
				          0.5 * lambda * sum((bHat - target)^2)
			}

			# assess convergence
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		if (ncol(X) >= nrow(X)){ 
			bHat <- target + crossprod(X, slh) / lambda
		}
	}

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (((max(abs(Dg)) != 0) && lambdaG != 0) && (ncol(U) > 0)){
		# generalized ridge penalization + unpenalized covariates

		# initiate
		loglikPrev <- -10^(10)
		if (ncol(X) >= nrow(X)){ 
			# evaluate subexpressions of the estimator
			Dg          <- qr.solve(diag(rep(lambda, ncol(X))) + 
			                        lambdaG * Dg, t(X), 
			                        tol=.Machine$double.eps)
			XDXT        <- crossprod(t(X), Dg)
			diagXDXTorg <- diag(XDXT)
		} else { 
			# evaluate subexpressions of the estimator
			Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
		}
		Xtarget <- t(tcrossprod(target, X))
		lp      <- rep(0, length(Y))
		for (iter in 1:maxIter){
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lp))
			W0    <- as.numeric(Ypred * (1 - Ypred))
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- 
				              .Machine$double.eps
			}

			# distinguish between the high- and low-dimensional cases
			if (ncol(X) >= nrow(X)){ 
				# obtain part one of the estimator
				Z <- lp + (Y - Ypred)/W0

				# now obtain the IRWLS update efficiently
				diag(XDXT) <- diag(XDXT) + 1/W0
				slh        <- qr.solve(XDXT, Z-Xtarget, 
				                       tol=.Machine$double.eps)
				gHat       <- solve(crossprod(U, solve(XDXT, U)),
					            crossprod(U, slh))
				slh        <- qr.solve(XDXT, Z-Xtarget - U %*% gHat,
				                       tol=.Machine$double.eps)
				diag(XDXT) <- diagXDXTorg
				lp         <- crossprod(XDXT, slh)
				penalty    <- sum(lp * slh) / 2
				lp         <- Xtarget + lp + U %*% gHat
				loglik    <- .loglikBLMlp(Y, lp) - penalty
			}
			if (ncol(X) < nrow(X)){ 
				# adjusted response
				Z <- W0*lp + (Y - Ypred) - W0*Xtarget

				# evaluate subexpressions of the estimator
				XTXpD <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
				tUTX  <- crossprod(X, sweep(U, 1, W0, FUN="*"))
				XTZ   <- crossprod(X, Z)
				UTZ   <- crossprod(U, Z)
				UTU   <- crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
				XTXpD <- XTXpD + Dg
	
				# evaluate unpenalized low-dim regression estimator
				gHat <- as.numeric(
				        solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
				             (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
				                    crossprod(X, sweep(U, 1, W0, FUN="*")))))))

				# obtain part one of the estimator
				XWZpT <- crossprod(Dg, target) + 
				         crossprod(X, W0*lp + (Y - Ypred) - 
				         W0*tcrossprod(gHat, U)[1,])[,1]

				# now obtain the IRWLS update efficiently
				bHat    <- qr.solve(XTXpD, XWZpT, 
				                    tol=.Machine$double.eps)
				penalty <- sum(crossprod(Dg, bHat - target) * 
				                            (bHat - target))/2
				lp      <- as.numeric(tcrossprod(X, t(bHat)) + 
				           as.numeric(tcrossprod(gHat, U)))

				# evaluate the loglikelihood
				loglik    <- .loglikBLMlp(Y, lp) - penalty
			}

			# assess convergence
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		if (ncol(X) >= nrow(X)){ 
			bHat <- target + as.numeric(tcrossprod(as.numeric(slh), Dg))
		}
	}
	return(c(gHat, bHat))
}



ridgeGLM <- function(Y, X, U=matrix(ncol=0, nrow=length(Y)), 
                     lambda, lambdaG=0, 
                     Dg=matrix(0, ncol=ncol(X), nrow=ncol(X)), 
                     target=rep(0, ncol(X)), model="linear", 
                     minSuccDiff=10^(-10), maxIter=100){

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
		return(.ridgeLM(Y, X, U, lambda, lambdaG, Dg, target))

	}
	if (model == "logistic"){
		return(.ridgeBLM(Y, X, U, lambda, lambdaG, Dg, target, 
		                 minSuccDiff, maxIter))
	}
}


ridgeGLMmultiT <- function(Y, X, U=matrix(ncol=0, nrow=length(Y)), 
                     lambdas, targetMat, model="linear", 
                     minSuccDiff=10^(-10), maxIter=100){

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
	mixedTarget  <- rep(0, ncol(X))
	lambdaTotal  <- sum(lambdas)
	for (k in 1:length(lambdas)){
		mixedTarget  <- mixedTarget + 
		                (lambdas[k] / lambdaTotal) * targetMat[,k]
	}

	# evaluate the regression estimator
	if (model == "linear"){
		return(.ridgeLM(Y, X, U, lambda=lambdaTotal, lambdaG=0, 
		                Dg=matrix(0, ncol(X), ncol(X)), 
		                target=mixedTarget))
	}
	if (model == "logistic"){
		return(.ridgeBLM(Y, X, U, lambda=lambdaTotal, lambdaG=0, 
		                 Dg=matrix(0, ncol(X), ncol(X)), 
		                 target=mixedTarget, 
		                 minSuccDiff, maxIter))
	}
}


##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Auxillary functions for cross-validation
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

makeFoldsGLMcv <- function(fold, Y, stratified=TRUE, model="linear"){

	## --------------------------------------------------------------------
	## Function that generates the folds for cross-validation
	## --------------------------------------------------------------------- 
	## Arguments
	## fold        : number of cross-validation folds
	## stratified  : boolean, should the data by stratified such that 
	##               range of the response is comparable among folds
	## Y           : response vector
	## model       : A character indicating which generalized linear model
	##               model instance is to be fitted.
	## ---------------------------------------------------------------------
	## Value:
	## A list with the folds for cross-validation
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# create folds
	if (!stratified){
	        fold    <- max(min(ceiling(fold), length(Y)), 2)
        	fold    <- rep(1:fold, ceiling(length(Y)/fold))[1:length(Y)]
        	shuffle <- sample(1:length(Y), length(Y))
        	folds   <- split(shuffle, as.factor(fold))
	}
	if (stratified){
		if (model == "linear"){	
			idSorted <- c(order(Y), rep(NA, ceiling(length(Y) / 
			                                fold) * fold - length(Y)))
			folds    <- matrix(idSorted, nrow=fold)
			folds    <- apply(folds, 2, function(Z){ 
			                            Z[sample(1:length(Z), 
			                              length(Z))] })
			folds    <- split(as.numeric(folds), 
			                  as.factor(matrix(rep(1:fold, 
			                            ncol(folds)), nrow=fold)))
			folds    <- lapply(folds, function(Z){ Z[!is.na(Z)] })
		}
		if (model == "logistic"){		
			idSorted <- c(which(Y == 0)[sample(1:sum(Y==0), sum(Y==0))], 
			              which(Y == 1)[sample(1:sum(Y==1), sum(Y==1))])
			idSorted <- c(idSorted, rep(NA, ceiling(length(Y) / 
			                                fold) * fold - length(Y)))
			folds    <- matrix(idSorted, nrow=fold)
			folds    <- split(as.numeric(folds), 
			                 as.factor(matrix(rep(1:fold, ncol(folds)), 
			                           nrow=fold)))
			folds    <- lapply(folds, function(Z){ Z[!is.na(Z)] })
		}		
	}
	return(folds)
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
	lp <- as.numeric(tcrossprod(betas, X))

	# evaluate the loglikelihood
	loglik1 <- exp(Y * lp) / (1 + exp(lp))
	loglik2 <- exp(((Y-1) * lp)) / (1 + exp(-lp))
	loglik1[!is.finite(log(loglik1))] <- NA
	loglik2[!is.finite(log(loglik2))] <- NA
	loglik  <- sum(log(apply(cbind(loglik1, loglik2), 1, mean, na.rm=TRUE)))
	return(loglik)
}



.loglikBLMlp <- function(Y, lp){

	## ---------------------------------------------------------------------
	## Function calculates the loglikelihood of the logistic regression model
	## --------------------------------------------------------------------- 
	## Arguments
	## Y       : response vector
	## lp      : linear predictor
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the loglikelihood of the logistic regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# evaluate the loglikelihood
	loglik1 <- exp(Y * lp) / (1 + exp(lp))
	loglik2 <- exp(((Y-1) * lp)) / (1 + exp(-lp))
	loglik1[!is.finite(log(loglik1))] <- NA
	loglik2[!is.finite(log(loglik2))] <- NA
	loglik  <- sum(log(apply(cbind(loglik1, loglik2), 1, mean, na.rm=TRUE)))
	return(loglik)
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

	return(-length(Y) * (log(2 * pi * sum((Y - as.numeric(tcrossprod(as.numeric(betas), X)))^2) / length(Y)) + 1) / 2)
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

	return(sum((Y - as.numeric(tcrossprod(as.numeric(betas), X)))^2))
}


##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Functions for cross-validation of the linear regression model
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

.kcvLMloss_PleqN_noUP_noGP <- function(lambda, Y, X, folds, target, loss="loglik"){

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
	                                          ridgeGLM(Y=Y[-folds[[k]]], 
		                                           X=X[-folds[[k]],,drop=FALSE], 
	                                                   lambda=lambda, 
	                                                   target=target, 
	                                                   model="linear"))
		}
		if (loss=="loglik"){
			cvLoss <- cvLoss - .loglikLM(Y[folds[[k]]], 
	                                             X[folds[[k]],,drop=FALSE], 
	                                             ridgeGLM(Y[-folds[[k]]], 
		                                              X[-folds[[k]],,drop=FALSE], 
	                                                      lambda=lambda, 
	                                                      target=target, 
	                                                      model="linear"))
		}
	}
	
	# average over the folds
	return(cvLoss / length(folds))
}


.kcvLMloss_PgeqN_noUP_noGP_org <- function(lambda, Y, X, folds, target, loss){

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
	XXT     <- tcrossprod(X) / lambda
	Xtarget <- t(tcrossprod(target, X))
	for (k in 1:length(folds)){
		Yhat   <- Xtarget[folds[[k]]] + 
		          tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
		                           XXT[-folds[[k]], -folds[[k]]], 
		                           Y[-folds[[k]]] - Xtarget[-folds[[k]]]), 
		                     XXT[folds[[k]], -folds[[k]], drop=FALSE])

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


.kcvLMloss_PgeqN_noUP_noGP <- function(lambda, Y, XXT, Xtarget, folds, loss){

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
	XXT    <- XXT / lambda
	for (k in 1:length(folds)){
		Yhat   <- Xtarget[folds[[k]]] + 
		          tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
		                           XXT[-folds[[k]], -folds[[k]]], 
		                           Y[-folds[[k]]] - Xtarget[-folds[[k]]]), 
		                     XXT[folds[[k]], -folds[[k]], drop=FALSE])

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


.optPenaltyLM_noUP_noGP.kCVauto <- function(Y, X, lambdaInit, folds, target, 
                                            loss, lambdaMin, implementation){

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
	## lambdaMin  : lower bound of search interval
	## lambdaMax  : upper bound of search interval
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the optimal cross-validated penalty parameter of the
	## ridge estimator of the linear regression model parameter
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# search by the Brent algorithm
	if (nrow(X) < ncol(X)){
		if (implementation=="alt"){
			XXT       <- tcrossprod(X)
			Xtarget   <- t(tcrossprod(target, X))
			optLambda <- optim(lambdaInit, 
			                   .kcvLMloss_PgeqN_noUP_noGP, Y=Y, 
			                   XXT=XXT, Xtarget=Xtarget, 
			                   folds=folds, loss=loss,
			                   lower=lambdaMin, upper=10^(10), 
			                   method="Brent")$par
		}
		if (implementation=="org"){
			optLambda <- optim(lambdaInit, 
			                   .kcvLMloss_PgeqN_noUP_noGP_org, Y=Y,
			                   X=X, target=target, folds=folds, 
			                   loss=loss, method="Brent", upper=10^(10), 
			                   lower=lambdaMin)$par
		}
	}
	if (nrow(X) >= ncol(X)){
		optLambda <- optim(lambdaInit, 
		                   .kcvLMloss_PleqN_noUP_noGP, Y=Y, X=X, 
		                   folds=folds, target=target, loss=loss, 
		                   method="Brent", lower=lambdaMin,
		                   upper=10^(10))$par
	}
	return(optLambda)
}


.kcvLMloss_PgeqN_UP_noGP <- function(lambda, Y, X, U, target, folds, loss){

	## ---------------------------------------------------------------------
	## Internal function yields the cross-validated loglikelihood of the 
	## ridge regression estimator with a two-dimensional covariate layout. 
	## --------------------------------------------------------------------- 
	## Arguments
	## lambdas : penalty parameter vector
	## Y       : response vector
	## X       : design matrix multiplied by the eigenvector matrix of the
	##           nonnegative definite matrix that specifies the structure 
	##           of spatial fused ridge penalty.
	## U       : design matrix of the unpenalized covariates
	## Xtarget : shrinkage target for regression parameter (multiplied by
	##           the eigenvector matrix of the fusion matrix)
	## Ds      : nonnegative eigenvalues of the nonnegative definite matrix 
	##           that specifies the structure of spatial fused ridge penalty.
	## folds   : list with folds
	## loss    : character, either 'loglik' of 'sos', specifying the loss
	##           criterion to be used in the cross-validation
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the linear regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# evaluate loss per fold
	cvLoss <- 0

	# calculation prior to cv-loop
	XXT     <- tcrossprod(X) / lambda
	Xtarget <- tcrossprod(X, matrix(target, nrow=1))

	for (k in 1:length(folds)){
		# evaluate the estimator of the unpenalized low-dimensional 
		# regression parameter on the training samples
		if (all(dim(U) > 0)){
			gHat <- solve(crossprod(U[-folds[[k]],,drop=FALSE], 
			                       solve(XXT[-folds[[k]], -folds[[k]]] + 
					             diag(nrow(XXT)-length(folds[[k]])), 
			                             U[-folds[[k]],,drop=FALSE])),
		                      crossprod(U[-folds[[k]],,drop=FALSE], 
			                       solve(XXT[-folds[[k]], -folds[[k]]] + 
					             diag(nrow(XXT)-length(folds[[k]])), 
			                       Y[-folds[[k]]] - Xtarget[-folds[[k]],])))
			gHat <- as.numeric(gHat)
		}

		# evaluate the linear predictor on the left-out samples
		if (all(dim(U) > 0)){
			Yhat <- Xtarget[folds[[k]]] + 
		        	tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
		        	                 XXT[-folds[[k]], -folds[[k]]], 
		        	                 Y[-folds[[k]]] - Xtarget[-folds[[k]]] - 
			                         as.numeric(crossprod(t(U[-folds[[k]],,drop=FALSE]), gHat))), 
		        	           XXT[folds[[k]], -folds[[k]], drop=FALSE])
		} else {
			Yhat <- Xtarget[folds[[k]]] +
			        tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
	        	                         XXT[-folds[[k]], -folds[[k]]],
	        	                         Y[-folds[[k]]] - Xtarget[-folds[[k]]]), 
			        	    XXT[folds[[k]], -folds[[k]], drop=FALSE])
		}
	
		if (loss == "loglik"){
			# evaluate the loglikelihood on the left-out samples
			cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
		                           (log(2 * pi * sum((Y[folds[[k]]] - Yhat)^2) / 
			                   length(folds[[k]])) + 1) / 2
		}
		if (loss == "sos"){
			# evaluate the sum-of-sqaures on the left-out samples
			cvLoss <- cvLoss + sum((Y[folds[[k]]] - Yhat)^2)
		}
	}
		
	# average over the folds
	return(cvLoss / length(folds))
}



.kcvLMloss_PleqN_UP_GPandNoGP <- function(lambdas, Y, X, U, Dg, folds, target, loss){

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
		# evaluate regression estimate
		bHat <- ridgeGLM(Y[-folds[[k]]], X[-folds[[k]],,drop=FALSE],
		                 U[-folds[[k]],,drop=FALSE], Dg=Dg,
		                 lambda =lambdas[1], lambdaG=lambdas[2],
		                 target=target, model="linear")

		# evaluate linear predictor		                 	
		Yhat <- X[folds[[k]],,drop=FALSE] %*% bHat[-c(1:ncol(U))] +
		        U[folds[[k]],,drop=FALSE] %*% bHat[c(1:ncol(U))]

		if (loss == "loglik"){
			# evaluate the loglikelihood on the left-out samples
			cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
		                           (log(2 * pi * sum((Y[folds[[k]]] - as.numeric(Yhat))^2) / 
			                   length(folds[[k]])) + 1) / 2
		}
		if (loss == "sos"){
			# evaluate the sum-of-sqaures on the left-out samples
			cvLoss <- cvLoss + sum((Y[folds[[k]]] - as.numeric(Yhat))^2)
		}
	}
	
	# average over the folds
	return(cvLoss / length(folds))
}


.kcvLMloss_PgeqN_UP_GP <- function(lambdas, Y, X, U, Xtarget, Ds, folds, loss){

	## ---------------------------------------------------------------------
	## Internal function yields the cross-validated loglikelihood of the 
	## ridge regression estimator with a two-dimensional covariate layout. 
	## --------------------------------------------------------------------- 
	## Arguments
	## lambdas : penalty parameter vector
	## Y       : response vector
	## X       : design matrix multiplied by the eigenvector matrix of the
	##           nonnegative definite matrix that specifies the structure 
	##           of spatial fused ridge penalty.
	## U       : design matrix of the unpenalized covariates
	## Xtarget : shrinkage target for regression parameter (multiplied by
	##           the eigenvector matrix of the fusion matrix)
	## Ds      : nonnegative eigenvalues of the nonnegative definite matrix 
	##           that specifies the structure of spatial fused ridge penalty.
	## folds   : list with folds
	## loss    : character, either 'loglik' of 'sos', specifying the loss
	##           criterion to be used in the cross-validation
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the cross-validated loss of the linear regression model
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	# evaluate loss per fold
	cvLoss <- 0

	# calculation prior to cv-loop
	X       <- sweep(X, 2, sqrt(lambdas[1] + lambdas[2] * Ds), FUN="/")
	XXT     <- tcrossprod(X)

	for (k in 1:length(folds)){
		# evaluate the estimator of the unpenalized low-dimensional 
		# regression parameter on the training samples
		if (all(dim(U) > 0)){
			gHat <- solve(crossprod(U[-folds[[k]],,drop=FALSE], 
			                       solve(XXT[-folds[[k]], -folds[[k]]] + 
					             diag(nrow(XXT)-length(folds[[k]])), 
			                             U[-folds[[k]],,drop=FALSE])),
		                      crossprod(U[-folds[[k]],,drop=FALSE], 
			                       solve(XXT[-folds[[k]], -folds[[k]]] + 
					             diag(nrow(XXT)-length(folds[[k]])), 
			                       Y[-folds[[k]]] - Xtarget[-folds[[k]],])))
			gHat <- as.numeric(gHat)
		}

		# evaluate the linear predictor on the left-out samples
		if (all(dim(U) > 0)){
			Yhat   <- Xtarget[folds[[k]]] + 
		        	  tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
		        	                   XXT[-folds[[k]], -folds[[k]]], 
		        	                   Y[-folds[[k]]] - Xtarget[-folds[[k]]] - 
			                           as.numeric(crossprod(t(U[-folds[[k]],,drop=FALSE]), gHat))), 
		        	             XXT[folds[[k]], -folds[[k]], drop=FALSE])
		} else {
			Yhat   <- Xtarget[folds[[k]]] + 
		        	  tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
		        	                   XXT[-folds[[k]], -folds[[k]]], 
		        	                   Y[-folds[[k]]] - Xtarget[-folds[[k]]]), 
		        	             XXT[folds[[k]], -folds[[k]], drop=FALSE])
		}


		if (loss == "loglik"){
			# evaluate the loglikelihood on the left-out samples
			cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
		                           (log(2 * pi * sum((Y[folds[[k]]] - as.numeric(Yhat))^2) / 
			                   length(folds[[k]])) + 1) / 2
		}
		if (loss == "sos"){
			# evaluate the sum-of-sqaures on the left-out samples
			cvLoss <- cvLoss + sum((Y[folds[[k]]] - as.numeric(Yhat))^2)
		}
	}
		
	# average over the folds
	return(cvLoss / length(folds))
}


##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Functions for cross-validation of the logistic regression model
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

.kcvBLMloss_PgeqN_noUP_noGP_alt <- function(lambda, Y, XXT, Xtarget, folds, 
                                            minSuccDiff, maxIter){

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
	XXT    <- XXT / lambda
	for (k in 1:length(folds)){
		# convergence criterion
		loglikPrev <- -10^(10)

		# evaluate the initial linear predictor
		lpOld <- rep(0, nrow(XXT))[-folds[[k]]]

		for (iter in 1:maxIter){
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lpOld))
			W0    <- Ypred * (1 - Ypred)
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- 
				              .Machine$double.eps
			}

			# obtain the adjusted response
			WZ <- W0 * lpOld + Y[-folds[[k]]] - Ypred  

			# obtain IRWLS update of the linear predictors efficiently
			XXTZpT <- solve(XXT[-folds[[k]], -folds[[k]]] + diag(1/W0), 
			                 (1/W0) * (WZ - W0 * Xtarget[-folds[[k]]]))
			lpAll  <- Xtarget + 
			          tcrossprod(XXT[,-folds[[k]],drop=FALSE], XXTZpT)
			lpOld  <- lpAll[-folds[[k]]]
			lpNew  <- lpAll[ folds[[k]]]

			# assess convergence
			loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld)
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
	}

	# average over the folds
	return(-cvLoss / length(folds))
}

.kcvBLMloss_PleqN_noUP_noGP_alt <- function(lambda, Y, X, folds, target, 
                                           minSuccDiff, maxIter){

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
		lpOld <- rep(0, nrow(X))[-folds[[k]]]

		for (iter in 1:maxIter){	
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lpOld))
			W0 <- Ypred * (1 - Ypred)
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- .Machine$double.eps
			}

			# obtain the adjusted response
			WZ <- W0 * lpOld + Y[-folds[[k]]] - Ypred  

			# evaluate X^t W X and add lambda to its diagonal
			XWZpT     <- lambda * target + crossprod(X[-folds[[k]], ], WZ)  
			XTX       <- crossprod(sweep(X[-folds[[k]], ], 1, sqrt(W0), FUN="*"))
			diag(XTX) <- diag(XTX) + lambda 
		
			# now obtain the IRWLS update of the linear predictor efficiently
			lpAll  <- tcrossprod(X, t(solve(XTX, XWZpT)))
			lpOld  <- lpAll[-folds[[k]]]
			lpNew  <- lpAll[ folds[[k]]]

			# assess convergence
			loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld)
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
	}

	# average over the folds
	return(-cvLoss / length(folds))
}


.kcvBLMloss_noUP_noGP <- function(lambda, Y, X, folds, target, minSuccDiff, maxIter){

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
	if (ncol(X) >= (nrow(X) - max(unlist(lapply(folds, length))))){ 
		XXT     <- tcrossprod(X) / lambda
		Xtarget <- t(tcrossprod(target, X))
	}

	for (k in 1:length(folds)){
		# convergence criterion
		loglikPrev <- -10^(10)

		# evaluate the initial linear predictor
		lpOld <- rep(0, nrow(X))[-folds[[k]]]

		for (iter in 1:maxIter){
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lpOld))
			W0 <- Ypred * (1 - Ypred)
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- .Machine$double.eps
			}

			# obtain the adjusted response
			WZ <- W0 * lpOld + Y[-folds[[k]]] - Ypred  

			# distinguish between the high- and low-dimensional cases
			if (ncol(X) >= (nrow(X) - max(unlist(lapply(folds, length))))){ 
				# obtain the IRWLS update of the linear predictors efficiently
				XXTZpT <- tcrossprod((1/W0) * (WZ - W0 * Xtarget[-folds[[k]]]), 
				                     solve(XXT[-folds[[k]], -folds[[k]]] + diag(1/W0)))
				lpAll  <- Xtarget + tcrossprod(XXT[,-folds[[k]],drop=FALSE], XXTZpT)
			} else {
				# obtain the IRWLS update of the linear predictors efficiently
				XWZpT     <- lambda * target + crossprod(X[-folds[[k]], ], WZ)  
				XTX       <- crossprod(sweep(X[-folds[[k]], ], 1, sqrt(W0), FUN="*"))
				diag(XTX) <- diag(XTX) + lambda 
				lpAll     <- tcrossprod(X, t(solve(XTX, XWZpT)))
			}
			lpOld  <- lpAll[-folds[[k]]]
			lpNew  <- lpAll[ folds[[k]]]

			# assess convergence
			loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld)
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
	}

	# average over the folds
	return(-cvLoss / length(folds))
}

.optPenaltyBLM_noUP_noGP.kCVauto <- function(Y, X, lambdaInit, folds, 
                                             target, lambdaMin, 
                                             minSuccDiff, maxIter, implementation){

	## ---------------------------------------------------------------------
	## Function finds the optimal penal11 parameter of the targeted ridge 
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
	if (implementation == "org"){
		optLambda <- optim(lambdaInit, .kcvBLMloss_noUP_noGP, Y=Y, X=X, 
		                   folds=folds, target=target, method="Brent",
		                   minSuccDiff=minSuccDiff, maxIter=maxIter,
		                   lower=lambdaMin, upper=10^(10))$par
	}
	if (implementation == "alt"){
		if (nrow(X) < ncol(X)){
			XXT       <- tcrossprod(X)
			Xtarget   <- t(tcrossprod(target, X))
			optLambda <- optim(lambdaInit, 
			                   .kcvBLMloss_PgeqN_noUP_noGP_alt, 
			                   Y=Y, XXT=XXT, Xtarget=Xtarget, 
			                   folds=folds, method="Brent",
			                   minSuccDiff=minSuccDiff, maxIter=maxIter,
			                   lower=lambdaMin, upper=10^(10))$par
		}
		if (nrow(X) >= ncol(X)){
			optLambda <- optim(lambdaInit, 
			                   .kcvBLMloss_PleqN_noUP_noGP_alt, 
			                   Y=Y, X=X, folds=folds, target=target, 
			                   method="Brent", 
			                   minSuccDiff=minSuccDiff, maxIter=maxIter,
			                   lower=lambdaMin, upper=10^(10))$par
		}
	}
	return(optLambda)
}



.kcvBLMloss_UP_GP <- function(lambdas, Y, X, U, Dg, target, folds, 
                              minSuccDiff, maxIter){

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
	cvLoss  <- 0
	lambda  <- lambdas[1]
	lambdaG <- lambdas[2]
			
	# evaluate subexpressions of the estimator
	if (ncol(X) >= nrow(X)){ 
		Dg          <- qr.solve(diag(rep(lambda, ncol(X))) + 
		                        lambdaG * Dg, t(X), 
		                        tol=.Machine$double.eps)
		XDXT        <- crossprod(t(X), Dg)
		diagXDXTorg <- diag(XDXT)
	}
	if (ncol(X) < nrow(X)){ 
		Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
	}
	Xtarget <- t(tcrossprod(target, X))

	for (k in 1:length(folds)){
		# convergence criterion
		loglikPrev <- -10^(10)
		penalty    <-  0

		# evaluate the initial linear predictor
		lpOld <- rep(0, length(Y))[-folds[[k]]]
		# lpOld <- Xtarget[-folds[[k]]]

		for (iter in 1:maxIter){
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lpOld))
			W0    <- Ypred * (1 - Ypred)
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- 
				              .Machine$double.eps
			}
			if (ncol(X) >= nrow(X)){ 
				# adjusted response
				Z <- lpOld + (Y[-folds[[k]]] - Ypred)/W0

				# evaluate unpenalized low-dim regression estimator
				diag(XDXT)[-folds[[k]]] <- diag(XDXT)[-folds[[k]]] + 1/W0
				slh        <- solve(XDXT[-folds[[k]], -folds[[k]]], 
				                       Z - Xtarget[-folds[[k]]])
				gHat       <- solve(crossprod(U[-folds[[k]],], 
				                    solve(XDXT[-folds[[k]], -folds[[k]]], 
				                          U[-folds[[k]],])),
					            crossprod(U[-folds[[k]],], slh))
				Ug         <- tcrossprod(U, t(gHat))

				# evaluate linear predictor without evaluating the estimator
				slh        <- solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE], 
				                       Z - Xtarget[-folds[[k]]] - Ug[-folds[[k]],])
				diag(XDXT)[-folds[[k]]] <- diagXDXTorg[-folds[[k]]]			
				lpAll      <- crossprod(XDXT[-folds[[k]], ,drop=FALSE], slh)
				penalty    <- 0.5 * sum(lpAll[-folds[[k]]] * slh)
				lpAll      <- Xtarget + lpAll + Ug
			}

			if (ncol(X) < nrow(X)){ 		
				# adjusted response
				Z <- W0*lpOld + (Y[-folds[[k]]] - Ypred) - 
				     W0*Xtarget[-folds[[k]]]

				# evaluate subexpressions of the estimator
				XTXpD <- crossprod(sweep(X[-folds[[k]], , drop=FALSE], 
				                         1, sqrt(W0), FUN="*")) 
				tUTX  <- crossprod(X[-folds[[k]], , drop=FALSE], 
				                   sweep(U[-folds[[k]], , drop=FALSE], 
				                         1, W0, FUN="*"))
				XTZ   <- crossprod(X[-folds[[k]], , drop=FALSE], Z)
				UTZ   <- crossprod(U[-folds[[k]], , drop=FALSE], Z)
				UTU   <- crossprod(sweep(U[-folds[[k]], , drop=FALSE], 
				                         1, sqrt(W0), FUN="*"))
				XTXpD <- XTXpD + Dg
	
				# evaluate unpenalized low-dim regression estimator
				gHat <- as.numeric(
				        solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
				             (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
				                    crossprod(X[-folds[[k]], , drop=FALSE], 
				                    sweep(U[-folds[[k]], , drop=FALSE], 
				                          1, W0, FUN="*")))))))
				Ug   <- as.numeric(tcrossprod(gHat, U))

				# evaluate linear predictor
				XWZpT <- crossprod(Dg, target) + 
				         as.numeric(
				         crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld + 
				         (Y[-folds[[k]]] - Ypred) - W0*Ug[-folds[[k]]]))
				bHat    <- qr.solve(XTXpD, XWZpT, 
				                    tol=.Machine$double.eps)
				lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug
				penalty <- 0.5 * sum(crossprod(Dg, bHat- target) * 
				                                  (bHat - target))
			}
			
			# split linear predictor by fold
			lpOld      <- lpAll[-folds[[k]]]
			lpNew      <- lpAll[ folds[[k]]]

			# assess convergence
			loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld) - penalty
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
	}

	# average over the folds
	return(-cvLoss / length(folds))
}


.kcvBLMloss_UP_noGP <- function(lambda, Y, X, U, target, folds, minSuccDiff, maxIter){

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
	cvLoss  <- 0
			
	# evaluate subexpressions of the estimator
	if (ncol(X) >= nrow(X)){ 
		XXT     <- tcrossprod(X)/lambda
	}
	Xtarget <- t(tcrossprod(target, X))

	for (k in 1:length(folds)){
		# convergence criterion
		loglikPrev <- -10^(10)

		# evaluate the initial linear predictor
		lpOld <- rep(0, length(Y))[-folds[[k]]]

		for (iter in 1:maxIter){
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lpOld))
			W0    <- Ypred * (1 - Ypred)
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- 
				              .Machine$double.eps
			}
			if (ncol(X) >= nrow(X)){ 
				# adjusted response
				Z <- lpOld + (Y[-folds[[k]]] - Ypred)/W0

				# evaluate unpenalized low-dim regression estimator
				diag(XXT)[-folds[[k]]] <- diag(XXT)[-folds[[k]]] + 1/W0
				slh        <- solve(XXT[-folds[[k]], -folds[[k]]], 
				                    Z - Xtarget[-folds[[k]]])
				gHat       <- solve(crossprod(U[-folds[[k]],], 
				                    solve(XXT[-folds[[k]], -folds[[k]]], 
				                          U[-folds[[k]],])),
					            crossprod(U[-folds[[k]],], slh))
				Ug         <- tcrossprod(U, t(gHat))

				# evaluate linear predictor without evaluating the estimator
				slh        <- solve(XXT[-folds[[k]], -folds[[k]], drop=FALSE], 
				                    Z - Xtarget[-folds[[k]]] - Ug[-folds[[k]],])
				diag(XXT)[-folds[[k]]] <- diag(XXT)[-folds[[k]]] - 1/W0				
				lpAll      <- crossprod(XXT[-folds[[k]], ,drop=FALSE], slh)
				penalty    <- sum(lpAll[-folds[[k]]] * slh) / 2
				lpAll      <- Xtarget + lpAll + Ug
			}
			if (ncol(X) < nrow(X)){ 
				# adjusted response
				Z <- W0*lpOld + (Y[-folds[[k]]] - Ypred) - W0*Xtarget[-folds[[k]]]

				# evaluate subexpressions of the estimator
				XTXpD <- crossprod(sweep(X[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*")) 
				tUTX  <- crossprod(X[-folds[[k]], , drop=FALSE], 
				sweep(U[-folds[[k]], , drop=FALSE], 1, W0, FUN="*"))
				XTZ   <- crossprod(X[-folds[[k]], , drop=FALSE], Z)
				UTZ   <- crossprod(U[-folds[[k]], , drop=FALSE], Z)
				UTU   <- crossprod(sweep(U[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*"))
				diag(XTXpD) <- diag(XTXpD) + lambda
	
				# evaluate unpenalized low-dim regression estimator
				gHat <- as.numeric(
				        solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
				             (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
				                    crossprod(X[-folds[[k]], , drop=FALSE], 
				                    sweep(U[-folds[[k]], , drop=FALSE], 
				                          1, W0, FUN="*")))))))
				Ug <- as.numeric(tcrossprod(gHat, U))

				# evaluate linear predictor
				XWZpT <- lambda * target + 
				         as.numeric(
				         crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld + 
				         (Y[-folds[[k]]] - Ypred) - W0*Ug[-folds[[k]]]))
  				bHat    <- qr.solve(XTXpD, XWZpT, 
 		                                   tol=.Machine$double.eps)
				lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug				
				penalty <- 0.5 * lambda * sum((bHat - target)^2) 
			}

			# split linear predictor by fold
			lpOld      <- lpAll[-folds[[k]]]
			lpNew      <- lpAll[ folds[[k]]]

			# assess convergence
			loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld) - penalty
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
	}

	# average over the folds
	return(-cvLoss / length(folds))
}


.kcvBLMloss_noUP_GP <- function(lambdas, Y, X, Dg, target, folds, minSuccDiff, maxIter){

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
	cvLoss  <- 0
	lambda  <- lambdas[1]
	lambdaG <- lambdas[2]
			
	# evaluate subexpressions of the estimator
	if (ncol(X) >= nrow(X)){ 
		Dg      <- qr.solve(diag(rep(lambda, ncol(X))) + 
		                    lambdaG * Dg, t(X), 
		                    tol=.Machine$double.eps)
		XDXT    <- crossprod(t(X), Dg)
		diagXDXTorg <- diag(XDXT)
	}
	if (ncol(X) < nrow(X)){ 
		Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
	}
	Xtarget <- t(tcrossprod(target, X))


	for (k in 1:length(folds)){
		# convergence criterion
		loglikPrev <- -10^(10)
		penalty    <- 0

		# evaluate the initial linear predictor
		lpOld <- rep(0, length(Y))[-folds[[k]]]
		if (.loglikBLMlp(Y[-folds[[k]]], lpOld) < 
		    .loglikBLMlp(Y[-folds[[k]]], Xtarget[-folds[[k]]])){
			lpOld <- Xtarget[-folds[[k]]]
		}

		for (iter in 1:maxIter){
			# calculate the weights
			Ypred <- 1 / (1 + exp(-lpOld))
			W0    <- Ypred * (1 - Ypred)
			if (min(W0) <= .Machine$double.eps){
				W0[which(W0 < .Machine$double.eps)] <- 
				              .Machine$double.eps
			}
			if (ncol(X) >= nrow(X)){ 
				# adjusted response
				Z <- lpOld + (Y[-folds[[k]]] - Ypred)/W0

				# evaluate linear predictor without evaluating the estimator
				diag(XDXT)[-folds[[k]]] <- diag(XDXT)[-folds[[k]]] + 1/W0
				slh        <- solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE], 
				                     Z - Xtarget[-folds[[k]]])
				# slh        <- qr.solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE],
		               #                        Z - Xtarget[-folds[[k]]], 
				#                        tol=sqrt(.Machine$double.eps))


				diag(XDXT)[-folds[[k]]] <- diagXDXTorg[-folds[[k]]]
				lpAll      <- crossprod(XDXT[-folds[[k]],,drop=FALSE], slh)
				penalty    <- sum(lpAll[-folds[[k]]] * slh) / 2
				lpAll      <- Xtarget + lpAll
			}
			if (ncol(X) < nrow(X)){ 
				# obtain part one of the estimator
				XWZpT <- crossprod(Dg, target) + 
				         as.numeric(
				         crossprod(X[-folds[[k]], , drop=FALSE], 
				                   W0*lpOld + 
				                   (Y[-folds[[k]]] - Ypred)))

				# evaluate subexpressions of the estimator
				XTXpD <- crossprod(sweep(X[-folds[[k]], , drop=FALSE], 
				                         1, sqrt(W0), FUN="*")) 
				XTXpD <- XTXpD + Dg
	
				# evaluate linear predictor
				bHat    <- qr.solve(XTXpD, XWZpT, 
 		                                   tol=.Machine$double.eps)
				lpAll   <- as.numeric(tcrossprod(X, t(bHat)))
				penalty <- 0.5 * sum(crossprod(Dg, bHat - target) * 
				                                  (bHat - target)) 
			}

			# split linear predictor by fold
			lpOld      <- lpAll[-folds[[k]]]
			lpNew      <- lpAll[ folds[[k]]]

			# assess convergence
			loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld) - penalty
			if (abs(loglik - loglikPrev) < minSuccDiff){ 
				break 
			} else {
				loglikPrev <- loglik
			}
		}
		cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
	}

	# average over the folds
	return(-cvLoss / length(folds))
}




################################################################################
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Functions for cross-validation of the GLM                                  ##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
################################################################################

optPenaltyGLM.kCVauto <- function(Y, 
                                  X, 
                                  U=matrix(ncol=0, nrow=length(Y)), 
                                  lambdaInit,
                                  lambdaGinit=0,
                                  Dg=matrix(0, ncol=ncol(X), nrow=ncol(X)),
                                  model="linear", 
                                  target=rep(0, ncol(X)), 
				   folds=makeFoldsGLMcv(min(10, length(X)), 
				                        Y, model=model),
                                  loss="loglik", 
                                  lambdaMin=10^(-5), 
                                  lambdaGmin=10^(-5),
                                  minSuccDiff=10^(-5), 
                                  maxIter=100, 
                                  implementation="org"){

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
	## folds       : list with folds
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
	## implementation : A character specifying the implementation to be 
	##                  used, "org" or "alt" (original or alternative)
	## ---------------------------------------------------------------------
	## Value:
	## A numeric, the optimal cross-validated penalty parameter of the
	## ridge estimator of the generalized linear regression model parameter
	## ---------------------------------------------------------------------
	## Authors : Wessel N. van Wieringen
	## ---------------------------------------------------------------------

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###

	if (((max(abs(Dg)) == 0) | lambdaGinit == 0) && (ncol(U) == 0)){
		# search by the Brent algorithm
		if (model == "linear"){
			return(.optPenaltyLM_noUP_noGP.kCVauto(Y, X, 
			                             lambdaInit, folds, target, 
			                             loss, lambdaMin, 
			                             implementation))
		}
		if (model == "logistic"){
			return(.optPenaltyBLM_noUP_noGP.kCVauto(Y, X, 
			                      lambdaInit, folds, target,
			                      lambdaMin, minSuccDiff, 
			                      maxIter, implementation))
		}
	}
	

	####----------------------------------------------------------------###
	####----------------------------------------------------------------###


	if (((max(abs(Dg)) == 0) | lambdaGinit == 0) && (ncol(U) != 0)){
		if (model == "linear"){
			if (ncol(X) >= nrow(X)){ 
				lambdasOpt <- optim(lambdaInit, 
			        	            .kcvLMloss_PgeqN_UP_noGP,
			        	            Y=Y,
			        	            X=X,
			        	            U=U,
			        	            target=target,
			        	            folds=folds,
			        	            method="Brent",
			        	            lower=lambdaMin,
			        	            upper=10^(10),
			        	            loss=loss)$par
				return(lambdasOpt)
			}
			if (ncol(X) < nrow(X)){ 
				lambdasOpt <- optim(lambdaInit, 
			        	            .kcvLMloss_PleqN_UP_GPandNoGP,
			        	            Y=Y,
			        	            X=X,
			        	            U=U,
			        	            Dg=lambdaGinit*Dg,
			        	            target=target,
			        	            folds=folds,
			        	            method="Brent",
			        	            lower=lambdaMin,
			        	            upper=10^(10),
			        	            loss=loss)$par
				return(lambdasOpt)
			}
		}
		if (model == "logistic"){
			lambdasOpt <- optim(lambdaInit,
			                    .kcvBLMloss_UP_noGP,
			                    Y=Y, 
			                    X=X, 
			                    U=U,
			                    target=target,
			                    folds=folds, 
			                    method="Brent",
			                    lower=lambdaMin,
			                    upper=10^(10),
			                    maxIter=maxIter,
			                    minSuccDiff=minSuccDiff)$par
			return(lambdasOpt)
		}
	}


	####----------------------------------------------------------------###
	####----------------------------------------------------------------###


	if (((max(abs(Dg)) != 0) && lambdaGinit != 0) && (ncol(U) == 0)){
		if (model == "linear"){
			if (ncol(X) >= nrow(X)){ 
				Dg <- eigen(Dg)
				lambdasOpt <- constrOptim(c(lambdaInit, lambdaGinit),
				                          .kcvLMloss_PgeqN_UP_GP,
				                          grad=NULL,
				                          Y=Y,
				                          X=X %*% Dg$vectors,
				                          U=U,
				                          Ds=Dg$values,
				                          Xtarget=tcrossprod(X, 
				                                  matrix(target, nrow=1)),
				                          folds=folds,
				                          ui=diag(2),
				                          ci=c(lambdaMin, lambdaGmin),
				                          loss=loss)$par
				return(lambdasOpt)
			}
			if (ncol(X) < nrow(X)){ 
				lambdasOpt <- constrOptim(c(lambdaInit, lambdaGinit),
				                          .kcvLMloss_PleqN_UP_GPandNoGP,
				                          grad=NULL,
				                          Y=Y,
				                          X=X,
				                          U=U,
				                          Dg=Dg,
				                          Xtarget=tcrossprod(X, 
				                                  matrix(target, nrow=1)),
				                          folds=folds,
				                          ui=diag(2),
				                          ci=c(lambdaMin, lambdaGmin),
				                          loss=loss)$par
				return(lambdasOpt)
			}
		}
		if (model == "logistic"){
			lambdasOpt <- constrOptim(c(lambdaInit, lambdaGinit),
			                          .kcvBLMloss_noUP_GP,
			                          grad=NULL,
			                          Y=Y,
			                          X=X,
			                          Dg=Dg,
			                          target=target,
			                          folds=folds,
			                          ui=diag(2),
			                          ci=c(lambdaMin, lambdaGmin),
			                          maxIter=maxIter,
			                          minSuccDiff=minSuccDiff)$par
			return(lambdasOpt)
		}
	}


	####----------------------------------------------------------------###
	####----------------------------------------------------------------###


	if (((max(abs(Dg)) != 0) && lambdaGinit != 0) && (ncol(U) != 0)){
		if (model == "linear"){
			if (ncol(X) >= nrow(X)){ 
				Dg         <- eigen(Dg)
				lambdasOpt <- constrOptim(c(lambdaInit, lambdaGinit),
				                          .kcvLMloss_PgeqN_UP_GP,
				                          grad=NULL,
				                          Y=Y,
				                          X=X %*% Dg$vectors,
				                          U=U,
				                          Ds=Dg$values,
				                          Xtarget=tcrossprod(X, 
				                                  matrix(target, nrow=1)),
				                          folds=folds,
				                          ui=diag(2),
				                          ci=c(lambdaMin, lambdaGmin),
				                          loss=loss)$par
				return(lambdasOpt)
			}
			if (ncol(X) < nrow(X)){ 
				lambdasOpt <- constrOptim(c(lambdaInit, lambdaGinit),
				                          .kcvLMloss_PleqN_UP_GPandNoGP,
				                          grad=NULL,
				                          Y=Y,
				                          X=X,
				                          U=U,
				                          Dg=Dg,
				                          target=target,
				                          folds=folds,
				                          ui=diag(2),
				                          ci=c(lambdaMin, lambdaGmin),
				                          loss=loss)$par
				return(lambdasOpt)
			}
		}
		if (model == "logistic"){
			lambdasOpt <- constrOptim(c(lambdaInit, lambdaGinit),
			                          .kcvBLMloss_UP_GP,
			                          grad=NULL,
			                          Y=Y,
			                          X=X,
			                          U=U,
			                          Dg=Dg,
			                          target=target,
			                          folds=folds,
			                          ui=diag(2),
			                          ci=c(lambdaMin, lambdaGmin),
			                          maxIter=maxIter,
			                          minSuccDiff=minSuccDiff)$par
			return(lambdasOpt)
		}
	}

	
}


################################################################################
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Functions for cross-validation of multi-targeted ridge regression          ##
## penalty parameters                                                         ##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
################################################################################

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
	                                                  lambda=lambdaTotal, target=mixedTarget, model="linear"))
		}
		if (loss=="loglik"){
			cvLoss <- cvLoss - .loglikLM(Y[folds[[k]]], 
	                                             X[folds[[k]],,drop=FALSE], 
	                                             ridgeGLM(Y[-folds[[k]]], 
		                                             X[-folds[[k]],,drop=FALSE], 
	                                                     lambda=lambdaTotal, target=mixedTarget, model="linear"))
		}
	}
	
	# average over the folds
	return(cvLoss / length(folds))
}

.kcvLMlossPgeqNorg_multiT <- function(lambdas, Y, X, folds, targetMat, loss="loglik"){

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
	Xtarget <- t(tcrossprod(mixedTarget, X)) 
	XXT     <- tcrossprod(X)/lambdaTotal

	# evaluate loss per fold
	cvLoss <- 0
	for (k in 1:length(folds)){
		# Yhat   <- Xtarget[folds[[k]]] + 
		#          tcrossprod(XXT[folds[[k]], -folds[[k]],drop=FALSE],
		#                     t(crossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
		#                     XXT[-folds[[k]], -folds[[k]]]),
		#                     (matrix(Y[-folds[[k]]], ncol=1) - Xtarget[-folds[[k]]]))))

		Yhat   <- Xtarget[folds[[k]]] + 
		          tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
		                           XXT[-folds[[k]], -folds[[k]]], 
		                           Y[-folds[[k]]] - Xtarget[-folds[[k]]]), 
		                     XXT[folds[[k]], -folds[[k]], drop=FALSE])

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
	return(.kcvBLMloss_noUP_noGP(lambdaTotal, Y, X, folds, mixedTarget, minSuccDiff, maxIter))
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
					  .kcvLMlossPgeqNorg_multiT, 
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




optPenaltyGLMmultiT.kCVauto <- function(Y, 
                                        X, 
                                        lambdaInit, 
                                        model="linear", 
				         targetMat,                                        
                                        folds=makeFoldsGLMcv(min(10, length(X)), 
				                              Y, model=model),
				         loss="loglik", 
				         lambdaMin=10^(-5),
				         minSuccDiff=10^(-5), 
				         maxIter=100){    

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

	# search by the Brent algorithm
	if (model == "linear"){
		return(.optPenaltyLMmultT.kCVauto(Y, X, lambdaInit, folds, targetMat, loss))
	}
	if (model == "logistic"){
		return(.optPenaltyBLMmultT.kCVauto(Y, X, lambdaInit, folds, targetMat, 
		                                   minSuccDiff, maxIter))
	}
}



