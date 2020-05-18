/* ---------------------------------------------------------------------------
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Below are functions related to
Aflakparast, M., de Gunst, M.C.M., van Wieringen, W.N. (2018) "Reconstruction 
of molecular network evolution from cross-sectional omics data", Biometrical 
Journal, 60(3), 547-563.
------------------------------------------------------------------------------
------------------------------------------------------------------------------
--------------------------------------------------------------------------- */


inline arma::mat armaRidgePanyTarget(const arma::mat & S,
                              const arma::mat & target,
                              const double lambda,
                              int invert = 2) {
	/*---------------------------------------------------------------------------------------------------
	## Compute the ridge estimate for general/arbitrary targets.
	## Depending on the value of "invert"" using matrix inversion (via
	## diagonalization) or avoiding it.
	## --------------------------------------------------------------------------------------------------
	## Arguments: 
	## S      : A sample covariance matrix. Should not contain NAs, Infs, or NaNs!
	## target : The target matrix with the same size as S
	## lambda : The the ridge penalty
	## invert : integer. Should the estimate be compute using inversion?
	##          0 = "no", 1 = "yes", 2 = "automatic" (default).
	## --------------------------------------------------------------------------------------------------
	## Value
	## Ridge covariance estimate
	## --------------------------------------------------------------------------------------------------
	## Authors      : Anders Bilgrau & Carel F.W. Peeters, 
	##                modified for current purpose by Wessel N. van Wieringen
	--------------------------------------------------------------------------------------------------*/

	// eigendecomposition of sample covariance matrix minus target
	arma::vec eigvals;
	arma::mat eigvecs = S - lambda*target;
	if (!eigvecs.is_finite()) {
		return target;
	}
	arma::eig_sym(eigvals, eigvecs, eigvecs, "dc");
	eigvals = 0.5*eigvals;
	arma::vec sqroot = sqrt(lambda + pow(eigvals, 2.0));

	// return target if shrunken evals are infinite and lambda is "large"
	// usually happens for lambda >= 1e154
	if (lambda > 1e6 && (!eigvals.is_finite() || !sqroot.is_finite())) {
		return target;
	}

	// determine to invert or not
	arma::vec D_inv = 1.0/(sqroot + eigvals); 
	if (invert == 2) { 
		if (lambda > 1) {  
			// generally, don't use inversion for "large" lambda
			invert = 0;
		} else {
			if (!D_inv.is_finite()) {
				invert = 0;
			} else {
				invert = 1;
			}
		}
	}

	// invert in most efficient way
	if (invert == 1) {
		eigvecs.each_row() %= arma::trans(arma::sqrt(D_inv));
	} else {
		eigvecs.each_row() %= arma::trans(arma::sqrt((sqroot - eigvals))/lambda);
	}
	return eigvecs * arma::trans(eigvecs);
}

inline arma::mat armaRidgePscalarTarget(const arma::mat& S,
                                 const double     alpha,
                                 const double     lambda,
                                 int              invert=2) {
	/*---------------------------------------------------------------------------------------------------
	## Compute the ridge estimate for rotational equivariant targets.
	## Depending on the value of "invert"" using matrix inversion (via
	## diagonalization) or avoiding it.
	## --------------------------------------------------------------------------------------------------
	## Arguments: 
	## S      : A sample covariance matrix. Should not contain NAs, Infs, or NaNs!
	## alpha  : The scaling of the identity matrix. Shoud not contain NaNs, Infs,
        ##          or NA.s
	## lambda : The the ridge penalty
	## invert : integer. Should the estimate be compute using inversion?
	##          0 = "no", 1 = "yes", 2 = "automatic" (default).
	## --------------------------------------------------------------------------------------------------
	## Value
	## Ridge covariance estimate
	## --------------------------------------------------------------------------------------------------
	## Authors      : Anders Bilgrau & Carel F.W. Peeters, 
	##                modified for current purpose by Wessel N. van Wieringen
	--------------------------------------------------------------------------------------------------*/

	// eigen decomposition of sample covariance matrix
	arma::vec eigvals;
	arma::mat eigvecs;
	arma::eig_sym(eigvals, eigvecs, S, "dc");

	// spank eigenvalues
	eigvals          = 0.5*(eigvals - lambda*alpha);
	arma::vec sqroot = arma::sqrt(lambda + arma::pow(eigvals, 2.0));

	// Return target if shrunken evals are infinite and lambda is "large"
	// Usually happens for lambda >= 1e154
	if (lambda > 1e6 && (!eigvals.is_finite() || !sqroot.is_finite())) {
		const int p = S.n_rows;
		return alpha*arma::eye<arma::mat>(p, p);
	}

	// determine to invert or not
	arma::vec D_inv = 1.0/(sqroot + eigvals); 
	if (invert == 2) {   
		if (lambda > 1) {  
			// generally, don't use inversion for "large" lambda
			invert = 0;
		} else {
			if (!D_inv.is_finite()) {
				invert = 0;
			} else {
				invert = 1;
			}
		}
	}

	// invert in most efficient way
	if (invert == 1) {
 		eigvecs.each_row() %= arma::trans(arma::sqrt(D_inv));
	} else {
 		eigvecs.each_row() %= arma::trans(arma::sqrt((sqroot - eigvals))/lambda);
	}
	return eigvecs * arma::trans(eigvecs);
}



inline arma::mat armaRidgeP(const arma::mat& S,
                            const arma::mat& target,
                            const double     lambda,
                            int              invert=2){

	/*---------------------------------------------------------------------------------------------------
	## The ridge estimator in C++. Wrapper for the subroutines
	## --------------------------------------------------------------------------------------------------
	## Arguments: 
	## S      : The sample covariance matrix (a numeric matrix on the R side)
	## target : Target matrix (a numeric matrix on the R side, same size as S)
	## lambda : The penalty (a numeric of length one on the R side)
	## invert : Should the estimate be compute using inversion?
        ##          0 = "no", 1 = "yes", 2 = "automatic", (default).
	## --------------------------------------------------------------------------------------------------
	## Value
	## Ridge covariance estimate
	## --------------------------------------------------------------------------------------------------
	## Authors      : Anders Bilgrau & Carel F.W. Peeters, 
	##                modified for current purpose by Wessel N. van Wieringen
	--------------------------------------------------------------------------------------------------*/

	// catch nonpositive penalty parameters
	if (lambda <= 0) {
		Rcpp::stop("The penalty (lambda) must be strictly postive");
	}

	// infinite penalty parameter: return target
	if (lambda == arma::datum::inf) {
		return target;
	}

	// create scalar target
	const int p          = S.n_rows;
	arma::vec scalTdiag  = target(0,0) * arma::ones(p);

	// assess whether target is scalar and use corresponding subroutine
	if (arma::all(arma::all(target == arma::diagmat(scalTdiag)))) {
		return armaRidgePscalarTarget(S,  target(0,0), lambda, invert);
  	} else {
    		return armaRidgePanyTarget(S, target, lambda, invert);
	}
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export(".armaRidgeGGMmixture")]]
Rcpp::List ridgeGGMmixture(const arma::mat& Y, 
                           const int        K, 
                           const double     lambda,
                           const arma::mat& target, 
                           const arma::mat& iWeights,
                           const int&       nInit, 
                           const double&    minSuccDiff,
                           const double&    minMixProp){

	/* ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## K		: number of mixture components.
	## lambda	: tuning parameter.
	## target	: target precision matrix.
	## iWeights	: sample-specific positive component weights 
	##                (may be larger than one as they are rescaled to sum to one)
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized 
	##                loglikehood) to be achieved
	## minMixProp	: smallest mixing probability tolerated.
	## ---------------------------------------------------------------------
	## Value
	## A list-object with the following slots:
	## mu		: estimated mean vectors 
	## P		: estimated mixture precision matrices 
	## pi		: estimated mixing probabilities 
	## weights	: estimated component memberships.
	## penLL	: penalized loglikelihood
	## ---------------------------------------------------------------------
	## Authors      : Mehran Aflakparast & Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	// extract number of samples and variates  
  	int n = Y.n_rows;	
	int p = Y.n_cols;
  
	// ensure weights sum to one, row-wise					
	arma::mat cWeights = arma::normalise(arma::trans(iWeights), 1);					
  
	// declare variables used in penalized EM algorithm 
	arma::vec pis;
	arma::mat cPs    = arma::mat(K*p, p);
	arma::mat cMeans = arma::mat(K, p);
	arma::mat cS;
	double maxW;
	double LL;
	double penLL = -10e+10;
	double penLLprev;
	arma::mat Ytilde;
	double val;
	double sign;

	// start iterating	
	bool stayInLoop = TRUE; 
	for (int u = 0; u < nInit && stayInLoop; ++u){
		//
		penLLprev = penLL;

		// update mixing probabilites
 		pis = arma::mean(cWeights, 1);

		// ensure mixing probabilities 
		arma::uvec nonzeroPis; 
		if (pis.min() < minMixProp){
			pis(arma::index_min(pis)) = minMixProp;
			nonzeroPis = arma::find(pis > minMixProp);
			pis(nonzeroPis) = (1-minMixProp) * pis(nonzeroPis) / arma::sum(pis(nonzeroPis));
		}

		// M-step
		for (int k = 0; k < K; ++k){
			// update means, sample covariances and precisions of cluster k
	 	        cMeans.row(k)      = cWeights.row(k) * Y / (n * pis(k));
			Ytilde             = Y.each_row() - cMeans.row(k);
			Ytilde.each_col() %= arma::trans(arma::sqrt(cWeights.row(k)));
			cS                 = (arma::trans(Ytilde) * Ytilde) / (n * pis(k));
			cPs.submat(k*p, 0, (k+1)*p-1, p-1) = armaRidgeP(cS, 
			                                                target, 
			                                                lambda / (n * pis(k)));
		}

		// E-step
		for (int k = 0; k < K; ++k){
			// calculate individual contributions on log-scale
			log_det(val, sign, cPs.submat(k*p, 0, (k+1)*p-1, p-1));
			Ytilde          = Y.each_row() - cMeans.row(k);
			cWeights.row(k) = log(pis(k)) + val / 2 - 
					  p * log(2 * arma::datum::pi) / 2 -
					  arma::trans(arma::sum((Ytilde * 
			                  cPs.submat(k*p, 0, (k+1)*p-1, p-1)) % Ytilde, 1)) / 2;
		}
		// subtract maximum value (to avoid infinity when taking exponential)
		maxW     = arma::max(arma::max(cWeights));
		cWeights = cWeights - maxW;
		// LL       = arma::conv_to<double>::from(arma::sum(
		//            arma::trans(arma::log(arma::sum(arma::exp(cWeights))))) 
		//            + n * maxW);										
		LL       = arma::accu(arma::trans(arma::log(arma::sum(arma::exp(cWeights))))) + n * maxW;										
		cWeights = arma::normalise(arma::exp(cWeights), 1);
	
		// evaluate penalty
		penLL = LL;
		for (int k = 0; k < K; ++k){
			penLL = penLL - lambda * arma::accu(
			                         arma::square(cPs.submat(k*p,       0, 
			                                                (k+1)*p-1, p-1) - target)) / 2;
		}
	        if (std::abs(penLL - penLLprev) < minSuccDiff || penLLprev > penLL){ 
			stayInLoop = FALSE; 
		}
	}
  
	// return stuff as a list-object
	return Rcpp::List::create(Rcpp::Named("mu")      = cMeans, 
                                  Rcpp::Named("P")       = cPs, 
                                  Rcpp::Named("pi")      = pis, 
                                  Rcpp::Named("weights") = cWeights, 
                                  Rcpp::Named("penLL")   = penLL);
}



inline double llGGMmixture(const arma::mat& Y, 
                           Rcpp::List       mixGGMhat){

	/*
	## Loglikehood of a mixture of GGMs
	## --------------------------------------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## mixGGMhat	: parameter estimates of the mixture of GGMs.
	##                (as provided by the rigdeGGMmixture-function)
	## --------------------------------------------------------------------------------------------------
	## Value
	## LL		: negative loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## --------------------------------------------------------------------------------------------------
	## Authors      : Mehran Aflakparast & Wessel N. van Wieringen
	## --------------------------------------------------------------------------------------------------
	*/

	// extract mixture GGM model parameters
	arma::mat mus = Rcpp::as<arma::mat>(mixGGMhat(0));
	arma::mat Ps  = Rcpp::as<arma::mat>(mixGGMhat(1));
	arma::vec pis = Rcpp::as<arma::vec>(mixGGMhat(2));

	// extract number of samples and variates  
  	int p = Y.n_cols;	
  	int n = Y.n_rows;
	int K = pis.n_elem;

	// declare other variables
	double LL = 0;
	double val;
	double sign;
	arma::mat LLparts(K, n);
	arma::mat Ytilde;
	for (int k = 0; k < K; ++k){
		// calculate individual contributions on log-scale
		log_det(val, sign, Ps.submat(k*p, 0, (k+1)*p-1, p-1));
		Ytilde          = Y.each_row() - mus.row(k);
		LLparts.row(k)  = log(pis(k)) + val / 2 - 
					  p * log(2 * arma::datum::pi) / 2 -
					  arma::trans(arma::sum((Ytilde * 
			                  Ps.submat(k*p, 0, (k+1)*p-1, p-1)) % Ytilde, 1)) / 2;
	}

	// subtract maximum value (to avoid infinity when taking exponential)
	double maxW = arma::max(arma::max(LLparts));
	LLparts     = LLparts - maxW;
	LL          = arma::accu(arma::trans(arma::log(arma::sum(arma::exp(LLparts))))) + n * maxW;			
	return(-LL);
}



// [[Rcpp::export(".armaKcvlGGMmixture")]]
double kcvlGGMmixture(double           lambda,
                      const arma::mat& Y, 
                      const int        K, 
                      const arma::mat& target, 
                      const arma::mat& iWeights,
                      const int&       nInit, 
                      const double&    minSuccDiff,
                      Rcpp::List       folds,
                      const double&    minMixProp){

	/*
	## Cross-validated loglikehood of a mixture of GGMs
	## --------------------------------------------------------------------------------------------------
	## Arguments 
	## lambda	: tuning parameter.
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## K		: number of mixture components.
	## target	: target precision matrix.
	## iWeights	: sample-specific positive component weights 
	##                (may be larger than one as they are rescaled to sum to one)
	## nInit	: target precision matrix.
	## minSuccDiff	: target precision matrix.
	## folds        : cross-validation sample splits.
	## minMixProp	: smallest mixing probability tolerated.
	## --------------------------------------------------------------------------------------------------
	## Value
	## cvLL		: negative cross-validated loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## --------------------------------------------------------------------------------------------------
	## Authors      : Mehran Aflakparast & Wessel N. van Wieringen
	## --------------------------------------------------------------------------------------------------
	*/

	// extract number samples and folds
  	int n = Y.n_rows;
	int nFolds = folds.size();

	// declare variables used in the for-loop
	double cvLL = 0;
	arma::uvec fold, allButFold;
	Rcpp::List mixGGMhat;
	int slh;
	for (int f = 0; f < nFolds; ++f){
		// sort out the fold and its complement
		fold = Rcpp::as<arma::uvec>(folds(f)) - 1;
		allButFold = arma::uvec(n-fold.n_elem);
		slh = 0;
		for (int i = 0; i < n; ++i){
			if (all(fold != i)){
				allButFold(slh) = i;
				slh = slh + 1;
			}
		}
		// estimate model from all but left-out fold
		mixGGMhat  = ridgeGGMmixture(Y.rows(allButFold), 
		                       K, 
		                       lambda, 
		                       target, 
		                       iWeights.rows(allButFold),
		                       nInit,
		                       minSuccDiff,
		                       minMixProp);
		// evaluate loglikelihood on left-out fold
		cvLL = cvLL + llGGMmixture(Y.rows(fold), mixGGMhat);
	}
	return(cvLL / nFolds);
}


