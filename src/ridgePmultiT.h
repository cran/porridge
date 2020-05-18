/* ---------------------------------------------------------------------------
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Functions related to:
van Wieringen, W.N., Stam, K.A., Peeters, C.F.W., van de Wiel, M.A. (2020), 
"Updating of the Gaussian graphical model through targeted penalized estimation", 
Journal of Multivariate Analysis, Volume 178, Article 104621
------------------------------------------------------------------------------
------------------------------------------------------------------------------
--------------------------------------------------------------------------- */


// [[Rcpp::export(".armaMixTargets")]]
arma::mat mixTargets(const Rcpp::List targetList,
		     const arma::vec lambda){

	/* ---------------------------------------------------------------------------
	## Function that computes a weighted average of targets
	## ---------------------------------------------------------------------------
	## Arguments
	## Tlist  : A list of target matrices
	## lambda : The penalty matrix
	## ---------------------------------------------------------------------------
	## Value
	## mixedTarget	: matrix, a weighted average of input matrices.
	## ---------------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------------- */

	// extract number of targets
	const int nTargets = lambda.n_elem;
	arma::mat target   = targetList[1];
	const int p        = target.n_cols;

	// reformulate lambda vector to weights
	const double lambdaTotal = arma::sum(lambda);
	arma::vec weights        = lambda / lambdaTotal;
		
	// define 
	arma::mat mixedTarget = arma::zeros(p, p);
	for (int k = 0; k < nTargets; k++){
		arma::mat target = targetList[k];
		mixedTarget      = mixedTarget + weights[k] * target;
	}
	
	// return ridge precision estimate
	return mixedTarget;
}



// [[Rcpp::export(".armaRidgePmultiT")]]
arma::mat armaRidgePmultiT(const arma::mat&  S,
                      arma::vec         lambda,
                      const Rcpp::List& targetList){
	/* ---------------------------------------------------------------------------
	## Function that computes the ridge precision matrix estimate. Essentially, it 
	## first computes a weighted average of targets, followed by a call to 
	## 'armaRidgePanyTarget' (a modified version of .... from the rags2ridges-package.
	## ---------------------------------------------------------------------------
	## Arguments
	## S          : A sample covariance matrix. Should not contain missings.
	## lambda     : A vector of penalty parameters
	## targetList : A list of target matrices
	## ---------------------------------------------------------------------------
	## Value
	## A matrix, a ridge precision matrix estimate.
	## ---------------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------------- */


	// extract number of targets
	const int nTargets = lambda.n_elem;
	arma::mat target   = targetList[1];
	const int p        = target.n_cols;

	// reformulate lambda vector to weights
	const double lambdaTotal = arma::sum(lambda);
	arma::vec weights        = lambda / lambdaTotal;
		
	// define 
	arma::mat mixedTarget    = arma::zeros(p, p);
	for (int k = 0; k < nTargets; k++){
		arma::mat target = targetList[k];
		mixedTarget      = mixedTarget + weights[k] * target;
	}
	
	// return ridge precision estimate
	return armaRidgePanyTarget(S, mixedTarget, lambdaTotal, 2);
}

