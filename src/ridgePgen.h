/* ---------------------------------------------------------------------------
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Below are functions related to:
Van Wieringen, W.N. (2019), "The generalized ridge estimator of the inverse 
covariance matrix", Journal of Computational and Graphical Statistics, 28(4), 
932-942.
------------------------------------------------------------------------------
------------------------------------------------------------------------------
--------------------------------------------------------------------------- */


inline double estEqP11(double    X, 
                       double    B11, 
                       double    L11, 
                       arma::vec BVc, 
                       arma::vec Du){ 
	/* ---------------------------------------------------------------------
	## internal auxillary function:
	## estimating equation of diagonal precision elements
	## ---------------------------------------------------------------------
	## Arguments : not specified
	## ---------------------------------------------------------------------
	## Value     : not specified
	## ---------------------------------------------------------------------
	## Authors   : Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	return 1/X - B11 - L11 * X + arma::accu(BVc % (Du / arma::square(Du + X)) % BVc);  
}


// [[Rcpp::export(".armaRidgePgen")]]
arma::mat ridgePgen(const arma::mat& S, 
                    const arma::mat& lambda, 
                    const arma::mat& target, 
                    const int&       nInit, 
                    const double&    minSuccDiff){
	/* ---------------------------------------------------------------------
	## Generalized ridge estimation of the inverse covariance matrix
	## ---------------------------------------------------------------------
	## Arguments 
	## S		   : sample covariance matrix.
	## lambda	   : penalty parameter matrix.
	## target	   : target precision matrix.
	## nInit	   : maximum number of iterations
	## minSuccDiff : minimum successive difference (in the penalized 
	##               loglikehood) to be achieved
	## ---------------------------------------------------------------------
	## Value
	## Phat		   : penalized precision matrix estimate
	## ---------------------------------------------------------------------
	## Authors     : Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	// declare variables
	arma::vec  eigvalC; 
	arma::mat  eigvecC;
	arma::vec  L12tilde;
	arma::vec  BVc;
	double     lowerZ;
	double     midZ; 	
	double     upperZ;
	double     slh;
	double     lowerF;
	double     midF; 
	double     upperF;
	double     midFtemp;
	double     midZslh; 
	double     midZtemp;
	double     lowerZslh; 
	double     upperZslh;
	double     lowerFslh; 
	double     upperFslh;
	bool       flag; 
	bool       cond1; 
	bool       cond2; 
	bool       cond3; 
	bool       cond4; 
	bool       cond5; 
	arma::uvec indices;
	arma::uvec j2uvec;
	int        counter;

	// initialize precision estimate
	arma::mat Pprev = arma::diagmat(1/arma::diagvec(S)) + target; 
	arma::mat Phat  = Pprev;

	// assemble matrix A
	arma::mat A = S - lambda % target;

	// start iterating
	bool stayInLoop = TRUE;
	for (int k = 0; k < nInit && stayInLoop; ++k) {
		// store current estimate
		Pprev = Phat;

		// start updating columns and rows of precision
		for (unsigned int j = 0; j < S.n_rows; ++j) {
			// index vector 
			//indices = arma::regspace<arma::uvec>(0,  S.n_rows-1)
			indices = arma::find(lambda.col(1) >= 0);
			j2uvec  = find(indices == j);
			indices.shed_row(j);

			// create and manipulate penalty parameter vector
			L12tilde = lambda.col(j);
			L12tilde.shed_row(j);
			L12tilde = 1/ arma::sqrt(L12tilde);
					
			// eigendecomposition of U
			arma::eig_sym(eigvalC, eigvecC, (A.submat(indices, indices) + 
				lambda.submat(indices, indices) % Phat.submat(indices, indices)) % (L12tilde * L12tilde.t()), "dc");

			// construct intermediate object
			BVc = eigvecC.t() * (A.submat(indices, j2uvec) % L12tilde);

			// initiate root finder
			// slh    = arma::accu(arma::square(BVc) / eigvalC);
			slh    = arma::accu(arma::square(BVc) / (eigvalC + 0.01));						
			lowerZ = 1 / (sqrt(lambda(j,j) + (A(j,j) * A(j,j))/4 ) + A(j,j)/2);
			upperZ = 1 / (sqrt(lambda(j,j) + (A(j,j) - slh) * (A(j,j) - slh)/4 ) + (A(j,j)-slh)/2);
			
            		if (lowerZ < 10e-5){ 
		                lowerZ = 10e-5; 
            		}       			
            		if (upperZ > 10e5){ 
		                upperZ = 10e5; 
            		}       			
            		if (upperZ < lowerZ){ 
		                upperZ = 10e5; 
            		}       			
    			lowerF = estEqP11(lowerZ, A(j,j), lambda(j,j), BVc, eigvalC);            
			upperF = estEqP11(upperZ, A(j,j), lambda(j,j), BVc, eigvalC);
			midZ   = lowerZ; 
			midF   = lowerF;
			flag   = TRUE;

			// start root finding
			counter = 0;
			do {						
				// suggested update
				if (std::abs(upperF - midF) > 10e-10 && std::abs(lowerF - midF) > 10e-10){						
					// inverse quadratic update
					midZtemp = lowerZ * upperF * midF   / ((lowerF - upperF) * (lowerF - midF)) +
					           upperZ * lowerF * midF   / ((upperF - lowerF) * (upperF - midF)) + 
					           midZ   * upperF * lowerF / ((midF - upperF) * (midF - lowerF));
				} else {
					// secant update
					midZtemp = upperZ - upperF * (upperZ - lowerZ) / (upperF - lowerF);
				}

				// check conditions for bisection
				if (counter > 0){
				    cond1 = !(midZtemp > (3 * lowerZ + upperZ)/4 && midZtemp < upperZ);
				    cond2 = ( flag && (std::abs(midZtemp - upperZ) >= std::abs(upperZ - midZ)/2));
				    cond3 = (!flag && (std::abs(midZtemp - upperZ) >= std::abs(midZslh - midZ)/2));
				    cond4 = ( flag && (std::abs(midZ - upperZ)  < 10e-10));
				    cond5 = (!flag && (std::abs(midZ - midZslh) < 10e-10));
				    if (cond1 || cond2 || cond3 || cond4 || cond5){
					    // bisection
					    midZtemp = (lowerZ + upperZ) / 2;
					    flag = TRUE;
				    } else {
					    flag = FALSE;
				    }
				}
                		counter = 1;
     
				// update
				midFtemp = estEqP11(midZtemp, A(j,j), lambda(j,j), BVc, eigvalC);
				midZslh  = midZ;
				midZ     = upperZ;
				if (lowerF * midFtemp < 0){ 
				    upperZ = midZtemp; 
				} else { 
				    lowerZ = midZtemp; 
				}				
				if (std::abs(lowerF) < std::abs(upperF)){ 
					lowerZslh = lowerZ;    upperZslh = upperZ;
					lowerZ    = upperZslh; upperZ    = lowerZslh; 
					lowerFslh = lowerF;    upperFslh = upperF;
					lowerF    = upperFslh; upperF    = lowerFslh; 
				}
			} while (std::abs(lowerZ - upperZ) > 10e-10);

			// update diagonal element
			Phat(j,j) = (lowerZ + upperZ)/2;

			// update remainder of column/row
			eigvalC += (lowerZ + upperZ)/2;
			eigvalC /= (lowerZ + upperZ)/2;
			Phat.submat(indices, j2uvec) = - L12tilde % (eigvecC * (BVc / eigvalC));
			Phat.submat(j2uvec, indices) = (Phat.submat(indices, j2uvec)).t();
		}
		// stop if convergence reached
		if ((arma::abs(Phat - Pprev)).max() < minSuccDiff){ stayInLoop = FALSE; }
	}
	 
	// return results
	return Phat;
}


// [[Rcpp::export(".armaKcvlossR")]]
double kcvlossR(arma::mat&    lambda, 
                arma::mat&    Y, 
                arma::mat&    target, 
                Rcpp::List&   folds,
                const int&    nInit, 
                const double& minSuccDiff){
	/* ---------------------------------------------------------------------
	## k-fold cross-validated negative (!) log-likelihood score for a 
	## specified ridge penalty matrix
	## ---------------------------------------------------------------------
	## Arguments 
	## -> lambda      : penalty matrix
	## -> Y           : data matrix with rows and columns containing samples and 
	##                  variables, respectively
	## -> target      : target matrix for the precision 
	## -> folds       : cross-validation sample splits
	## -> nInit	      : maximum number of iterations. Passed on to the 
	##                  ridgePgen-function.
	## -> minSuccDiff : minimum successive difference (in the penalized 
	##                  loglikehood) to be achieved. Passed on to the 
	##                  ridgePgen-function.
	## ---------------------------------------------------------------------
	## Value
	## cvloss	  : negative (!) cross-validated loglikelihood
	## ---------------------------------------------------------------------
	## Authors        : Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	// declare variables
	arma::uvec leftin;         
	arma::uvec leftout;
	arma::mat  Phat;
	double     valPdet;
	double     signPdet;
	int        slh; 

	// extract number samples and folds
  	int n      = Y.n_rows;
	int nFolds = folds.size();

	// initiate the loss
	double cvLoss = 0;
    
	for (int f = 0; f < nFolds; ++f){
		// sort out the fold and its complement
		leftout = Rcpp::as<arma::uvec>(folds(f)) - 1;
		leftin  = arma::uvec(n-leftout.n_elem);
		slh     = 0;
		for (int i = 0; i < n; ++i){
			if (all(leftout != i)){
				leftin(slh) = i;
				slh         = slh + 1;
			}
		}

		// calculate loss
        	Phat = ridgePgen(arma::trans(Y.rows(leftin)) * Y.rows(leftin) / leftin.size(), 
        	                 lambda, target, nInit, minSuccDiff);
		log_det(valPdet, signPdet, Phat); 
		cvLoss = cvLoss + valPdet - 
                	          arma::accu((arma::trans(Y.rows(leftout)) * 
	                                      Y.rows(leftout)) % Phat) / leftout.size();
	}
	return -cvLoss / folds.size();
}


// [[Rcpp::export(".armaPenaltyPgen_banded")]]
arma::mat penaltyPgen_banded(double&        lambda, 
			     int            p,
                             arma::uvec&    zerosR, 
                             arma::uvec&    zerosC, 
                             bool           penalize_diag){  
	/* ---------------------------------------------------------------------
	## Construction of a penalty matrix with groups of variates sharing 
	## a common penalty parameter, possibly with unpenalized diagonal elements, 
	## and infinitely (i.e., very large) parameter values for zero elements.
	## ---------------------------------------------------------------------
	## Arguments 
	## -> lambda         : vector with penalty parameter values, one per group.
	##                     values should be specified in the same order as the 
	##                     first appearance of a group representative.  
	## -> groups         : vector indicating to which group a variate belongs.
	## -> zerosR         : row-index of zero precision elements.
	## -> zerosC         : column-index of zero precision elements.
	## -> penalized_diag : logical indicating whether the diagonal should 
	##                     be penalized
	## ---------------------------------------------------------------------
	## Value
	## cvloss	     : matrix with precision element-wise penalties
	## ---------------------------------------------------------------------
	## Authors           : Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	// declare variables
	arma::mat lambdaMat = arma::zeros(p, p);

	// assign each variate with a penalty value
	for (int j1 = 0; j1 < p; ++j1){
		for (int j2 = j1; j2 < p; ++j2){
			// lambdaMat(j2, j1) = std::pow(1 + lambda, -(p-(j2-j1)));
			lambdaMat(j2, j1) = lambda * (1+std::abs(j2-j1));
			lambdaMat(j1, j2) = lambdaMat(j2, j1);
		}
	}
	
	// set diagonal penalties to (virtually) zero if precision 
	// diagonal is not to be penalized
	if (!penalize_diag){
		lambdaMat.diag().zeros();
		lambdaMat.diag() += 10e-10;
	}
    
	// set penalty of known precision zeros to infinity
	if (zerosR.n_elem > 0){
		for (unsigned int z = 0; z < zerosR.n_elem; ++z){
			lambdaMat(zerosR[z], zerosC[z]) = 10e10;
		} 
	}

	// return penalty matrix
	return lambdaMat; 
}


// [[Rcpp::export(".armaPenaltyPgen_groups")]]
arma::mat penaltyPgen_groups(arma::vec&  lambda, 
                             arma::vec&  groups, 
                             arma::uvec& zerosR, 
                             arma::uvec& zerosC, 
                             bool        penalize_diag){  
	/* ---------------------------------------------------------------------
	## Construction of a penalty matrix with groups of variates sharing 
	## a common penalty parameter, possibly with unpenalized diagonal elements, 
	## and infinitely (i.e., very large) parameter values for zero elements.
	## ---------------------------------------------------------------------
	## Arguments 
	## -> lambda         : vector with penalty parameter values, one per group.
	##                     values should be specified in the same order as the 
	##                     first appearance of a group representative.  
	## -> groups         : vector indicating to which group a variate belongs.
	## -> zerosR         : row-index of zero precision elements.
	## -> zerosC         : column-index of zero precision elements.
	## -> penalized_diag : logical indicating whether the diagonal should 
	##                     be penalized
	## ---------------------------------------------------------------------
	## Value
	## cvloss	     : matrix with precision element-wise penalties
	## ---------------------------------------------------------------------
	## Authors           : Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	// declare variables
	int p = groups.n_elem;
	arma::vec lambdaVec = arma::zeros(p); 
	arma::mat lambdaMat = arma::zeros(p, p);
	arma::uvec groupIds;
    
	// assign each variate with a penalty value
	arma::uword counter = 0; 
	for (unsigned int g = 0; g < groups.n_elem; ++g){
		groupIds = arma::find(groups == groups(g));
		if (all(groupIds >= g) && (counter < lambda.n_elem)){
			lambdaVec.rows(groupIds) += lambda(counter);
			counter = counter + 1;
		}
	}
	
	// fill penalty matrix
	lambdaMat.each_col() += lambdaVec / 2;
	lambdaMat.each_row() += arma::trans(lambdaVec) / 2;

	// set diagonal penalties to (virtually) zero if precision 
	// diagonal is not to be penalized
	if (!penalize_diag){
		lambdaMat.diag().zeros();
		lambdaMat.diag() += 10e-10;
	}
    
	// set penalty of known precision zeros to infinity
	if (zerosR.n_elem > 0){
		for (unsigned int z = 0; z < zerosR.n_elem; ++z){
			lambdaMat(zerosR[z], zerosC[z]) = 10e10;
		} 
	}

	// return penalty matrix
	return lambdaMat; 
}


// [[Rcpp::export(".armaKCVlossR_groups")]]
double kcvlossR_groups(arma::vec&   lambdaGrps, 
                      arma::mat&    Y, 
                      arma::mat&    target, 
                      Rcpp::List&   folds, 
                      arma::vec&    groups,
                      arma::uvec&   zerosR, 
                      arma::uvec&   zerosC, 
                      bool          penalize_diag,    
                      const int&    nInit, 
                      const double& minSuccDiff){

	/* ---------------------------------------------------------------------
	## K-fold cross-validated negative (!) loglikelihood for a specified 
	## penalty matrix, which assumes that variates are grouped and penalized 
	## group-wise. Effectively, this is a wrapper around the 
	## penaltyPgen_groups and kcvloss-functions.
	## ---------------------------------------------------------------------
	## Arguments 
	## -> lambdaGrps     : vector with penalty parameter values, one per group.
	##                     values should be specified in the same order as the 
	##                     first appearance of a group representative.  
	## -> Y              : data matrix with rows and columns containing samples
	##                     and variables, respectively
	## -> target         : target matrix for the precision 
	## -> folds          : cross-validation sample splits
	## -> groups         : vector indicating to which group a variate belongs.
	## -> zerosR         : row-index of zero precision elements.
	## -> zerosC         : column-index of zero precision elements.
	## -> penalized_diag : logical indicating whether the diagonal should 
	##                     be penalized    
	## -> nInit	         : maximum number of iterations. Passed on to the 
	##                     ridgePgen-function.
	## -> minSuccDiff    : minimum successive difference (in the penalized 
	##                     loglikehood) to be achieved. Passed on to the 
	##                     ridgePgen-function.
	## ---------------------------------------------------------------------
	## Value
	## cvloss	      : negative (!) cross-validated loglikelihood
	## ---------------------------------------------------------------------
	## Authors        : Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	// build penalty matrix
	arma::mat lambda = penaltyPgen_groups(lambdaGrps, 
	                                      groups, 
                                              zerosR,
                                              zerosC,
                                              penalize_diag);
                                      
    return kcvlossR(lambda, Y, target, folds, nInit, minSuccDiff);
}


// [[Rcpp::export(".armaKCVlossR_banded")]]
double kcvlossR_banded(double&   lambda, 
                      arma::mat&    Y, 
                      arma::mat&    target, 
                      Rcpp::List&   folds, 
                      arma::uvec&   zerosR, 
                      arma::uvec&   zerosC, 
                      bool          penalize_diag,    
                      const int&    nInit, 
                      const double& minSuccDiff){
	/* ---------------------------------------------------------------------
	## K-fold cross-validated negative (!) loglikelihood for a specified 
	## penalty matrix, which assumes that variates are grouped and penalized 
	## group-wise. Effectively, this is a wrapper around the 
	## penaltyPgen_groups and kcvloss-functions.
	## ---------------------------------------------------------------------
	## Arguments 
	## -> lambdaGrps     : vector with penalty parameter values, one per group.
	##                     values should be specified in the same order as the 
	##                     first appearance of a group representative.  
	## -> Y              : data matrix with rows and columns containing samples
	##                     and variables, respectively
	## -> target         : target matrix for the precision 
	## -> folds          : cross-validation sample splits
	## -> groups         : vector indicating to which group a variate belongs.
	## -> zerosR         : row-index of zero precision elements.
	## -> zerosC         : column-index of zero precision elements.
	## -> penalized_diag : logical indicating whether the diagonal should 
	##                     be penalized    
	## -> nInit	         : maximum number of iterations. Passed on to the 
	##                     ridgePgen-function.
	## -> minSuccDiff    : minimum successive difference (in the penalized 
	##                     loglikehood) to be achieved. Passed on to the 
	##                     ridgePgen-function.
	## ---------------------------------------------------------------------
	## Value
	## cvloss	     : negative (!) cross-validated loglikelihood
	## ---------------------------------------------------------------------
	## Authors           : Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	// build penalty matrix
	arma::mat lambdaMat = penaltyPgen_banded(lambda, 
                                          	 Y.n_cols, 
                                          	 zerosR,
                                          	 zerosC,
                                          	 penalize_diag);
                                      
	return kcvlossR(lambdaMat, Y, target, folds, nInit, minSuccDiff);
}

