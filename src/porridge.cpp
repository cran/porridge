// Rcpp::Rcout << "wessel is gek 1" << std::endl;					

#include <RcppArmadillo.h>
#include <cmath>

// These are not needed:
// using namespace std;
// using namespace Rcpp;
// using namespace RcppArmadillo;
// [[Rcpp::depends("Rcpp")]]
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::interfaces(r, cpp)]]


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





///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


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



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


inline arma::mat armaRidgeSanyTarget(const arma::mat & S,
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
	arma::vec D_inv = (sqroot + eigvals); 
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
		eigvecs.each_row() %= arma::trans(lambda/arma::sqrt((sqroot - eigvals)));
	}
	return eigvecs * arma::trans(eigvecs);
}

inline arma::mat armaRidgeSscalarTarget(const arma::mat& S,
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
	arma::vec D_inv = (sqroot + eigvals); 
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
 		eigvecs.each_row() %= arma::trans(arma::sqrt(lambda/(sqroot - eigvals)));
	}
	return eigvecs * arma::trans(eigvecs);
}



inline arma::mat armaRidgeS(const arma::mat& S,
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
		return armaRidgeSscalarTarget(S,  target(0,0), lambda, invert);
  	} else {
    		return armaRidgeSanyTarget(S, target, lambda, invert);
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


inline double ridgePrepCompLL(const arma::mat&  Y, 
                              const arma::ivec& ids, 
                              arma::mat&        Sz,
                              arma::mat&        Se){

	/*---------------------------------------------------------------------------------------------------
	## Loglikehood of a mixture of GGMs
	## --------------------------------------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## repPhat	: parameter estimates of the mixture of GGMs.
	##                (as provided by the rigdeGGMmixture-function)
	## --------------------------------------------------------------------------------------------------
	## Value
	## LL		: negative loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## --------------------------------------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------------------------------------*/

	// extract number of samples and variates  
  	int        nk = Y.n_rows;
  	arma::ivec is = arma::unique(ids);		
	int        n  = is.n_elem;
	int        p  = Y.n_cols;

	// declare variables used in penalized EM algorithm 
	double       compLL;
	double       valSz;
	double       signSz;
	double       valSe;
	double       signSe;
	unsigned int Ki;
	arma::uvec   kis;
	arma::mat    Yi;
	arma::uvec   i2uvec;
	arma::mat    Zis = arma::mat(n, p);
	arma::mat    sampleSz = arma::mat(p, p);
	arma::mat    sampleSe = arma::mat(p, p);

	// evaluate the log-likelihood
	sampleSe = arma::zeros(p, p);
	for (int i = 0; i < n; ++i){
		// obtain number of and row ids of repetitions
		kis    = arma::find(ids == is(i));
		Ki     = kis.n_elem; 

		// E-step: conditional signal
		if (Ki > 0){
			// calculate expectation of the signal, given data and parameters
			i2uvec            = arma::find(is == is(i));
			Sz               *= Ki;
			Zis.rows(i2uvec)  = arma::sum(Y.rows(kis), 0) * arma::inv_sympd(Sz + Se) * Sz;
			Zis.rows(i2uvec) /= Ki;
			Sz               /= Ki;

			// calculate sample's contribution to sample error covariance matrix
			Yi             = Y.rows(kis);
			Yi.each_row() += -Zis.rows(i2uvec);
			sampleSe       = sampleSe + Yi.t() * Yi;
		}
	}
	sampleSz  = arma::trans(Zis) * Zis;

	// calculate complete loglikelihood
	log_det(valSz, signSz, Sz);
	log_det(valSe, signSe, Se);
	compLL = n * valSz + nk * valSe + arma::accu(arma::inv_sympd(Sz) % sampleSz) + arma::accu(arma::inv_sympd(Se) % sampleSe);

	// return the loglikelihood
	return -compLL;
}



inline double ridgePrepQloss(const arma::mat&  Y, 
                              const arma::ivec& ids, 
                              arma::mat&        Sz,
                              arma::mat&        Se){

	/*---------------------------------------------------------------------------------------------------
	## Loglikehood of a mixture of GGMs
	## --------------------------------------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## repPhat	: parameter estimates of the mixture of GGMs.
	##                (as provided by the rigdeGGMmixture-function)
	## --------------------------------------------------------------------------------------------------
	## Value
	## LL		: negative loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## --------------------------------------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------------------------------------*/

	// extract number of samples and variates  
  	int        nk = Y.n_rows;
  	arma::ivec is = arma::unique(ids);		
	int        n  = is.n_elem;
	int        p  = Y.n_cols;

	// declare variables used in penalized EM algorithm 
	double       quadLoss;
	unsigned int Ki;
	arma::uvec   kis;
	arma::mat    Yi;
	arma::uvec   i2uvec;
	arma::mat    Zis = arma::mat(n, p);
	arma::mat    sampleSz = arma::mat(p, p);
	arma::mat    sampleSe = arma::mat(p, p);

	// evaluate the log-likelihood
	sampleSe = arma::zeros(p, p);
	for (int i = 0; i < n; ++i){
		// obtain number of and row ids of repetitions
		kis    = arma::find(ids == is(i));
		Ki     = kis.n_elem; 

		// E-step: conditional signal
		if (Ki > 0){
			// calculate expectation of the signal, given data and parameters
			i2uvec            = arma::find(is == is(i));
			Sz               *= Ki;
			Zis.rows(i2uvec)  = arma::sum(Y.rows(kis), 0) * arma::inv_sympd(Sz + Se) * Sz;
			Zis.rows(i2uvec) /= Ki;
			Sz               /= Ki;

			// calculate sample's contribution to sample error covariance matrix
			Yi             = Y.rows(kis);
			Yi.each_row() += -Zis.rows(i2uvec);
			sampleSe       = sampleSe + Yi.t() * Yi;
		}
	}
	sampleSz  = arma::trans(Zis) * Zis;

	// calculate the quadratic loss
	sampleSz /= n;
	sampleSe /= nk;
 	quadLoss  = arma::accu(arma::square(arma::inv_sympd(Sz) * sampleSz - arma::eye(p,p))) + arma::accu(arma::square(arma::inv_sympd(Se) * sampleSe - arma::eye(p,p)));

	// return the negative quadratic loss (negative for the purpose of maximization in R)
	return -quadLoss;
}



inline double ridgePrepQloss2(const arma::mat&  Y, 
                              const arma::ivec& ids, 
                              arma::mat&        Sz,
                              arma::mat&        Se){

	/*---------------------------------------------------------------------------------------------------
	## Loglikehood of a mixture of GGMs
	## --------------------------------------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## repPhat	: parameter estimates of the mixture of GGMs.
	##                (as provided by the rigdeGGMmixture-function)
	## --------------------------------------------------------------------------------------------------
	## Value
	## LL		: negative loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## --------------------------------------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------------------------------------*/

	// extract number of samples and variates  
  	int        nk = Y.n_rows;
  	arma::ivec is = arma::unique(ids);		
	int        n  = is.n_elem;
	int        p  = Y.n_cols;

	// declare variables used in penalized EM algorithm 
	double       quadLoss;
	arma::uvec   kis;
	arma::mat    Yi;
	arma::uvec   i2uvec;
	arma::mat    sampleSz = arma::mat(p, p);
	arma::mat    sampleSe = arma::mat(p, p);
	arma::mat    Yav      = arma::mat(n, p);
	double       kweight  = 0;

	// initiate
	for (unsigned int i = 0; i < std::abs(n); ++i){
		kis        = arma::find(ids == is(i));
		kweight    = kweight + 1/kis.n_elem;
		Yav.row(i) = arma::mean(Y.rows(kis));
	}
	sampleSz  = arma::trans(Yav) * Yav;
	sampleSz /= n;
	sampleSe  = arma::trans(Y) * Y;
	sampleSe /= nk;
	sampleSe  = sampleSe - sampleSz;
	sampleSe /= (1 - kweight/n);
	sampleSz  = sampleSz - sampleSe * (kweight/n);

	// calculate the quadratic loss
	sampleSz /= n;
	sampleSe /= nk;
 	quadLoss  = arma::accu(arma::square(arma::inv_sympd(Sz) * sampleSz - arma::eye(p,p))) + arma::accu(arma::square(arma::inv_sympd(Se) * sampleSe - arma::eye(p,p)));

	// return the negative quadratic loss (negative for the purpose of maximization in R)
	return -quadLoss;
}



inline double ridgePrepLL(const arma::mat&  Y, 
                          const arma::ivec& ids, 
                          arma::mat&        Sz,
                          arma::mat&        Se){

	/*---------------------------------------------------------------------------------------------------
	## Loglikehood of a mixture of GGMs
	## --------------------------------------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size and p is dimension.
	## repPhat	: parameter estimates of the mixture of GGMs.
	##                (as provided by the rigdeGGMmixture-function)
	## --------------------------------------------------------------------------------------------------
	## Value
	## LL		: negative loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## --------------------------------------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------------------------------------*/

	// extract number of samples and variates  
  	arma::ivec is = arma::unique(ids);		
	int        n  = is.n_elem;
	int        p  = Y.n_cols;

	// declare variables used in penalized EM algorithm 
	double       LL;
	arma::mat    Ytilde;
	double       valDetSe;
	double       valDetSze;
	unsigned int Ki;
	arma::uvec   kis;
	arma::mat    Yi;
	arma::uvec   i2uvec;
	arma::mat    Syi;
	arma::mat    Sei;
	arma::vec    eigvalSze;
	arma::mat    eigvecSze;
	arma::mat    Seinvsqrt;

	// some linear algebra for efficient evaluation of the loglikelihood 
	// when the number of repetitions varies among samples
	// eigendecomposition of error covariance matrix
  	arma::eig_sym(eigvalSze, eigvecSze, Se, "dc");

	// determinant of error covariance matrix
	valDetSe = arma::accu(arma::log(eigvalSze));

	// evaluate matrix square root of inverse of error covariance matrix
	Seinvsqrt = eigvecSze * arma::diagmat(arma::pow(arma::sqrt(eigvalSze), -1)) * arma::trans(eigvecSze);

	// obtain eigen values of Se^(-1/2) Sz Se^(-1/2)
  	arma::eig_sym(eigvalSze, eigvecSze, Seinvsqrt * Sz * Seinvsqrt, "dc");

	// matrix multiplication that would otherwise be done repetitive inside the loop
	eigvecSze = Seinvsqrt * eigvecSze;

	// evaluate the log-likelihood
	LL = 0;
	for (int i = 0; i < n; ++i){
		// obtain number of and row ids of repetitions
		kis    = arma::find(ids == is(i));
		Ki     = kis.n_elem; 
		i2uvec = arma::find(is == is(i));

		// calculate the sample's loglikelihood contribution
		Yi     = Y.rows(kis);
		if (Ki > 1){
			// calculate auxillary matrices
			Syi = arma::zeros(p, p);
			Sei = arma::zeros(p, p);
			for (unsigned ki1 = 0; ki1 < Ki-1; ++ki1){
				for (unsigned ki2 = ki1+1; ki2 < Ki; ++ki2){
					Syi = Syi + arma::trans(Yi.row(ki1)) * Yi.row(ki2);
				}
			}
			Sei = arma::trans(Yi) * Yi;
			Syi = Syi + arma::trans(Syi) + Sei;

			// actual loglikelihood evaluation
			eigvalSze *= Ki;
			eigvalSze += 1;
			valDetSze  = arma::accu(arma::log(eigvalSze));
			LL         = LL - Ki * valDetSe - valDetSze -
                                          arma::accu((Seinvsqrt * Seinvsqrt) % Sei) +
      	                                  arma::accu((eigvecSze * arma::diagmat(1-arma::pow(eigvalSze, -1)) *
                                                         arma::trans(eigvecSze)) % Syi) / Ki;
			eigvalSze += -1;
			eigvalSze /= Ki;
		} 
		if (Ki == 1){
			// loglikelihood evaluation for samples without repetitions
			eigvalSze += 1;
			valDetSze  = arma::accu(arma::log(eigvalSze));
			LL         = LL - valDetSze - arma::accu(Yi * eigvecSze * arma::diagmat(arma::pow(eigvalSze, -1)) *
                                                            arma::trans(eigvecSze) * arma::trans(Yi));
			eigvalSze += -1;
		}
	}

	// return the loglikelihood
	return LL;
}


// [[Rcpp::export(".armaRidgePrepEM")]]
Rcpp::List ridgePrepEM(arma::mat         Y, 
                       arma::ivec        ids, 
                       const double      lambdaZ,
                       const double      lambdaE,
                       const arma::mat&  targetZ,
                       const arma::mat&  targetE,  
                       const int&        nInit, 
                       const double&     minSuccDiff){

	/* ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size (including the 
	##                repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## lambdaZ      : ridge penalty parameter for the signal precision matrix
        ## lambdaE      : ridge penalty parameter for the error precision matrix
        ## targetZ      : target matrix towards the signal precision matrix 
	##                is shrunken
        ## targetE      : target matrix towards the error precision matrix 
	##                is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized 
	##                loglikehood) to be achieved
	## ---------------------------------------------------------------------
	## Value
	## A list-object with the following slots:
	## Pz		: estimated signal precision matrices
	## Pe		: estimated error precision matrices  
	## penLL	: penalized loglikelihood
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------- */


	//----------- variable declaration ----------------------------------//
	// extract number of samples and variates  
  	int        nk = Y.n_rows;
  	arma::ivec is = arma::unique(ids);		
	int        n  = is.n_elem;
	int        p  = Y.n_cols;

	// declare variables used in penalized EM algorithm 
	double       penLL = -10e+10;
	double       penLLprev;
	arma::mat    Ytilde;
	arma::mat    Zis = arma::mat(n, p);
	unsigned int Ki;
	arma::uvec   kis;
	arma::mat    Yi;
	arma::uvec   i2uvec;
	arma::mat    Sz;
	arma::mat    Se;
	arma::mat    Yav = arma::mat(n,p);
	double       kweight = 0;
	arma::ivec   Ks = arma::ivec(n);
	arma::mat    VarZ  = arma::zeros(p,p);
	arma::mat    VarZk = arma::zeros(p,p);
	arma::mat    VarZslh;
	arma::mat    invSe;
	arma::mat    invSz;
	arma::uvec   isameK;
	arma::mat    Sze;
	//----------- end of variable declaration ---------------------------//



	//----------- estimator initialization ------------------------------//
	// use moment estimators followed by penalization
	// first, create nxp matrix with sample-wise averaged observations
	for (unsigned int i = 0; i < std::abs(n); ++i){
		kis        = arma::find(ids == is(i));
		Ks(i)      = kis.n_elem;
		kweight    = kweight + 1/Ks(i);
		Yav.row(i) = arma::mean(Y.rows(kis));
	}

	// obtain moment estimators from moment decomposition
	Sz  = arma::trans(Yav) * Yav;
	Sz /= n;
	Se  = arma::trans(Y) * Y;
	Se /= nk;
	Se  = Se - Sz;
	Se /= (1 - kweight/n);
	Sz  = Sz - Se * (kweight/n);

	// obtain initial estimates from moment estimators
	Sz  = armaRidgeS(Sz, targetZ, lambdaZ/n);
	Se  = armaRidgeS(Se, targetE, lambdaE/nk);
	//----------- end of initialization ---------------------------------//


	//----------- reformat data -----------------------------------------//
	// sort data by number of replications (for faster E-step calculation)
	arma::ivec uniqKs = arma::unique(Ks);
	if (Ks.min() != Ks.max()){
		if (uniqKs.n_elem > 1){
			arma::uvec shuffle = arma::uvec(nk);
			arma::uvec id2Ks;
			int counter = 0;
			for (unsigned int K = 0; K < uniqKs.n_elem; ++K){
				id2Ks = arma::find(Ks == uniqKs(K));
				for (unsigned int i = 0; i < id2Ks.n_elem; ++i){
					kis = arma::find(ids == is(id2Ks(i)));
					shuffle.subvec(counter, counter+kis.n_elem-1) = kis;
					counter = counter + kis.n_elem;
				}
			}

			// now reshuffle for ascending number of replications order
			Y   = Y.rows(shuffle);
			ids = ids.elem(shuffle);
		}
	}
	//----------- end of reformating ------------------------------------//



	//----------- EM algorithm ------------------------------------------//
	// start iterating	
	bool stayInLoop = TRUE;
	for (int u = 0; u < nInit && stayInLoop; ++u){
		// store previously achieved penalized loglikelihood
		penLLprev = penLL;


		//----------- E-step ----------------------------------------//
		// obtain the two sufficient statistics
		// the balanced and unbalanced (replication-wise) cases are dealt 
		// with seperately, as the former allows to speed up computations
		
		// inverse of current signal and error covariance matrices
		// need in both steps
		invSz    = arma::inv_sympd(Sz);
		invSe    = arma::inv_sympd(Se);

		// unbalanced (replication-wise) design
		if (Ks.min() != Ks.max()){
			VarZ   = arma::zeros(p,p);
			VarZk  = arma::zeros(p,p);

			for (unsigned int r = 0; r < uniqKs.n_elem; ++r){	
				isameK = arma::find(Ks == uniqKs(r));
				Ki     = uniqKs(r);

				// calculate expectation of the signal, given data and parameters
				Sz  *= Ki;
				Sze  = arma::inv_sympd(Sz + Se) * Sz;

				for (unsigned int i = 0; i < isameK.n_elem; ++i){
					// obtain number of and row ids of repetitions
					kis               = arma::find(ids == is(isameK(i)));
					i2uvec            = arma::find(is  == is(isameK(i)));

					// calculate expectation of the signal, given data and parameters
					Zis.rows(i2uvec)  = arma::sum(Y.rows(kis), 0) * Sze;
					Zis.rows(i2uvec) /= Ki;
				}
				Sz /= Ki;

				// calculate variance of the signal, given data and parameters
				invSe        *= Ki;
				VarZslh       = arma::inv_sympd(invSz + invSe);
				VarZslh      *= isameK.n_elem;
				VarZ          = VarZ  + VarZslh;
				VarZslh      *= Ki;
				VarZk         = VarZk + VarZslh;
				invSe        /= Ki;
				VarZslh      /= Ki * isameK.n_elem;
			}

			// average the variances
			VarZ  /= n;
			VarZk /= nk;
		}

		// balanced (replication-wise) design
		if (Ks.min() == Ks.max()){
			// calculation of sufficient statistics
			Sz *= Ks.min();
			Sz  = arma::inv_sympd(Sz + Se) * Sz;
			for (int i = 0; i < n; ++i){
				// obtain number of and row ids of repetitions
				kis               = arma::find(ids == is(i));
				Ki                = kis.n_elem; 

				// calculate expectation of the signal, given data and parameters
				i2uvec            = arma::find(is == is(i));
				Zis.rows(i2uvec)  = arma::sum(Y.rows(kis), 0) * Sz;
				Zis.rows(i2uvec) /= Ki;
			}

			// calculate the variance of the signal, given data and parameters
			invSe        *= Ks.min();
			VarZk         = arma::inv_sympd(invSz + invSe);
			VarZ          = VarZk;
		}
		//----------- end of E-step ---------------------------------//



		//----------- M-step ----------------------------------------//
		// M-step: estimation of signal precision matrix
		Sz  = arma::trans(Zis) * Zis;
		Sz /= n;
		Sz = Sz + VarZ;
		Sz  = armaRidgeS(Sz, targetZ, lambdaZ/n);

		// M-step: estimation of error precision matrix		
		Se = arma::zeros(p, p);
		for (int i = 0; i < n; ++i){
			// obtain number of and row ids of repetitions
			kis = arma::find(ids == is(i));
			Ki  = kis.n_elem; 
			if (Ki > 0){
				// calculate sample's contribution to sample error covariance matrix
				i2uvec         = arma::find(is == is(i));
				Yi             = Y.rows(kis);
				Yi.each_row() += -Zis.rows(i2uvec);
				Se             = Se + Yi.t() * Yi;
			}
		}
		Se /= nk;
		Se = Se + VarZk;
		Se  = armaRidgeS(Se, targetE, lambdaE/nk);
		//----------- end of M-step ---------------------------------//



		//----------- convergence assessment ------------------------//
		// evaluate the loglikelihood
		penLL = ridgePrepLL(Y, ids, Sz, Se);

		// add the penalty's contribution
		penLL = penLL - 0.5 * lambdaZ * arma::accu(arma::square(arma::inv_sympd(Sz) - targetZ)) - 
		                0.5 * lambdaE * arma::accu(arma::square(arma::inv_sympd(Se) - targetE));			

		// assess convergence		
	    	if (std::abs((penLLprev - penLL)/penLLprev) < minSuccDiff){ 
			stayInLoop = FALSE; 
		}
		//----------- end of convergence assessment -----------------//
	}
	//----------- end of EM algorithm -----------------------------------//


	// return stuff as a list-object
	return Rcpp::List::create(Rcpp::Named("Pz")    = arma::inv_sympd(Sz), 
                                  Rcpp::Named("Pe")    = arma::inv_sympd(Se), 
                                  Rcpp::Named("penLL") = penLL);
}


inline double ridgePrepEMinternal(arma::mat         Ytrain,
				  const arma::mat&  Ytest,  
	                          arma::ivec        idsTrain,
	                          const arma::ivec& idsTest,  
	                          const double      lambdaZ,
        	                  const double      lambdaE,
        	                  const arma::mat&  targetZ,
        	                  const arma::mat&  targetE,  
        	                  const int&        nInit, 
        	                  const double&     minSuccDiff,
				  std::string       CVcrit){

	/* ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size (including the 
	##                repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## lambdaZ      : ridge penalty parameter for the signal precision matrix
        ## lambdaE      : ridge penalty parameter for the error precision matrix
        ## targetZ      : target matrix towards the signal precision matrix 
	##                is shrunken
        ## targetE      : target matrix towards the error precision matrix 
	##                is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized 
	##                loglikehood) to be achieved
	## ---------------------------------------------------------------------
	## Value
	## The cross-validated loss for a single fold
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------- */

	//----------- variable declaration ----------------------------------//
	// extract number of samples and variates  
  	int          nk = Ytrain.n_rows;
  	arma::ivec   is = arma::unique(idsTrain);		
	unsigned int n  = is.n_elem;
	unsigned int p  = Ytrain.n_cols;

	// declare variables used in penalized EM algorithm 
	double       penLL = -10e+10;
	double       penLLprev;
	arma::mat    Ytilde;
	arma::mat    Zis = arma::mat(n, p);
	unsigned int Ki;
	arma::uvec   kis;
	arma::mat    Yi;
	arma::uvec   i2uvec;
	arma::mat    Sz;
	arma::mat    Se;
	arma::mat    Yav = arma::mat(n,p);
	double       kweight = 0;
	arma::ivec   Ks = arma::ivec(n);
	arma::mat    VarZ  = arma::zeros(p,p);
	arma::mat    VarZk = arma::zeros(p,p);
	arma::mat    VarZslh;
	arma::mat    invSe;
	arma::mat    invSz;
	arma::uvec   isameK;
	arma::mat    Sze;
	double       CVloss = 10e-10;
	//----------- end of variable declaration ---------------------------//



	//----------- estimator initialization ------------------------------//
	// use moment estimators followed by penalization
	// first, create nxp matrix with sample-wise averaged observations
	for (unsigned int i = 0; i < n; ++i){
		kis        = arma::find(idsTrain == is(i));
		Ks(i)      = kis.n_elem;
		kweight    = kweight + 1/Ks(i);
		Yav.row(i) = arma::mean(Ytrain.rows(kis));
	}

	// obtain moment estimators from moment decomposition
	Sz  = arma::trans(Yav) * Yav;
	Sz /= n;
	Se  = arma::trans(Ytrain) * Ytrain;
	Se /= nk;
	Se  = Se - Sz;
	Se /= (1 - kweight/n);
	Sz  = Sz - Se * (kweight/n);

	// obtain initial estimates from moment estimators
	Sz  = armaRidgeS(Sz, targetZ, lambdaZ);
	Se  = armaRidgeS(Se, targetE, lambdaE);
	//----------- end of initialization ---------------------------------//




	//----------- reformat data -----------------------------------------//
	// sort data by number of replications (for faster E-step calculation)
	arma::ivec uniqKs = arma::unique(Ks);
	if (uniqKs.n_elem > 1){
		arma::uvec shuffle = arma::uvec(nk);
		arma::uvec id2Ks;
		int counter = 0;
		for (unsigned int K = 0; K < uniqKs.n_elem; ++K){
			id2Ks = arma::find(Ks == uniqKs(K));
			for (unsigned int i = 0; i < id2Ks.n_elem; ++i){
				kis = arma::find(idsTrain == is(id2Ks(i)));
				shuffle.subvec(counter, counter+kis.n_elem-1) = kis;
				counter = counter + kis.n_elem;
			}
		}

		// now reshuffle for ascending number of replications order
		Ytrain   = Ytrain.rows(shuffle);
		idsTrain = idsTrain.elem(shuffle);
	}
	//----------- end of reformating ------------------------------------//




	//----------- EM algorithm ------------------------------------------//
	// start iterating	
	bool stayInLoop = TRUE;
	for (int u = 0; u < nInit && stayInLoop; ++u){
		// store previously achieved penalized loglikelihood
		penLLprev = penLL;


		//----------- E-step ----------------------------------------//
		// obtain the two sufficient statistics
		// the balanced and unbalanced (replication-wise) cases are dealt 
		// with seperately, as the former allows to speed up computations
		
		// inverse of current signal and error covariance matrices
		// need in both steps
		invSz    = arma::inv_sympd(Sz);
		invSe    = arma::inv_sympd(Se);

		// unbalanced (replication-wise) design
		if (Ks.min() != Ks.max()){
			VarZ   = arma::zeros(p,p);
			VarZk  = arma::zeros(p,p);

			for (unsigned int r = 0; r < uniqKs.n_elem; ++r){	
				isameK = arma::find(Ks == uniqKs(r));
				Ki     = uniqKs(r);

				// calculate expectation of the signal, given data and parameters
				Sz  *= Ki;
				Sze  = arma::inv_sympd(Sz + Se) * Sz;

				for (unsigned int i = 0; i < isameK.n_elem; ++i){
					// obtain number of and row ids of repetitions
					kis               = arma::find(idsTrain == is(isameK(i)));
					i2uvec            = arma::find(is       == is(isameK(i)));

					// calculate expectation of the signal, given data and parameters
					Zis.rows(i2uvec)  = arma::sum(Ytrain.rows(kis), 0) * Sze;
					Zis.rows(i2uvec) /= Ki;
				}
				Sz /= Ki;

				// calculate variance of the signal, given data and parameters
				invSe        *= Ki;
				VarZslh       = arma::inv_sympd(invSz + invSe);
				VarZslh      *= isameK.n_elem;
				VarZ          = VarZ  + VarZslh;
				VarZslh      *= Ki;
				VarZk         = VarZk + VarZslh;
				invSe        /= Ki;
				VarZslh      /= Ki * isameK.n_elem;
			}

			// average the variances
			VarZ  /= n;
			VarZk /= nk;
		}

		// balanced (replication-wise) design
		if (Ks.min() == Ks.max()){
			// calculation of sufficient statistics
			Sz *= Ks.min();
			Sz  = arma::inv_sympd(Sz + Se) * Sz;
			for (unsigned int i = 0; i < n; ++i){
				// obtain number of and row ids of repetitions
				kis               = arma::find(idsTrain == is(i));
				Ki                = kis.n_elem; 

				// calculate expectation of the signal, given data and parameters
				i2uvec            = arma::find(is == is(i));
				Zis.rows(i2uvec)  = arma::sum(Ytrain.rows(kis), 0) * Sz;
				Zis.rows(i2uvec) /= Ki;
			}

			// calculate the variance of the signal, given data and parameters
			invSe        *= Ks.min();
			VarZk         = arma::inv_sympd(invSz + invSe);
			VarZ          = VarZk;
		}
		//----------- end of E-step ---------------------------------//

	

		//----------- M-step ----------------------------------------//
		// M-step: estimation of signal precision matrix
		Sz  = arma::trans(Zis) * Zis;
		Sz /= n;
		Sz = Sz + VarZ;
		Sz  = armaRidgeS(Sz, targetZ, lambdaZ/n);

		// M-step: estimation of error precision matrix		
		Se = arma::zeros(p, p);
		for (unsigned int i = 0; i < n; ++i){
			// obtain number of and row ids of repetitions
			kis = arma::find(idsTrain == is(i));
			Ki  = kis.n_elem; 
			if (Ki > 0){
				// calculate sample's contribution to sample error covariance matrix
				i2uvec         =  arma::find(is == is(i));
				Yi             =  Ytrain.rows(kis);
				Yi.each_row() += -Zis.rows(i2uvec);
				Se             =  Se + Yi.t() * Yi;
			}
		}
		Se /= nk;
		Se = Se + VarZk;
		Se  = armaRidgeS(Se, targetE, lambdaE/nk);
		//----------- end of M-step ---------------------------------//



		//----------- convergence assessment ------------------------//
		// add the penalty's loglikelihood contribution
		penLL = ridgePrepLL(Ytrain, idsTrain, Sz, Se);
		penLL = penLL - 0.5 * lambdaZ * arma::accu(arma::square(arma::inv_sympd(Sz) - targetZ)) -
                                0.5 * lambdaE * arma::accu(arma::square(arma::inv_sympd(Se) - targetE));

		// assess convergence		
	    	if (std::abs((penLLprev - penLL)/penLLprev) < minSuccDiff){ stayInLoop = FALSE; }
		//----------- end of convergence assessment -----------------//

	}
	//----------- end of EM algorithm -----------------------------------//


	// select cross-validation criterion
	if (CVcrit == "LL"){ 
		CVloss = ridgePrepLL(Ytest, idsTest, Sz, Se);
	}
	if (CVcrit == "compLL"){ 
		CVloss = ridgePrepCompLL(Ytest, idsTest, Sz, Se);
	}
	if (CVcrit == "Qloss"){ 
		CVloss = ridgePrepQloss(Ytest, idsTest, Sz, Se);
	}
	if (CVcrit == "Qloss2"){ 
		CVloss = ridgePrepQloss2(Ytest, idsTest, Sz, Se);
	}
	return CVloss;
}



// [[Rcpp::export(".armaRidgePrepKcvLL")]]
double ridgePrepKcvLL(const arma::vec   lambdaZE,
                      const arma::mat&  Y, 
                      const arma::ivec& ids, 
                      const arma::mat&  targetZ,
                      const arma::mat&  targetE,  
                      const int&        nInit, 
                      const double&     minSuccDiff,
                      Rcpp::List        folds,
		      std::string       CVcrit){

	/* --------------------------------------------------------------------------------------------------
	## Cross-validated loglikehood of a mixture of GGMs
	## --------------------------------------------------------------------------------------------------
	## Arguments 
        ## lambdaZE     : ridge penalty parameters for the signal and error precision matrix
	## Y		: n*p data matrix where n is sample size (including the repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## targetZ      : target matrix towards the signal precision matrix is shrunken
        ## targetE      : target matrix towards the error precision matrix is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized 
	##                loglikehood) to be achieved
	## folds        : cross-validation sample splits.
	## --------------------------------------------------------------------------------------------------
	## Value
	## cvLL		: negative cross-validated loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## --------------------------------------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	-------------------------------------------------------------------------------------------------- */

	// extract number samples and folds
  	int n      = Y.n_rows;
	int nFolds = folds.size();

	// declare variables used in the CV-loop
	double     cvLL = 0;
	arma::uvec fold, allButFold;
	Rcpp::List ridgePrepHat;
	int        slh;

	// start of cross-validation
	for (int f = 0; f < nFolds; ++f){
		// sort out the fold and its complement
		fold       = Rcpp::as<arma::uvec>(folds(f)) - 1;
		allButFold = arma::uvec(n-fold.n_elem);
		slh        = 0;
		for (int i = 0; i < n; ++i){
			if (all(fold != i)){
				allButFold(slh) = i;
				slh = slh + 1;
			}
		}

		// estimate model from all but left-out fold
		cvLL = cvLL + ridgePrepEMinternal(Y.rows(allButFold), 
						  Y.rows(fold),
						  ids.rows(allButFold),
						  ids.rows(fold),  
						  lambdaZE(0),
						  lambdaZE(1),  
						  targetZ,
						  targetE,  
						  nInit,
						  minSuccDiff, 
						  CVcrit);
	}
	
	// average over folds and return negative CV-loglikelihood
	return -cvLL / nFolds;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

// [[Rcpp::export(".armaRidgePrepEMdiag")]]
Rcpp::List ridgePrepEMdiag(arma::mat            Y, 
                           arma::ivec           ids, 
                           const double         lambdaZ,
                           const arma::mat&     targetZ,
                           const unsigned int&  nInit, 
                           const double&        minSuccDiff){

	/* ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size (including the 
	##                repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## lambdaZ      : ridge penalty parameter for the signal precision matrix
        ## lambdaE      : ridge penalty parameter for the error precision matrix
        ## targetZ      : target matrix towards the signal precision matrix 
	##                is shrunken
        ## targetE      : target matrix towards the error precision matrix 
	##                is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized 
	##                loglikehood) to be achieved
	## ---------------------------------------------------------------------
	## Value
	## A list-object with the following slots:
	## Pz		: estimated signal precision matrices
	## Pe		: estimated error precision matrices  
	## penLL	: penalized loglikelihood
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------- */


	//----------- variable declaration ----------------------------------//
	// extract number of samples and variates  
  	int          nk = Y.n_rows;
  	arma::ivec   is = arma::unique(ids);		
	unsigned int n  = is.n_elem;
	unsigned int p  = Y.n_cols;

	// declare variables used in penalized EM algorithm 
	double       penLL = -10e+10;
	double       penLLprev;
	arma::mat    Ytilde;
	arma::mat    Zis = arma::mat(n, p);
	unsigned int Ki;
	arma::uvec   kis;
	arma::mat    Yi;
	arma::uvec   i2uvec;
	arma::mat    Sz;
	arma::mat    Se;
	arma::mat    Yav = arma::mat(n,p);
	double       kweight = 0;
	arma::ivec   Ks = arma::ivec(n);
	arma::mat    VarZ  = arma::zeros(p,p);
	arma::mat    VarZk = arma::zeros(p,p);
	arma::mat    VarZslh;
	arma::mat    invSe;
	arma::mat    invSz;
	arma::vec    diagSe;
	arma::uvec   isameK;
	arma::mat    Sze;
	//----------- end of variable declaration ---------------------------//

	//----------- estimator initialization ------------------------------//
	// uses moment estimators followed by penalization
	// first, create nxp matrix with sample-wise averaged observations
	for (unsigned int i = 0; i < n; ++i){
		kis        = arma::find(ids == is(i));
		Ks(i)      = kis.n_elem;
		kweight    = kweight + 1/Ks(i);
		Yav.row(i) = arma::mean(Y.rows(kis));
	}
	
	// obtain moment estimators from moment decomposition
	Sz  = arma::trans(Yav) * Yav;
	Sz /= n;
	Se  = arma::trans(Y) * Y;
	Se /= nk;
	Se  = Se - Sz;
	Se /= (1 - kweight/n);
	Sz  = Sz - Se * (kweight/n);

	// obtain initial estimates from moment estimators
	Sz  = armaRidgeS(Sz, targetZ, lambdaZ/n);
	Se  = arma::diagmat(arma::diagvec(Se));
	//----------- end of initialization ---------------------------------//



	//----------- reformat data -----------------------------------------//
	// sort data by number of replications (for faster E-step calculation)
	arma::ivec uniqKs = arma::unique(Ks);
	if (Ks.min() != Ks.max()){
		if (uniqKs.n_elem > 1){
			arma::uvec shuffle = arma::uvec(nk);
			arma::uvec id2Ks;
			int counter = 0;
			for (unsigned int K = 0; K < uniqKs.n_elem; ++K){
				id2Ks = arma::find(Ks == uniqKs(K));
				for (unsigned int i = 0; i < id2Ks.n_elem; ++i){
					kis = arma::find(ids == is(id2Ks(i)));
					shuffle.subvec(counter, counter+kis.n_elem-1) = kis;
					counter = counter + kis.n_elem;
				}
			}

			// now reshuffle for ascending number of replications order
			Y   = Y.rows(shuffle);
			ids = ids.elem(shuffle);
		}
	}
	//----------- end of reformating ------------------------------------//


	//----------- EM algorithm ------------------------------------------//
	// start iterating
	bool stayInLoop = TRUE;
	for (unsigned int u = 0; u < nInit && stayInLoop; ++u){
		// store previously achieved penalized loglikelihood
		penLLprev = penLL;


		//----------- E-step ----------------------------------------//
		// obtain the two sufficient statistics
		// the balanced and unbalanced (replication-wise) cases are dealt 
		// with seperately, as the former allows to speed up computations
		
		// inverse of current signal and error covariance matrices
		// need in both steps
		invSz    = arma::inv_sympd(Sz);
		invSe    = arma::inv_sympd(Se);

		// unbalanced (replication-wise) design
		if (Ks.min() != Ks.max()){
			VarZ   = arma::zeros(p,p);
			VarZk  = arma::zeros(p,p);

			for (unsigned int r = 0; r < uniqKs.n_elem; ++r){	
				isameK = arma::find(Ks == uniqKs(r));
				Ki     = uniqKs(r);

				// calculate expectation of the signal, given data and parameters
				Sz  *= Ki;
				Sze  = arma::inv_sympd(Sz + Se) * Sz;

				for (unsigned int i = 0; i < isameK.n_elem; ++i){
					// obtain number of and row ids of repetitions
					kis               = arma::find(ids == is(isameK(i)));
					i2uvec            = arma::find(is  == is(isameK(i)));

					// calculate expectation of the signal, given data and parameters
					Zis.rows(i2uvec)  = arma::sum(Y.rows(kis), 0) * Sze;
					Zis.rows(i2uvec) /= Ki;
				}
				Sz /= Ki;

				// calculate variance of the signal, given data and parameters
				invSe        *= Ki;
				VarZslh       = arma::inv_sympd(invSz + invSe);
				VarZslh      *= isameK.n_elem;
				VarZ          = VarZ  + VarZslh;
				VarZslh      *= Ki;
				VarZk         = VarZk + VarZslh;
				invSe        /= Ki;
				VarZslh      /= Ki * isameK.n_elem;
			}

			// average the variances
			VarZ  /= n;
			VarZk /= nk;
		}

		// balanced (replication-wise) design
		if (Ks.min() == Ks.max()){
			// calculation of sufficient statistics
			Sz *= Ks.min();
			Sz  = arma::inv_sympd(Sz + Se) * Sz;
			for (unsigned int i = 0; i < n; ++i){
				// obtain number of and row ids of repetitions
				kis               = arma::find(ids == is(i));
				Ki                = kis.n_elem; 

				// calculate expectation of the signal, given data and parameters
				i2uvec            = arma::find(is == is(i));
				Zis.rows(i2uvec)  = arma::sum(Y.rows(kis), 0) * Sz;
				Zis.rows(i2uvec) /= Ki;
			}

			// calculate the variance of the signal, given data and parameters
			invSe        *= Ks.min();
			VarZk         = arma::inv_sympd(invSz + invSe);
			VarZ          = VarZk;
		}
		//----------- end of E-step ---------------------------------//



		//----------- M-step ----------------------------------------//
		// estimate the signal precision matrix
		Sz  = arma::trans(Zis) * Zis;
		Sz /= n;
		Sz  = Sz + VarZ;
		Sz  = armaRidgeS(Sz, targetZ, lambdaZ/n);

		// estimate the error precision matrix		
		diagSe  = arma::zeros(p);
		for (unsigned int i = 0; i < n; ++i){
			// obtain number of and row ids of repetitions
			kis            = arma::find(ids == is(i));
			Ki             = kis.n_elem; 
			if (Ki > 0){
				// calculate sample's contribution to sample error covariance matrix
				i2uvec         =  arma::find(is == is(i));
				Yi             =  Y.rows(kis);
				Yi.each_row() += -Zis.rows(i2uvec);
				diagSe        +=  arma::trans(arma::sum(arma::square(Yi)));
			}
		}
		diagSe /= nk;
		Se      = arma::diagmat(diagSe + arma::diagvec(VarZk));
		//----------- end of M-step ---------------------------------//



		//----------- convergence assessment ------------------------//
		// add the penalty's loglikelihood contribution
		// evaluate the loglikelihood
		penLL = ridgePrepLL(Y, ids, Sz, Se);

		// add the penalty
		penLL = penLL - 0.5 * lambdaZ * arma::accu(arma::square(arma::inv_sympd(Sz) - targetZ));				

		// assess convergence	
	    	if (std::abs((penLLprev - penLL)/penLLprev) < minSuccDiff){ 
			stayInLoop = FALSE; 
		}
		//----------- end of convergence assessment -----------------//	        
	}
	//----------- end of EM algorithm -----------------------------------//

	// return stuff as a list-object
	return Rcpp::List::create(Rcpp::Named("Pz")      = arma::inv_sympd(Sz), 
                                  Rcpp::Named("Pe")      = arma::diagmat(arma::pow(arma::diagvec(Se), -1)), 
                                  Rcpp::Named("penLL")   = penLL);
}


inline double ridgePrepEMdiagInternal(arma::mat         Ytrain,
				      const arma::mat&  Ytest,  
	                              arma::ivec        idsTrain,
	                              const arma::ivec& idsTest,  
	                              const double      lambdaZ,
        	                      const arma::mat&  targetZ,
        	                      const int&        nInit, 
        	                      const double&     minSuccDiff,
				      std::string       CVcrit){

	/* ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a mixture of GMMs
	## ---------------------------------------------------------------------
	## Arguments 
	## Y		: n*p data matrix where n is sample size (including the 
	##                repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## lambdaZ      : ridge penalty parameter for the signal precision matrix
        ## targetZ      : target matrix towards the signal precision matrix 
	##                is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized 
	##                loglikehood) to be achieved
	## ---------------------------------------------------------------------
	## Value
	## The (negative?) cross-validated loss
	## ---------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	--------------------------------------------------------------------- */



	//----------- variable declaration ----------------------------------//
	// extract number of samples and variates  
  	int          nk = Ytrain.n_rows;
  	arma::ivec   is = arma::unique(idsTrain);		
	unsigned int n  = is.n_elem;
	unsigned int p  = Ytrain.n_cols;

	// declare variables used in penalized EM algorithm 
	double       penLL = -10e+10;
	double       penLLprev;
	arma::mat    Ytilde;
	arma::mat    Zis = arma::mat(n, p);
	unsigned int Ki;
	arma::uvec   kis;
	arma::mat    Yi;
	arma::uvec   i2uvec;
	arma::mat    Sz;
	arma::mat    Se;
	arma::mat    Yav = arma::mat(n,p);
	double       kweight = 0;
	arma::ivec   Ks = arma::ivec(n);
	arma::mat    VarZ  = arma::zeros(p,p);
	arma::mat    VarZk = arma::zeros(p,p);
	arma::mat    VarZslh;
	arma::mat    invSe;
	arma::mat    invSz;
	arma::vec    diagSe;
	arma::uvec   isameK;
	arma::mat    Sze;
	double       CVloss = 10e-10;
	//----------- end of variable declaration ---------------------------//



	//----------- estimator initialization ------------------------------//
	// uses moment estimators followed by penalization
	// first, create nxp matrix with sample-wise averaged observations
	for (unsigned int i = 0; i < n; ++i){
		kis        = arma::find(idsTrain == is(i));
		Ks(i)      = kis.n_elem;
		kweight    = kweight + 1/Ks(i);
		Yav.row(i) = arma::mean(Ytrain.rows(kis));
	}
	
	// obtain moment estimators from moment decomposition
	Sz  = arma::trans(Yav) * Yav;
	Sz /= n;
	Se  = arma::trans(Ytrain) * Ytrain;
	Se /= nk;
	Se  = Se - Sz;
	Se /= (1 - kweight/n);
	Sz  = Sz - Se * (kweight/n);

	// obtain initial estimates from moment estimators
	Sz  = armaRidgeS(Sz, targetZ, lambdaZ/n);
	Se  = arma::diagmat(arma::diagvec(Se));
	//----------- end of initialization ---------------------------------//



	//----------- reformat data -----------------------------------------//
	// sort data by number of replications (for faster E-step calculation)
	arma::ivec uniqKs = arma::unique(Ks);
	if (uniqKs.n_elem > 1){
		arma::uvec shuffle = arma::uvec(nk);
		arma::uvec id2Ks;
		int counter = 0;
		for (unsigned int K = 0; K < uniqKs.n_elem; ++K){
			id2Ks = arma::find(Ks == uniqKs(K));
			for (unsigned int i = 0; i < id2Ks.n_elem; ++i){
				kis = arma::find(idsTrain == is(id2Ks(i)));
				shuffle.subvec(counter, counter+kis.n_elem-1) = kis;
				counter = counter + kis.n_elem;
			}
		}

		// now reshuffle for ascending number of replications order
		Ytrain   = Ytrain.rows(shuffle);
		idsTrain = idsTrain.elem(shuffle);
	}
	//----------- end of reformating ------------------------------------//



	//----------- EM algorithm ------------------------------------------//
	// start iterating
	bool stayInLoop = TRUE;
	for (unsigned int u = 0; u < std::abs(nInit) && stayInLoop; ++u){
		// store previously achieved penalized loglikelihood
		penLLprev = penLL;


		//----------- E-step ----------------------------------------//
		// obtain the two sufficient statistics
		// the balanced and unbalanced (replication-wise) cases are dealt 
		// with seperately, as the former allows to speed up computations
		
		// inverse of current signal and error covariance matrices
		// need in both steps
		invSz    = arma::inv_sympd(Sz);
		invSe    = arma::inv_sympd(Se);

		// unbalanced (replication-wise) design
		if (Ks.min() != Ks.max()){
			VarZ   = arma::zeros(p,p);
			VarZk  = arma::zeros(p,p);

			for (unsigned int r = 0; r < uniqKs.n_elem; ++r){	
				isameK = arma::find(Ks == uniqKs(r));
				Ki     = uniqKs(r);

				// calculate expectation of the signal, given data and parameters
				Sz  *= Ki;
				Sze  = arma::inv_sympd(Sz + Se) * Sz;

				for (unsigned int i = 0; i < isameK.n_elem; ++i){
					// obtain number of and row ids of repetitions
					kis               = arma::find(idsTrain == is(isameK(i)));
					i2uvec            = arma::find(is       == is(isameK(i)));

					// calculate expectation of the signal, given data and parameters
					Zis.rows(i2uvec)  = arma::sum(Ytrain.rows(kis), 0) * Sze;
					Zis.rows(i2uvec) /= Ki;
				}
				Sz /= Ki;

				// calculate variance of the signal, given data and parameters
				invSe        *= Ki;
				VarZslh       = arma::inv_sympd(invSz + invSe);
				VarZslh      *= isameK.n_elem;
				VarZ          = VarZ  + VarZslh;
				VarZslh      *= Ki;
				VarZk         = VarZk + VarZslh;
				invSe        /= Ki;
				VarZslh      /= Ki * isameK.n_elem;
			}

			// average the variances
			VarZ  /= n;
			VarZk /= nk;
		}

		// balanced (replication-wise) design
		if (Ks.min() == Ks.max()){
			// calculation of sufficient statistics
			Sz *= Ks.min();
			Sz  = arma::inv_sympd(Sz + Se) * Sz;
			for (unsigned int i = 0; i < n; ++i){
				// obtain number of and row ids of repetitions
				kis               = arma::find(idsTrain == is(i));
				Ki                = kis.n_elem; 

				// calculate expectation of the signal, given data and parameters
				i2uvec            = arma::find(is == is(i));
				Zis.rows(i2uvec)  = arma::sum(Ytrain.rows(kis), 0) * Sz;
				Zis.rows(i2uvec) /= Ki;
			}

			// calculate the variance of the signal, given data and parameters
			invSe        *= Ks.min();
			VarZk         = arma::inv_sympd(invSz + invSe);
			VarZ          = VarZk;
		}
		//----------- end of E-step ---------------------------------//



		//----------- M-step ----------------------------------------//
		// estimate the signal precision matrix
		Sz  = arma::trans(Zis) * Zis;
		Sz /= n;
		Sz  = Sz + VarZ;
		Sz  = armaRidgeS(Sz, targetZ, lambdaZ/n);

		// estimate the error precision matrix		
		diagSe  = arma::zeros(p);
		for (unsigned int i = 0; i < n; ++i){
			// obtain number of and row ids of repetitions
			kis            = arma::find(idsTrain == is(i));
			Ki             = kis.n_elem; 
			if (Ki > 0){
				// calculate sample's contribution to sample error covariance matrix
				i2uvec         =  arma::find(is == is(i));
				Yi             =  Ytrain.rows(kis);
				Yi.each_row() += -Zis.rows(i2uvec);
				diagSe        +=  arma::trans(arma::sum(arma::square(Yi)));
			}
		}
		diagSe /= nk;
		Se      = arma::diagmat(diagSe + arma::diagvec(VarZk));
		//----------- end of M-step ---------------------------------//



		//----------- convergence assessment ------------------------//
		// add the penalty's loglikelihood contribution
		// evaluate the loglikelihood
		penLL = ridgePrepLL(Ytrain, idsTrain, Sz, Se);
		// add the penalty
		penLL = penLL - 0.5 * lambdaZ * arma::accu(arma::square(arma::inv_sympd(Sz) - targetZ));

		// assess convergence	
	    	if (std::abs((penLLprev - penLL)/penLLprev) < minSuccDiff){ 
			stayInLoop = FALSE; 
		}
		//----------- end of convergence assessment ------------------------//	     
	}
	//----------- end of EM algorithm -----------------------------------//

	// select cross-validation criterion
	if (CVcrit == "LL"){ 
		CVloss = ridgePrepLL(Ytest, idsTest, Sz, Se);
	}
	if (CVcrit == "compLL"){ 
		CVloss = ridgePrepCompLL(Ytest, idsTest, Sz, Se);
	}
	if (CVcrit == "Qloss"){ 
		CVloss = ridgePrepQloss(Ytest, idsTest, Sz, Se);
	}
	if (CVcrit == "Qloss2"){ 
		CVloss = ridgePrepQloss2(Ytest, idsTest, Sz, Se);
	}
	return CVloss;
}



// [[Rcpp::export(".armaRidgePrepKcvLLdiag")]]
double ridgePrepKcvLLdiag(const double      lambdaZ,
                          const arma::mat&  Y, 
                          const arma::ivec& ids, 
                          const arma::mat&  targetZ,
                          const int&        nInit, 
                          const double&     minSuccDiff,
                          Rcpp::List        folds,
		          std::string       CVcrit){

	/* --------------------------------------------------------------------------------------------------
	## Cross-validated loglikehood of a mixture of GGMs
	## --------------------------------------------------------------------------------------------------
	## Arguments 
        ## lambdaZ      : ridge penalty parameter for the signal precision matrix
	## Y		: n*p data matrix where n is sample size (including the repetitions) and p is dimension.
	## ids		: numeric indicating which rows of Y belong to the same 
	##                individal
        ## targetZ      : target matrix towards the signal precision matrix is shrunken
	## nInit	: maximum number of iterations
	## minSuccDiff	: minimum successive difference (in the penalized 
	##                loglikehood) to be achieved
	## folds        : cross-validation sample splits.
	## --------------------------------------------------------------------------------------------------
	## Value
	## cvLL		: negative cross-validated loglikelihood
	##                (negation as function to be minimized in cross-validation)
	## --------------------------------------------------------------------------------------------------
	## Authors      : Wessel N. van Wieringen
	-------------------------------------------------------------------------------------------------- */

	// extract number samples and folds
  	unsigned int n      = Y.n_rows;
	unsigned int nFolds = folds.size();

	// declare variables used in the CV-loop
	double     cvLL = 0;
	arma::uvec fold, allButFold;
	Rcpp::List ridgePrepHat;
	int        slh;

	// start of cross-validation
	for (unsigned int f = 0; f < nFolds; ++f){
		// sort out the fold and its complement
		fold       = Rcpp::as<arma::uvec>(folds(f)) - 1;
		allButFold = arma::uvec(n-fold.n_elem);
		slh        = 0;
		for (unsigned int i = 0; i < n; ++i){
			if (all(fold != i)){
				allButFold(slh) = i;
				slh = slh + 1;
			}
		}

		// estimate model from all but left-out fold
		cvLL = cvLL + ridgePrepEMdiagInternal(Y.rows(allButFold), 
						  Y.rows(fold),
						  ids.rows(allButFold),
						  ids.rows(fold),  
						  lambdaZ,
						  targetZ,
						  nInit,
						  minSuccDiff, 
						  CVcrit);
	}
	
	// average over folds and return negative CV-loglikelihood
	return -cvLL / nFolds;
}


