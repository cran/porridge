/* ---------------------------------------------------------------------------
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Below are some functions related to
van Wieringen, W.N., Chen, Y. (2019), "Penalized estimation of the
Gaussian graphical model from data with replicates", submitted
------------------------------------------------------------------------------
------------------------------------------------------------------------------
--------------------------------------------------------------------------- */


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
	## Loglikehood of a simple multivariate random effects model
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
	## Ridge penalized EM algorithm for a simple multivariate random effects 
	## model
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
	## Ridge penalized EM algorithm for a simple multivariate random effects 
	## model
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
	## Cross-validated loglikehood a simple multivariate random effects model
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



// [[Rcpp::export(".armaRidgePrepEMdiag")]]
Rcpp::List ridgePrepEMdiag(arma::mat            Y, 
                           arma::ivec           ids, 
                           const double         lambdaZ,
                           const arma::mat&     targetZ,
                           const unsigned int&  nInit, 
                           const double&        minSuccDiff){

	/* ---------------------------------------------------------------------
	## Ridge penalized EM algorithm for a simple multivariate random effects 
	## model
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
	// Se  = arma::diagmat(arma::diagvec(Se));
	diagSe = arma::diagvec(Se);
	double diagSeMax = std::abs(diagSe.max());
	arma::uvec diagSeIDneg = arma::find(diagSe <= 0);
	if (diagSeIDneg.n_elem > 0){
		for (unsigned int K = 0; K < diagSeIDneg.n_elem; ++K){
			diagSe(diagSeIDneg(K)) = diagSeMax;
		}
	}
	Se  = arma::diagmat(diagSe);
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
	## Ridge penalized EM algorithm for a simple multivariate random 
	## effects model
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
	// Se  = arma::diagmat(arma::diagvec(Se));
	diagSe = arma::diagvec(Se);
	double diagSeMax = std::abs(diagSe.max())+0.1;
	arma::uvec diagSeIDneg = arma::find(diagSe <= 0);
	if (diagSeIDneg.n_elem > 0){
		for (unsigned int K = 0; K < diagSeIDneg.n_elem; ++K){
			diagSe(diagSeIDneg(K)) = diagSeMax;
		}
	}
	Se  = arma::diagmat(diagSe);
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
	## Cross-validated loglikehood of a simple multivariate random effects model
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



