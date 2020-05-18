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

// pull in functions 
// (the order below matters if there is cross-reference)
#include "ridgePgen.h"
#include "ridgeGGMmix.h"
#include "ridgePrep.h"
#include "ridgePmultiT.h"


