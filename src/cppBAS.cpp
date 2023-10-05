// cppBAS.cpp

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

#include <Rcpp.h>
#include <RcppThread.h>

//#include "cppBASMastersample.h"

using namespace Rcpp;

#define RCPPTHREAD_OVERRIDE_COUT 1   // std::cout override
#define RCPPTHREAD_OVERRIDE_THREAD 1 // std::thread override


