// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SolveCongruence
NumericVector SolveCongruence(NumericMatrix& A, NumericVector& base, NumericVector J);
RcppExport SEXP _uc511_SolveCongruence(SEXP ASEXP, SEXP baseSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type base(baseSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(SolveCongruence(A, base, J));
    return rcpp_result_gen;
END_RCPP
}
// GetBoxIndices
NumericMatrix GetBoxIndices(NumericMatrix& lxy, IntegerVector& base, IntegerVector J);
RcppExport SEXP _uc511_GetBoxIndices(SEXP lxySEXP, SEXP baseSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type lxy(lxySEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type base(baseSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(GetBoxIndices(lxy, base, J));
    return rcpp_result_gen;
END_RCPP
}
// cppHaltonSeq
NumericVector cppHaltonSeq(const int& k, double& base, int& n);
RcppExport SEXP _uc511_cppHaltonSeq(SEXP kSEXP, SEXP baseSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< double& >::type base(baseSEXP);
    Rcpp::traits::input_parameter< int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(cppHaltonSeq(k, base, n));
    return rcpp_result_gen;
END_RCPP
}
// compareBoxesBoxInit
Rcpp::NumericVector compareBoxesBoxInit(Rcpp::NumericVector boxes, Rcpp::NumericVector boxInit, int intB);
RcppExport SEXP _uc511_compareBoxesBoxInit(SEXP boxesSEXP, SEXP boxInitSEXP, SEXP intBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type boxes(boxesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type boxInit(boxInitSEXP);
    Rcpp::traits::input_parameter< int >::type intB(intBSEXP);
    rcpp_result_gen = Rcpp::wrap(compareBoxesBoxInit(boxes, boxInit, intB));
    return rcpp_result_gen;
END_RCPP
}
// cppProductPoweredElements
int cppProductPoweredElements(NumericVector& J, NumericVector& bases, int numElements);
RcppExport SEXP _uc511_cppProductPoweredElements(SEXP JSEXP, SEXP basesSEXP, SEXP numElementsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type J(JSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type bases(basesSEXP);
    Rcpp::traits::input_parameter< int >::type numElements(numElementsSEXP);
    rcpp_result_gen = Rcpp::wrap(cppProductPoweredElements(J, bases, numElements));
    return rcpp_result_gen;
END_RCPP
}
// cppWhere2Start
NumericVector cppWhere2Start(NumericVector& J, IntegerVector& seeds, NumericVector& bases, NumericVector& boxes);
RcppExport SEXP _uc511_cppWhere2Start(SEXP JSEXP, SEXP seedsSEXP, SEXP basesSEXP, SEXP boxesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type J(JSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type seeds(seedsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type bases(basesSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type boxes(boxesSEXP);
    rcpp_result_gen = Rcpp::wrap(cppWhere2Start(J, seeds, bases, boxes));
    return rcpp_result_gen;
END_RCPP
}
// log_a_to_base_b
double log_a_to_base_b(long long a, int b);
RcppExport SEXP _uc511_log_a_to_base_b(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< long long >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(log_a_to_base_b(a, b));
    return rcpp_result_gen;
END_RCPP
}
// cppRSHalton
NumericVector cppRSHalton(int n, IntegerVector seeds, NumericVector bases, NumericVector boxes, NumericVector J);
RcppExport SEXP _uc511_cppRSHalton(SEXP nSEXP, SEXP seedsSEXP, SEXP basesSEXP, SEXP boxesSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type seeds(seedsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bases(basesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type boxes(boxesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(cppRSHalton(n, seeds, bases, boxes, J));
    return rcpp_result_gen;
END_RCPP
}
// sample_int
Rcpp::IntegerVector sample_int(int n, int min, int max);
RcppExport SEXP _uc511_sample_int(SEXP nSEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type min(minSEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int(n, min, max));
    return rcpp_result_gen;
END_RCPP
}
// removeDuplicates
Rcpp::NumericVector removeDuplicates(Rcpp::NumericVector vec);
RcppExport SEXP _uc511_removeDuplicates(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(removeDuplicates(vec));
    return rcpp_result_gen;
END_RCPP
}
// cppBASpts
Rcpp::List cppBASpts(int n, IntegerVector seeds, NumericVector bases);
RcppExport SEXP _uc511_cppBASpts(SEXP nSEXP, SEXP seedsSEXP, SEXP basesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type seeds(seedsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bases(basesSEXP);
    rcpp_result_gen = Rcpp::wrap(cppBASpts(n, seeds, bases));
    return rcpp_result_gen;
END_RCPP
}
// cppRSHalton_br
Rcpp::List cppRSHalton_br(int n, NumericVector bases, NumericVector seeds);
RcppExport SEXP _uc511_cppRSHalton_br(SEXP nSEXP, SEXP basesSEXP, SEXP seedsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bases(basesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type seeds(seedsSEXP);
    rcpp_result_gen = Rcpp::wrap(cppRSHalton_br(n, bases, seeds));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_uc511_SolveCongruence", (DL_FUNC) &_uc511_SolveCongruence, 3},
    {"_uc511_GetBoxIndices", (DL_FUNC) &_uc511_GetBoxIndices, 3},
    {"_uc511_cppHaltonSeq", (DL_FUNC) &_uc511_cppHaltonSeq, 3},
    {"_uc511_compareBoxesBoxInit", (DL_FUNC) &_uc511_compareBoxesBoxInit, 3},
    {"_uc511_cppProductPoweredElements", (DL_FUNC) &_uc511_cppProductPoweredElements, 3},
    {"_uc511_cppWhere2Start", (DL_FUNC) &_uc511_cppWhere2Start, 4},
    {"_uc511_log_a_to_base_b", (DL_FUNC) &_uc511_log_a_to_base_b, 2},
    {"_uc511_cppRSHalton", (DL_FUNC) &_uc511_cppRSHalton, 5},
    {"_uc511_sample_int", (DL_FUNC) &_uc511_sample_int, 3},
    {"_uc511_removeDuplicates", (DL_FUNC) &_uc511_removeDuplicates, 1},
    {"_uc511_cppBASpts", (DL_FUNC) &_uc511_cppBASpts, 3},
    {"_uc511_cppRSHalton_br", (DL_FUNC) &_uc511_cppRSHalton_br, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_uc511(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
