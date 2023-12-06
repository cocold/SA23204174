// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// EM
NumericMatrix EM(int nc, int ni, int nt);
RcppExport SEXP _SA23204174_EM(SEXP ncSEXP, SEXP niSEXP, SEXP ntSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type ni(niSEXP);
    Rcpp::traits::input_parameter< int >::type nt(ntSEXP);
    rcpp_result_gen = Rcpp::wrap(EM(nc, ni, nt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA23204174_EM", (DL_FUNC) &_SA23204174_EM, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA23204174(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
