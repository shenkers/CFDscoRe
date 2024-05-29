// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// private_optimal_alignment
Rcpp::List private_optimal_alignment(Rcpp::List activity_scores, Rcpp::CharacterVector query, Rcpp::CharacterVector genome);
RcppExport SEXP _CFDscoRe_private_optimal_alignment(SEXP activity_scoresSEXP, SEXP querySEXP, SEXP genomeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type activity_scores(activity_scoresSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type query(querySEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type genome(genomeSEXP);
    rcpp_result_gen = Rcpp::wrap(private_optimal_alignment(activity_scores, query, genome));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CFDscoRe_private_optimal_alignment", (DL_FUNC) &_CFDscoRe_private_optimal_alignment, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CFDscoRe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}