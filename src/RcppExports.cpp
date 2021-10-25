// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// colSumSq_dge
NumericVector colSumSq_dge(NumericVector x, IntegerVector dim);
RcppExport SEXP _bravo_colSumSq_dge(SEXP xSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(colSumSq_dge(x, dim));
    return rcpp_result_gen;
END_RCPP
}
// colSumSq_matrix
NumericVector colSumSq_matrix(NumericMatrix x);
RcppExport SEXP _bravo_colSumSq_matrix(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(colSumSq_matrix(x));
    return rcpp_result_gen;
END_RCPP
}
// colMSD_dgc
NumericVector colMSD_dgc(S4 mat, NumericVector m);
RcppExport SEXP _bravo_colMSD_dgc(SEXP matSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(colMSD_dgc(mat, m));
    return rcpp_result_gen;
END_RCPP
}
// colSUMIDX_dgc
NumericVector colSUMIDX_dgc(S4 mat);
RcppExport SEXP _bravo_colSUMIDX_dgc(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(colSUMIDX_dgc(mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bravo_colSumSq_dge", (DL_FUNC) &_bravo_colSumSq_dge, 2},
    {"_bravo_colSumSq_matrix", (DL_FUNC) &_bravo_colSumSq_matrix, 1},
    {"_bravo_colMSD_dgc", (DL_FUNC) &_bravo_colMSD_dgc, 2},
    {"_bravo_colSUMIDX_dgc", (DL_FUNC) &_bravo_colSUMIDX_dgc, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_bravo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
