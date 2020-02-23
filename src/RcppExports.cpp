// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace Eigen;

// KRP
MatrixXd KRP(MatrixXd A, MatrixXd B);
RcppExport SEXP _tensorApp_KRP(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(KRP(A, B));
    return rcpp_result_gen;
END_RCPP
}
// TuckerALS
MatrixXd TuckerALS(MatrixXd T1, int d0, VectorXi dims, VectorXi rs, List D0, List optsList);
RcppExport SEXP _tensorApp_TuckerALS(SEXP T1SEXP, SEXP d0SEXP, SEXP dimsSEXP, SEXP rsSEXP, SEXP D0SEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T1(T1SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type rs(rsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(TuckerALS(T1, d0, dims, rs, D0, optsList));
    return rcpp_result_gen;
END_RCPP
}
// CPTPM
MatrixXd CPTPM(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList);
RcppExport SEXP _tensorApp_CPTPM(SEXP T0SEXP, SEXP d0SEXP, SEXP dSEXP, SEXP dimsSEXP, SEXP D0SEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(CPTPM(T0, d0, d, dims, D0, optsList));
    return rcpp_result_gen;
END_RCPP
}
// CPTPMorthogon
MatrixXd CPTPMorthogon(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList);
RcppExport SEXP _tensorApp_CPTPMorthogon(SEXP T0SEXP, SEXP d0SEXP, SEXP dSEXP, SEXP dimsSEXP, SEXP D0SEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(CPTPMorthogon(T0, d0, d, dims, D0, optsList));
    return rcpp_result_gen;
END_RCPP
}
// CPTPM_dr
MatrixXd CPTPM_dr(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList);
RcppExport SEXP _tensorApp_CPTPM_dr(SEXP T0SEXP, SEXP d0SEXP, SEXP dSEXP, SEXP dimsSEXP, SEXP D0SEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(CPTPM_dr(T0, d0, d, dims, D0, optsList));
    return rcpp_result_gen;
END_RCPP
}
// CPTPMorthogon_dr
MatrixXd CPTPMorthogon_dr(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList);
RcppExport SEXP _tensorApp_CPTPMorthogon_dr(SEXP T0SEXP, SEXP d0SEXP, SEXP dSEXP, SEXP dimsSEXP, SEXP D0SEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(CPTPMorthogon_dr(T0, d0, d, dims, D0, optsList));
    return rcpp_result_gen;
END_RCPP
}
// CPALS
MatrixXd CPALS(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList);
RcppExport SEXP _tensorApp_CPALS(SEXP T0SEXP, SEXP d0SEXP, SEXP dSEXP, SEXP dimsSEXP, SEXP D0SEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(CPALS(T0, d0, d, dims, D0, optsList));
    return rcpp_result_gen;
END_RCPP
}
// CPTPMsym2
MatrixXd CPTPMsym2(MatrixXd T0, int d, int k1, int k2, VectorXi dims, List D0, List optsList);
RcppExport SEXP _tensorApp_CPTPMsym2(SEXP T0SEXP, SEXP dSEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP dimsSEXP, SEXP D0SEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(CPTPMsym2(T0, d, k1, k2, dims, D0, optsList));
    return rcpp_result_gen;
END_RCPP
}
// CPTPMsym2Orth
MatrixXd CPTPMsym2Orth(MatrixXd T0, int d, int k1, int k2, VectorXi dims, List D0, List optsList, List optsList_pen);
RcppExport SEXP _tensorApp_CPTPMsym2Orth(SEXP T0SEXP, SEXP dSEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP dimsSEXP, SEXP D0SEXP, SEXP optsListSEXP, SEXP optsList_penSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    Rcpp::traits::input_parameter< List >::type optsList_pen(optsList_penSEXP);
    rcpp_result_gen = Rcpp::wrap(CPTPMsym2Orth(T0, d, k1, k2, dims, D0, optsList, optsList_pen));
    return rcpp_result_gen;
END_RCPP
}
// TuckerALSsym2
MatrixXd TuckerALSsym2(MatrixXd T0, int k1, int k2, VectorXi dims, VectorXi rs, List D0, List optsList);
RcppExport SEXP _tensorApp_TuckerALSsym2(SEXP T0SEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP dimsSEXP, SEXP rsSEXP, SEXP D0SEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type rs(rsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(TuckerALSsym2(T0, k1, k2, dims, rs, D0, optsList));
    return rcpp_result_gen;
END_RCPP
}
// SCPTPM
List SCPTPM(MatrixXd T0, int d0, int d, VectorXi dims, List D1, VectorXd lambda, List optsList, List optsList_pen);
RcppExport SEXP _tensorApp_SCPTPM(SEXP T0SEXP, SEXP d0SEXP, SEXP dSEXP, SEXP dimsSEXP, SEXP D1SEXP, SEXP lambdaSEXP, SEXP optsListSEXP, SEXP optsList_penSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D1(D1SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    Rcpp::traits::input_parameter< List >::type optsList_pen(optsList_penSEXP);
    rcpp_result_gen = Rcpp::wrap(SCPTPM(T0, d0, d, dims, D1, lambda, optsList, optsList_pen));
    return rcpp_result_gen;
END_RCPP
}
// SCPTPM_part
List SCPTPM_part(MatrixXd T0, int d0, int d, VectorXi dims, VectorXi actives, List D1, VectorXd lambda, List optsList, List optsList_pen);
RcppExport SEXP _tensorApp_SCPTPM_part(SEXP T0SEXP, SEXP d0SEXP, SEXP dSEXP, SEXP dimsSEXP, SEXP activesSEXP, SEXP D1SEXP, SEXP lambdaSEXP, SEXP optsListSEXP, SEXP optsList_penSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type actives(activesSEXP);
    Rcpp::traits::input_parameter< List >::type D1(D1SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    Rcpp::traits::input_parameter< List >::type optsList_pen(optsList_penSEXP);
    rcpp_result_gen = Rcpp::wrap(SCPTPM_part(T0, d0, d, dims, actives, D1, lambda, optsList, optsList_pen));
    return rcpp_result_gen;
END_RCPP
}
// TransferModalUnfoldingsT
MatrixXd TransferModalUnfoldingsT(MatrixXd T, int d1, int d2, VectorXi dims);
RcppExport SEXP _tensorApp_TransferModalUnfoldingsT(SEXP TSEXP, SEXP d1SEXP, SEXP d2SEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< int >::type d2(d2SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(TransferModalUnfoldingsT(T, d1, d2, dims));
    return rcpp_result_gen;
END_RCPP
}
// gtsem0
MatrixXd gtsem0(MatrixXd S, int r1, int r2, VectorXi dims);
RcppExport SEXP _tensorApp_gtsem0(SEXP SSEXP, SEXP r1SEXP, SEXP r2SEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< int >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(gtsem0(S, r1, r2, dims));
    return rcpp_result_gen;
END_RCPP
}
// setuplambdaPC
VectorXd setuplambdaPC(MatrixXd T0, int d0, VectorXi dims, List D0, int nlam, VectorXd setlam);
RcppExport SEXP _tensorApp_setuplambdaPC(SEXP T0SEXP, SEXP d0SEXP, SEXP dimsSEXP, SEXP D0SEXP, SEXP nlamSEXP, SEXP setlamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< List >::type D0(D0SEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type setlam(setlamSEXP);
    rcpp_result_gen = Rcpp::wrap(setuplambdaPC(T0, d0, dims, D0, nlam, setlam));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tensorApp_KRP", (DL_FUNC) &_tensorApp_KRP, 2},
    {"_tensorApp_TuckerALS", (DL_FUNC) &_tensorApp_TuckerALS, 6},
    {"_tensorApp_CPTPM", (DL_FUNC) &_tensorApp_CPTPM, 6},
    {"_tensorApp_CPTPMorthogon", (DL_FUNC) &_tensorApp_CPTPMorthogon, 6},
    {"_tensorApp_CPTPM_dr", (DL_FUNC) &_tensorApp_CPTPM_dr, 6},
    {"_tensorApp_CPTPMorthogon_dr", (DL_FUNC) &_tensorApp_CPTPMorthogon_dr, 6},
    {"_tensorApp_CPALS", (DL_FUNC) &_tensorApp_CPALS, 6},
    {"_tensorApp_CPTPMsym2", (DL_FUNC) &_tensorApp_CPTPMsym2, 7},
    {"_tensorApp_CPTPMsym2Orth", (DL_FUNC) &_tensorApp_CPTPMsym2Orth, 8},
    {"_tensorApp_TuckerALSsym2", (DL_FUNC) &_tensorApp_TuckerALSsym2, 7},
    {"_tensorApp_SCPTPM", (DL_FUNC) &_tensorApp_SCPTPM, 8},
    {"_tensorApp_SCPTPM_part", (DL_FUNC) &_tensorApp_SCPTPM_part, 9},
    {"_tensorApp_TransferModalUnfoldingsT", (DL_FUNC) &_tensorApp_TransferModalUnfoldingsT, 4},
    {"_tensorApp_gtsem0", (DL_FUNC) &_tensorApp_gtsem0, 4},
    {"_tensorApp_setuplambdaPC", (DL_FUNC) &_tensorApp_setuplambdaPC, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_tensorApp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
