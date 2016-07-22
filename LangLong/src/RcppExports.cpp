// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// estimateLong_cpp
List estimateLong_cpp(Rcpp::List in_list);
RcppExport SEXP LANGlong_estimateLong_cpp(SEXP in_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List >::type in_list(in_listSEXP);
    __result = Rcpp::wrap(estimateLong_cpp(in_list));
    return __result;
END_RCPP
}
// predictLong_cpp
List predictLong_cpp(Rcpp::List in_list);
RcppExport SEXP LANGlong_predictLong_cpp(SEXP in_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List >::type in_list(in_listSEXP);
    __result = Rcpp::wrap(predictLong_cpp(in_list));
    return __result;
END_RCPP
}
// simulateLong_cpp
List simulateLong_cpp(List obs_, List operator_, List theta_);
RcppExport SEXP LANGlong_simulateLong_cpp(SEXP obs_SEXP, SEXP operator_SEXP, SEXP theta_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type obs_(obs_SEXP);
    Rcpp::traits::input_parameter< List >::type operator_(operator_SEXP);
    Rcpp::traits::input_parameter< List >::type theta_(theta_SEXP);
    __result = Rcpp::wrap(simulateLong_cpp(obs_, operator_, theta_));
    return __result;
END_RCPP
}
// testSimulateX_cpp
List testSimulateX_cpp(List operator_, List theta_);
RcppExport SEXP LANGlong_testSimulateX_cpp(SEXP operator_SEXP, SEXP theta_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type operator_(operator_SEXP);
    Rcpp::traits::input_parameter< List >::type theta_(theta_SEXP);
    __result = Rcpp::wrap(testSimulateX_cpp(operator_, theta_));
    return __result;
END_RCPP
}
// testSimulateX2_cpp
List testSimulateX2_cpp(List operator_, List theta_);
RcppExport SEXP LANGlong_testSimulateX2_cpp(SEXP operator_SEXP, SEXP theta_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type operator_(operator_SEXP);
    Rcpp::traits::input_parameter< List >::type theta_(theta_SEXP);
    __result = Rcpp::wrap(testSimulateX2_cpp(operator_, theta_));
    return __result;
END_RCPP
}
// simulateLongGH_cpp
List simulateLongGH_cpp(Rcpp::List in_list);
RcppExport SEXP LANGlong_simulateLongGH_cpp(SEXP in_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List >::type in_list(in_listSEXP);
    __result = Rcpp::wrap(simulateLongGH_cpp(in_list));
    return __result;
END_RCPP
}
// estimateProcess_cpp
List estimateProcess_cpp(Rcpp::List in_list);
RcppExport SEXP LANGlong_estimateProcess_cpp(SEXP in_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List >::type in_list(in_listSEXP);
    __result = Rcpp::wrap(estimateProcess_cpp(in_list));
    return __result;
END_RCPP
}
