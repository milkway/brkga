// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// tsp_brkga
Rcpp::List tsp_brkga(std::string instanceFile);
RcppExport SEXP _brkga_tsp_brkga(SEXP instanceFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type instanceFile(instanceFileSEXP);
    rcpp_result_gen = Rcpp::wrap(tsp_brkga(instanceFile));
    return rcpp_result_gen;
END_RCPP
}
// mdp_brkga
Rcpp::List mdp_brkga(const arma::mat DistanceMatrix, const unsigned m, const unsigned method, const unsigned MAX_TIME, const unsigned p, const double pe, const double pm, const double rhoe, const unsigned K, const unsigned MAXT, const unsigned X_INTVL, const unsigned X_NUMBER, const unsigned MAX_GENS, const unsigned RESET_AFTER, const bool verbose, const long unsigned rngSeed);
RcppExport SEXP _brkga_mdp_brkga(SEXP DistanceMatrixSEXP, SEXP mSEXP, SEXP methodSEXP, SEXP MAX_TIMESEXP, SEXP pSEXP, SEXP peSEXP, SEXP pmSEXP, SEXP rhoeSEXP, SEXP KSEXP, SEXP MAXTSEXP, SEXP X_INTVLSEXP, SEXP X_NUMBERSEXP, SEXP MAX_GENSSEXP, SEXP RESET_AFTERSEXP, SEXP verboseSEXP, SEXP rngSeedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type DistanceMatrix(DistanceMatrixSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type m(mSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type MAX_TIME(MAX_TIMESEXP);
    Rcpp::traits::input_parameter< const unsigned >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type pe(peSEXP);
    Rcpp::traits::input_parameter< const double >::type pm(pmSEXP);
    Rcpp::traits::input_parameter< const double >::type rhoe(rhoeSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type K(KSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type MAXT(MAXTSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type X_INTVL(X_INTVLSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type X_NUMBER(X_NUMBERSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type MAX_GENS(MAX_GENSSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type RESET_AFTER(RESET_AFTERSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const long unsigned >::type rngSeed(rngSeedSEXP);
    rcpp_result_gen = Rcpp::wrap(mdp_brkga(DistanceMatrix, m, method, MAX_TIME, p, pe, pm, rhoe, K, MAXT, X_INTVL, X_NUMBER, MAX_GENS, RESET_AFTER, verbose, rngSeed));
    return rcpp_result_gen;
END_RCPP
}
// nl_brkga
Rcpp::List nl_brkga(SEXP func_, arma::vec lowerLimit, arma::vec upperLimit, const unsigned p, const double pe, const double pm, const double rhoe, const unsigned K, const unsigned MAXT, const bool verbose, const long unsigned rngSeed, const unsigned X_INTVL, const unsigned X_NUMBER, const unsigned MAX_GENS);
RcppExport SEXP _brkga_nl_brkga(SEXP func_SEXP, SEXP lowerLimitSEXP, SEXP upperLimitSEXP, SEXP pSEXP, SEXP peSEXP, SEXP pmSEXP, SEXP rhoeSEXP, SEXP KSEXP, SEXP MAXTSEXP, SEXP verboseSEXP, SEXP rngSeedSEXP, SEXP X_INTVLSEXP, SEXP X_NUMBERSEXP, SEXP MAX_GENSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type func_(func_SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lowerLimit(lowerLimitSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type upperLimit(upperLimitSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type pe(peSEXP);
    Rcpp::traits::input_parameter< const double >::type pm(pmSEXP);
    Rcpp::traits::input_parameter< const double >::type rhoe(rhoeSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type K(KSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type MAXT(MAXTSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const long unsigned >::type rngSeed(rngSeedSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type X_INTVL(X_INTVLSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type X_NUMBER(X_NUMBERSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type MAX_GENS(MAX_GENSSEXP);
    rcpp_result_gen = Rcpp::wrap(nl_brkga(func_, lowerLimit, upperLimit, p, pe, pm, rhoe, K, MAXT, verbose, rngSeed, X_INTVL, X_NUMBER, MAX_GENS));
    return rcpp_result_gen;
END_RCPP
}
// Ackleys
SEXP Ackleys(const std::vector< double >& X);
RcppExport SEXP _brkga_Ackleys(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Ackleys(X));
    return rcpp_result_gen;
END_RCPP
}
// AluffiPentini
SEXP AluffiPentini(const std::vector< double >& X);
RcppExport SEXP _brkga_AluffiPentini(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(AluffiPentini(X));
    return rcpp_result_gen;
END_RCPP
}
// BeckerLago
SEXP BeckerLago(const std::vector< double >& X);
RcppExport SEXP _brkga_BeckerLago(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(BeckerLago(X));
    return rcpp_result_gen;
END_RCPP
}
// Bohachevsky1
SEXP Bohachevsky1(const std::vector< double >& X);
RcppExport SEXP _brkga_Bohachevsky1(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Bohachevsky1(X));
    return rcpp_result_gen;
END_RCPP
}
// Bohachevsky2
SEXP Bohachevsky2(const std::vector< double >& X);
RcppExport SEXP _brkga_Bohachevsky2(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Bohachevsky2(X));
    return rcpp_result_gen;
END_RCPP
}
// Branin
SEXP Branin(const std::vector< double >& X, const double a, const double b, const double c, const double d, const double e, const double f);
RcppExport SEXP _brkga_Branin(SEXP XSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP eSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type d(dSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(Branin(X, a, b, c, d, e, f));
    return rcpp_result_gen;
END_RCPP
}
// Camel3
SEXP Camel3(const std::vector< double >& X);
RcppExport SEXP _brkga_Camel3(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Camel3(X));
    return rcpp_result_gen;
END_RCPP
}
// Camel6
SEXP Camel6(const std::vector< double >& X);
RcppExport SEXP _brkga_Camel6(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Camel6(X));
    return rcpp_result_gen;
END_RCPP
}
// CosMixN
SEXP CosMixN(const std::vector< double >& X);
RcppExport SEXP _brkga_CosMixN(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(CosMixN(X));
    return rcpp_result_gen;
END_RCPP
}
// DekkersAarts
SEXP DekkersAarts(const std::vector< double >& X);
RcppExport SEXP _brkga_DekkersAarts(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(DekkersAarts(X));
    return rcpp_result_gen;
END_RCPP
}
// Easom
SEXP Easom(const std::vector< double >& X);
RcppExport SEXP _brkga_Easom(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Easom(X));
    return rcpp_result_gen;
END_RCPP
}
// EMichalewicz
SEXP EMichalewicz(const std::vector< double >& X, const double m);
RcppExport SEXP _brkga_EMichalewicz(SEXP XSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(EMichalewicz(X, m));
    return rcpp_result_gen;
END_RCPP
}
// Expo
SEXP Expo(const std::vector< double >& X);
RcppExport SEXP _brkga_Expo(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Expo(X));
    return rcpp_result_gen;
END_RCPP
}
// GoldPrice
SEXP GoldPrice(const std::vector< double >& X);
RcppExport SEXP _brkga_GoldPrice(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(GoldPrice(X));
    return rcpp_result_gen;
END_RCPP
}
// Griewank
SEXP Griewank(const std::vector< double >& X);
RcppExport SEXP _brkga_Griewank(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Griewank(X));
    return rcpp_result_gen;
END_RCPP
}
// Gulf
SEXP Gulf(const std::vector< double >& X);
RcppExport SEXP _brkga_Gulf(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Gulf(X));
    return rcpp_result_gen;
END_RCPP
}
// Hartman3
SEXP Hartman3(const std::vector< double >& X);
RcppExport SEXP _brkga_Hartman3(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Hartman3(X));
    return rcpp_result_gen;
END_RCPP
}
// Hartman6
SEXP Hartman6(const std::vector< double >& X);
RcppExport SEXP _brkga_Hartman6(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Hartman6(X));
    return rcpp_result_gen;
END_RCPP
}
// Hosaki
SEXP Hosaki(const std::vector< double >& X);
RcppExport SEXP _brkga_Hosaki(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Hosaki(X));
    return rcpp_result_gen;
END_RCPP
}
// Kowalik
SEXP Kowalik(const std::vector< double >& X);
RcppExport SEXP _brkga_Kowalik(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Kowalik(X));
    return rcpp_result_gen;
END_RCPP
}
// LM1
SEXP LM1(const std::vector< double >& X);
RcppExport SEXP _brkga_LM1(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(LM1(X));
    return rcpp_result_gen;
END_RCPP
}
// McCormic
SEXP McCormic(const std::vector< double >& X);
RcppExport SEXP _brkga_McCormic(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(McCormic(X));
    return rcpp_result_gen;
END_RCPP
}
// MeyerRoth
SEXP MeyerRoth(const std::vector< double >& X);
RcppExport SEXP _brkga_MeyerRoth(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MeyerRoth(X));
    return rcpp_result_gen;
END_RCPP
}
// MieleCantrell
SEXP MieleCantrell(const std::vector< double >& X);
RcppExport SEXP _brkga_MieleCantrell(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MieleCantrell(X));
    return rcpp_result_gen;
END_RCPP
}
// Modlangerman
SEXP Modlangerman(const std::vector< double >& X);
RcppExport SEXP _brkga_Modlangerman(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Modlangerman(X));
    return rcpp_result_gen;
END_RCPP
}
// ModRosenbrock
SEXP ModRosenbrock(const std::vector< double >& X);
RcppExport SEXP _brkga_ModRosenbrock(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(ModRosenbrock(X));
    return rcpp_result_gen;
END_RCPP
}
// MultiGauss
SEXP MultiGauss(const std::vector< double >& X);
RcppExport SEXP _brkga_MultiGauss(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MultiGauss(X));
    return rcpp_result_gen;
END_RCPP
}
// Neumaier2
SEXP Neumaier2(const std::vector< double >& X);
RcppExport SEXP _brkga_Neumaier2(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Neumaier2(X));
    return rcpp_result_gen;
END_RCPP
}
// Neumaier3
SEXP Neumaier3(const std::vector< double >& X);
RcppExport SEXP _brkga_Neumaier3(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Neumaier3(X));
    return rcpp_result_gen;
END_RCPP
}
// Paviani
SEXP Paviani(const std::vector< double >& X);
RcppExport SEXP _brkga_Paviani(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Paviani(X));
    return rcpp_result_gen;
END_RCPP
}
// Periodic
SEXP Periodic(const std::vector< double >& X);
RcppExport SEXP _brkga_Periodic(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Periodic(X));
    return rcpp_result_gen;
END_RCPP
}
// PowellQ
SEXP PowellQ(const std::vector< double >& X);
RcppExport SEXP _brkga_PowellQ(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(PowellQ(X));
    return rcpp_result_gen;
END_RCPP
}
// PriceTransistor
SEXP PriceTransistor(const std::vector< double >& X);
RcppExport SEXP _brkga_PriceTransistor(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(PriceTransistor(X));
    return rcpp_result_gen;
END_RCPP
}
// Rastrigin
SEXP Rastrigin(const std::vector< double >& X);
RcppExport SEXP _brkga_Rastrigin(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Rastrigin(X));
    return rcpp_result_gen;
END_RCPP
}
// Rosenbrock
SEXP Rosenbrock(const std::vector< double >& X);
RcppExport SEXP _brkga_Rosenbrock(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Rosenbrock(X));
    return rcpp_result_gen;
END_RCPP
}
// Salomon
SEXP Salomon(const std::vector< double >& X);
RcppExport SEXP _brkga_Salomon(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Salomon(X));
    return rcpp_result_gen;
END_RCPP
}
// Schaffer1
SEXP Schaffer1(const std::vector< double >& X);
RcppExport SEXP _brkga_Schaffer1(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Schaffer1(X));
    return rcpp_result_gen;
END_RCPP
}
// Schaffer2
SEXP Schaffer2(const std::vector< double >& X);
RcppExport SEXP _brkga_Schaffer2(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Schaffer2(X));
    return rcpp_result_gen;
END_RCPP
}
// Schubert
SEXP Schubert(const std::vector< double >& X);
RcppExport SEXP _brkga_Schubert(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Schubert(X));
    return rcpp_result_gen;
END_RCPP
}
// Schwefel
SEXP Schwefel(const std::vector< double >& X);
RcppExport SEXP _brkga_Schwefel(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Schwefel(X));
    return rcpp_result_gen;
END_RCPP
}
// Shekel10
SEXP Shekel10(const std::vector< double >& X);
RcppExport SEXP _brkga_Shekel10(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Shekel10(X));
    return rcpp_result_gen;
END_RCPP
}
// Shekel5
SEXP Shekel5(const std::vector< double >& X);
RcppExport SEXP _brkga_Shekel5(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Shekel5(X));
    return rcpp_result_gen;
END_RCPP
}
// Shekel7
SEXP Shekel7(const std::vector< double >& X);
RcppExport SEXP _brkga_Shekel7(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Shekel7(X));
    return rcpp_result_gen;
END_RCPP
}
// Shekelfox5
SEXP Shekelfox5(const std::vector< double >& X);
RcppExport SEXP _brkga_Shekelfox5(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Shekelfox5(X));
    return rcpp_result_gen;
END_RCPP
}
// Wood
SEXP Wood(const std::vector< double >& X);
RcppExport SEXP _brkga_Wood(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Wood(X));
    return rcpp_result_gen;
END_RCPP
}
// Zeldasine
SEXP Zeldasine(const std::vector< double >& X, const double A, const double B);
RcppExport SEXP _brkga_Zeldasine(SEXP XSEXP, SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type A(ASEXP);
    Rcpp::traits::input_parameter< const double >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(Zeldasine(X, A, B));
    return rcpp_result_gen;
END_RCPP
}
// reg_beta
SEXP reg_beta(const std::vector< double >& X, const arma::mat& data);
RcppExport SEXP _brkga_reg_beta(SEXP XSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< double >& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(reg_beta(X, data));
    return rcpp_result_gen;
END_RCPP
}
// reg_brkga
Rcpp::List reg_brkga(SEXP func_, arma::vec lowerLimit, arma::vec upperLimit, const unsigned p, const double pe, const double pm, const double rhoe, const unsigned K, const unsigned MAXT, const bool verbose, Rcpp::Nullable<Rcpp::NumericMatrix> data, const long unsigned rngSeed, const unsigned X_INTVL, const unsigned X_NUMBER, const unsigned MAX_GENS);
RcppExport SEXP _brkga_reg_brkga(SEXP func_SEXP, SEXP lowerLimitSEXP, SEXP upperLimitSEXP, SEXP pSEXP, SEXP peSEXP, SEXP pmSEXP, SEXP rhoeSEXP, SEXP KSEXP, SEXP MAXTSEXP, SEXP verboseSEXP, SEXP dataSEXP, SEXP rngSeedSEXP, SEXP X_INTVLSEXP, SEXP X_NUMBERSEXP, SEXP MAX_GENSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type func_(func_SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lowerLimit(lowerLimitSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type upperLimit(upperLimitSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type pe(peSEXP);
    Rcpp::traits::input_parameter< const double >::type pm(pmSEXP);
    Rcpp::traits::input_parameter< const double >::type rhoe(rhoeSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type K(KSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type MAXT(MAXTSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const long unsigned >::type rngSeed(rngSeedSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type X_INTVL(X_INTVLSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type X_NUMBER(X_NUMBERSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type MAX_GENS(MAX_GENSSEXP);
    rcpp_result_gen = Rcpp::wrap(reg_brkga(func_, lowerLimit, upperLimit, p, pe, pm, rhoe, K, MAXT, verbose, data, rngSeed, X_INTVL, X_NUMBER, MAX_GENS));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_brkga_tsp_brkga", (DL_FUNC) &_brkga_tsp_brkga, 1},
    {"_brkga_mdp_brkga", (DL_FUNC) &_brkga_mdp_brkga, 16},
    {"_brkga_nl_brkga", (DL_FUNC) &_brkga_nl_brkga, 14},
    {"_brkga_Ackleys", (DL_FUNC) &_brkga_Ackleys, 1},
    {"_brkga_AluffiPentini", (DL_FUNC) &_brkga_AluffiPentini, 1},
    {"_brkga_BeckerLago", (DL_FUNC) &_brkga_BeckerLago, 1},
    {"_brkga_Bohachevsky1", (DL_FUNC) &_brkga_Bohachevsky1, 1},
    {"_brkga_Bohachevsky2", (DL_FUNC) &_brkga_Bohachevsky2, 1},
    {"_brkga_Branin", (DL_FUNC) &_brkga_Branin, 7},
    {"_brkga_Camel3", (DL_FUNC) &_brkga_Camel3, 1},
    {"_brkga_Camel6", (DL_FUNC) &_brkga_Camel6, 1},
    {"_brkga_CosMixN", (DL_FUNC) &_brkga_CosMixN, 1},
    {"_brkga_DekkersAarts", (DL_FUNC) &_brkga_DekkersAarts, 1},
    {"_brkga_Easom", (DL_FUNC) &_brkga_Easom, 1},
    {"_brkga_EMichalewicz", (DL_FUNC) &_brkga_EMichalewicz, 2},
    {"_brkga_Expo", (DL_FUNC) &_brkga_Expo, 1},
    {"_brkga_GoldPrice", (DL_FUNC) &_brkga_GoldPrice, 1},
    {"_brkga_Griewank", (DL_FUNC) &_brkga_Griewank, 1},
    {"_brkga_Gulf", (DL_FUNC) &_brkga_Gulf, 1},
    {"_brkga_Hartman3", (DL_FUNC) &_brkga_Hartman3, 1},
    {"_brkga_Hartman6", (DL_FUNC) &_brkga_Hartman6, 1},
    {"_brkga_Hosaki", (DL_FUNC) &_brkga_Hosaki, 1},
    {"_brkga_Kowalik", (DL_FUNC) &_brkga_Kowalik, 1},
    {"_brkga_LM1", (DL_FUNC) &_brkga_LM1, 1},
    {"_brkga_McCormic", (DL_FUNC) &_brkga_McCormic, 1},
    {"_brkga_MeyerRoth", (DL_FUNC) &_brkga_MeyerRoth, 1},
    {"_brkga_MieleCantrell", (DL_FUNC) &_brkga_MieleCantrell, 1},
    {"_brkga_Modlangerman", (DL_FUNC) &_brkga_Modlangerman, 1},
    {"_brkga_ModRosenbrock", (DL_FUNC) &_brkga_ModRosenbrock, 1},
    {"_brkga_MultiGauss", (DL_FUNC) &_brkga_MultiGauss, 1},
    {"_brkga_Neumaier2", (DL_FUNC) &_brkga_Neumaier2, 1},
    {"_brkga_Neumaier3", (DL_FUNC) &_brkga_Neumaier3, 1},
    {"_brkga_Paviani", (DL_FUNC) &_brkga_Paviani, 1},
    {"_brkga_Periodic", (DL_FUNC) &_brkga_Periodic, 1},
    {"_brkga_PowellQ", (DL_FUNC) &_brkga_PowellQ, 1},
    {"_brkga_PriceTransistor", (DL_FUNC) &_brkga_PriceTransistor, 1},
    {"_brkga_Rastrigin", (DL_FUNC) &_brkga_Rastrigin, 1},
    {"_brkga_Rosenbrock", (DL_FUNC) &_brkga_Rosenbrock, 1},
    {"_brkga_Salomon", (DL_FUNC) &_brkga_Salomon, 1},
    {"_brkga_Schaffer1", (DL_FUNC) &_brkga_Schaffer1, 1},
    {"_brkga_Schaffer2", (DL_FUNC) &_brkga_Schaffer2, 1},
    {"_brkga_Schubert", (DL_FUNC) &_brkga_Schubert, 1},
    {"_brkga_Schwefel", (DL_FUNC) &_brkga_Schwefel, 1},
    {"_brkga_Shekel10", (DL_FUNC) &_brkga_Shekel10, 1},
    {"_brkga_Shekel5", (DL_FUNC) &_brkga_Shekel5, 1},
    {"_brkga_Shekel7", (DL_FUNC) &_brkga_Shekel7, 1},
    {"_brkga_Shekelfox5", (DL_FUNC) &_brkga_Shekelfox5, 1},
    {"_brkga_Wood", (DL_FUNC) &_brkga_Wood, 1},
    {"_brkga_Zeldasine", (DL_FUNC) &_brkga_Zeldasine, 3},
    {"_brkga_reg_beta", (DL_FUNC) &_brkga_reg_beta, 2},
    {"_brkga_reg_brkga", (DL_FUNC) &_brkga_reg_brkga, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_brkga(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}