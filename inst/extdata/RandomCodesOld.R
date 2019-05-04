

MDGa1n500m50 <- MDGa1n500m50[1:500,1:500]

write_rds(MDGa1n500m50, path = "inst/extdata/MDGa1n500m50.rds")

library(RcppXPtrUtils)

cpp <- "SEXP obj(const std::vector< double >& X) { 
        double myFitness;  
      	for(unsigned i = 0; i < X.size(); ++i) 
		        myFitness += (double(i + 1) * X[i]);
        return wrap(myFitness); 
        }"
func_cpp <- cppXPtr(cpp)
checkXPtr(func_cpp, "SEXP", c("const std::vector< double >&"))
brkga::nl_brkga(func_ = func_cpp, lowerLimit = rep(0, 10), upperLimit = rep(1, 10), K = 50)


cpp <- "SEXP Ackleys(const std::vector< double >& X){
  const int N = X.size();
  const arma::vec P = arma::conv_to< arma::vec >::from(X);
  const double subsum1 = arma::sum(P%P);
  const double subsum2 = arma::sum(arma::cos(2*arma::datum::pi*P));
  const double f = -20.*exp(-.2*sqrt(subsum1/(double)N))-exp(subsum2/(double)N)+20.+exp(1); 
  /* note that GAMS uses factor of -0.02 not 0.2  */
  return(Rcpp::wrap(f));}"

func_cpp <- cppXPtr(cpp, depends = "RcppArmadillo")
func_cpp <- cppXPtr(.beta_rega, depends = "RcppArmadillo")

checkXPtr(func_cpp, "SEXP", c("const std::vector< double >&"))
checkXPtr(func_cpp, "SEXP", c("const std::vector< double >&", "const arma::mat&"))

dados <- matrix(runif(100), ncol = 4)
brkga::nl_brkga(func_ = func_cpp, 
                lowerLimit = rep(-2, 5), 
                upperLimit = rep(2, 5), 
                K = 3, 
                p = 100, 
                data = dados, 
                X_INTVL = 5, 
                MAX_GENS = 100)

brkga::nl_brkga(func_ = func_cpp, lowerLimit = rep(-5, 5), upperLimit = rep(5, 5), K = 1, p = 1000)

c("Ackleys", "AluffiPentini", "BeckerLago",
  "Bohachevsky1", "Bohachevsky2", "Branin",
  "Camel3", "Camel6", "CosMix2", "CosMix4",
  "DekkersAarts", "Easom", "EMichalewicz",
  "Expo", "GoldPrice", "Griewank", "Gulf",
  "Hartman3", "Hartman6", "Hosaki", "Kowalik",
  "LM1", "LM2n10", "LM2n5", "McCormic",
  "MeyerRoth", "MieleCantrell", "Modlangerman",
  "ModRosenbrock", "MultiGauss", "Neumaier2",
  "Neumaier3", "Paviani", "Periodic",
  "PowellQ", "PriceTransistor", "Rastrigin",
  "Rosenbrock", "Salomon", "Schaffer1",
  "Schaffer2", "Schubert", "Schwefel",
  "Shekel10", "Shekel5", "Shekel7",
  "Shekelfox5", "Wood", "Zeldasine10",
  "Zeldasine20").