


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
  "Zeldasine20")



library(tidyverse)
result <- readr::read_rds("~/Dropbox/resultado.rds")
result %>%  filter(LSEr >= 0)
tour <- result %>% filter(Seed == 2651) %>% select(Tour) %>% unnest() %>% unlist(use.names = FALSE)

write.table(result$Tour[3] %>% unlist(),  sep = " ", file = "~/Dropbox/Trabalhos/pdm/PDM/calcularBeneficio/tour_mdg_a2.txt", row.names = FALSE, col.names = FALSE)

write.table(M %>% unlist(),  sep = " ", file = "~/Dropbox/Trabalhos/pdm/PDM/calcularBeneficio/tour_mdg_a3.txt", row.names = FALSE, col.names = FALSE)

write.table(res$PDMSolution,  sep = " ", file = "~/Dropbox/Trabalhos/pdm/PDM/calcularBeneficio/tour_mdg_a1.txt", row.names = FALSE, col.names = FALSE)

write.table(M,  sep = " ", file = "~/Dropbox/Trabalhos/pdm/PDM/calcularBeneficio/tour_mdg_a3.txt", row.names = FALSE, col.names = FALSE)


textfile <- system.file("extdata", package = "brkga",  "conv_MDG-a_5_n500_m50.txt")
teste <- read.table(file = textfile, header = FALSE, sep = " ", skip = 1L) %>% filter(row_number() <= 500) %>% select(1:500)

for(i in 1:nrow(lista)){
  cat("\nProcessing file ", i, " of ", nrow(lista))
  textfile <- paste0("~/Downloads/MDG-c/", lista$Name[i]) #system.file("~/Downloads/MDG-c", lista$Name[i])
  write_rds(path = paste0("inst/extdata/", lista$Type[i], ".", i, ".", lista$SubType[i], ".n", lista$n[i], "m", lista$m[i]), read.table(file = textfile, header = FALSE, sep = " ", skip = 1L) %>% filter(row_number() <= lista$n[i]) %>% select(1:lista$n[i]) %>% as.matrix())
}

#assign(paste0(".mdg", ".a.n", 500, "m", 50), read.table(file = textfile, header = FALSE, sep = " ", skip = 1L) %>% filter(row_number() <= 500) %>% select(1:500) %>% as.matrix())
#dist_matrix <- read_rds(paste0("inst/extdata/", lista$Type[i], ".", i, ".", lista$SubType[i], ".n", lista$n[i], "m", lista$m[i])) 
library(brkga)
lista <- read_rds(system.file("extdata", package = "brkga",  "mdplib.rds")) %>% 
  mutate(Name = paste0("conv_", Instance), 
         Type = str_sub(Name, 6,8), 
         SubType = str_sub(Name, 10,10), 
         n = as.integer(str_remove(str_extract(Name,pattern = "n\\d+"), pattern = "n")),
         m = as.integer(str_remove(str_extract(Name,pattern = "m\\d+"), pattern = "m"))) %>% 
  filter(Type == "MDG", SubType == "a")

index = 3
dist_matrix <- read_rds(system.file("extdata", package = "brkga",  paste0(lista$Type[index], ".", index, ".", lista$SubType[index], ".n", lista$n[index], "m", lista$m[index])))

res <- brkga::mdp_brkga(DistanceMatrix = dist_matrix,
                        m = lista$m[index],
                        MAX_TIME = 10, 
                        p = 150, 
                        pe = .2, 
                        pm = .2,
                        rhoe = .75,
                        MAXT=8, 
                        K=3, 
                        MAX_GENS = 1500,
                        RESET_AFTER = 200, 
                        rngSeed = 265)
