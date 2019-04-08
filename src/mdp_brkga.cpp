#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


#include <vector>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <time.h>

#include "MdpDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"

//#ifndef _OPENMP
//#define _OPENMP 1

std::vector< double >  dist;	// Vector of distances ROW*COL elements
unsigned m;


// Runs in \Theta(n \log n):
Rcpp::List localSearch(const std::vector< double >& chromosome){
  
  // Fitness of choromosome
  double myFitness = 0.0;
  
  std::vector< std::pair< double, unsigned > > ranking(chromosome.size());
  
  // ranking is the chromosome and vector of indices [0, 1, ..., n-1]
  for(unsigned i = 0; i < chromosome.size(); ++i) {
    ranking[i] = std::pair< double, unsigned >(chromosome[i], i);
  }
  
  // Here we sort 'permutation', which will then produce a permutation of [n] in pair::second:
  std::sort(ranking.begin(), ranking.end());
  
  // permutation[i].second is in {0, ..., n - 1}; a permutation can be obtained as follows
  std::vector< unsigned > permutation;
  for(auto i = ranking.begin();
      i != ranking.begin()+m; ++i) {
    permutation.push_back(i->second);
  }
  
  // Here we compute the fitness of chromosome by permutation obtained
  unsigned elem1, elem2;
  for (auto i = permutation.begin(); i != (permutation.end()); ++i) {
    elem1 = *i;
    for (auto j = i;j != permutation.end(); ++j) {
      elem2 = *j;
      myFitness += dist[elem1*chromosome.size()+elem2];
    }
  }
  //Rcpp::Rcout << "Best solution found has objective value = " << -algorithm.getBestFitness() << "em " << difftime(end_t,start_t)/CLOCKS_PER_SEC << " segundos" << std::endl;
  
  // Print the PDM solution
  /*
  for(std::list< unsigned >::const_iterator pt = permutation.begin(); pt != permutation.end(); ++pt) {
  std::cout << " " << *pt;
  }
  std::cout << std::endl;
  */
  // Execute the local search
  
  // Set N-M
  std::vector <unsigned> N_M; //Non selected elements of N
  std::vector <unsigned> M(permutation); //selected elements of N
  permutation.push_back(-1);
  
  for(unsigned i=0; i < chromosome.size(); i++){
    auto it = std::find(permutation.begin(), permutation.end(), i);
    if (it == permutation.end()) N_M.push_back(i);
  }
  
  /*
  for(std::list< unsigned >::const_iterator pt = N_M.begin(); pt != N_M.end(); ++pt) {
  std::cout << " " << *pt;
  }
  std::cout << std::endl;*/
  double delta_Z, max_delta_Z=0, LocalSearchFitness=myFitness;
  //std::list< unsigned >::const_iterator max_it_m, max_it_n_m;
  unsigned max_it_m, max_it_n_m;
  
  //First Improvement
  for(auto it_m = M.begin(); it_m != M.end(); ++it_m) {
    for(auto it_n_m = N_M.begin(); it_n_m != N_M.end(); ++it_n_m) {
      // calculando o delta z (vizinho - melhor_Solucao)
      delta_Z=0;
      Rcpp::checkUserInterrupt();
      for(auto item = M.begin(); item != M.end(); item++) {
        if (*item != *it_m) 
          delta_Z+=-dist[*it_m*chromosome.size()+*item]+dist[*it_n_m*chromosome.size()+*item];
      }
      if (delta_Z > 0) {
        //Rcpp::Rcout << "Atualização da solução : "<< delta_Z << " M "<< it_m - M.begin() <<  " N " << it_n_m - N_M.begin() << std::endl; 
        std::iter_swap(it_m, it_n_m);
        LocalSearchFitness+=delta_Z;
        //it_m = M.begin();
        //it_n_m = N_M.begin();
        std::prev(it_m, it_m - M.begin() - 1);
        std::prev(it_n_m, it_n_m - N_M.begin() - 1);
      } 
    }
  }
  //myFitness=myFitness+max_delta_Z;
  
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("Chromosome") = chromosome,
    Rcpp::Named("Solution") = M,
    Rcpp::Named("Fitness") = myFitness,
    Rcpp::Named("LocalSearchFitness") = LocalSearchFitness);
  
  return rst;
}
  
  
//' Execute a MDP solution search using BRKGA
//' @details Escrever um detalhamento bem legal e super bacana
//'  aqui.
//' @param instanceFile Path to file containing the instance 
//' @param MAX_TIME  Max time of execution in seconds
//' @param p  Size of population
//' @param pe Fraction of population to be the elite-set
//' @param pm Fraction of population to be replaced by mutants
//' @param rhoe Frobability that offspring inherit an allele from elite parent
//' @param K number of independent populations
//' @param MAXT number of threads for parallel decoding
//' @param X_INTVL Exchange best individuals at every 100 generations
//' @param X_NUMBER Exchange top 2 best
//' @param MAX_GENS Max number of generations
//' @param SEM_MELHORA  Pergunte para Geiza
//' @param verbose Level of output information
//' @param rngSeed Seed to the random number generator
//' @return A numeric vector of random values
//' @seealso \code{\link{api-usage}} and \url{https://github.com/milkway/brkga}
//' @export 
// [[Rcpp::export]]
Rcpp::List mdp_brkga(std::string instanceFile, 
                     const unsigned     MAX_TIME = 10,	 // run for 10 seconds
                     const unsigned            p = 500,	 // size of population
                     const double             pe = 0.20, // fraction of population to be the elite-set
                     const double             pm = 0.10, // fraction of population to be replaced by mutants
                     const double           rhoe = 0.70, // probability that offspring inherit an allele from elite parent
                     const unsigned            K = 1,		 // number of independent populations
                     const unsigned         MAXT = 8,    // number of threads for parallel decoding
                     const unsigned      X_INTVL = 100,	 // exchange best individuals at every 100 generations
                     const unsigned     X_NUMBER = 2,	   // exchange top 2 best
                     const unsigned     MAX_GENS = 1000, // run for 1000 gens
                     const unsigned  SEM_MELHORA = 0,
                     const bool          verbose = false, 
                     const long unsigned rngSeed = 0	   // seed to the random number generator
) {
  
  unsigned n;				// size of chromosomes
  MdpDecoder decoder;			// initialize the decoder
  MTRand rng(rngSeed);				// initialize the random number generator
  
  // reading input file with adjacency matrix
  // Formato: n n dist(i,j)
  std::ifstream myReadFile;
  myReadFile.open(instanceFile);
  if (!myReadFile) {
    Rcpp::Rcout << "Nao consigo abrir arquivo " << std::endl;
    return(-1);
  }
  
  // initialize the adjacency matrix
  double t;
  unsigned primeironumeroarquivo, segundonumeroarquivo;
  
  myReadFile >> primeironumeroarquivo;	
  myReadFile >> segundonumeroarquivo;
  
  n = primeironumeroarquivo;     // size of chromosomes - initialize n
  
  while (myReadFile >> t){
    dist.push_back (t);
  }
  
  // initialize the BRKGA-based heuristic
  BRKGA< MdpDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
  
  unsigned generation = 0;	// current generation
  const unsigned fitness = 0;
  m = n*0.10;			// número de elementos do subconjunto de máxima diversidade
  clock_t start_t, end_t;
  
  //verificar se esta estagnado, se estiver por X_INTVL iteracoes, reestart.
  start_t = clock();
  double BestLocalSearchFitness=0;
  std::list< unsigned > BestSolution;
  do {
    algorithm.evolve();	// evolve the population for one generation
    Rcpp::Rcout << "Best solution generation " << generation << " found has objective value = " << -algorithm.getBestFitness() << std::endl;
    
    // Local Search
    Rcpp::List search = localSearch(algorithm.getBestChromosome());
    double LocalSearchFitness = search["LocalSearchFitness"];
    std::list< unsigned > solution = search["Solution"];
    if (LocalSearchFitness > BestLocalSearchFitness) {
      BestLocalSearchFitness = LocalSearchFitness;
      BestSolution = solution;
    }

    Rcpp::Rcout << "Local Search  of generation :" << LocalSearchFitness << std::endl;
    Rcpp::Rcout << "--------------" << std::endl;
    
    
    if((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
    end_t = clock();
    //} while (generation < MAX_GENS);
  } while ((double)(difftime(end_t,start_t)/CLOCKS_PER_SEC) <= MAX_TIME);
  
  Rcpp::Rcout << "Best solution found has objective value = " << -algorithm.getBestFitness() << "em " << difftime(end_t,start_t)/CLOCKS_PER_SEC << " segundos" << std::endl;
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("PDMSolution") = BestSolution,
    Rcpp::Named("FitnessBRKGA") = algorithm.getBestFitness(),
    Rcpp::Named("FitnessLS") = BestLocalSearchFitness);
  
  return rst;
  }
