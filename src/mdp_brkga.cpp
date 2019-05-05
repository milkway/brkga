#include "RcppArmadillo.h"
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include <vector>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <time.h>
#include <functional>

#include "MdpDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"



// Used to sort ranking by second argument
bool compare(const std::pair<double, int>&i, const std::pair<double, int>&j)
{
  return i.second < j.second;
}

//Get M set
std::set< unsigned > get_M(const std::vector< double >& chromosome, unsigned m){
  std::vector< std::pair< double, unsigned > > ranking(chromosome.size());
  
  // ranking is the chromosome and vector of indices [0, 1, ..., n-1]
  for(unsigned i = 0; i < chromosome.size(); ++i) {
    ranking[i] = std::pair< double, unsigned >(chromosome[i], i);
  }
  
  // Here we sort 'permutation', which will then produce a permutation of [n] in pair::second:
  std::sort(ranking.begin(), ranking.end());
  std::set< unsigned> M;
  for(auto it = ranking.begin(); it != ranking.begin() + m; ++it) 
    M.insert((*it).second);
  return(M);
}

//Get N-M set
std::set< unsigned > get_N(unsigned n, std::set< unsigned > M){
  std::set< unsigned > N;
  for(unsigned i = 0; i < n; i++){
    auto it = M.find(i);
    if (it == M.end()) N.insert(i);
  }
  return(N);
}


// Runs in \Theta(n \log n):
Rcpp::List localSearch(const std::vector< double >& chromosome, double myFitness, const arma::mat& DistMat, unsigned m){
  
  std::vector< std::pair< double, unsigned > > ranking(chromosome.size());
  
  // ranking is the chromosome and vector of indices [0, 1, ..., n-1]
  for(unsigned i = 0; i < chromosome.size(); ++i) {
    ranking[i] = std::pair< double, unsigned >(chromosome[i], i);
  }
  
  // Here we sort 'permutation', which will then produce a permutation of [n] in pair::second:
  std::sort(ranking.begin(), ranking.end());
  std::set< unsigned> M;
  for(auto it = ranking.begin(); it != ranking.begin() + m; ++it) 
    M.insert((*it).second);
  std::sort(ranking.begin(),ranking.end(), compare);
  
  // Set N-M
  std::set< unsigned > N_M;
  for(unsigned i = 0; i < chromosome.size(); i++){
    auto it = M.find(i);
    if (it == M.end()) N_M.insert(i);
  }

  double LocalSearchFitness = myFitness;

  //First Improvement
  double delta_Z = 0;
  for(auto it_m = M.begin(); it_m != M.end(); ++it_m) {
    for(auto it_n_m = N_M.begin(); it_n_m != N_M.end(); ++it_n_m) {
      // calculando o delta z (vizinho - melhor_Solucao)
      delta_Z = 0;
      Rcpp::checkUserInterrupt();
      for(auto item = M.begin(); item != M.end(); item++) {
        if (*item != *it_m) 
          delta_Z += -DistMat(*it_m, *item) + DistMat(*it_n_m, *item);
      }
      if (delta_Z > 0) {
        M.insert(*it_n_m);
        double aux = ranking[*it_m].first;
        ranking[*it_m].first = ranking[*it_n_m].first;
        ranking[*it_n_m].first = aux;
        N_M.insert(*it_m);
        M.erase(*it_m);
        N_M.erase(*it_n_m);
        
        // PENSAR: Modificar o cromossomo com a sol da busca
        LocalSearchFitness -= delta_Z;
        it_m = M.begin();
        it_n_m = N_M.begin();
      } 
    }
  }
  
  std::vector<double> chrm(chromosome.size());
  for(auto it = ranking.begin(); it != ranking.end(); ++it) 
    chrm[it - ranking.begin()] = (*it).first;
  
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("Chromosome") = chrm,
    Rcpp::Named("Solution") = M,
    Rcpp::Named("lsFitness") = LocalSearchFitness);
  
  return rst;
}
  
  
//' Execute a MDP solution search using BRKGA
//' @details Escrever um detalhamento bem legal e super bacana
//'  aqui.
//' @param \code{DistanceMatrix} distances matrix of the instance 
//' @param \code{m} Number of elements selected (m <= n)
//' @param \code{MAX_TIME}  Max time of execution in seconds
//' @param \code{p}  Size of population
//' @param \code{pe} Fraction of population to be the elite-set
//' @param \code{pm} Fraction of population to be replaced by mutants
//' @param \code{rhoe} Frobability that offspring inherit an allele from elite parent
//' @param \code{lambda} probability update population with local search
//' @param \code{K} number of independent populations
//' @param \code{THREADS} number of threads for parallel decoding
//' @param \code{X_INTVL} Exchange best individuals at every 100 generations
//' @param \code{X_NUMBER} Exchange top 2 best
//' @param \code{MAX_GENS} Max number of generations
//' @param \code{verbose} Level of output information
//' @param \code{rngSeed} Seed to the random number generator
//' @return A numeric vector of random values
//' @seealso \code{\link{api-usage}} and \url{https://github.com/milkway/brkga}
//' @export 
// [[Rcpp::export]]
Rcpp::List mdp_brkga(const arma::mat   DistanceMatrix,  
                     const unsigned                 m,   // Number of elements selected (m <= n)
                     const unsigned       method =  0,   // 0: BRKGA; 1: BRKGA + LS; 2: BRKGA(LS) 
                     const unsigned  LS_INTERVAL = 10,  // Generations between local searches 
                     const unsigned     MAX_TIME = 10,	 // run for 10 seconds
                     const unsigned            p = 500,	 // size of population
                     const double             pe = 0.20, // fraction of population to be the elite-set
                     const double             pm = 0.10, // fraction of population to be replaced by mutants
                     const double           rhoe = 0.70, // probability that offspring inherit an allele from elite parent
                     const double         lambda = 0.01, // probability update population with local search
                     const unsigned            K = 3,		 // number of independent populations
                     const unsigned      THREADS = 8,    // number of threads for parallel decoding
                     const unsigned      X_INTVL = 100,	 // exchange best individuals at every 100 generations
                     const unsigned     X_NUMBER = 2,	   // exchange top 2 best
                     const unsigned     MAX_GENS = 100,  // Max generations without improvement 
                     const unsigned RESET_AFTER = 200, 
                     const bool          verbose = false, 
                     const long unsigned rngSeed = 0	   // seed to the random number generator
) {
  
  if (DistanceMatrix.n_cols != DistanceMatrix.n_rows) Rcpp::stop("Distance matrix must be square!");
  unsigned n = DistanceMatrix.n_cols;				// size of chromosomes
  
  MTRand rng(rngSeed);				// initialize the random number generator
 
  // Rand generator for Pop and Allele
  std::mt19937::result_type seed_K = rngSeed + 1;
  std::mt19937::result_type seed_n = rngSeed + 2;
  std::mt19937::result_type seed_u = rngSeed + 3;
  auto K_rand = std::bind(std::uniform_int_distribution<unsigned>(0,K-1),
                          std::mt19937(seed_K));
  auto n_rand = std::bind(std::uniform_int_distribution<unsigned>(0,n-1),
                          std::mt19937(seed_n));
  auto u_rand = std::bind(std::uniform_real_distribution<double>(0,1),
                             std::mt19937(seed_u));

  MdpDecoder decoder(DistanceMatrix, m);			// initialize the decoder

  // initialize the BRKGA-based heuristic
  BRKGA< MdpDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, THREADS);
  
  unsigned generation = 0;	// current generation
  //double lastBKfitness = 0;
  double lastLSfitness = 0;
  int gensLoosing = 0;
  
  //m = n*0.10;			// número de elementos do subconjunto de máxima diversidade
  clock_t start_t, end_t;
  
  //verificar se esta estagnado, se estiver por X_INTVL iteracoes, reestart.
  start_t = clock();
  double BestLocalSearchFitness = 0;
  unsigned relevantGeneration = 0;	// last relevant generation: best updated or reset called
  unsigned ImprovedSol = 0;
  std::set< unsigned > BestSolution;
  std::vector< double > BestChromosome(n);
  Rcpp::List search;
  Rprintf("+--------------+--------------+------------+---------+----------+\n");
  Rprintf("| Best Fitness | Local Search | Generation | Loosing | Duration |\n");
  double loopTime = 0;
  do {
    gensLoosing++;
    algorithm.evolve();	// evolve the population for one generation

    
    // New Local Search
    if (generation % LS_INTERVAL == 0){
      std::vector<double> chromosome = algorithm.getBestChromosome();
      double LocalSearchFitness = (double)algorithm.getBestFitness();
      std::vector< std::pair< double, unsigned > > ranking(n);
      // ranking is the chromosome and vector of indices [0, 1, ..., n-1]
      for(unsigned i = 0; i < n; ++i) {
        ranking[i] = std::pair< double, unsigned >(chromosome[i], i);
      }
      
      std::set< unsigned > M = get_M(algorithm.getBestChromosome(), m);
      std::set< unsigned > N = get_N(n, M);
      
      //Look for improvements and update a population
      //instance <- read_rds("inst/extdata/MDG.1.a.n500m50.rds")
      double delta_Z = 0;
      for(auto it_m = M.begin(); it_m != M.end(); ++it_m) {
        for(auto it_n_m = N.begin(); it_n_m != N.end(); ++it_n_m) {
          // calculando o delta z (vizinho - melhor_Solucao)
          delta_Z = 0;
          Rcpp::checkUserInterrupt();
          for(auto item = M.begin(); item != M.end(); item++) {
            if (*item != *it_m) 
              delta_Z += -DistanceMatrix(*it_m, *item) + DistanceMatrix(*it_n_m, *item);
          }
          if (delta_Z > 0) {
            M.insert(*it_n_m);
            double aux = chromosome[*it_m];
            chromosome[*it_m] = chromosome[*it_n_m];
            chromosome[*it_n_m] = aux;
            N.insert(*it_m);
            M.erase(*it_m);
            N.erase(*it_n_m);
            //Update Fitness
            LocalSearchFitness -= delta_Z;
            it_m = M.begin();
            it_n_m = N.begin();
            // Update best chromossome
            if (LocalSearchFitness < BestLocalSearchFitness){
              BestChromosome = chromosome;
              BestLocalSearchFitness = LocalSearchFitness;
              BestSolution = M;
              gensLoosing = 0;
            }
            //result <- mdp_brkga(DistanceMatrix = instance, m = 50, K = 3, MAX_TIME = 10, rngSeed = as.integer(Sys.time()))
            if (u_rand() < lambda) 
              algorithm.exchangeAlleles(chromosome, K_rand(), n_rand(), LocalSearchFitness);
            Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %8.0f |", \
                    algorithm.getBestFitness(),                   \
                    BestLocalSearchFitness,                       \
                    generation,                                   \
                    gensLoosing,                                  \
                    (double)(difftime(clock(),start_t)/CLOCKS_PER_SEC));
            Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
          } 
        }
      }
    }
    if (lastLSfitness < BestLocalSearchFitness){      
      ImprovedSol++;
    } 

    //lastBKfitness = algorithm.getBestFitness();
    lastLSfitness = BestLocalSearchFitness;
    
    //  Evolution strategy: restart
    if(generation - relevantGeneration > RESET_AFTER) {
      algorithm.reset();	// restart the algorithm with random keys
      relevantGeneration = generation;
      }
    
    if((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
    end_t = clock();
    loopTime = (double)(difftime(end_t,start_t)/CLOCKS_PER_SEC);
    Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %8.0f |", \
            algorithm.getBestFitness(),                   \
            BestLocalSearchFitness,                       \
            generation,                                   \
            gensLoosing,                                  \
            (double)(difftime(clock(),start_t)/CLOCKS_PER_SEC));
    Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
  } while ((loopTime <= MAX_TIME) && (gensLoosing <= MAX_GENS));
  Rprintf("\n\nThis is the end. The Doors.\n");
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("LSChr") = BestChromosome,
    Rcpp::Named("BKChr") = algorithm.getBestChromosome(),
    Rcpp::Named("Tour") = BestSolution,
    Rcpp::Named("LSFitness") = algorithm.getBestFitness(),
    Rcpp::Named("BKFitness") = BestLocalSearchFitness,
    Rcpp::Named("Generations Number") = generation,
    Rcpp::Named("Improvement Number") = ImprovedSol,
    Rcpp::Named("Duration") = loopTime);
  return rst;
  }
