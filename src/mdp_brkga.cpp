#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include <vector>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <time.h>
#include <chrono>
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
std::unordered_set< unsigned > get_M(const std::vector< double >& chromosome, unsigned m){
  std::vector< std::pair< double, unsigned > > ranking(chromosome.size());
  
  // ranking is the chromosome and vector of indices [0, 1, ..., n-1]
  for(unsigned i = 0; i < chromosome.size(); ++i) {
    ranking[i] = std::pair< double, unsigned >(chromosome[i], i);
  }
  
  // Here we sort 'permutation', which will then produce a permutation of [n] in pair::second:
  std::sort(ranking.begin(), ranking.end());
  std::unordered_set< unsigned> M;
  for(auto it = ranking.begin(); it != ranking.begin() + m; ++it) 
    M.insert((*it).second);
  return(M);
}

//Get N-M set
std::unordered_set< unsigned > get_N(unsigned n, std::unordered_set< unsigned > M){
  std::unordered_set< unsigned > N;
  for(unsigned i = 0; i < n; i++){
    auto it = M.find(i);
    if (it == M.end()) N.insert(i);
  }
  return(N);
}

///' Get Alpha vector
///' @details Get the distance from nodes to every node in M
///' @param \code{Tour} Set of nodes;
///' @param \code{Distances} Distance matrix
///' @export
/// [[Rcpp::export]]
arma::vec getAlphaVector(const std::unordered_set<unsigned>& Tour, const arma::mat& Distances){
  unsigned n = Distances.n_cols;
  arma::vec alpha = arma::zeros<arma::vec>(n);
  for(unsigned i = 0; i < n; i++){
    auto it = Tour.find(i);
    if (it != Tour.end())
      alpha[i] = 1.0;
  }
  return(Distances*alpha);
}

//' Get fitness from tour
//' @details Get fitness using the tour, m and distance matrix
//' @param \code{Tour} Set of tour's nodes.
//' @param \code{Distances} Distance matrix
//' @return A double value representing the chromosome fitness
//' @export 
// [[Rcpp::export]]
double getTourFitness(const std::vector<unsigned>& Tour, const arma::mat& Distances){
  double Fitness = 0;
  for (auto i = Tour.begin(); i != Tour.end(); ++i) {
    for (auto j = i; j != Tour.end(); ++j) {
      Fitness += Distances(*i, *j);
    }
  }
  return(Fitness);
}

//' Get fitness from chromosome
//' @details Get fitness using chromossome, m and distance matrix
//' @param \code{chromosome} Chromosome
//' @param \code{Distances} Distance matrix
//' @param \code{m} Size of tour
//' @return A double value representing the chromosome fitness
//' @export 
// [[Rcpp::export]]
double getChromosomeFitness(const std::vector< double >& chromosome, const arma::mat& Distances, unsigned m){
  std::unordered_set< unsigned > M = get_M(chromosome, m);
  std::vector<unsigned> Tour(M.begin(),M.end());
  double Fitness = getTourFitness(Tour, Distances);
  return(Fitness);
}
  
  
//' Execute a MDP solution search using BRKGA
//' @details Escrever um detalhamento bem legal e super bacana
//'  aqui.
//' @param \code{DistanceMatrix} distances matrix of the instance 
//' @param \code{m} Number of elements selected (m <= n)
//' @param \code{LS_INTVL} Generations between local searches
//' @param \code{GEN_INTVL} Generations in each evolution;
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
                     const unsigned     LS_INTVL =  5,  // Generations between local searches 
                     const unsigned    GEN_INTVL =  5,  // Generations in each evolution;                      
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
  auto n_rand = std::bind(std::uniform_int_distribution<unsigned>(0,p-1),
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
  //clock_t start_t, end_t;
  
  // Timing using chrono library
  auto start = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff;
  
  //verificar se esta estagnado, se estiver por X_INTVL iteracoes, reestart.
  double time_elapsed = 0;
  double BestLocalSearchFitness = 0;
  unsigned relevantGeneration = 0;	// last relevant generation: best updated or reset called
  unsigned ImprovedSol = 0;
  std::unordered_set< unsigned > BestSolution;
  std::vector< double > BestChromosome(n);
  
  Rprintf("+--------------+--------------+------------+---------+----------+\n");
  Rprintf("| Best Fitness | Local Search | Generation | Loosing | Duration |\n");
  double loopTime = 0;
  do {
    gensLoosing++;
    algorithm.evolve(GEN_INTVL);	// evolve the population for one generation

    // New Local Search
    if ((generation % LS_INTVL == 0) && (method == 1)){
      std::vector<double> chromosome(n);
      std::copy(algorithm.getBestChromosome().begin(), algorithm.getBestChromosome().end(), chromosome.begin());
      double LocalSearchFitness = (double)algorithm.getBestFitness();
      std::vector< std::pair< double, unsigned > > ranking(n);
      // ranking is the chromosome and vector of indices [0, 1, ..., n-1]
      for(unsigned i = 0; i < n; ++i) {
        ranking[i] = std::pair< double, unsigned >(chromosome[i], i);
      }
      
      std::unordered_set< unsigned > M = get_M(chromosome, m);
      std::unordered_set< unsigned > N = get_N(n, M);
      
      //Look for improvements and update a population
      double delta_Z = 0;
      for(auto it_m = M.begin(); it_m != M.end(); ++it_m) {
        for(auto it_n = N.begin(); it_n != N.end(); ++it_n) {
          // calculando o delta z (vizinho - melhor_Solucao)
          delta_Z = 0;
          Rcpp::checkUserInterrupt();
          for(auto item = M.begin(); item != M.end(); item++) {
            if (*item != *it_m) 
              delta_Z += -DistanceMatrix(*it_m, *item) + DistanceMatrix(*it_n, *item);
          }
          if (delta_Z > 0) {
            M.insert(*it_n);
            double aux = chromosome[*it_m];
            chromosome[*it_m] = chromosome[*it_n];
            chromosome[*it_n] = aux;
            N.insert(*it_m);
            M.erase(*it_m);
            N.erase(*it_n);
            //Update Fitness
            LocalSearchFitness -= delta_Z;
            it_m = M.begin();
            it_n = N.begin();
            // Update best chromossome
            if (LocalSearchFitness < BestLocalSearchFitness){
              BestChromosome = chromosome;
              BestLocalSearchFitness = LocalSearchFitness;
              BestSolution = M;
              gensLoosing = 0;
            }
            if (u_rand() < lambda){
              //Rprintf("\nK rand: %i, n rand: % i\n", K_rand(), n_rand());
              algorithm.exchangeAlleles(chromosome, K_rand(), n_rand(), LocalSearchFitness);
            } 
            auto time = std::chrono::steady_clock::now();
            diff = time - start;
            time_elapsed = std::chrono::duration <double, std::milli> (diff).count()/1000;
            Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |", \
                    -algorithm.getBestFitness(),                   \
                    -BestLocalSearchFitness,                       \
                    generation,                                   \
                    gensLoosing,                                  \
                    time_elapsed);
            Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
          } 
          if (time_elapsed > MAX_TIME) break;
          }
        if (time_elapsed > MAX_TIME) break;
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
    generation += GEN_INTVL; 
    if((generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
    auto time = std::chrono::steady_clock::now();
    diff = time - start;
    loopTime = std::chrono::duration <double, std::milli> (diff).count()/1000;
    Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |", \
            -algorithm.getBestFitness(),                   \
            -BestLocalSearchFitness,                       \
            generation,                                   \
            gensLoosing,                                  \
            loopTime);
    Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
  } while ((loopTime <= MAX_TIME) && (gensLoosing <= MAX_GENS));
  Rprintf("\n\nThis is the end. The Doors.\n");
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("LSChr") = BestChromosome,
    Rcpp::Named("BKChr") = algorithm.getBestChromosome(),
    Rcpp::Named("Tour") = BestSolution,
    Rcpp::Named("LSFitness") = -BestLocalSearchFitness,
    Rcpp::Named("BKFitness") = -algorithm.getBestFitness(),
    Rcpp::Named("Generations Number") = generation,
    Rcpp::Named("Improvement Number") = ImprovedSol,
    Rcpp::Named("Duration") = loopTime);
  return rst;
  }

///////////////////////////////////////////////////////////////////////////////////////////////


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
Rcpp::List mdp_brkgals(const arma::mat   DistanceMatrix,  
                     const unsigned                 m,   // Number of elements selected (m <= n)
                     const unsigned       method =  0,   // 0: BRKGA; 1: BRKGA + LS; 2: BRKGA(LS) 
                     const unsigned     LS_INTVL = 10,   // Generations between local searches 
                     const unsigned    GEN_INTVL = 5,    // Interval between Generations                      
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
                     const bool          verbose = true, 
                     const long unsigned rngSeed = 0	   // seed to the random number generator
) {

  if ( THREADS > 0 )
    omp_set_num_threads( THREADS );
  Rprintf("\nNumber of threads=%i\n", omp_get_max_threads());


  if (DistanceMatrix.n_cols != DistanceMatrix.n_rows) Rcpp::stop("Distance matrix must be square!");
  unsigned n = DistanceMatrix.n_cols;				// size of chromosomes
  
  MTRand rng(rngSeed);				// initialize the random number generator
  
  // Rand generator for Pop and Allele
  std::mt19937::result_type seed_K = rngSeed + 1;
  std::mt19937::result_type seed_n = rngSeed + 2;
  std::mt19937::result_type seed_u = rngSeed + 3;
  auto K_rand = std::bind(std::uniform_int_distribution<unsigned>(0,K-1),
                          std::mt19937(seed_K));
  auto n_rand = std::bind(std::uniform_int_distribution<unsigned>(0,p-1),
                          std::mt19937(seed_n));
  auto u_rand = std::bind(std::uniform_real_distribution<double>(0,1),
                          std::mt19937(seed_u));
  
  MdpDecoder decoder(DistanceMatrix, m);			// initialize the decoder
  
  // initialize the BRKGA-based heuristic
  BRKGA< MdpDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, THREADS);
  unsigned generation = 0;	// current generation
  double BestLSFitness = 0;
  arma::uvec LS_tour(m, arma::fill::zeros);
  double lastLSfitness = 0;
  int gensLoosing = 0;
  
  // Timing using chrono library
  auto start = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff;
  
  //verificar se esta estagnado, se estiver por X_INTVL iteracoes, reestart.
  unsigned relevantGeneration = 0;	// last relevant generation: best updated or reset called
  unsigned ImprovedSol = 0;
  if (verbose) {
    Rprintf("+--------------+--------------+------------+---------+----------+\n");
    Rprintf("| Best Fitness | Local Search | Generation | Loosing | Duration |\n");
  }
  

  double loopTime = 0;
  do {
    
    gensLoosing++;
    algorithm.evolve(GEN_INTVL);	// evolve the population for one generation
    
    if ((generation % LS_INTVL == 0) && (method == 1)){
      Progress p(1, false); // create the progress monitor
      arma::umat M(m, K, arma::fill::zeros);
      arma::umat BestTour(m, K, arma::fill::zeros);
      arma::umat N(n - m, K, arma::fill::zeros);
      arma::vec K_fitness(K, arma::fill::zeros);
      
      for (unsigned i_ = 0; i_ < K; i_++){
        std::vector<double> chromosome =  algorithm.getBestPopulationChromosome(i_);
        std::unordered_set< unsigned > M_ = get_M(chromosome, m);
        std::unordered_set< unsigned > N_ = get_N(n, M_);
        std::vector<unsigned> Tour(M_.begin(),M_.end());
        K_fitness(i_) = -getTourFitness(Tour, DistanceMatrix);
        auto it = M_.begin();
        for (unsigned j_ = 0; j_ < m; j_++){
          M(j_,i_) = *it;
          std::advance(it,1);
        }
        it = N_.begin();
        for (unsigned j_ = 0; j_ < (n - m); j_++){
          N(j_,i_) = *it;
          std::advance(it,1);
        }
      }
// candidate loop.
      #pragma omp parallel for schedule(dynamic)
      for(unsigned k_ = 0; k_ < K; k_++){
          double LocalSearchFitness = K_fitness(k_);
          double BestLocalSearchFitness = 0;
          
          //Look for improvements and update a population
          double delta = 0;
          for(unsigned m_ = 0; m_ < m; m_++) {
            for(unsigned n_ = 0; n_ < (n-m); n_++) {
              // calculando o delta z (vizinho - melhor_Solucao)
              delta = 0;
              for(unsigned item = 0; item < m; item++) {
                if (M(item, k_) != M(m_, k_)) 
                  delta += -DistanceMatrix(M(m_, k_), M(item, k_)) + DistanceMatrix(N(n_,k_), M(item, k_));
              }
              if (delta > 0) {
                unsigned aux = M(m_,k_);
                M(m_,k_) = N(n_,k_);
                N(n_,k_) = aux;
                //Update Fitness
                LocalSearchFitness -= delta;
                m_ = 0;
                n_ = 0;
                // Update best
                if (LocalSearchFitness < BestLocalSearchFitness){
                  BestLocalSearchFitness = LocalSearchFitness;
                  BestTour.col(k_) = M.col(k_);
                }
              }
              auto time = std::chrono::steady_clock::now();
              diff = time - start;
              if (std::chrono::duration <double, std::milli> (diff).count()/1000 > MAX_TIME)
                break; // Get out of here!
              if (p.is_aborted()) break; // Get out of here!
            }
            if (std::chrono::duration <double, std::milli> (diff).count()/1000 > MAX_TIME)
              break; // Get out of here!
            if (p.is_aborted()) break; // Get out of here!
          }
          K_fitness(k_) = BestLocalSearchFitness;          
      }
      //Rcpp::checkUserInterrupt();
      
       for(unsigned kk = 0; kk < K; kk++){
        //Rprintf("\nBestLocalSearchFitness %f\n\n", K_fitness(kk));
        if(K_fitness(kk) < BestLSFitness){
          BestLSFitness = K_fitness(kk);
          LS_tour = BestTour.col(kk);
        }
      }
       p.cleanup();   
    }

    // for(unsigned i = 0; i < m; i++ )
    //   Rcpp::Rcout << " " << LS_tour(i);
    // Rcpp::Rcout << std::endl << std::endl;
    
    if (lastLSfitness < BestLSFitness){
      ImprovedSol++;
    } 
    
    //lastBKfitness = algorithm.getBestFitness();
    lastLSfitness = BestLSFitness;
    
    //  Evolution strategy: restart
    if(generation - relevantGeneration > RESET_AFTER) {
      algorithm.reset();	// restart the algorithm with random keys
      relevantGeneration = generation;
    }
    
    if((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
    
    // end of timing
    auto end = std::chrono::steady_clock::now();
    // Store the time difference between start and end
    diff = end - start;
    loopTime = std::chrono::duration <double, std::milli> (diff).count()/1000;
    if (verbose) {
      Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |", \
              -algorithm.getBestFitness(),                   \
              -BestLSFitness,                                \
              generation,                                    \
              gensLoosing,                                   \
              loopTime);
      Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
    }
  } while ((loopTime <= MAX_TIME) && (gensLoosing <= MAX_GENS));
  Rprintf("\n\nThis is the end. The Doors.\n");
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("Tour") = LS_tour,
    Rcpp::Named("LSFitness") = -BestLSFitness,
    Rcpp::Named("BKFitness") = -algorithm.getBestFitness(),
    Rcpp::Named("Generations Number") = generation,
    Rcpp::Named("Improvement Number") = ImprovedSol,
    Rcpp::Named("Duration") = loopTime);
  return rst;
}

/////////////////////////////////

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
Rcpp::List mdp_brkgaus(const arma::mat   DistanceMatrix,  
                       const unsigned                 m,   // Number of elements selected (m <= n)
                       const unsigned       method =  0,   // 0: BRKGA; 1: BRKGA + LS; 2: BRKGA(LS) 
                       const unsigned     LS_INTVL = 10,   // Generations between local searches 
                       const unsigned    GEN_INTVL = 5,    // Interval between Generations                      
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
                       const bool          verbose = true, 
                       const long unsigned rngSeed = 0	   // seed to the random number generator
) {
  
  if ( THREADS > 0 )
    omp_set_num_threads( THREADS );
  Rprintf("\nNumber of threads=%i\n", omp_get_max_threads());
  
  if (DistanceMatrix.n_cols != DistanceMatrix.n_rows) Rcpp::stop("Distance matrix must be square!");
  unsigned n = DistanceMatrix.n_cols;				// size of chromosomes
  
  MTRand rng(rngSeed);				// initialize the random number generator
  MdpDecoder decoder(DistanceMatrix, m);			// initialize the decoder
  
  // initialize the BRKGA-based heuristic
  BRKGA< MdpDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, THREADS);
  std::unordered_set< unsigned > BestTour;
  double Best_LS_Fitness = 0;
  double last_LS_fitness = 0;
  double Best_BK_Fitness = 0;
  double Best_Fitness = 0;
  double last_BK_fitness = 0;  
  // Timing using chrono library
  auto start = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff;
  double time_elapsed = 0;
  //verificar se esta estagnado, se estiver por X_INTVL iteracoes, reestart.
  unsigned generation = 0;	// current generation
  unsigned bestGeneration = 0;
  unsigned gensLoosing = 0;
  unsigned relevantGeneration = 0;	// last relevant generation: best updated or reset called
  unsigned ImprovedSol = 0;
  // if (verbose) {
  //   Rprintf("+--------------+--------------+------------+---------+----------+\n");
  //   Rprintf("| Best Fitness | Local Search | Generation |  Best   | Duration |\n");
  // }
  
  
  double loopTime = 0;
  do {
    
    gensLoosing++;
    algorithm.evolve(GEN_INTVL);	// evolve the population for one generation
    
    if ((generation % LS_INTVL == 0) && (method == 1)){
      Progress p(1, false); // create the progress monitor
      std::list< std::unordered_set<unsigned> > tours;
      std::unordered_set< unsigned > BestLSTour;
      #pragma omp parallel for schedule(dynamic)
      for(unsigned k_ = 0; k_ < K; k_++){
        std::vector<double> chromosome = algorithm.getBestPopulationChromosome(k_);
        std::unordered_set< unsigned > M = get_M(chromosome, m);
        std::unordered_set< unsigned > N = get_N(n, M);
        
        //Look for improvements and update a population
        arma::vec alpha  = getAlphaVector(M, DistanceMatrix);
        for(auto it_m = M.begin(); it_m != M.end(); ++it_m) {
          //p.increment();
          for(auto it_n = N.begin(); it_n != N.end(); ++it_n) {
            // calculando o delta z (vizinho - melhor_Solucao)
            double delta = alpha(*it_n) - alpha(*it_m) - DistanceMatrix(*it_m,*it_n);
            // double delta = 0;
            // for(auto item = M.begin(); item != M.end(); item++) {
            //   if (*item != *it_m) 
            //     delta += -DistanceMatrix(*it_m, *item) + DistanceMatrix(*it_n, *item);
            // }
            if (delta > 0) {
              alpha = alpha - DistanceMatrix.col(*it_m) + DistanceMatrix.col(*it_n);
              M.insert(*it_n);
              N.insert(*it_m);
              M.erase(*it_m);
              N.erase(*it_n);
              it_m = M.begin();
              it_n = N.begin();
              auto time = std::chrono::steady_clock::now();
              diff = time - start;
              time_elapsed = std::chrono::duration <double, std::milli> (diff).count()/1000;
            } 
            if (p.is_aborted()) break; // Get out of here!
            if (time_elapsed > MAX_TIME) break; // Get out of here!
          }
          if (p.is_aborted()) break; // Get out of here!
          if (time_elapsed > MAX_TIME) break; // Get out of here!
        }
        tours.push_front(M);
      }

      for(auto it = tours.begin(); it != tours.end(); it++){
        std::vector<unsigned> Tour((*it).begin(),(*it).end());
        double fitness = -getTourFitness(Tour, DistanceMatrix);
        if(fitness < Best_LS_Fitness){
          Best_LS_Fitness = fitness;
          BestTour = *it;
        }
      }
      p.cleanup();   
    }
    
    if (algorithm.getBestFitness() < Best_BK_Fitness) {
      Best_BK_Fitness = algorithm.getBestFitness();
      if (Best_BK_Fitness < Best_Fitness){
        bestGeneration = generation;
        Best_Fitness = Best_BK_Fitness;
      } 
    }

    if ((last_LS_fitness > Best_LS_Fitness) || (last_BK_fitness > Best_BK_Fitness)){
      ImprovedSol++;
      if (Best_LS_Fitness < Best_Fitness) {
        Best_Fitness = Best_LS_Fitness;
        bestGeneration = generation;
      }
    } 
    
    last_BK_fitness = Best_BK_Fitness;
    last_LS_fitness = Best_LS_Fitness;
    
    //  Evolution strategy: restart
    if(generation - relevantGeneration > RESET_AFTER) {
      algorithm.reset();	// restart the algorithm with random keys
      relevantGeneration = generation;
    }
    
    if((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
    
    // end of timing
    auto end = std::chrono::steady_clock::now();
    // Store the time difference between start and end
    diff = end - start;
    loopTime = std::chrono::duration <double, std::milli> (diff).count()/1000;
    
    // if (verbose) {
    //   Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
    //   Rprintf("| Best Fitness | Local Search | Generation |  Best   | Duration |\n");
    //   Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");      
    //   Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |", \
    //           -Best_BK_Fitness,                   \
    //           -Best_LS_Fitness,                                \
    //           generation,                                    \
    //           bestGeneration,                                   \
    //           loopTime);
    //   Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
    // }
  } while ((loopTime <= MAX_TIME) && (gensLoosing <= MAX_GENS));
  
  if (verbose) {
    Rprintf("\n+--------------+--------------+------------+---------+----------+");
    Rprintf("\n| Best Fitness | Local Search | Generation |  Best   | Duration |");
    Rprintf("\n+--------------+--------------+------------+---------+----------+");      
    Rprintf("\n| %12.2f | %12.2f | %10i | %7i | %7.1fs |",    \
            -Best_BK_Fitness,                                 \
            -Best_LS_Fitness,                                 \
            generation,                                       \
            bestGeneration,                                   \
            loopTime);
    Rprintf("\n+--------------+--------------+------------+---------+----------+");
  }
  
  Rprintf("\n\nThis is the end. The Doors.\n");
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("Tour") = BestTour,
    Rcpp::Named("Teste") =  getAlphaVector(BestTour, DistanceMatrix),    
    Rcpp::Named("LSFitness") = -Best_LS_Fitness,
    Rcpp::Named("BKFitness") = -Best_BK_Fitness,
    Rcpp::Named("Generations Number") = generation,
    Rcpp::Named("Best Generation") = bestGeneration,
    Rcpp::Named("Improvement Number") = ImprovedSol,
    Rcpp::Named("Duration") = loopTime);
  return rst;
}

