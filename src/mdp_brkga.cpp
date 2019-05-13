#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]
// // [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
// #include <progress_bar.hpp>

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


//' Get M 
//' @details Get the distance from nodes to every node in M
//' @param \code{Chromosome} Chromosome;
//' @param \code{m} Tour size
//' @export
// [[Rcpp::export(getM)]]
arma::uvec Arma_get_M(const arma::vec& chromosome, unsigned m){
  arma::uvec M = arma::sort_index(chromosome);
  return(M.head(m));
}




//' Get N 
//' @details Get the distance from nodes to every node in M
//' @param \code{M}
//' @param \code{n}
//' @export
// [[Rcpp::export("getN")]]
arma::uvec Arma_get_N(const arma::vec& chromosome, unsigned m){
  arma::uvec M = arma::sort_index(chromosome);
  unsigned n = chromosome.size();
  return(M.tail(n-m));
}


//' Get Alpha vector
//' @details Get the distance from nodes to every node in M
//' @param \code{Tour} Set of nodes;
//' @param \code{Distances} Distance matrix
//' @export
// [[Rcpp::export("getAlphaVector")]]
arma::vec Arma_getAlphaVector(const arma::uvec& Tour, const arma::mat& Distances){
  unsigned n = Distances.n_cols;
  unsigned m = Tour.size();
  arma::vec alpha = arma::zeros<arma::vec>(n);
  for(unsigned i = 0; i < m; i++){
      alpha(Tour(i)) = 1.0;
  }
  return(Distances*alpha/2);
}



//' Get Binet tour
//' @details Get the binary representation of the tour
//' @param \code{Tour} Set of nodes;
//' @param \code{n} Number of nodes
//' @export
// [[Rcpp::export("getTourBinary")]]
arma::uvec Arma_getBinaryTour(const arma::uvec& Tour, unsigned n){
  arma::uvec alpha = arma::zeros<arma::uvec>(n);
  unsigned m = Tour.size();
  for(unsigned i = 0; i < m; i++){
    alpha(Tour(i)) = 1;
  }
  return(alpha);
}


//' Get fitness from tour
//' @details Get fitness using the tour, m and distance matrix
//' @param \code{Tour} Set of tour's nodes.
//' @param \code{Distances} Distance matrix
//' @return A double value representing the chromosome fitness
//' @export 
// [[Rcpp::export("getTourFitness")]]
double Arma_getTourFitness(const arma::uvec& Tour, const arma::mat& Distances){
  double Fitness = 0;
  unsigned m = Tour.size();
  for (unsigned i = 0; i < m; i++) {
    for (auto j = i; j <  m; j++) {
      Fitness += Distances(Tour(i), Tour(j));
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
// [[Rcpp::export("getChromosomeFitness")]]
double Arma_getChromosomeFitness(const arma::vec& chromosome, const arma::mat& Distances, unsigned m){
  arma::uvec Ranking = arma::sort_index(chromosome);
  arma::uvec M = Ranking.head(m);
  double Fitness = Arma_getTourFitness(M, Distances);
  return(Fitness);
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
// [[Rcpp::export("mdp_brkga")]]
Rcpp::List mdp_brkgaArma(const arma::mat   DistanceMatrix,  
                      const unsigned                 m,   // Number of elements selected (m <= n)
                      const unsigned     LS_INTVL = 10,   // Generations between local searches 
                      const unsigned    GEN_INTVL = 5,    // Interval between Generations                      
                      const unsigned     MAX_TIME = 10,	 // run for 10 seconds
                      const unsigned            p = 500,	 // size of population
                      const double             pe = 0.20, // fraction of population to be the elite-set
                      const double             pm = 0.10, // fraction of population to be replaced by mutants
                      const double           rhoe = 0.70, // probability that offspring inherit an allele from elite parent
                      const unsigned            K = 4,		 // number of independent populations
                      const unsigned      THREADS = 8,    // number of threads for parallel decoding
                      const unsigned      X_INTVL = 100,	 // exchange best individuals at every 100 generations
                      const unsigned     X_NUMBER = 2,	   // exchange top 2 best
                      const unsigned     MAX_GENS = 100,  // Max generations without improvement 
                      const unsigned RESET_AFTER = 200, 
                      const unsigned     verbose = 2, 
                      const long unsigned rngSeed = 0	   // seed to the random number generator
) {
  
  if (DistanceMatrix.n_cols != DistanceMatrix.n_rows) Rcpp::stop("Distance matrix must be square!");
  unsigned n = DistanceMatrix.n_cols;				// size of chromosomes
  
  MTRand rng(rngSeed);				// initialize the random number generator
  MdpDecoder decoder(DistanceMatrix, m);			// initialize the decoder
  
  // initialize the BRKGA-based heuristic
  BRKGA< MdpDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, THREADS);
  //std::unordered_set< unsigned > BestTour;
  arma::uvec BestTour(m, arma::fill::zeros); 
  double Best_LS_Fitness = -arma::datum::eps;
  double last_LS_fitness = -arma::datum::eps;
  double Best_BK_Fitness = -arma::datum::eps;
  double Best_Fitness = -arma::datum::eps;
  double last_BK_fitness = -arma::datum::eps;  
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
  if (verbose >= 1) {
    Rprintf("+--------------+--------------+------------+---------+----------+\n");
    Rprintf("| Best BRKGA   | Local Search | Generation | Bst Gen | Duration |\n");
    Rprintf("+--------------+--------------+------------+---------+----------+\n");      
  }
  
  double loopTime = 0;
  double delta = 0;
  //Progress mp(1, false); // create the progress monitor
  do {
    gensLoosing++;
    algorithm.evolve(GEN_INTVL);	// evolve the population for one generation
    
    if (algorithm.getBestFitness() < Best_BK_Fitness) {
      Best_BK_Fitness = algorithm.getBestFitness();
      if (Best_BK_Fitness < Best_Fitness){
        bestGeneration = generation;
        Best_Fitness = Best_BK_Fitness;
        gensLoosing = 0;
      } 
    }
    
    if (verbose == 2) {
      Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |",    \
              -Best_BK_Fitness,                                 \
              -Best_LS_Fitness,                                 \
              generation,                                       \
              bestGeneration,                                   \
              loopTime);
      Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
    }
    
    
    if ((generation % LS_INTVL == 0)){
      arma::umat tours(m, K, arma::fill::zeros);
#pragma omp parallel for  schedule(static)    
      for(unsigned k = 0; k < K; k++){
        std::vector<double> chromosome = algorithm.getBestPopulationChromosome(k);
        arma::uvec Ranking = arma::sort_index(arma::conv_to< arma::vec >::from(chromosome));
        arma::uvec M = Ranking.head(m);
        arma::uvec N = Ranking.tail(n-m);
        //Look for improvements and update a population
        arma::vec alpha  = Arma_getAlphaVector(M, DistanceMatrix);
        bool flag = true;
        while (flag) {
          //Rcpp::checkUserInterrupt();
          M = arma::shuffle(M);
          N = arma::shuffle(N);
          flag = false;
          for(unsigned i = 0; i < m; i++) {
            for(unsigned j = 0; j < (n-m); j++) {
              // calculando o delta z (vizinho - melhor_Solucao)
              delta = alpha(N(j)) - alpha(M(i)) - DistanceMatrix(M(i),N(j));
              if (std::isgreaterequal(delta, 0.0001)) {
                flag = true;
                alpha = alpha - DistanceMatrix.col(M(i)) + DistanceMatrix.col(N(j));
                unsigned aux = M[i];
                M[i] = N[j];
                N[j] = aux;
              } 
            }
            auto time = std::chrono::steady_clock::now();
            diff = time - start;
            time_elapsed = std::chrono::duration <double, std::milli> (diff).count()/1000;
            // if (verbose == 2) {
            //   Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |",    \
            //           -Best_BK_Fitness,                                 \
            //           -Best_LS_Fitness,                                 \
            //           generation,                                       \
            //           bestGeneration,                                   \
            //           time_elapsed);
            //   Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
            // }
          } 
          if (time_elapsed > MAX_TIME) break; // Get out of here!
        }
        // double fitness = -Arma_getTourFitness(M, DistanceMatrix);
        // if(fitness < Best_LS_Fitness){
        //   Best_LS_Fitness = fitness;
        //   BestTour = M;
        // }
        tours.col(k) = M;
      }
      // End of parallel loop
      
      for(unsigned k = 0; k < K; k++){
        double fitness = -Arma_getTourFitness(tours.col(k), DistanceMatrix);
        if(fitness < Best_LS_Fitness){
          Best_LS_Fitness = fitness;
          BestTour = tours.col(k);
        }
      }
    }
    
    
    
    if ((last_LS_fitness > Best_LS_Fitness) || (last_BK_fitness > Best_BK_Fitness)){
      ImprovedSol++;
      if (Best_LS_Fitness < Best_Fitness) {
        Best_Fitness = Best_LS_Fitness;
        bestGeneration = generation;
        gensLoosing = 0;
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
  } while ((loopTime <= MAX_TIME) && (gensLoosing <= MAX_GENS));
  if (verbose >= 1 ) {
    Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |",    \
            -Best_BK_Fitness,                                 \
            -Best_LS_Fitness,                                 \
            generation,                                       \
            bestGeneration,                                   \
            loopTime);
    Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
  }
  
  if (verbose) Rprintf("\n\nThis is the end. The Doors.\n");
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("Tour") = BestTour,
    Rcpp::Named("LSFitness") = -Best_LS_Fitness,
    Rcpp::Named("BKFitness") = -Best_BK_Fitness,
    Rcpp::Named("Generations Number") = generation,
    Rcpp::Named("Best Generation") = bestGeneration,
    Rcpp::Named("Improvement Number") = ImprovedSol,
    Rcpp::Named("Duration") = loopTime);
  return rst;
}


/////////////////////////////////////////////////////////////////



///' Get fitness from tour
///' @details Get fitness using the tour, m and distance matrix
///' @param \code{Tour} Set of tour's nodes.
///' @param \code{Distances} Distance matrix
///' @return A double value representing the chromosome fitness
///' @export 
/// [[Rcpp::export]]
double getTourFitness(const std::unordered_set<unsigned>& Tour, const arma::mat& Distances){
  double Fitness = 0;
  for (auto i = Tour.begin(); i != Tour.end(); ++i) {
    for (auto j = i; j != Tour.end(); ++j) {
      Fitness += Distances(*i, *j);
    }
  }
  return(Fitness);
}

//////////////////////////////////////////////////////////

///' Get Binet tour
///' @details Get the binary representation of the tour
///' @param \code{Tour} Set of nodes;
///' @param \code{n} Number of nodes
///' @export
/// [[Rcpp::export]]
arma::uvec getBinaryTour(const std::unordered_set<unsigned>& Tour, unsigned n){
  arma::uvec alpha = arma::zeros<arma::uvec>(n);
  for(unsigned i = 0; i < n; i++){
    auto it = Tour.find(i);
    if (it != Tour.end())
      alpha(i) = 1;
  }
  return(alpha);
}
////////////////////////////////////////////////////////


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
      alpha(i) = 1.0;
  }
  return(Distances*alpha/2);
}


//////////////////////////////////////////////////



//Get N-M set
std::unordered_set< unsigned > get_N(unsigned n, std::unordered_set< unsigned > M){
  std::unordered_set< unsigned > N;
  for(unsigned i = 0; i < n; i++){
    auto it = M.find(i);
    if (it == M.end()) N.insert(i);
  }
  return(N);
}

///////////////////////////

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

/////////////////////////////////////////////////////


// Used to sort ranking by second argument
bool compare(const std::pair<double, int>&i, const std::pair<double, int>&j)
{
  return i.second < j.second;
}


//////////////////////////////////////////////////////////////

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
Rcpp::List mdp_brkga_uset(const arma::mat   DistanceMatrix,  
                       const unsigned                 m,   // Number of elements selected (m <= n)
                       const unsigned     LS_INTVL = 10,   // Generations between local searches 
                       const unsigned    GEN_INTVL = 5,    // Interval between Generations                      
                       const unsigned     MAX_TIME = 10,	 // run for 10 seconds
                       const unsigned            p = 500,	 // size of population
                       const double             pe = 0.20, // fraction of population to be the elite-set
                       const double             pm = 0.10, // fraction of population to be replaced by mutants
                       const double           rhoe = 0.70, // probability that offspring inherit an allele from elite parent
                       const unsigned            K = 4,		 // number of independent populations
                       const unsigned      THREADS = 8,    // number of threads for parallel decoding
                       const unsigned      X_INTVL = 100,	 // exchange best individuals at every 100 generations
                       const unsigned     X_NUMBER = 2,	   // exchange top 2 best
                       const unsigned     MAX_GENS = 100,  // Max generations without improvement 
                       const unsigned RESET_AFTER = 200, 
                       const unsigned     verbose = 2, 
                       const long unsigned rngSeed = 0	   // seed to the random number generator
) {
  
  if ( THREADS > 0 )
    omp_set_num_threads( THREADS );
  //Rprintf("\nNumber of threads=%i\n", omp_get_max_threads());
  
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
  if (verbose >= 1) {
    Rprintf("+--------------+--------------+------------+---------+----------+\n");
    Rprintf("| Best BRKGA   | Local Search | Generation | Bst Gen | Duration |\n");
    Rprintf("+--------------+--------------+------------+---------+----------+\n");      
  }
  
  double loopTime = 0;
  do {
    gensLoosing++;
    algorithm.evolve(GEN_INTVL);	// evolve the population for one generation
    
    if (algorithm.getBestFitness() < Best_BK_Fitness) {
      Best_BK_Fitness = algorithm.getBestFitness();
      if (Best_BK_Fitness < Best_Fitness){
        bestGeneration = generation;
        Best_Fitness = Best_BK_Fitness;
        gensLoosing = 0;
      } 
    }
    
    if (verbose == 2) {
      Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |",    \
              -Best_BK_Fitness,                                 \
              -Best_LS_Fitness,                                 \
              generation,                                       \
              bestGeneration,                                   \
              loopTime);
      Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
    }
    
    if ((generation % LS_INTVL == 0)){
      std::list< std::unordered_set<unsigned> > tours;
#pragma omp parallel for  schedule(static)
      for(unsigned k_ = 0; k_ < K; k_++){
        std::vector<double> chromosome = algorithm.getBestPopulationChromosome(k_);
        std::unordered_set< unsigned > M = get_M(chromosome, m);
        std::unordered_set< unsigned > N = get_N(n, M);
        
        //Look for improvements and update a population
        arma::vec alpha  = getAlphaVector(M, DistanceMatrix);
          for(auto it_m = M.begin(); it_m != M.end(); ++it_m) {
            for(auto it_n = N.begin(); it_n != N.end(); ++it_n) {
              // calculando o delta z (vizinho - melhor_Solucao)
              double delta = alpha(*it_n) - alpha(*it_m) - DistanceMatrix(*it_m,*it_n);
              if (std::isgreaterequal(delta, 0.0001)) {
                alpha = alpha - DistanceMatrix.col(*it_m) + DistanceMatrix.col(*it_n);
                M.insert(*it_n);
                N.insert(*it_m);
                M.erase(*it_m);
                N.erase(*it_n);
                it_m = M.begin();
                it_n = N.begin();
              } 
            }
            auto time = std::chrono::steady_clock::now();
            diff = time - start;
            time_elapsed = std::chrono::duration <double, std::milli> (diff).count()/1000;
            if (time_elapsed > MAX_TIME) break; // Get out of here!
            // if (verbose == 2) {
            //   Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |",    \
            //           -Best_BK_Fitness,                                 \
            //           -Best_LS_Fitness,                                 \
            //           generation,                                       \
            //           bestGeneration,                                   \
            //           time_elapsed);
            //   Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
            // }
          }
            tours.push_front(M);
          }
          
          for(auto it = tours.begin(); it != tours.end(); it++){
            double fitness = -getTourFitness(*it, DistanceMatrix);
            if(fitness < Best_LS_Fitness){
              Best_LS_Fitness = fitness;
              BestTour = *it;
            }
          }
      }            
    //       }
    //       double fitness = -getTourFitness(M, DistanceMatrix);
    //       if(fitness < Best_LS_Fitness){
    //         Best_LS_Fitness = fitness;
    //         BestTour = M;
    //       }
    //   }
    // }
    
    if ((last_LS_fitness > Best_LS_Fitness) || (last_BK_fitness > Best_BK_Fitness)){
      ImprovedSol++;
      if (Best_LS_Fitness < Best_Fitness) {
        Best_Fitness = Best_LS_Fitness;
        bestGeneration = generation;
        gensLoosing = 0;
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
  } while ((loopTime <= MAX_TIME) && (gensLoosing <= MAX_GENS));
  if (verbose >= 1) {
    Rprintf("\r| %12.2f | %12.2f | %10i | %7i | %7.1fs |",    \
            -Best_BK_Fitness,                                 \
            -Best_LS_Fitness,                                 \
            generation,                                       \
            bestGeneration,                                   \
            loopTime);
    Rprintf("\n+--------------+--------------+------------+---------+----------+\r\b\r");
  }
  
  if (verbose) Rprintf("\n\nThis is the end. The Doors.\n");
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("Tour") = BestTour,
    Rcpp::Named("LSFitness") = -Best_LS_Fitness,
    Rcpp::Named("BKFitness") = -Best_BK_Fitness,
    Rcpp::Named("Generations Number") = generation,
    Rcpp::Named("Best Generation") = bestGeneration,
    Rcpp::Named("Improvement Number") = ImprovedSol,
    Rcpp::Named("Duration") = loopTime);
  return rst;
}




