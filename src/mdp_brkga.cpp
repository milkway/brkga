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
// Used to sort ranking by second argument
bool compare(const std::pair<double, int>&i, const std::pair<double, int>&j)
{
  return i.second < j.second;
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

  // std::random_device rd;
  // std::mt19937 g(rd());
  // 
  // std::shuffle(M.begin(), M.end(),g);
  // std::shuffle(N_M.begin(), N_M.end(),g);
  
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
//' @param distances Distance matrix of the instance 
//' @param m Number of elements selected (m <= n)
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
Rcpp::List mdp_brkga(const arma::mat DistanceMatrix, 
                     const unsigned m, // Number of elements selected (m <= n)
                     const unsigned method =  0, // 0: BRKGA; 1: BRKGA + LS; 2: BRKGA(LS) 
                     const unsigned     MAX_TIME = 10,	 // run for 10 seconds
                     const unsigned            p = 500,	 // size of population
                     const double             pe = 0.20, // fraction of population to be the elite-set
                     const double             pm = 0.10, // fraction of population to be replaced by mutants
                     const double           rhoe = 0.70, // probability that offspring inherit an allele from elite parent
                     const unsigned            K = 1,		 // number of independent populations
                     const unsigned         MAXT = 8,    // number of threads for parallel decoding
                     const unsigned      X_INTVL = 100,	 // exchange best individuals at every 100 generations
                     const unsigned     X_NUMBER = 2,	   // exchange top 2 best
                     const unsigned     MAX_GENS = 100, // Max generations without improvement 
                     const unsigned RESET_AFTER = 200, 
                     const bool          verbose = false, 
                     const long unsigned rngSeed = 0	   // seed to the random number generator
) {
  
  if (DistanceMatrix.n_cols != DistanceMatrix.n_rows) Rcpp::stop("Distance matrix must be square!");
  unsigned n = DistanceMatrix.n_cols;				// size of chromosomes
  
  MTRand rng(rngSeed);				// initialize the random number generator
  
  MdpDecoder decoder(DistanceMatrix, m);			// initialize the decoder
  
  
  // initialize the BRKGA-based heuristic
  BRKGA< MdpDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
  
  unsigned generation = 0;	// current generation
  double lastfitness = 0;
  double gensLoosing = 0;
  
  //m = n*0.10;			// número de elementos do subconjunto de máxima diversidade
  clock_t start_t, end_t;
  
  //verificar se esta estagnado, se estiver por X_INTVL iteracoes, reestart.
  start_t = clock();
  double BestLocalSearchFitness = 0;
  unsigned relevantGeneration = 0;	// last relevant generation: best updated or reset called
  unsigned ImprovedSol = 0;
  std::vector< unsigned > BestSolution;
  std::vector< double > BestChromosome(n);
  Rcpp::List search;
  
  double loopTime = 0;
  do {
    gensLoosing++;
    algorithm.evolve();	// evolve the population for one generation
    Rcpp::Rcout << "Best solution generation " << generation << " found has objective value = " << algorithm.getBestFitness() << std::endl;
    
    // Local Search
    if (algorithm.getBestFitness() != lastfitness){
      search = localSearch(algorithm.getBestChromosome(), (double)algorithm.getBestFitness(), DistanceMatrix, m);
      std::vector<double> cr = Rcpp::as<std::vector<double>>(search["Chromosome"]);
      if (method == 2) 
        algorithm.exchangeAlleles(search["Chromosome"], 0, n-1, search["lsFitness"]);
      //std::list< unsigned > solution = search["Solution"];
      if ((double)search["lsFitness"] < BestLocalSearchFitness) {
        BestLocalSearchFitness = search["lsFitness"];
        BestSolution = Rcpp::as<std::vector<unsigned>>(search["Solution"]);
        BestChromosome = Rcpp::as<std::vector<double>>(search["Chromosome"]);
        ImprovedSol++;
        gensLoosing = 0;
        relevantGeneration = generation;
        //Rcpp::Rcout << "\n Inserting best of local search in population\n\n";
      }
      
    }
    
    Rcpp::Rcout << "Local Search  of generation : " << BestLocalSearchFitness << "*********" << std::endl;
    Rcpp::Rcout << "--------------" << std::endl;
    
    lastfitness = algorithm.getBestFitness();
    
    //  Evolution strategy: restart
    if(generation - relevantGeneration > RESET_AFTER) {
      algorithm.reset();	// restart the algorithm with random keys
      relevantGeneration = generation;
      
      Rcpp::Rcout << "\t" << generation << ") Reset at generation "
                  << generation << std::endl;
    }
    
    if((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
    }
    end_t = clock();
    loopTime = (double)(difftime(end_t,start_t)/CLOCKS_PER_SEC);
    //} while (generation < MAX_GENS);
  } while ((loopTime <= MAX_TIME) && (gensLoosing <= MAX_GENS));
  
  Rcpp::Rcout << "Best solution found has objective value = " << algorithm.getBestFitness() << "em " << difftime(end_t,start_t)/CLOCKS_PER_SEC << " segundos" << std::endl << std::endl;
  Rcpp::List rst = Rcpp::List::create(
    Rcpp::Named("LSChr") = BestChromosome,
    Rcpp::Named("BKChr") = algorithm.getBestChromosome(),
    Rcpp::Named("Tour") = BestSolution,
    Rcpp::Named("Fitness") = algorithm.getBestFitness(),
    Rcpp::Named("lsFitness") = BestLocalSearchFitness,
    Rcpp::Named("Generations Number") = generation,
    Rcpp::Named("Improvement Number") = ImprovedSol,
    Rcpp::Named("Duration") = loopTime);
  
  return rst;
  }
