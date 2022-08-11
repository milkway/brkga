#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <algorithm>
#include "BRKGA.h"
#include "MTRand.h"

#include "caladoSolver.h"
#include "caladoDecoder.h"
#include "caladoInstance.h"
//#include "caladoInstance.cpp"
//caladoDL *caladodl = create_caladoDL();



//' Execute a TSPDL solution search using BRKGA
//' @details Escrever um detalhamento bem legal e super bacana
//'  aqui.
//' @param \code{DistanceVector} distances matrix of the instance 
//' @param \code{DemandVector} distances matrix of the instance 
//' @param \code{DraftVector} distances matrix of the instance 
//' @param \code{N} Number of nodes
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
// [[Rcpp::export(" calado_brkga")]]
Rcpp::List calado_brkgaArma(const std::vector<double>    DistanceVector,  
                         const std::vector<double>    DemandVector,  
                         const std::vector<double>    DraftVector,  
                         const unsigned                 N,   // Number of nodes
                         const unsigned     LS_INTVL = 10,   // Generations between local searches 
                         const unsigned    GEN_INTVL = 5,    // Interval between Generations                      
                         const unsigned     MAX_TIME = 10,	  // run for 10 seconds
                         const unsigned            p = 500,	// size of population
                         const double             pe = 0.20, // fraction of population to be the elite-set
                         const double             pm = 0.10, // fraction of population to be replaced by mutants
                         const double           rhoe = 0.70, // probability that offspring inherit an allele from elite parent
                         const unsigned            K = 4,		// number of independent populations
                         const unsigned      THREADS = 8,    // number of threads for parallel decoding
                         const unsigned      X_INTVL = 100,	// exchange best individuals at every 100 generations
                         const unsigned     X_NUMBER = 2,	  // exchange top 2 best
                         const unsigned     MAX_GENS = 1000,  // Max generations without improvement 
                         const unsigned RESET_AFTER = 200, 
                         const unsigned     verbose = 2, 
                         long unsigned rngSeed = 0	    // seed to the random number generator
) {
  
  
  if ( THREADS > 0 )
    omp_set_num_threads( THREADS );
	//const clock_t begin = clock();

	if (0 == rngSeed) 
	  rngSeed = time(0);
	
  // Read the instance:
	caladoInstance instance(DistanceVector, DemandVector, DraftVector, N); 	// initialize the instance
	// std::cout << "Instance read; here's the info:"
	// 		<< "\n\tN: " << instance.num_nodes << endl;
	// 		cout<< "demand: "<< endl;
	// 		for(unsigned i=0; i < instance.num_nodes; i++)
	// 			cout<< instance.demand[i] << " ";
	// 		cout<<endl;
	// 		cout<< "draft: "<<endl;
	// 		for(unsigned i=0; i < instance.num_nodes; i++)
	// 			cout << instance.draft[i] << " ";
	// 		cout<<endl;

	caladoDecoder decoder(instance);		// initialize the decoder

	//const long unsigned rngSeed = time(0);	// seed to the random number generator
	MTRand rng(rngSeed);					// initialize the random number generator

	// const unsigned n = instance.num_nodes;		// size of chromosomes
	// const unsigned p = 256;		// size of population
	// const double pe = 0.10;		// fraction of population to be the elite-set
	// const double pm = 0.10;		// fraction of population to be replaced by mutants
	// const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
	// const unsigned K = 3;		// number of independent populations
	// const unsigned MAX_TIME = 2;	// number of threads for parallel decoding

	// initialize the BRKGA-based heuristic
	BRKGA< caladoDecoder, MTRand > algorithm(N, p, pe, pm, rhoe, decoder, rng, K, MAX_TIME);

	// BRKGA inner loop (evolution) configuration: Exchange top individuals
	//const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
	//const unsigned X_NUMBER = 2;	// exchange top 2 best
	//const unsigned MAX_GENS = 1000;	// run for 1000 gens

	// BRKGA evolution configuration: restart strategy
	unsigned relevantGeneration = 1;	// last relevant generation: best updated or reset called
	//const unsigned RESET_AFTER = 200;
	std::vector< double > bestChromosome;
	double bestFitness = std::numeric_limits< double >::max();

	auto start = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff;
	double time_elapsed = 0;
	double loopTime = 0;
	
	if (verbose >= 1) {
	  Rprintf("+--------------+------------+---------+----------+\n");
	  Rprintf("| Best BRKGA   | Generation | Bst Gen | Duration |\n");
	  Rprintf("+--------------+------------+---------+----------+\n");      
	}
	
	// Run the evolution loop:
	unsigned generation = 1;		// current generation
	do {
		algorithm.evolve();	// evolve the population for one generation

		// Bookeeping: has the best solution thus far improved?
		if(algorithm.getBestFitness() < bestFitness) {
			// Save the best solution to be used after the evolution chain:
			relevantGeneration = generation;
			bestFitness = algorithm.getBestFitness();
			bestChromosome = algorithm.getBestChromosome();
			
		//	std::cout << "\t" << generation
		//			<< ") Improved best solution thus far: "
		//			<< bestFitness << std::endl;
		}

		//  Evolution strategy: restart
		if(generation - relevantGeneration > RESET_AFTER) {
			algorithm.reset();	// restart the algorithm with random keys
			relevantGeneration = generation;
			
		//	std::cout << "\t" << generation << ") Reset at generation "
		//			<< generation << std::endl;
		}

		// Evolution strategy: exchange top individuals among the populations
		if(generation % X_INTVL == 0 && relevantGeneration != generation) {
			algorithm.exchangeElite(X_NUMBER);
			
		//	std::cout << "\t" << generation
		//			<< ") Exchanged top individuals." << std::endl;
		}
		if (verbose == 2) {
		  Rprintf("\r| %12.2f | %10i | %7i | %7.1fs |",    \
            bestFitness,                              \
            generation,                                       \
            relevantGeneration,                                   \
            loopTime);
		  Rprintf("\n+--------------+------------+---------+----------+\r\b\r");
		}
		

		// Next generation?
		++generation;
		
		// end of timing
		auto end = std::chrono::steady_clock::now();
		// Store the time difference between start and end
		diff = end - start;
		loopTime = std::chrono::duration <double, std::milli> (diff).count()/1000;
	} while ((loopTime <= MAX_TIME) && (generation < MAX_GENS));
	if (verbose >= 1 ) {
	  Rprintf("\r| %12.2f | %10i | %7i | %7.1fs |",    \
           bestFitness,                                 \
           generation,                                       \
           relevantGeneration,                                   \
           loopTime);
	  Rprintf("\n+--------------+------------+---------+----------+\r\b\r");
	}
	// print the fitness of the top 10 individuals of each population:
	// std::cout << "Fitness of the top 10 individuals of each population:" << std::endl;
	//const unsigned bound = std::min(p, unsigned(10));	// makes sure we have 10 individuals
	// for(unsigned i = 0; i < K; ++i) {
	// 	std::cout << "Population #" << i << ":" << std::endl;
	// 	for(unsigned j = 0; j < bound; ++j) {
	// 		std::cout << "\t" << j << ") "
	// 				<< algorithm.getPopulation(i).getFitness(j) << std::endl;
	// 	}
	// }

	// rebuild the best solution:
	caladoSolver bestSolution(instance, bestChromosome);

	//bestsolution.bestTour()


	// print its distance:
	//std::cout << "Best solution found has objective value = "
	 //		<< bestSolution.getTourDistance() << std::endl;

	//print its best tour:
	//std::cout << "Best tour:";
	//const std::list< unsigned > bestTour = bestSolution.getTour();
	//for(std::list< unsigned >::const_iterator pt = bestTour.begin(); pt != bestTour.end(); ++pt) {
	//std::cout << " " << *pt;
	//}
	//std::cout << std::endl;

	//std::cout << "My tour:";
	//const std::list< unsigned > myTour = bestSolution.myTour();
	//for(std::list< unsigned >::const_iterator tp = myTour.begin(); tp != myTour.end(); ++tp) {
	//std::cout << " " << *tp;
	//}
	//std::cout << std::endl;

	/*std::cout << "it:";
	const std::list< unsigned > it = bestSolution.ite();
	for(std::list< unsigned >::const_iterator tp = ite.begin(); tp != ite.end(); ++tp) {
		std::cout << " " << *tp;
	}
	std::cout << std::endl;*/

	//const clock_t end = clock();
	//std::cout << "BRKGA run finished in " << (end - begin) / double(CLOCKS_PER_SEC) << " s." << std::endl;

	if (verbose) Rprintf("\n\nThis is the end. The Doors.\n");
	Rcpp::List rst = Rcpp::List::create(
	  Rcpp::Named("Tour") = bestSolution.getTour(),
	  Rcpp::Named("MyTour") = bestSolution.myTour(),
	  Rcpp::Named("BKFitness") = bestSolution.getTourDistance() ,
	  Rcpp::Named("Generations Number") = generation,
	  Rcpp::Named("Best Generation") = relevantGeneration,
	  Rcpp::Named("Duration") = loopTime);
	return rst;
}
