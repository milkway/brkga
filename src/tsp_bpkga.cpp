/**
 * brkga-tsp.cpp
 *
 * Driver class with a simple example of how to instantiate and use the BRKGA API to find solutions
 * to the symmetric traveling salesman problem (TSP) on TSPLIB instances.
 *
 * See TSPDecoder.h for details on the decoder's implementation.
 *
 * Created on : Nov 17, 2011 by rtoso
 * Authors    : Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 *              Mauricio G.C. Resende <mgcr@research.att.com>
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2018
 * Rodrigo Franco Toso (rfrancotoso@gmail.com) and
 * Mauricio G.C. Resende
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


#include <algorithm>
#include <sstream>
#include "BRKGA.h"
#include "MTRand.h"

#include "TSPSolver.h"
#include "TSPDecoder.h"
#include "TSPInstance.h"

//' Execute a MDP solution search using BRKGA
//' @details Escrever um detalhamento bem legal e super bacana
//'  aqui.
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
Rcpp::List tsp_brkga(std::string instanceFile) {

	Rcpp::Rcout << "Welcome to the BRKGA API sample driver.\nFinding a (heuristic) minimizer for "
			<< " the TSP." << std::endl;

	const clock_t begin = clock();

	Rcpp::Rcout << "Instance file: " << instanceFile << std::endl;

	// Read the instance:
	TSPInstance instance(instanceFile); 	// initialize the instance
	Rcpp::Rcout << "Instance read; here's the info:"
			<< "\n\tName: " << instance.getName()
			<< "\n\tComment: " << instance.getComment()
			<< "\n\tDimension: " << instance.getNumNodes()
			<< "\n\tEdge type: " << instance.getProblemType()
			<< "\n\tEdge Weight Type: " << instance.getEdgeWeightType() << std::endl;

	TSPDecoder decoder(instance);		// initialize the decoder

	const long unsigned rngSeed = time(0);	// seed to the random number generator
	MTRand rng(rngSeed);					// initialize the random number generator

	const unsigned n = instance.getNumNodes();		// size of chromosomes
	const unsigned p = 256;		// size of population
	const double pe = 0.10;		// fraction of population to be the elite-set
	const double pm = 0.10;		// fraction of population to be replaced by mutants
	const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
	const unsigned K = 3;		// number of independent populations
	const unsigned MAXT = 2;	// number of threads for parallel decoding

	// initialize the BRKGA-based heuristic
	BRKGA< TSPDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);

	// BRKGA inner loop (evolution) configuration: Exchange top individuals
	const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
	const unsigned X_NUMBER = 2;	// exchange top 2 best
	const unsigned MAX_GENS = 1000;	// run for 1000 gens

	// BRKGA evolution configuration: restart strategy
	unsigned relevantGeneration = 1;	// last relevant generation: best updated or reset called
	const unsigned RESET_AFTER = 200;
	std::vector< double > bestChromosome;
	double bestFitness = std::numeric_limits< double >::max();

	// Print info about multi-threading:
	#ifdef _OPENMP
	Rcpp::Rcout << "Running for " << MAX_GENS << " generations using " << MAXT
				<< " out of " << omp_get_max_threads()
				<< " available thread units..." << std::endl;
	#endif
	#ifndef _OPENMP
	Rcpp::Rcout << "Running for " << MAX_GENS
				<< " generations without multi-threading..." << std::endl;
	#endif
	
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
			
			Rcpp::Rcout << "\t" << generation
					<< ") Improved best solution thus far: "
					<< bestFitness << std::endl;
		}

		//  Evolution strategy: restart
		if(generation - relevantGeneration > RESET_AFTER) {
			algorithm.reset();	// restart the algorithm with random keys
			relevantGeneration = generation;
			
			Rcpp::Rcout << "\t" << generation << ") Reset at generation "
					<< generation << std::endl;
		}

		// Evolution strategy: exchange top individuals among the populations
		if(generation % X_INTVL == 0 && relevantGeneration != generation) {
			algorithm.exchangeElite(X_NUMBER);
			
			Rcpp::Rcout << "\t" << generation
					<< ") Exchanged top individuals." << std::endl;
		}

		// Next generation?
		++generation;
	} while (generation < MAX_GENS);

	// print the fitness of the top 10 individuals of each population:
	Rcpp::Rcout << "Fitness of the top 10 individuals of each population:" << std::endl;
	const unsigned bound = std::min(p, unsigned(10));	// makes sure we have 10 individuals
	for(unsigned i = 0; i < K; ++i) {
	  Rcpp::Rcout << "Population #" << i << ":" << std::endl;
		for(unsigned j = 0; j < bound; ++j) {
		  Rcpp::Rcout << "\t" << j << ") "
					<< algorithm.getPopulation(i).getFitness(j) << std::endl;
		}
	}

	// rebuild the best solution:
	TSPSolver bestSolution(instance, bestChromosome);

	// print its distance:
	Rcpp::Rcout << "Best solution found has objective value = "
	 		<< bestSolution.getTourDistance() << std::endl;

	// print its best tour:
	Rcpp::Rcout << "Best tour:";
	const std::list< unsigned > bestTour = bestSolution.getTour();
	for(std::list< unsigned >::const_iterator pt = bestTour.begin(); pt != bestTour.end(); ++pt) {
	  Rcpp::Rcout << " " << *pt;
	}
	Rcpp::Rcout << std::endl;

	const clock_t end = clock();
	Rcpp::Rcout << "BRKGA run finished in " << (end - begin) / double(CLOCKS_PER_SEC) << " s." << std::endl;

	return Rcpp::List::create(Rcpp::Named("Chromosome") = algorithm.getBestChromosome(),
                           Rcpp::Named("BestFitness") = algorithm.getBestFitness(),
                           Rcpp::Named("Tour") = bestTour);
}
