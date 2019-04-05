/**
 * reg-brkga.cpp
 *
 * Driver class with a simple example of how to instantiate and use the BRKGA API.
 * See SampleDecoder.h for details on the decoder's implementation.
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

//#include <iostream>
#include <algorithm>
#include "BRKGA.h"
#include "MTRand.h"
#include "regDecoder.h"

// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Rcpp::List reg_brkga(SEXP func_,
              arma::vec lowerLimit, // lower limit
              arma::vec upperLimit, // upper limit
              const unsigned p = 100,		// size of population
              const double pe = 0.10,		// fraction of population to be the elite-set
              const double pm = 0.10,		// fraction of population to be replaced by mutants
              const double rhoe = 0.70,	// probability that offspring inherit an allele from elite parent
              const unsigned K = 3,	  	// number of independent populations
              const unsigned MAXT = 1,  	// number of threads for parallel decoding
              const bool verbose = false,
              Rcpp::Nullable<Rcpp::NumericMatrix> data  = R_NilValue,
              const long unsigned rngSeed = 0,	// seed to the random number generator
              const unsigned X_INTVL = 100,	// exchange best individuals at every X_INTVL generations
              const unsigned X_NUMBER = 2,	// exchange top X_NUMBER best
              const unsigned MAX_GENS = 1000	// run for MAX_GENS gens
              ) { // 
  
  //if (lowerLimit.n_cols > 1 || upperLimit.n_cols  > 1) Rcpp::stop("Limits vectors must have only une column");
  if (lowerLimit.n_rows != upperLimit.n_rows) Rcpp::stop("Limit vector must have the same lenght");
  
  const unsigned n = lowerLimit.n_rows; // size of chromosomes

  
  
  //if (data.isNotNull()) {
  Rcpp::NumericMatrix nm_data;
  if (data.isNotNull()) nm_data = Rcpp::as<Rcpp::NumericMatrix>(data);
  arma::mat data_arma = Rcpp::as<arma::mat>(nm_data);
  

	regDecoder decoder(func_, lowerLimit, upperLimit, data_arma);				// initialize the decoder
	
	//const long unsigned rngSeed = 0;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator
	
	
	// initialize the BRKGA-based heuristic
	BRKGA< regDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
	
	unsigned generation = 0;		// current generation
	// const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
	// const unsigned X_NUMBER = 2;	// exchange top 2 best
	// const unsigned MAX_GENS = 1000;	// run for 1000 gens
	std::cout << "Running for " << MAX_GENS << " generations..." << std::endl;
	do {
		algorithm.evolve();	// evolve the population for one generation
		
		if((++generation) % X_INTVL == 0) {
			algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
		}
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
	
	Rcpp::Rcout << "Best solution found has objective value = "
	 		<< algorithm.getBestFitness() << std::endl;
	
	Rcpp::NumericVector Point(n);
	Rcpp::CharacterVector rowNames(n);
	Rcpp::NumericVector Chromosome(n);
	auto BestChromosome = algorithm.getBestChromosome();
	for (int i = 0; i < n; i++){
	  std::ostringstream varname;
	  varname << "x_" << std::to_string(i+1);
	  rowNames[i] = varname.str();
	  Chromosome[i] = BestChromosome[i];
	  Point[i] = BestChromosome[i]*(upperLimit[i] - lowerLimit[i]) + lowerLimit[i];
	}
	Chromosome.attr("names") = rowNames;
	Point.attr("names") = rowNames;
	Rcpp::List rst = Rcpp::List::create(
	  Rcpp::Named("Chromosome") =  Chromosome,
	  Rcpp::Named("Point") =  Point,
	  Rcpp::Named("Optimum") = algorithm.getBestFitness()
	);
	return rst;
}