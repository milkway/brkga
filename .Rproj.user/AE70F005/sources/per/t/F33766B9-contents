/*
 * regDecoder.cpp
 *
 * For more information, see SampleDecoder.h
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
#include "regDecoder.h"

regDecoder::regDecoder(SEXP _func_, arma::vec lowerLimit, arma::vec upperLimit, arma::mat data) :  func_(_func_), lLim(lowerLimit), uLim(upperLimit), data_(data) { }

regDecoder::~regDecoder() { }

// Runs in O(n \log n):
double regDecoder::decode(const std::vector< double >& chromosome) const {
  funcPtr func = *Rcpp::XPtr<funcPtr>(func_);
	typedef std::pair< double, unsigned > ValueKeyPair;
	std::vector< ValueKeyPair > rank(chromosome.size());
	int N = chromosome.size();
	std::vector< double > Point(N);
	for(unsigned i = 0; i < N; ++i) {
	  rank[i] = ValueKeyPair(chromosome[i], i);
	  Point[i] = chromosome[i]*(uLim[i] - lLim[i]) + lLim[i];
	}
	double myFitness = Rcpp::as<double>(func(Point, data_));
	// Here we sort 'permutation', which will then produce a permutation of [n]
	// stored in ValueKeyPair::second:
	std::sort(rank.begin(), rank.end());

	// permutation[i].second is in {0, ..., n - 1}; a permutation can be obtained as follows
	std::list< unsigned > permutation;
	for(std::vector< ValueKeyPair >::const_iterator i = rank.begin(); i != rank.end(); ++i) {
		permutation.push_back(i->second);
	}
	
	//Rcpp::Rcout << "Data Matrix " << data_.n_rows << " : " << data_.n_cols << " : "<< arma::size(data_) << std::endl;

	// Return the fitness:
	return myFitness;
}
