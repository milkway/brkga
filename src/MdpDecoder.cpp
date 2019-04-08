/*
 * MdpDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */
//#include <omp.h>
#include "MdpDecoder.h"

using namespace std;

extern 	std::vector< double >  dist;	// Vector of distances ROW*COL elements
extern unsigned m;	// number of elements set M

MdpDecoder::MdpDecoder() { }

MdpDecoder::~MdpDecoder() { }


// Runs in \Theta(n \log n):
double MdpDecoder::decode(const std::vector< double >& chromosome) const {


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
	std::list< unsigned > permutation;
	for(std::vector< std::pair< double, unsigned > >::const_iterator i = ranking.begin();
			i != ranking.begin()+m; ++i) {
		permutation.push_back(i->second);
	}

	// Here we compute the fitness of chromosome by permutation obtained
	
	unsigned elem1, elem2;
  	
	for (std::list< unsigned >::const_iterator i = permutation.begin(); i != (permutation.end()); ++i) {
		elem1 = *i;
		for (std::list< unsigned >::const_iterator j = i;j != permutation.end(); ++j) {
			elem2 = *j;
			myFitness += dist[elem1*chromosome.size()+elem2];
		}
	}

	return -myFitness;
}


