/*
 * caladoSolver.h
 *
 *  Created on: Mar 16, 2013
 *      Author: Rodrigo
 */

#ifndef CALADOSOLVER_H
#define CALADOSOLVER_H

#include <list>
#include <limits>
#include <vector>
#include <algorithm>
#include "caladoInstance.h"

class caladoSolver {
public:
	// The constructor 'solves' the problem in O(n log n) by transforming the chromosome into
	// a tour (by getting a permutation out of the chromosome):
	caladoSolver(const caladoInstance& instance, const std::vector< double >& chromosome);
	virtual ~caladoSolver();

	unsigned getTourDistance() const;		// Returns the tour distance
	std::list< unsigned > getTour() const;	// Returns the tour (first node not copied in the end)
	list<unsigned>myTour() const;
	list<unsigned> it() const;
	//list<unsigned>ite;

private:
	typedef std::pair< double, unsigned > ValueKeyPair;

	int distance;
	std::vector< ValueKeyPair > tour;
	list<unsigned>iteration;

};

#endif
