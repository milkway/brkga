/*
 * caladoDecoder.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: Rodrigo
 */

#include "caladoDecoder.h"

caladoDecoder::caladoDecoder(const caladoInstance& _instance) : instance(_instance) {
}

caladoDecoder::~caladoDecoder() {
}

double caladoDecoder::decode(const std::vector< double >& chromosome) const {
	// 1) Solve the problem (i.e., create a tour out of this chromosome):
	// Avoids race conditions by making sure we have a single caladoSolver for each thread calling
	// ::decode (as long as caladoSolver does not make use of 'static' or other gimmicks):
	caladoSolver solver(instance, chromosome);

	// 2) Extract the fitness (tour distance):
	const unsigned fitness = solver.getTourDistance();

	// 3) Return:
	return double(fitness);
}

