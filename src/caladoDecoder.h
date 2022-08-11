/*
 * caladoDecoder.h
 *
 *  Created on: Mar 16, 2013
 *      Author: Rodrigo
 */

#ifndef CALADODECODER_H
#define CALADODECODER_H

#include "caladoSolver.h"
#include "caladoInstance.h"

class caladoDecoder {
public:
	caladoDecoder(const caladoInstance& instance);
	virtual ~caladoDecoder();

	// Decodes a chromosome into a solution to the calado:
	double decode(const std::vector< double >& chromosome) const;

private:
	const caladoInstance& instance;
};

#endif
