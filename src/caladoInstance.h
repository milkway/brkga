/*
 * caladoInstance.h
 *
 * Reads an instance from caladoLIB (Symmetric calado).
 *
 * Here's the URL: http://www.iwr.uni-heidelberg.de/groups/comopt/software/caladoLIB95/calado/
 *
 * Here's the format:
 *
 * NAME : a280
 * COMMENT : drilling problem (Ludwig)
 * TYPE : calado
 * DIMENSION: 280
 * EDGE_WEIGHT_TYPE : EUC_2D
 * NODE_COORD_SECTION
 * 1 288 149
 * 2 288 129
 * 3 270 133
 * 4 256 141
 * ...
 * EOF
 *
 *  Created on: Mar 16, 2013
 *      Author: Rodrigo
 */

#ifndef CALADOINSTANCE_H
#define CALADOINSTANCE_H



#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
using namespace std;

class caladoInstance {
public:
    /// Default Constructor.
	typedef std::runtime_error Error;
    caladoInstance(std::vector<double> _DistanceMatrix_, std::vector<double> _DemandVector_, std::vector<double> _DraftVector_, unsigned _N_);

	//caladoInstance(const std::string& instanceFile) throw(Error);
	//virtual ~caladoInstance();

    /// Return the distance between nodes i and j.
    double menordistance(unsigned Li) const;

    /// Number of nodes.
    unsigned num_nodes;

    /// Distances between the nodes.
    std::vector<double> distances;

    /// Demand
    std::vector<double> demand;

    /// Draft
    std::vector<double> draft;
};

#endif
