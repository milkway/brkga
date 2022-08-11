/*
 * caladoInstance.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: Rodrigo
 */

#include "caladoInstance.h"
using namespace std;

//-----------------------------[ Constructor ]--------------------------------//

caladoInstance::caladoInstance(std::vector<double> _DistanceVector_, std::vector<double> _DemandVector_, std::vector<double> _DraftVector_, unsigned _N_):
    num_nodes(0), distances(), demand(), draft()
{
    //ifstream fin;
    //fin.open(instanceFile/*.c_str()*/);
    //if(!fin)
    //    throw runtime_error("Cannot open instance file");



    //fin.exceptions(ifstream::failbit | ifstream::badbit);
    try {
      num_nodes = _N_;

    	distances.resize(num_nodes * num_nodes);

    	for(unsigned i = 0; i < distances.size(); ++i)
        	distances[i] = _DistanceVector_[i];

    	demand.resize(num_nodes);
    	for(unsigned i = 0; i < demand.size() ; ++i)
    		demand[i] = _DemandVector_[i];
    	draft.resize(num_nodes);
    	for(unsigned i = 0; i < draft.size() ; ++i)
    		draft[i] = _DraftVector_[i];
    }
    catch(std::ifstream::failure& e) {
        throw fstream::failure("Error reading the instance file");
    }
}

//-------------------------------[ Distance ]---------------------------------//

double caladoInstance::menordistance(unsigned Li) const {
    unsigned i, j;
    //unsigned menor=999;

    for(i = 1; i<num_nodes;++i){
    	if(/*distances[destination*num_nodes+ i] < menor &&*/ distances[(i-1)*num_nodes+i] != 0 && Li<=draft[i]){
    		//menor = distances[destination*num_nodes+i];
    		j = i;
    	}
    }
    return j;
}
