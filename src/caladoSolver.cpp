/*
 * TSPSolver.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: Rodrigo
 */

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "caladoSolver.h"
#include <iostream>
using namespace std;


caladoSolver::caladoSolver(const caladoInstance& instance, const std::vector< double >& chromosome) :
		distance(0), tour(instance.num_nodes), iteration(){

	 unsigned Li=0;
	// unsigned menor=9999;
	 //list<unsigned> iteration1;


	 /* Li vai ser a soma das demandas,
	   cada vez que o navio passar por um porto vai ser diminuido essa demanda para compara com o limite de calado*/
	 for(unsigned i =0; i<instance.demand.size();i++){
		 Li += instance.demand[i];

	 }
	 //cout<< Li << endl;
	 //int cont=0;
	 //int menor;

	// 1) Obtain a permutation out of the chromosome -- this will be the tour:
	for(unsigned i = 0; i < chromosome.size(); ++i) { tour[i] = ValueKeyPair(chromosome[i], i); }

	// Here we sort 'rank', which will produce a permutation of [n] stored in ValueKeyPair::second:
	std::sort(tour.begin(), tour.end());
	int i=1;
	// 2) Compute the distance of the tour given by the permutation:
		//while(Li>=0){
		for(unsigned i = 1;i<tour.size() ; ++i) {
			// Compute distance(i-1, i) in the permutation:
			const unsigned& source = tour[i-1].second;
			const unsigned&  destination = tour[i].second;


				if(instance.draft[destination ] >= Li ){
						distance += instance.distances[source*tour.size()+destination];
						Li -= instance.demand[destination];
						iteration.push_back(destination);

				}

				/*else{
					i++;
					distance += instance.distances[source*tour.size()+destination];
					Li -= instance.demand[destination];
					iteration.push_back(destination);
				}*/


			}


	const unsigned& last = tour.back().second;
	const unsigned& first = tour.front().second;
	distance += instance.distances[last*tour.size()+first];
	iteration.push_back(first);
}

caladoSolver::~caladoSolver() {
}

unsigned caladoSolver::getTourDistance() const { return distance; }

std::list< unsigned > caladoSolver::getTour() const {
	std::list< unsigned > tourSequence;

	for(unsigned i = 0; i < tour.size(); ++i) { tourSequence.push_back(tour[i].second); }

	return tourSequence;
}

list<unsigned> caladoSolver::myTour() const{return iteration;}

//list<unsigned> TSPSolver::it() const{return ite;}
