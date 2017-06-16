/*
 * Population.h
 *
 * Encapsulates a population of chromosomes represented by a vector of doubles. We don't decode
 * nor deal with random numbers here; instead, we provide private support methods to set the
 * fitness of a specific chromosome as well as access methods to each allele. Note that the BRKGA
 * class must have access to such methods and thus is a friend.
 *
 *  Created on : Jun 21, 2010 by rtoso
 *  Last update: Nov 15, 2010 by rtoso
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 */

#ifndef POPULATION_H
#define POPULATION_H

#include <vector>
#include <algorithm>
#include <exception>
#include <stdexcept>

class Population {
	template< class Decoder, class RNG >
	friend class BRKGA;

public:
	unsigned getN() const;	// Size of each chromosome
	unsigned getP() const;	// Size of population

	//double operator()(unsigned i, unsigned j) const;	// Direct access to allele j of chromosome i

	// These methods REQUIRE fitness to be sorted, and thus a call to sortFitness() beforehand
	// (this is done by BRKGA, so rest assured: everything will work just fine with BRKGA).
	double getBestFitness() const;			// Returns the best fitness in this population
	double getFitness(unsigned i) const;	// Returns the fitness of chromosome i
	const std::vector< double >& getChromosome(unsigned i) const;	// Returns i-th best chromosome

private:
	Population(const Population& other);
	Population(unsigned n, unsigned p);
	~Population();

	std::vector< std::vector< double > > population;		// Population as vectors of prob.
	std::vector< std::pair< double, unsigned > > fitness;	// Fitness (double) of a each chromosome

	void sortFitness();									// Sorts 'fitness' by its first parameter
	void setFitness(unsigned i, double f);				// Sets the fitness of chromosome i
	std::vector< double >& getChromosome(unsigned i);	// Returns a chromosome

	double& operator()(unsigned i, unsigned j);		// Direct access to allele j of chromosome i
	std::vector< double >& operator()(unsigned i);	// Direct access to chromosome i
};

Population::Population(const Population& pop) :
		population(pop.population),
		fitness(pop.fitness) {
}

Population::Population(const unsigned n, const unsigned p) :
		population(p, std::vector< double >(n, 0.0)), fitness(p) {
	if(p == 0) { throw std::range_error("Population size p cannot be zero."); }
	if(n == 0) { throw std::range_error("Chromosome size n cannot be zero."); }
}

Population::~Population() {
}

unsigned Population::getN() const {
	return population[0].size();
}

unsigned Population::getP() const {
	return population.size();
}

double Population::getBestFitness() const {
	return getFitness(0);
}

double Population::getFitness(unsigned i) const {
	return fitness[i].first;
}

const std::vector< double >& Population::getChromosome(unsigned i) const {
	return population[ fitness[i].second ];
}

std::vector< double >& Population::getChromosome(unsigned i) {
	return population[ fitness[i].second ];
}

void Population::setFitness(unsigned i, double f) {
	fitness[i].first = f;
	fitness[i].second = i;
}

void Population::sortFitness() {
	sort(fitness.begin(), fitness.end());
}

//double Population::operator()(unsigned chromosome, unsigned allele) const {
//	return population[chromosome][allele];
//}

double& Population::operator()(unsigned chromosome, unsigned allele) {
	return population[chromosome][allele];
}

std::vector< double >& Population::operator()(unsigned chromosome) {
	return population[chromosome];
}

#endif
