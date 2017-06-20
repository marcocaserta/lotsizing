#include <iostream>
#include "SampleDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) 
{
   cout << "starting with brkga ... " << endl;
   ofstream fsol("bestBRKGA.csv", ios::out);
   fsol << "Initializing file " << endl;
   cout << "--------------> Generation[0] " << endl;

   const unsigned n = 4;		// size of chromosomes
   const unsigned p = 50;	// size of population
   const double pe = 0.20;		// fraction of population to be the elite-set
   const double pm = 0.10;		// fraction of population to be replaced by mutants
   const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
   const unsigned K = 1;		// number of independent populations
   const unsigned MAXT = 2;	// number of threads for parallel decoding
	
   SampleDecoder decoder;			// initialize the decoder
	
   const long unsigned rngSeed = time(0);	// seed to the random number generator
   MTRand rng(rngSeed);				// initialize the random number generator

   // initialize the BRKGA-based heuristic
   BRKGA< SampleDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
	
   unsigned generation = 0;		// current generation
   cout << "--------------> Generation[0] " << endl;

   const unsigned X_INTVL = 5;	// exchange best individuals at every 100 generations
   const unsigned X_NUMBER = 2;	// exchange top 2 best
   const unsigned MAX_GENS = 20;	// run for 1000 gens
   cout << "--------------> Generation[" << generation << "] " << endl;

   do {

      algorithm.evolve();	// evolve the population for one generation
      cout << "--------------> Generation[" << generation << "] Current best is " << algorithm.getBestFitness() << endl;

      fsol << generation << "\t" << algorithm.getBestFitness();
      std::vector <double> bestX =algorithm.getBestChromosome();
      for (unsigned j = 0; j < bestX.size(); j++)
	 fsol << " \t " << bestX[j];
      fsol << endl;

      if((++generation) % X_INTVL == 0) {
	 algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
      }
   } while (generation < MAX_GENS);
	
   std::cout << "Best solution found has objective value = "
	     << algorithm.getBestFitness() << std::endl;
   std::vector <double> bestX =algorithm.getBestChromosome();
   for (unsigned j = 0; j < bestX.size(); j++)
      cout << "    " << bestX[j];
   cout << endl;
   for (unsigned j = 0; j < bestX.size(); j++)
      fsol << "    " << bestX[j];
   fsol << endl;
      
   fsol.close();
   return 0;

}
