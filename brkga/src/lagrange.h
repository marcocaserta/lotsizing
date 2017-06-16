#ifndef LAGRANGE_H
#define LAGRANGE_H

#include <limits>
#include <math.h>
#include <map>

void lambda_initialization2(double * lambda);
void lambda_initialization(double * lambda);
void lambda_random_perturbation(double * lambda, double * lambdaBest);
double lagrange_step(double * lambda, int * xL, int iter, double ** rc);
bool lambda_update(double * lambda, double & delta, int * xL, double & bestLagr, 
		     double & worstLagr, double zL, double lb, int iter, 
		   double & best300, double & start300);
#endif
