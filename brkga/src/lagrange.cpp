/***************************************************************************
 *   copyright (C) 2015 by Marco Caserta                                   *
 *   marco dot caserta at ie dot edu                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/**! \file lagrange.cpp
  \brief Lagrangean framework for the Muldichoice Multidimension Knapsack Problem.

Author: Marco Caserta (marco dot caserta at ie dot edu)\n
Started: 08.02.15\n
Ended:   \n 

*/

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cassert>

#include "lagrange.h"

using namespace std;

extern ofstream flagr;

struct INSTANCE {
    int nR;			// number of resources (constraints)
    int nC;			// number of classes
    int * ri;			// number of items in each class i

    double  ** c;
    int *** w;
    int   * R;
};

extern INSTANCE inp;
const double EPSI      = 0.00000001;
const long _MAXRANDOM  = 2147483647;
const double ZERO      = 0.0e0;
extern double INFTY;

void lambda_initialization2(double * lambda)
{
    int minW;
    for (int k = 0; k < inp.nR; k++)
    {
        lambda[k] = ZERO;
        for (int i = 0; i < inp.nC; i++)
        {
            minW = _MAXRANDOM;
            for (int j = 0; j < inp.ri[i]; j++)
                if (inp.w[i][j][k] < minW)
                    minW = inp.w[i][j][k];

            lambda[k] += (double)minW;
        }
        lambda[k] /= (double)inp.R[k];
    }
}

void lambda_initialization(double * lambda)
{
    double  score, maxW;
    for (int k = 0; k < inp.nR; k++)
    {
        lambda[k] = ZERO;
        for (int i = 0; i < inp.nC; i++)
        {
            maxW = -INFTY;
            for (int j = 0; j < inp.ri[i]; j++)
            {
                if (inp.w[i][j][k] > 0)
                    score = (double)inp.c[i][j]/(double)inp.w[i][j][k];
                if (score > maxW)
                    maxW = score;
            }
            lambda[k] += (double)maxW;
        }
        lambda[k] = (double)inp.R[k]/lambda[k];
    }
}

void lambda_random_perturbation(double * lambda, double * lambdaBest)
{
    double rr;
    for (int k = 0; k < inp.nR; k++)
    {
        //cout << "old lambdaBest( " << k << ") = " << lambdaBest[k] << endl;
        // random nr between 0.0 and 0.2
        rr = -0.1 + ((double)(rand()+1)/((double)(RAND_MAX)+1.0))/5.0;
        //cout << "rr is " << rr << endl;

        lambda[k] = lambdaBest[k]*(1.0 + rr);
        //cout << "new lambda( " << k << ") = " << lambda[k] << endl;
    }
}


/// Get lagrangean solution and compute obj function value of relaxation
double lagrange_step(double * lambda, int * xL, int iter, double ** rc)
{
    double rcBest;
    double zL       = ZERO;		// lagrangean function

    // find best item in each  class and build Lagrangean solution
    for (int i = 0; i < inp.nC; i++)
    {
        //cout << "COMPONENT " << i << endl;

        rcBest = -INFTY;
        //rcBest = ZERO;
        xL[i]  = -1;
        for (int j = 0; j < inp.ri[i]; j++)
        {
            rc[i][j] = inp.c[i][j];

            for (int k = 0; k < inp.nR; k++)
                rc[i][j] -= lambda[k]*inp.w[i][j][k];

            if (rc[i][j] > rcBest)
            {
                rcBest = rc[i][j];
                xL[i]  = j;
            }
            //cout << " ... j = " << j << " " << rc[i][j] << " vs " << rcBest << endl;
        }
        assert(xL[i] != -1);
        // add term to lagrangean function
        zL += rcBest;
    }

    // complete lagrangean function
    for (int k = 0; k < inp.nR; k++)
        zL += lambda[k]*inp.R[k];

#ifdef M_DEBUG
    cout << "ub(" << iter << ") = " << zL << endl;

    flagr << iter << "\t" << zL << endl;
#endif
    return zL;
}

/** This function is written for a max problem, i.e., we want to minimize the value
  of the lagrangean function. Every lagrangean value zL is an upper bound of the
  optimal solution.
  */
bool lambda_update(double * lambda, double & delta, int * xL, double & bestLagr, 
        double & worstLagr, double zL, double lb, int iter, double & best300, double & start300)
{
    if ((iter % 300) == 0)
    {
        if (((start300 - best300)/start300) < 0.001) return true;

        best300  = zL;
        start300 = zL;
    }

    if (zL < best300)
        best300 = zL;

    // update step size
    if ((iter % 20) == 0)
    {
        //cout << "z* = " << lb << " with DELTA = " << delta << " (best and worst Lagr " << bestLagr 
        //     << ", " << worstLagr << ")" << endl;
        if (((worstLagr - bestLagr)/worstLagr) > 0.01)
        {
            delta *= 0.5;
            //cout << "delta is halved " << endl;
        }
        else if (((worstLagr - bestLagr)/worstLagr) < 0.001)
        {
            delta *= 1.5;
            //cout << "delta multipled by 1.5 " << endl;
        }

        // now reset best and worst lagrangean values
        bestLagr  = zL;
        worstLagr = zL;
    }
    else
        if (zL < bestLagr)
            bestLagr = zL;
        else if (zL > worstLagr)
            worstLagr = zL;


    long double normS = ZERO;
    double * s        = new double[inp.nR];
    // compute subgradient components
    for (int k = 0; k < inp.nR; k++)
    {
        s[k] = ZERO;
        for (int i = 0; i < inp.nC; i++)
            s[k] += (double)inp.w[i][xL[i]][k];

        s[k] -= (double)inp.R[k];

        normS += (long double)s[k]*(long double)s[k];
    }
    normS = sqrt(normS);

    if (normS < EPSI)
    {
        cout << "Norm of subgradient is zero " << endl;
        exit(123);
    }

    // update lambda_k
    for (int k = 0; k < inp.nR; k++)
    {
        lambda[k] +=  delta*(zL-lb)*s[k]/normS;
        if (lambda[k] < ZERO) lambda[k] = ZERO;
    }

    return false;		// return not finished
}

