/********************************en*******************************************
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

/*! \file robust.cpp

  \brief Robust Lagrangean for the robust MMKP
  \author Marco Caserta 2017 (c) 
  \version v. 1.0.0
  \date Start Date : 26.01.17
  \date Last Update: 
 **/
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN


#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cassert>
#include <map>

#include "timer.h"

using namespace std;

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

typedef IloArray < IloNumVarArray > TwoD;

extern double zBest;
extern int * xBest;
extern int max_iter;             // max number lagrangean iterations

extern double time_limit;        // !< wall-clock time limit
extern double corridorWidthBase; // base value for corridor width
extern int nSolBase;             // nr of sols to be found by cplex within corridor
extern double ubStar;        // best upper bound (lagr solution)
extern double propFixed0;        // % of vars to be fixed to zero
extern double propFixed1;        // % of vars to be fixed to one
extern int add_oldSol_cut;       // boolean: cut out old feasible solution
extern int  * F;             // F[i] = 3 ==> var 3 in class "i" is set to 1
extern int **F0;             // F0[i][j] = 1 ==> var x_{ij} is set to 0
extern double best_time;     // time to best solution
extern int bestIter;         // lagr iteration in which best was found

extern timer tTime;            //!< Ojbect clock to measure REAL and VIRTUAL (cpu) time

// ============================================================================
double defineRobustModel(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        IloObjective & obj, IloNumVar & Q_ilo, double Omega, double * sigma2);
bool lambda_update(double * lambda, double & delta, int * xL, double & bestLagr, 
        double & worstLagr, double zL, double lb, int iter, double & best300, 
        double & start300, double * sigma2);
double get_first_lb(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, 
        int nSol, IloNumVar & Q_ilo, double Omega, double * sigma2);
double lagrange_step(double * lambda, int * xL, int iter, double ** rc, 
        double Omega, double * sigma2);
double lagrange_step_SOCP(INSTANCE & inp, double Omega, double * sigma2, 
        int * xL, double * lambda, int iter, double ** rc);

void lambda_initialization(double * lambda);
int stopping_criteria(int iter, int max_iter, double time_limit);
void update_best(int * xLBest, int * xL, double & ubStar, double zL, double * lambda, double * lambdaBest);
void lambda_random_perturbation(double * lambda, double * lambdaBest);
double solve_KNAP(IloModel model, IloCplex cplex, int solLimit, int displayLimit, int timeLim);
void find_worst(int ** F0, double ** rc, double propFixed0);
void find_best(int * F, double ** rc);
double get_cplex_sol(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xIlo);
void add_z_cut(IloModel & model, IloCplex & cplex, TwoD & x_ilo, double zBest);
double refine_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        IloObjective & obj, int * xL, double & cWidth, int & nSol, 
        double & zBest, int * xIlo, int & iterImproved, int iter, 
        double ** rc, bool fix2zero, bool fix2one, IloRangeArray & cutSol, 
        int lagrIter);
// ============================================================================

void robust_lagrangean_phase(INSTANCE & inp, IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        IloObjective & obj, IloNumVar & Q_ilo, double Omega, double * sigma2)
{

    IloEnv env = model.getEnv();
    IloRangeArray cutSol(env);

    map<int,int> zList;
    // initialize data structure for Lagrangean phase
    double * lambda     = new double[inp.nR];
    double * lambdaBest = new double[inp.nR];
    int * xL            = new int[inp.nC]; // current lagrangean solution
    int * xLBest        = new int[inp.nC]; // best lagrangean solution (infeasible)
    int * xIlo          = new int[inp.nC]; // best feasible solution from cplex
    double lb           = ZERO;
    
    // nr feasible solutions to be found by cplex
    zBest = get_first_lb(model, cplex, x_ilo, obj, 1, Q_ilo, Omega, sigma2);
    cout << "First feasible solution = " << zBest << endl;

    lambda_initialization(lambda);

#ifdef M_DEBUG
    for (int k = 0; k < inp.nR; k++)
        cout << "l(" << k << ") = " << lambda[k] << endl;
#endif

    double ** rc = new double*[inp.nC]; // lagrangean reduced costs
    for (int i = 0; i < inp.nC; i++)
        rc[i] = new double[inp.ri[i]];


    int lagrIter = 0;
    int iter     = 0;
    while (lagrIter < 3 && !stopping_criteria(iter, max_iter, time_limit))
    {
        cout <<"Start Lagrangean Cycle Nr. " << lagrIter << endl;
        double cWidth     = corridorWidthBase;   // <---------- used to be 0.8
        int iterImproved  = 200;
        int nSol          = nSolBase;		// nr of solutions for cplex
        int Freq          = 20;

        double zL, bestLagr, worstLagr;

        double start300  = INFTY;
        double best300   = INFTY;
        // double delta     = 0.1;
        double delta     = 0.01;
        bool stopping    = false;
        bool fix2zero    = false;
        bool fix2one     = false;

        while (!stopping_criteria(iter, max_iter, time_limit) && !stopping)
        {
            // note: zL in an "upper bound" of the optimal value
            zL = lagrange_step(lambda, xL, iter, rc, Omega, sigma2);
            // zL = lagrange_step_SOCP(inp, Omega, sigma2, xL, lambda, iter, rc);
            if ((iter % (2*Freq)) == 0)
                cout << "Lagr(" << iter << ") = " << zL << endl;
            if (zL < ubStar)
                update_best(xLBest, xL, ubStar, zL, lambda, lambdaBest);

            stopping = lambda_update(lambda, delta, xL, bestLagr, worstLagr, 
                    zL, zBest, iter, best300, start300, sigma2); 
            if (iter > 199) 	// start fixing schemes
            {
                fix2one = true;
                fix2zero = true;
            }

            // refine lagrange solution to obtain a valid lower bound
            // if (iter > 100 && (iter % Freq) == 0)
            // if (iter >= 300 )
            if (iter >= 100 && (iter % Freq)== 0)
            {
                lb = refine_solution(model, cplex, x_ilo, obj, xL, cWidth, 
                        nSol, zBest, xIlo, iterImproved, iter, rc, fix2zero, 
                        fix2one, cutSol, lagrIter);
            }

            if ((iter > 250) && (iter - iterImproved) > 100)
            {
                nSol++;
                cWidth       *= 0.9;
                iterImproved  = iter;
                Freq         /= 2;
                if (Freq < 5) Freq = 5;

                //cout << "Enhancing corridor to : " << cWidth << " and nsol to " << nSol << endl;
            }

            iter++;
        } // end lagrangean cycle

        cout << "** ** ** SUMMARY ** ** ** " << endl;
        cout << "** ** ** Lagrangean Iteration Nr. " << lagrIter << " :: Best LB = " << zBest << endl;
        cout << setw(49) << " :: Best UB = " << ubStar << endl;
#ifdef M_DEBUG
        cout << "Before removing cuts " << cplex.getNrows() << " constraints " << endl;
#endif

        model.remove(cutSol);

#ifdef M_DEBUG
        cout << "After removing cuts " << cplex.getNrows() << " constraints " << endl;
#endif

        lagrIter++;
        iter = 0;		// reset counter

        if (lagrIter < 3)
            lambda_random_perturbation(lambda, lambdaBest);

    } // END OF LAGRANGEAN CYCLE -- repeat 3 times

    // copy best feasible solution into xBest
    for (int i = 0; i < inp.nC; i++) 
        xBest[i] = xIlo[i];

}

/// Get first feasible solution
double get_first_lb(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, 
        int nSol, IloNumVar & Q_ilo, double Omega, double * sigma2)
{

    defineRobustModel(model, cplex, x_ilo, obj, Q_ilo, Omega, sigma2);


    double statusBin = solve_KNAP(model, cplex, 3, 0, 10000); // call cplex
    // double statusBin = solve_KNAP(model, cplex, 9999, 4, 10000); // call cplex
    cout << "Status :: " << statusBin << endl;
    cout << "CPLEX = " << cplex.getStatus() << endl;

    return statusBin;
}

// Quadratic Model SOCP
double defineRobustModel(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
       IloObjective & obj, IloNumVar & Q_ilo, double Omega, double * sigma2)
{
    IloEnv env = model.getEnv();

    // robust constraint
    IloExpr lhsRobust(env);
    double rhs = 0.0;
    lhsRobust = -Q_ilo*Q_ilo;
    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
            lhsRobust += x_ilo[i][j]*x_ilo[i][j]*sigma2[i];

    model.add(lhsRobust <= rhs);

    for (int k = 0; k < inp.nR; k++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nC; i++)
            for (int j = 0; j < inp.ri[i]; j++)
                sum += (inp.w[i][j][k]*x_ilo[i][j]);

        model.add(sum + Omega*Q_ilo <= inp.R[k]);
        sum.end();
    }
    // multi-choice constraint    
    for (int i = 0; i < inp.nC; i++)
    {
        IloExpr sum(env);
        for (int j = 0; j < inp.ri[i]; j++)
            sum += x_ilo[i][j];

        model.add(sum == 1.0);
    }
    // add objective function
    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
            obj.setLinearCoef(x_ilo[i][j], inp.c[i][j]);
    model.add(obj);
}


// build Lagrangean solution (approximate Lagrange)
double lagrange_step(double * lambda, int * xL, int iter, double ** rc, 
        double Omega, double * sigma2)
{
    double rcBest;
    double zL       = ZERO;		// lagrangean function

    /* cout << "... lambdas :: ";
     * for (int k= 0; k < inp.nR; k++)
     *     cout << " " << lambda[k];
     * cout << endl; */

    // find best item in each  class and build Lagrangean solution
    for (int i = 0; i < inp.nC; i++)
    {
        // cout << "COMPONENT " << i << endl;
        rcBest = -INFTY;
        //rcBest = ZERO;
        xL[i]  = -1;
        for (int j = 0; j < inp.ri[i]; j++)
        {
            // cout << "c[ " << i <<"," << j << " ] " << inp.c[i][j] << endl;
            rc[i][j] = inp.c[i][j];

            double ss = 0.0;
            for (int k = 0; k < inp.nR; k++)
            {
                rc[i][j] -= lambda[k]*(inp.w[i][j][k] 
                        + Omega*sqrt(sigma2[i])/sqrt(inp.nC));
                // rc[i][j] -= lambda[k]*(inp.w[i][j][k] - Omega*sigma2[i]);
                // rc[i][j] -= lambda[k]*inp.w[i][j][k];

                // ss += lambda[k]*Omega*sqrt(sigma2[i]);
            }
            
            // ss /= sqrt(inp.nC);
            // ss /= (inp.nC);
            // rc[i][j] -= ss;

            if (rc[i][j] > rcBest)
            {
                rcBest = rc[i][j];
                xL[i]  = j;
            }
            // cout << " ... j = " << j << " " << rc[i][j] << " vs " << rcBest << endl;
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
    for (int i = 0; i < inp.nC; i++)
        cout << " " << xL[i];
    cout << endl;

    //flagr << iter << "\t" << zL << endl;
    int aka;
    cin >> aka;
#endif
    return zL;
}

/// Solve the Lagrangean problem to optimality using SOCP and cplex.
// This function provides a real upper bound to the original problem.
double lagrange_step_SOCP(INSTANCE & inp, double Omega, double * sigma2, 
        int * xL, double * lambda, int iter, double ** rc)
{
    /* cout << "... lambdas :: ";
     * for (int k= 0; k < inp.nR; k++)
     *     cout << " " << lambda[k];
     * cout << endl; */


    // computing rc for fixing scheme
    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
        {
            rc[i][j] = inp.c[i][j];
    
            for (int k = 0; k < inp.nR; k++)
                rc[i][j] -= lambda[k]*inp.w[i][j][k];
        }


    IloEnv envL;
    IloModel modelL(envL);
    IloCplex cplexL(modelL);
    TwoD xL_ilo(envL, inp.nC);
    for (int i = 0; i < inp.nC; i++)
        xL_ilo[i] = IloNumVarArray(envL, inp.ri[i], 0, 1, ILOINT);

    IloObjective obj = IloAdd(modelL, IloMaximize(envL, 0));

    IloNumVar W_ilo(envL, 0.0, IloInfinity, IloNumVar::Float);

    IloExpr lhsRobust(envL);
    double rhs = 0.0;
    lhsRobust = -W_ilo*W_ilo;
    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
            lhsRobust += xL_ilo[i][j]*xL_ilo[i][j]*sigma2[i];

    modelL.add(lhsRobust <= rhs);

    for (int i = 0; i < inp.nC; i++)
    {
        IloExpr sum(envL);
        for (int j = 0; j < inp.ri[i]; j++)
            sum += xL_ilo[i][j];

        modelL.add(sum == 1.0);
    }

    // add objective function
    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
        {
            double coeff = inp.c[i][j];
            for (int k = 0; k < inp.nR; k++)
                coeff -= lambda[k]*inp.w[i][j][k];

            obj.setLinearCoef(xL_ilo[i][j], coeff);
        }

    double totLambda = 0.0;
    for (int k = 0; k < inp.nR; k++)
        totLambda += lambda[k];

    obj.setLinearCoef(W_ilo, -Omega*totLambda);
    modelL.add(obj);

    double zL = -1;
    try 
    {
        zL = solve_KNAP(modelL, cplexL, 99998, 0, 1); // call cplex
        // cout << "zL(cplex)= " << zL;
        if (zL != -1)
        {
            // add final term to lagragean function
            for (int k = 0; k < inp.nR; k++)
                zL += lambda[k]*inp.R[k];

            // cout << " Total = " << zL << endl;
            // get lagrangean solution

            for (int i = 0; i < inp.nC; i++)
                for (int j = 0; j < inp.ri[i]; j++)
                    if (cplexL.getValue(xL_ilo[i][j]) >= 1.0 - EPSI)
                        xL[i] = j;

        }

    }
    catch(IloException & e)
    {
        cout << "CLEX exception caught " << e << endl;
        return -1.0;
    }
    catch (...)
    {
        cout << "CPLEX unknown exception caught " << endl;
    }

    return zL;
}


/** This function is written for a max problem, i.e., we want to minimize the value
  of the lagrangean function. Every lagrangean value zL is an upper bound of the
  optimal solution.

Note: This is no longer true, since we are using an approximate Lagrangean.

*/
bool lambda_update(double * lambda, double & delta, int * xL, double & bestLagr, 
        double & worstLagr, double zL, double lb, int iter, double & best300, 
        double & start300, double * sigma2)
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
        double SSsigma = ZERO;
        for (int i = 0; i < inp.nC; i++)
        {
            s[k]    += (double)inp.w[i][xL[i]][k];
            SSsigma += sigma2[i];
        }
        s[k] += sqrt(SSsigma);
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
        // cout << "Initial lambda [ " << k << "] = " << lambda[k];
        lambda[k] +=  delta*(zL-lb)*s[k]/normS;
        if (lambda[k] < ZERO) lambda[k] = ZERO;
        // cout << " vs new = " << lambda[k] << endl;
    }

    return false;		// return not finished
}


