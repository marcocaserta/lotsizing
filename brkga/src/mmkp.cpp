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
/**! \file mmkp.cpp
  \brief Algorithm for the Multiple-choice Multidimensional Knapsack Problem

Author: Marco Caserta (marco dot caserta at ie dot edu)\n
Started: 05.02.15\n
Ended:   19.03.15\n 
Updated: January 2017 \n

Compile with: make
Execute with: ./bin/mmkp -f "input_data_file" 

To have a look at the flags that can be defined via command line, type:
> ./bin/mmkp -h

The options for the command line are defined in the file "options.cpp".

PROJECT STRUCTURE
=================

The project is composed of the following files:
- mmkp.cpp: The main file, which reads the instance data, calls the 
lagrangean phase, and implements the refinement procedure.
- lagrange.cpp: The implementation of the lagrangean relaxation with
subgradient optimization. 
- cplex.cpp : It manages the call to cplex (model creation+call to solver)
- options.cpp : Command-line parameters management.
*/

//#define M_DEBUG  // activate this to see the output for debugging
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cassert>


#include "timer.h"
#include "options.h"

#include "SampleDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"

#include "lagrange.h"
#include "cplex.h"

using namespace std;
char* _FILENAME;          // !< name of the instance file
double corridorWidthBase; // base value for corridor width
double time_limit;        // !< wall-clock time limit
int nSolBase;             // nr of sols to be found by cplex within corridor
double propFixed0;        // % of vars to be fixed to zero
double propFixed1;        // % of vars to be fixed to one
int add_oldSol_cut;       // boolean: cut out old feasible solution
int max_iter;             // max number lagrangean iterations
double Omega;             // parameter of robust formulation
double sigmaSq;           // parameter of individual item

ofstream flagr("lagrange.txt", ios::out);
ofstream fsol("solution.txt", ios::out);
const char * solutionFile = "optSolution.txt";
const char * solFileName;

/// Implementation of the basic data structure to hold the instance data
struct INSTANCE {
    int nR;       // number of resources (constraints)
    int nC;       // number of classes
    int * ri;     // number of items in each class i

    double  ** c; // obj function coefficients
    int *** w;    // knapsack constraints coefficients
    int   * R;    // r.h.s. values of the knapsack constraints
};

INSTANCE inp;     // INSTANCE DATA

double * sigma2;  // sigma^2 for robust formulation
int  * F;         // F[i] = 3 ==> var 3 in class "i" is set to 1
int **F0;         // F0[i][j] = 1 ==> var x_{ij} is set to 0
int **fF0;        // permanently fixed to 0 (based on gap)
double best_time; // time to best solution
int bestIter;     // lagr iteration in which best was found
double zBest;     // best feasible solution (lower bound)
double ubStar;    // best upper bound (lagr solution)
int * xBest;      // best solution overall

// definition of constant values
double INFTY = std::numeric_limits<double>::infinity();
const long _MAXRANDOM = 2147483647;
const double ZERO     = 0.0e0;
const double EPSI     = 0.00001;

double bestLagrHeur   = 0.0;
double zBestCE        = 0.0;
int * xBestCE;
timer tTime;            //!< Object clock to measure REAL and VIRTUAL (cpu) time

/*****************************************************************************/
/*                           FUNCTIONS                                       */
/*****************************************************************************/
void read_problem_data(INSTANCE & inp);
void print_options(INSTANCE inp);
void lagrangean_phase(INSTANCE & inp, IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        IloObjective & obj, IloNumVar & Q_ilo, double Omega, double sigma);
double get_first_lb(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, 
        int nSol, IloNumVar & Q_ilo, double Omega, double sigma);
void update_best(int * xLBest, int * xL, double & ubStar, double zL, double * lambda, double * lambdaBest);
double refine_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, int * xL, 
        double & cWidth, int & nSol, double & zBest, int * xIlo, int & iterImproved, 
        int iter, double ** rc, bool fix2zero, bool fix2one, IloRangeArray & cutSol,
        int lagrIter);
void find_best(int * F, double ** rc);
void find_best(int * F, double * delta);
void find_worst(int ** F0, double ** rc, double propFixed0);
void find_kworst(int ** F0, double ** rc, double propFixed0);
void cut_out_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xIlo, IloRangeArray & cutSol);
void add_z_cut(IloModel & model, IloCplex & cplex, TwoD & x_ilo, double zBest);

double compute_fitness_value(double corridorWidthBase, int nSolBase, double propFixed0, 
        double propFixed1, int add_z_cut);
int stopping_criteria(int iter, int max_iter, double time_limit);

double solve_robust_problem(INSTANCE & inp, IloModel & model, IloCplex & cplex, 
        TwoD & x_ilo, IloObjective & obj, IloNumVar & Q_ilo, 
        double Omega, double * sigma2, int * xBest, double &zBest);

double solve_nominal_problem(INSTANCE & inp, IloModel & model, IloCplex & cplex, 
        TwoD & x_ilo, IloObjective & obj, int * xNominal, double & zNominal,
        int max_time);

double computeSigmaSq(INSTANCE inp, double * sigma2);
void robust_lagrangean_phase(INSTANCE & inp, IloModel & model, IloCplex & cplex,
        TwoD & x_ilo, IloObjective & obj, IloNumVar & Q_ilo, double Omega, 
        double * sigma2);

// void dp_scheme(INSTANCE inp);
void quickSort(int* rc, int* x, int low, int high);
void quickSort(double * rc, int* x, int low, int high);
int is_feasible(const char * solutionFile, INSTANCE inp, int * xNominal, 
        double &zNominal);
int is_feasible(INSTANCE inp, int *x);
void robust_simulation(INSTANCE inp, int * xBest, int * xNominal, double zNominal,
        double zBest, int nRuns, double SSsigma);

double lagrangean_heuristic(INSTANCE inp, double ** rc);
void CE(INSTANCE inp, double ** rc, bool fix2zero);
double compute_fitness(INSTANCE inp, int * x);
void add_cover_inequality(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        int * x, double Omega, double * sigma2);
void add_cover_inequality2(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        int * x, double Omega, double * sigma2); 
void uplift_cover(INSTANCE inp, IloExpr & uplifted, TwoD & x_ilo, int isCover, 
        int * x, int group, int el);
int is_minimal_cover(INSTANCE & inp, int * x, double excess, int group, int el,
        int constr); 
void solve_bertsimas(INSTANCE & inp, IloModel &model, IloCplex & cplex, 
        TwoD & x_ilo, IloObjective & obj, double Omega, double sigma);
void solve_relaxed(IloModel model, TwoD x_ilo, double ** xLP, int * xL);
void compute_stability(INSTANCE & inp, double ** rc, double * delta, int * xL);
/*****************************************************************************/
/*                      END FUNCTIONS                                        */
/*****************************************************************************/





/************************ MAIN PROGRAM ******************************/
// [Lagrangean + Corridor Method + Fixing Schemes]
/************************ MAIN PROGRAM ******************************/
int main(int argc, char *argv[])
{

    srand(time(0));

    //freopen("debug.txt", "w", stdout); //!< redirect output to a file
    int err = parseOptions(argc, argv);
    if ( err != 0)
    { if (err != -1) cout << "Error argument " << err+1 << endl; exit(1); }

    tTime.resetTime();		      // start clock

    read_problem_data(inp); // read instance data
    print_options(inp);


    //=========================================================================
    // read solution from disk file (solutionFile)
    // activate this when robustness study is carried out (to compare the
    // number of infeasible sols reached by nominal vs. robust solution)
    int feasibilityCheck = 1;
    double zNominal      = 0.0;
    int * xNominal       = new int[inp.nC];
    xBestCE              = new int[inp.nC];
    if (feasibilityCheck)
    {
        cout << "Is Solution Feasible = " << is_feasible(solutionFile, inp, 
                xNominal, zNominal ) << endl;
        // exit(123);
    }
    //=========================================================================

    sigma2         = new double[inp.nC];
    double SSsigma = computeSigmaSq(inp, sigma2);

    tTime.resetTime();

    // vectors used to defined variables fixed to 0 and 1
    F   = new int[inp.nC];
    F0  = new int*[inp.nC];
    fF0 = new int*[inp.nC];
    for (int i = 0; i < inp.nC; i++)
    {
        F0[i]  = new int[inp.ri[i]];
        fF0[i] = new int[inp.ri[i]];
        for (int j = 0; j < inp.ri[i]; j++) 
            fF0[i][j] = 0;
    }
    zBest  = ZERO;
    ubStar = INFTY;
    xBest  = new int[inp.nC]; // store best solution

    // cplex stuff
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    TwoD x_ilo(env, inp.nC);
    for (int i = 0; i < inp.nC; i++)
        x_ilo[i] = IloNumVarArray(env, inp.ri[i], 0, 1, ILOINT);
    IloObjective obj = IloAdd(model, IloMaximize(env, 0));

    IloNumVar Q_ilo(env, 0.0, IloInfinity, IloNumVar::Float);

    // solve_bertsimas(inp, model, cplex, x_ilo, obj, Omega, sigma);

    // NOTE: Activate this to compare Nominal vs Robust
    // To avoid solving the nominal problem all the time, we store the solution
    // in a disk file and simply read the nominal solution from disk.
    // solve_nominal_problem(inp, model, cplex, x_ilo, obj, xNominal, zNominal,
            // 3600);
    // best_time    = tTime.elapsedTime(timer::REAL); // measure wall-clock time

    // solve_robust_problem(inp, model, cplex, x_ilo, obj, Q_ilo, Omega,
    // sigma2, xBest, zBest);
    // call lagrangean phase (subgradient optimization)
    // lagrangean_phase(inp, model, cplex, x_ilo, obj, Q_ilo, Omega, sigmaSq);
    // cout << "BEST LAGRANGEAN FOUND = " << bestLagrHeur << endl;
    robust_lagrangean_phase(inp, model, cplex, x_ilo, obj, Q_ilo, Omega,
    sigma2);

    fsol << _FILENAME << "\t" << zBest << "\t" << best_time << "\t" 
        << ubStar << "\t" << bestIter << "\t" << corridorWidthBase 
        << "\t" << nSolBase << "\t" << propFixed0 << "\t" 
        << propFixed1 << "\t" << add_oldSol_cut << endl;

    flagr.close();
    fsol.close();
    
    // compute hamming distrance between robust and determ solutions
    int count = 0;
    for (int i = 0; i < inp.nC; i++) 
        if (xNominal[i] == xBest[i])
            count++;

    cout << "Nr. common components " << count << endl;
    ofstream fCommon("common.txt", ios::app);
    fCommon << _FILENAME << "\t" << Omega << "\t" << count << endl;
    fCommon.close();

    // save solution in a named file (contains the name of the instance)
    string sBase = "solCplex-";
    int ll = strlen(_FILENAME);
    for (int i = 0; i < ll; i++) 
        if (_FILENAME[ll-1-i] == '/')
        {
            for (int j = 0; j < ll; j++) 
                sBase += _FILENAME[ll-i+j];
            break;
        }           

    solFileName = sBase.c_str();
    ofstream fName(solFileName, ios::out);
    fName << _FILENAME << "\t" << zNominal << "\t" << best_time << endl;
    fName.close();
    cout << "Solution saved in file " << solFileName << "." << endl;


    // irace minimizes (change sign since we are maximizing)
    ofstream fIrace ("solution4irace.txt", ios::out);
    fIrace << -zBest << endl;
    env.end();

#ifdef W_OPT
    cout << "Do you want to rewrite the opt solution on disk? (1 --> Yes) ";
    int answer;
    cin >> answer;
    if (answer == 1)
    {
        // save best solution to disk
        ofstream ffsol(solutionFile, ios::out);
        for (int i = 0; i < inp.nC; i++) 
            ffsol << xBest[i] << "\t";
        fifsol << endl;
    }
#endif
#ifdef W_ANALYSIS
    // analysis of robust solution for I07 --> zBest (nominal) = 24595
    int nRuns       = 1000;
    robust_simulation(inp, xBest, xNominal, zNominal, zBest, nRuns, SSsigma);
#endif
    return zBest;		// return value to brkGA for 

}
/************************ MAIN PROGRAM ******************************/
// [Lagrangean + Corridor Method + Fixing Schemes]
/************************ MAIN PROGRAM ******************************/



/************************ FUNCTIONS ******************************/
/// Functions
/************************ FUNCTIONS ******************************/
void solve_bertsimas(INSTANCE & inp, IloModel &model, IloCplex & cplex, 
        TwoD & x_ilo, IloObjective & obj, double Omega, double sigma)
{

    IloEnv env = model.getEnv();
    int totEl = inp.nC*inp.ri[0];
    for (int u = 1; u < totEl; u++) 
    {
        IloModel mm(env);
        IloCplex cplex(mm);
        cout << "U = " << u << endl;
        cout << "Nr constraints is " << cplex.getNrows() << endl;
        double statusBin = solve_KNAP(mm, cplex, 99999, 0, 10000); // call cplex

        cout << "... z = " << statusBin << endl;
        for (int i = 0; i < inp.nC; i++) 
            for (int j = 0; j < inp.ri[i]; j++) 
                if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
                    cout << " " << "[" << i << "," << j << "] ";
        cout << endl;



    }


}



void read_problem_data(INSTANCE & inp)
{

    int temp;

    ifstream fdata(_FILENAME, ios::in);
    if (!fdata)
    {
        cerr << "Cannot open file " << _FILENAME << endl; exit(1);
    }

    fdata >> inp.nC >> temp >> inp.nR;

    // all the classes have the same nr of elements
    inp.ri = new int[inp.nC];
    inp.c  = new double*[inp.nC];
    inp.w  = new int**[inp.nC];
    for (int i = 0; i < inp.nC; i++)
    {
        inp.ri[i] = temp;
        inp.c[i] = new double[inp.ri[i]];
        inp.w[i] = new int*[inp.ri[i]];
        for (int j = 0; j < inp.ri[i]; j++)
            inp.w[i][j] = new int[inp.nR];
    }

    // read rhs of knapsack constraints
    inp.R = new int[inp.nR];
    for (int k = 0; k < inp.nR; k++)
        fdata >> inp.R[k];

    // read data for each class
    for (int i = 0; i < inp.nC; i++)
    {
        fdata >> temp;
        assert(temp == (i+1));

        for (int j = 0; j < inp.ri[i]; j++)
        {
            fdata >> inp.c[i][j];
            for (int k = 0; k < inp.nR; k++)
                fdata >> inp.w[i][j][k];
        }
    }
#ifdef M_DEBUG
    // analysis of first constraint 
    for (int i = 0; i < inp.nC; i++) 
    {
        cout << "CLASS " << i << " :: " ;
        for (int j = 0; j < inp.ri[j]; j++) 
        {
            cout << "\t" << inp.c[i][j]/(double)inp.w[i][j][0];
        }
        cout << endl;
    }
#endif
}

int is_feasible(INSTANCE inp, int * x)
{

    for (int k = 0; k < inp.nR; k++) 
    {
        int totR = 0;
        for (int i = 0; i < inp.nC; i++) 
            totR += inp.w[i][x[i]][k];
        if (totR > inp.R[k]) 
            return 0;
    }
    return 1; 
}

// Check whether a solution saved on disk is feasible
int is_feasible(const char * solutionFile, INSTANCE inp, int * xNominal, 
        double & zNominal)
{
    ifstream fsol(solutionFile, ios::in);
    // int * x = new int[inp.nC];
    for (int i = 0; i < inp.nC; i++) 
        fsol >> xNominal[i];
    fsol.close();

    int isFeasible = 1;
    // double z = 0.0;
    zNominal = 0.0;
    // check obj function value
    for (int i = 0; i < inp.nC; i++) 
        zNominal += inp.c[i][xNominal[i]];
    cout << "z*_nom = " << zNominal  << endl; 

    // check feasible solution
    for (int k = 0; k < inp.nR; k++) 
    {
        int totR = 0;
        for (int i = 0; i < inp.nC; i++) 
            totR += inp.w[i][xNominal[i]][k];

        cout << "lhs( " << k << ") = " << totR << " vs rhs = "<< inp.R[k];
        if (totR > inp.R[k])
        {
            isFeasible = 0;
            cout << " *** " << endl;
        }
        else
            cout << "  v  " << endl;
    }

    return isFeasible;
}


double computeSigmaSq(INSTANCE inp, double * sigma2)
{
    // some statistics
    double * avg = new double[inp.nC];
    double tot = 0.0;

    double SSsigma = 0.0;
    for (int i = 0; i < inp.nC; i++)
    {
        // avg profit per group
        avg[i]= 0.0;
        for (int j = 0; j < inp.ri[i]; j++)
        {
            avg[i] += inp.c[i][j];
        }
        avg[i] /= (double)inp.ri[i];
        tot += avg[i];

        // cout << "Avg Profit Class " << i << " = " << avg[i] << endl;
    }

    for (int i = 0; i < inp.nC; i++)
    {
        sigma2[i] =  0.8 + avg[i]/tot;
        sigma2[i] = 1.0;
        SSsigma += sigma2[i];
        // cout << "SigmaSq[" << i << "] = " << sigma2[i] << endl;
    }
    return SSsigma;
}

/// Print instance info and algorithmic parameters.
void print_options(INSTANCE inp)
{
    cout << "-------------------------------------" << endl;
    cout << "- INSTANCE DATA : " << endl;
    cout << "-------------------------------------" << endl;
    cout << " Nr. Classes\t \t :: " << inp.nC << endl;
    cout << " Nr. Resources\t \t :: " << inp.nR << endl;
    cout << " Nr. Items per class\t :: " << inp.ri[0] << endl;
    cout << "-------------------------------------" <<  endl << endl;   
    cout << "-------------------------------------" << endl;
    cout << "- ALGORITHMIC PARAMETERS : " << endl;
    cout << "-------------------------------------" << endl;
    cout << " Time Limit \t \t :: " << time_limit << endl;
    cout << " Iter Limit \t \t :: " << max_iter << endl;
    cout << " Corridor Base\t \t :: " << corridorWidthBase << endl;
    cout << " Base Sol. Nr.\t \t :: " << nSolBase << endl;
    cout << " % Fix to Zero \t \t :: " << propFixed0 << endl;
    cout << " % Fix to One \t \t :: " << propFixed1 << endl;
    cout << " With Cut Sol \t \t :: " << add_oldSol_cut << endl;
    cout << " Omega \t \t \t :: " << Omega << endl;
    cout << " Sigma^2 \t \t :: " << sigmaSq << endl;
    cout << "-------------------------------------" << endl;
    cout << "-------------------------------------" << endl << endl;;


}

/// Define stopping criteria of the algorithm
int stopping_criteria(int iter, int max_iter, double time_limit)
{
    return ( (tTime.elapsedTime(timer::REAL) >= time_limit) ||
            (iter >= max_iter) );
}

double lagrangean_heuristic(INSTANCE inp, double ** rc)
{
    int ** index = new int*[inp.nC];
    std::vector < std::vector <int> > heap;

    double z = 0.0;
    int * x = new int[inp.nC];
    for (int i = 0; i < inp.nC; i++) 
    {
        index[i] = new int[inp.ri[i]];
        double maxC = 0.0;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            index[i][j] = j;
            if (inp.c[i][j] > maxC) 
            {
                maxC = inp.c[i][j];
                x[i] = j;
            }
        }  

        z += inp.c[i][x[i]];
    }
    // check feasibility
    if (is_feasible(inp, x))
    {
        cout << "OPT SOL FOUND " << endl;
        exit(124);
    }


    for (int i = 0; i < inp.nC; i++) 
        quickSort(rc[i], index[i], 0, inp.ri[i]-1);

    // copy index into heap 
    for (int i = 0; i < inp.nC; i++) 
    {
        std::vector <int> aux;
        for (int j = 0; j < inp.ri[i]; j++) 
            if (index[i][j] != x[i])    
                aux.push_back(index[i][j]);
        heap.push_back(aux);
    }
#ifdef M_DEBUG
    cout << "The Heap is :: " << endl;
    for (int i = 0; i < inp.nC; i++) 
    {
        cout << "class [" << i << "]  ";
        for (int j = 0; j < inp.ri[i]; j++) 
            cout << " "  << index[i][j];
        cout << endl; 
        for (std::vector<int>::iterator it = heap[i].begin(); it != heap[i].end(); ++it) 
        {
            cout << " "  << *it;
        }

        cout << endl;
    } 
#endif

    // repeat cycle until feasibility is reached
    while (!is_feasible(inp, x))
    {
        double delta = -INFTY;
        int istar = -1;
        for (int i = 0; i < inp.nC; i++) 
        {
            int cc = heap[i].back();
            if (rc[i][cc] > delta)
            {
                delta = rc[i][cc];
                istar = i;
            }
        }
        // remove candidate from the heap and update solution
        z -= inp.c[istar][x[istar]];
        x[istar] = heap[istar].back();
        z += inp.c[istar][x[istar]];
        heap[istar].pop_back();

    }
    cout << "Lagragean Heuristic Solution is " << z << endl;
    if (z > bestLagrHeur)
        bestLagrHeur = z;
}


void CE(INSTANCE inp, double ** rc, bool fix2zero)
{

    double prop0CE = 0.30;
    int N          = 1000;
    double rho     = 0.2;
    double alpha   = 0.8;
    int maxCE      = 500;
    int iterCE     = 0;
    int nFixed     = 0;

    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
            F0[i][j] = 0;

    fix2zero = 0;
    if (fix2zero)
    {

        find_worst(F0, rc, prop0CE);

#ifdef M_DEBUG
        cout << "List of columns fixed to zero (" << prop0CE << ") :: " << endl;
        for (int i = 0; i < inp.nC; i++) 
        {
            cout << "class (" << i << ") :: ";
            for (int j = 0; j < inp.ri[i]; j++) 
                cout << " " << rc[i][j];
            cout << endl;
            for (int j = 0; j < inp.ri[i]; j++) 
            {
                if (F0[i][j] == 1)
                    cout << " " << j;
            }
            cout << endl;
        }
#endif
    }

    // Note: rc can be both positive and negative
    // We shift everything by the minimum value and, this way,
    // we effectively prevent from selecting the item with min rc
    int ** xCE  = new int*[N];
    for (int i  = 0; i < N; i++)
        xCE[i]  = new int[inp.nC];

    double * z  = new double[N];
    int * index  = new int[N];
    double ** p = new double*[inp.nC];
    for (int i  = 0; i < inp.nC; i++)
    {
        if (fix2zero)
            nFixed = ceil(prop0CE*(double)inp.ri[i]);	
        p[i] = new double[inp.ri[i]];
        double totL = 0.0;
        double minL = INFTY;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            if (F0[i][j] == 0)
            {
                if (rc[i][j] < minL)
                    minL = rc[i][j];
                totL += rc[i][j];
            }
        }
        // shift rc if the min is negative
        if (minL < 0.0)
            totL += (inp.ri[i]-nFixed)*minL;
        else
            minL = 0.0;

        for (int j = 0; j < inp.ri[i]; j++) 
        {
            if (F0[i][j] == 1)
                p[i][j] = 0.0;
            else
                p[i][j] = (minL + rc[i][j]) / totL;

        }
    }
#ifdef M_DEBUG
    // check probabilities
    for (int i  = 0; i < inp.nC; i++)
    {
        cout << "c( " << i << ") == ";
        double totP = 0.0;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            cout << " " << p[i][j];
            totP += p[i][j];
        }
        cout << " \t TOT = " << totP << endl;
    }
#endif


    for (int kk = 0; kk < maxCE; kk++) 
    {

        int generate = 0;
        for (int k = 0; k < N; k++) 
        {
            do
            {
                generate = 0;
                for (int i = 0; i < inp.nC; i++) 
                {

                    // select item for this class
                    double rr   = (double)(rand()+1)/((double)(RAND_MAX)+1.0);
                    int j       = 0;
                    double totP = 0.0;
                    while (totP + p[i][j] < rr)
                    {
                        totP += p[i][j];
                        j++;
                    }
                    xCE[k][i] = j;
                }
                // if feasible, compute fitness value
                if (is_feasible(inp, xCE[k]))
                    z[k] = compute_fitness(inp, xCE[k]);
                else
                {
                    generate += 1;
                    if (generate > 100)
                    {
                        cout << "PROBLEMS with CE here ... " << endl;
                        int abc;
                        cin >> abc;
                    }

                }
            }
            while (generate > 0);
        }    

        // sort and find quantiles (ascending order)
        for (int i = 0; i < N; i++) 
            index[i]= i;
        quickSort(z, index, 0,N-1); 

        int limit = (int)(rho*(double)N) ;

        if (z[index[N-1]] > zBestCE)
        {
            zBestCE = z[index[N-1]];
            for (int i = 0; i < inp.nC; i++) 
                xBestCE[i] = xCE[index[N-1]][i];

            cout << "... best CE[" << iterCE << "] = " << zBestCE << endl;
        }

        // update probs
        for (int i = 0; i < inp.nC; i++) 
        {
            int * y = new int[inp.ri[i]];
            for (int j = 0; j < inp.ri[i]; j++) 
                y[j] = 0;

            for (int k = N-limit; k < N; k++) 
                y[xCE[index[k]][i]]++;

            // update probs for this class
            for (int j = 0; j < inp.ri[i]; j++) 
                p[i][j] = alpha*p[i][j] + (1.0-alpha)*(double)y[j]/(double)limit;

            double totP = 0.0;
            for (int j = 0; j < inp.ri[i]; j++) 
                totP += p[i][j];

            assert(totP > 1.0-EPSI & totP < 1.0+EPSI);
        }
        iterCE++;
    } // end CE cycle
}


double compute_fitness(INSTANCE inp, int * x)
{
    double z = 0.0;
    for (int i = 0; i < inp.nC; i++) 
        z += inp.c[i][x[i]];
    return z;
}


//****************************************************************************/
/// LAGRANGEAN PHASE                                                         */
//****************************************************************************/ 
/*  The implementation of the lagrangean relaxation coupled with subgradient
    optimization is defined in the file "lagrange.cpp." 

    The Lagrangean relaxation method produces both UPPER and LOWER bounds:
    - at each lagrange iteration, i.e., for each set of lagrange multipliers, we
    solve to optimality the corresponding relaxed problem (we obtain a valid
    UPPER BOUND zL)
    - starting from the infeasible lagrangean solution, we apply the corridor
    method to look for feasible solutions in the neighborhood of the lagrangean
    solution. Any feasible solution obtained during this refining procedure is
    a valid LOWER BOUND.

    The key idea is:
    - we repeat the lagrangean cycle 3 times (each time with a different
    starting set of lagrange multipliers)
    - each cycle is divided into two main parts:
    i.  the first, e.g., 100 iterations are used to quickly get near-optimal
    lagrange multipliers
    ii. once near-optimal multipliers are available, we periodically call a 
    refining procedure (based on the corridor method) to achieve better 
    lower bounds.
 **/
void lagrangean_phase(INSTANCE & inp, IloModel & model, IloCplex & cplex, 
        TwoD & x_ilo, IloObjective & obj, IloNumVar & Q_ilo, double Omega, 
        double sigma)
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

    zBest = get_first_lb(model, cplex, x_ilo, obj, 1, Q_ilo, Omega, sigma);// nr feasible solutions to be found by cplex
    cout << "First feasible solution = " << zBest << endl;

    lambda_initialization(lambda);

    double ** rc = new double*[inp.nC]; // lagrangean reduced costs
    for (int i = 0; i < inp.nC; i++)
        rc[i] = new double[inp.ri[i]];

    int lagrIter = 0;
    int iter     = 0;
    while (lagrIter < 3 && !stopping_criteria(iter, max_iter, time_limit))
    {

        double zL, bestLagr, worstLagr;

        double cWidth    = corridorWidthBase; // <---------- used to be 0.8
        int nSol         = nSolBase;          // nr of solutions for cplex
        double start300  = INFTY;
        double best300   = INFTY;
        int iterImproved = 200;
        int Freq         = 40;
        double delta     = 0.1;
        bool stopping    = false;
        bool fix2zero    = false;
        bool fix2one     = false;

        while (!stopping_criteria(iter, max_iter, time_limit) && !stopping)
        {
            // note: zL in an "upper bound" of the optimal value
            zL = lagrange_step(lambda, xL, iter, rc);
            if ((iter % (2*Freq)) == 0)
                cout << "LB( " << iter << " ) = " << zL << endl;
            if (zL < ubStar)
                update_best(xLBest, xL, ubStar, zL, lambda, lambdaBest);

            // if (iter > 51)
            // CE(inp, rc, fix2zero);
#ifdef AAA
            int count0 = 0;
            for (int i = 0; i < inp.nC; i++) 
                for (int j = 0; j < inp.ri[i]; j++) 
                    if (fF0[i][j] == 0)
                        if (zL - rc[i][xL[i]] + rc[i][j] <= zBest)
                        {
                            fF0[i][j] = 1;
                            count0++;
                            model.add(x_ilo[i][j] == ZERO);
                        }
#endif
            stopping = lambda_update(lambda, delta, xL, bestLagr, worstLagr, zL, 
                    zBest, iter, best300, start300);


            // if (iter > 199) 	// start fixing schemes
            if (iter > 50) 	// start fixing schemes
            {
                fix2one = true;
                fix2zero = true;
            }


            // refine lagrange solution to obtain a valid lower bound
            if (iter > 100 && (iter % Freq) == 0)
            {
                // lagrangean_heuristic(inp, rc);
                lb = refine_solution(model, cplex, x_ilo, obj, xL, cWidth, nSol,
                        zBest, xIlo, iterImproved, iter, rc, fix2zero, fix2one,
                        cutSol, lagrIter);
            }

            if ((iter > 250) && (iter - iterImproved) > 100)
            {
                nSol++;
                cWidth       *= 0.9;
                iterImproved  = iter;
                Freq         /= 2;
                if (Freq < 5) 
                    Freq = 5;
            }

            iter++;
        } // end lagrangean cycle

        cout << "** ** ** SUMMARY ** ** ** " << endl;
        cout << "** ** ** Lagrangean Iteration Nr. " << lagrIter 
            << " :: Best LB = " << zBest << endl;
        cout << setw(49) << " :: Best UB = " << ubStar << endl;

        model.remove(cutSol);


        lagrIter++;
        iter = 0;		// reset counter

        if (lagrIter < 3)
            lambda_random_perturbation(lambda, lambdaBest);

    } // END OF LAGRANGEAN CYCLE -- repeat 3 times

    // store best solution found (save to disk to check feasibility)
    for (int i = 0; i < inp.nC; i++) 
        xBest[i] = xIlo[i];

}



/// Update best solution
void update_best(int * xLBest, int * xL, double & ubStar, double zL, 
        double * lambda, double * lambdaBest)
{
    ubStar = zL;
    for (int i = 0; i < inp.nC; i++)
        xLBest[i] = xL[i];

    for (int k = 0; k < inp.nR; k++)
        lambdaBest[k] = lambda[k];
}

/// Get first feasible solution
double get_first_lb(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, 
        int nSol, IloNumVar & Q_ilo, double Omega, double sigma)
{

    defineModel(model, cplex, x_ilo, obj); // get the first cplex solution

    //defineRobustModel(model, cplex, x_ilo, obj, Q_ilo, Omega, sigma);

    // cycle for all the possible values of u in W	
    /*     for (int cc = 0; cc < inp.nC; cc++)
     *     {
     *
     *
     *         double u = (double)(cc+1)*sigma;
     *         cout << "** ** ** U[" << cc+1 <<"] = " << u << " ** ** ** " << endl;
     *         defineRobustDet(model, cplex, x_ilo, obj, Omega, sigma, u);
     *         double statusBin = solve_KNAP(model, cplex, 99999, 4, 10000); // call cplex
     *         int abc;
     *         cin >> abc;
     *     } */


    double statusBin = solve_KNAP(model, cplex, nSol, 0, 10); // call cplex

    // double statusBin = solve_KNAP(model, cplex, 9999, 4, 10000); call cplex

    // get solution
    /*     int * xIlo = new int[inp.nC];
     *
     *     for (int i = 0; i < inp.nC; i++)
     *         for (int j = 0; j < inp.ri[i]; j++)
     *             if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
     *                 xIlo[i] = j;
     *     for (int i = 0; i < inp.nC; i++)
     *         cout << "x[" << i <<  "] = " << xIlo[i] << endl;
     *
     *     exit(145); */
    return statusBin;
}

// Solve LP relaxation, with the addition of a "diversity constraint"
void solve_relaxed(IloModel model, TwoD x_ilo, double ** xLP, int * xL)
{
    double rhs = (double)inp.nC*propFixed1;


    IloEnv env = model.getEnv();
    IloModel relax(env);
    relax.add(model);

    // now add constraint to force a different LP solution
    IloExpr lhs(env);
    for (int i = 0; i < inp.nC; i++) 
        lhs += x_ilo[i][xL[i]]; 

    relax.add(lhs <= rhs);

    for (int i = 0; i < inp.nC; i++) 
        relax.add(IloConversion(env, x_ilo[i], ILOFLOAT));

    IloCplex cplex(relax);

    cplex.setParam(IloCplex::SimDisplay, 0); 
    cplex.solve();
    cout << "LP Solution Status " << cplex.getStatus() << " with z = " 
        << cplex.getObjValue() << endl;
    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
            xLP[i][j] = cplex.getValue(x_ilo[i][j]);

    // for (int i = 0; i < inp.nC; i++)
    // {
    // for (int j = 0; j < inp.ri[i]; j++)
    // cout << " x["<<i<<","<<j<<"] = " << xLP[i][j];
    // cout << endl;
    // }
}



/// Refine lagrange solution to obtain a valid lower bound (i.e., a feasible solution)
/*  Use the Corridor Method, with two fixing scheme, to solve the constrained
 *  MMKP. The constrained formulation makes use of the following:
 *  - corridor constraint, i.e., a distance around the lagrangean solution
 *  - a fixing-to-zero scheme, based on the variables lagrangean costs
 *  - a fixing-to-one scheme. Since this is quite delicate, we explored
 *    different variants, based on the lagrangean costs, the LP relaxation,
 *    etc.
 * */
double refine_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        IloObjective & obj, int * xL, double & cWidth, int & nSol, 
        double & zBest, int * xIlo, int & iterImproved, int iter, 
        double ** rc, bool fix2zero, bool fix2one, IloRangeArray & cutSol, 
        int lagrIter)
{
    IloEnv env = model.getEnv();

    int    nFixed1;
    double width1;
    std::vector <int> fixed2OneLP;

    // FIXING SCHEMES (keep them separate, to activate them at different times)
    //=========================================================================
    // if fix-to-zero is activated, select strategy
    //=========================================================================
    if (fix2zero)
    {

        for (int i = 0; i < inp.nC; i++)
            for (int j = 0; j < inp.ri[i]; j++)
                F0[i][j] = 0;

        // find_worst(F0, rc, propFixed0);
        find_kworst(F0, rc, propFixed0);    
    }

    int fixLC     = 1; // activate fix-to-one using Lagrangean costs
    int fixLP     = 0; // activate fix-to-one using LP relaxation
    int fixStable = 0; // activate fix-to-one using stability criterion
    // if fix-to-one scheme is activated, select strategy
    if (fix2one)
    {
        nFixed1 = ceil(propFixed1*(double)inp.nC); // used to be 0.25
        width1  = ceil(1.25*cWidth*(double)nFixed1);
        if (width1 > nFixed1)
            width1 = nFixed1;

        if (fixLC == 1)
        {
            for (int i = 0; i < inp.nC; i++)
                F[i] = -1;

            for (int k = 0; k < nFixed1; k++)
                find_best(F, rc);
        }
        if (fixStable == 1)
        {
            double * delta = new double[inp.nC];
            compute_stability(inp, rc, delta, xL);
            for (int i = 0; i < inp.nC; i++)
                F[i] = -1;

            for (int k = 0; k < nFixed1; k++)
                find_best(F, delta);

        }

        if (fixLP == 1)
        {
            // SOLVE LP to develop fixing scheme
            double ** xLP = new double *[inp.nC];
            for (int i = 0; i < inp.nC; i++) 
                xLP[i] = new double[inp.ri[i]];

            solve_relaxed(model, x_ilo, xLP, xL);
            int * ones = new int[inp.nC];
            for (int i = 0; i < inp.nC; i++) 
            {
                for (int j = 0; j < inp.ri[i]; j++) 
                {
                    if (j == xL[i]) // examine variable
                    {
                        if (xLP[i][j] >= 1.0 - EPSI)
                            ones[i] = j;
                        else
                            ones[i] = -1;
                        // cout << " * xLP[" << i << "," << j << "] = " << xLP[i][j];
                    }
                }
                // cout << endl;
            }
            int totOnes = 0;
            for (int i = 0; i < inp.nC; i++) 
            {
                if (ones[i] != -1)
                {
                    fixed2OneLP.push_back(i);
                    // cout << " " << i;
                    totOnes++;
                }
            }
            // cout << " ... total = " << totOnes << endl;
        }
    }
    //=================== end fixing schemes ==================================

    //=========================================================================
    // add corridor constraint
    //========================================================================= 
    IloExpr lhsCorridor(env);
    IloRangeArray neighborhood(env);
    double rhs = cWidth*(double)(inp.nC);

    for (int i = 0; i < inp.nC; i++)
        lhsCorridor += x_ilo[i][xL[i]];

    neighborhood.add(lhsCorridor >= rhs);
    model.add(neighborhood);
    //=========================================================================
    // END corridor constraint
    //========================================================================= 



    IloExpr         lhs0(env);
    IloExpr         lhs1(env);
    IloRangeArray   constraint0(env);
    IloRangeArray   constraint1(env);

    //=========================================================================
    // Fix-to-zero :: The matrix F0 has a flag (1/0) to indicate whether a var
    // x_{ij} should be fixed (1) or not (0)
    //=========================================================================
    int count0 = 0;
    if (fix2zero)
    {
        for (int i = 0; i < inp.nC; i++)
            for (int j = 0; j < inp.ri[i]; j++)
                if (F0[i][j] == 1 && fF0[i][j] == 0)
                {
                    // lhs0 += x_ilo[i][j];
                    constraint0.add(x_ilo[i][j] == ZERO);
                    count0++;
                }
        //fix_to_zero.add(lhs0 == ZERO); // an old way of fixing to zero
        model.add(constraint0);
        cout << " ... fixed to zero via LC = " << count0 << endl;
    }
    //=========================================================================


    //=========================================================================
    // Fix-to-one :: The matrix F has a flag (1/0) to indicate whether a var
    // x_{ij} should be fixed (1) or not (0)
    //=========================================================================
    if (fix2one)
    {
        if (fixStable == 1)
        {
            int count1 = 0;
            for (int i = 0; i < inp.nC; i++) 
            {
                if (F[i] != -1)
                {
                    count1++;
                    constraint1.add(x_ilo[i][xL[i]] == 1.0);
                }
            }
            cout << " ... fixed to one via Stability = " << count1 << endl;
        }

        if (fixLC == 1)
        {
            int count1 = 0;
            for (int i = 0; i < inp.nC; i++)
                if (F[i] != -1)
                {
                    count1++; 
                    lhs1 += x_ilo[i][F[i]];
                }

            constraint1.add(lhs1 >= width1);
            model.add(constraint1);
            cout << " ... fixed to one via LC = " << count1 << endl;
        }
        if (fixLP == 1)
        {
            std::vector<int>::iterator it; 
            for (it = fixed2OneLP.begin(); it != fixed2OneLP.end(); ++it) 
            {
                // cout << "Fixing to 1 class " << *it << endl;
                constraint1.add(x_ilo[*it][xL[*it]] == 1.0);
            }
            model.add(constraint1);
            cout << " ... fixed to one via LP = " << fixed2OneLP.size() << endl;
        }
    }

    //cout << "Nr constraints after adding ALL CUTS is " << cplex.getNrows() << endl;

    //=========================================================================
    // Solve CORRIDOR MMKP [corridor + fix2zero + fix2one]
    //=========================================================================
    double statusBin = solve_KNAP(model, cplex, nSol, 0, 10);	// call cplex

    if (statusBin >= -1.0+EPSI) 
    {
#ifdef M_DEBUG
        cout << " ... repaired solution is " << statusBin << "(z* = " 
            <<  zBest << ") with cplex status :: " << cplex.getStatus() 
            << ". Params: [ zeros = " << count0 << ", ones = " << width1 << "/" 
            << nFixed1 << "]" << endl;
#endif

        int * xRepaired = new int[inp.nC];
        get_cplex_sol(model, cplex, x_ilo, xRepaired);
        //cut_out_solution(model, cplex, x_ilo, xRepaired, cutSol);
        // cout << "Nr. Constraints after cutting out solutions is : " 
        // << cplex.getNrows() << endl;

        // Update when a new LB has been found
        if (statusBin > (zBest+EPSI))
        {
            zBest        = statusBin;
            best_time    = tTime.elapsedTime(timer::REAL); // measure wall-clock time
            bestIter     = lagrIter;
            iterImproved = iter; // save iteration of last improvement
            cWidth       = corridorWidthBase;
            nSol         = nSolBase;
            for (int i = 0; i < inp.nC; i++) 
                xIlo[i] = xRepaired[i];

            if(add_oldSol_cut)
                add_z_cut(model, cplex, x_ilo, statusBin);

            cout << "... improved best lb(" << iter << ") = " << zBest 
                << " found after " << best_time << " seconds." << endl;
#ifdef M_DEBUG
            cout << "Nr constraints before adding CUTS is " << cplex.getNrows() << endl;
            add_cover_inequality2(model, cplex, x_ilo, xRepaired, Omega, sigma2);
            cout << "Nr constraints after adding CUTS is " << cplex.getNrows() << endl;
#endif

        }
    }

    //========================================================================= 
    // Eliminate fixing schemes and corridor
    //========================================================================= 
    if (fix2one) // remove fix to one constraints
    {
        model.remove(constraint1);
        constraint1.end();
        lhs1.end();
    }
    if (fix2zero) // remove fix to zero constraints
    {
        model.remove(constraint0);
        constraint0.end();
        lhs0.end();
    }

    // remove neighborhood constraint
    model.remove(neighborhood);
    lhsCorridor.end();
    //cout << "Nr constraints after removing ALL CUTS is " << cplex.getNrows() 
    //<< endl;
    //========================================================================= 

    return statusBin;
}   // END OF refine_solution() -- corridor method


void find_best(int * F, double * delta)
{
    double bestRC = -INFTY;
    int group = -1;
    for (int i = 0; i < inp.nC; i++)
    {
        if (F[i] != -1) continue;
        if (delta[i] > bestRC)
        {
            bestRC = delta[i];
            group  = i;
        }
    }
    if (group != -1)
        F[group] = 1;
    else
    {
        cout << "Unable to fix variable " << endl;
        exit(144);
    }
}


void find_best(int * F, double ** rc)
{
    double bestRC = -INFTY;
    int group = -1;
    int item  = -1;
    for (int i = 0; i < inp.nC; i++)
    {
        if (F[i] != -1) continue;
        for (int j = 0; j < inp.ri[i]; j++)
            if (rc[i][j] > bestRC)
            {
                bestRC = rc[i][j];
                group  = i;
                item   = j;
            }
    }
    if (group != -1)
        F[group] = item;
    else
    {
        cout << "Unable to fix variable " << endl;
        exit(144);
    }
}

void find_kworst(int ** F0, double ** rc, double propFixed0)
{

    for (int i = 0; i < inp.nC; i++) 
    {
        int nFixed  = ceil(propFixed0*(double)inp.ri[i]);	
        int * index = new int[inp.ri[i]];
        for (int j = 0; j < inp.ri[i]; j++) 
            index[j] = j;
        
        quickSort(rc[i], index, 0, inp.ri[i]-1);
        int counter = 0;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            int pos = index[j];
            if (fF0[i][pos] == 1 || F0[i][pos] == 1) continue;
            F0[i][pos] = 1;
            counter++;
            if (counter == nFixed) break;
        }
#ifdef M_DEBUG
        cout << "Fixing to 0 for class " << i << endl;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            cout << "rc["<<j<<"] = " << rc[i][j] << " .. fF = " << fF0[i][j] << " F = " << F0[i][j] << endl;
        }
#endif        
    } 
}


/// Find variables with worst lagrangean costs
void find_worst(int ** F0, double ** rc, double propFixed0)
{

    for (int i = 0; i < inp.nC; i++)
    {
        int nFixed = ceil(propFixed0*(double)inp.ri[i]);	
        //cout << "for module " << i << " we fix " << nFixed << " elements " << endl;

        for (int counter = 0; counter < nFixed; counter++)
        {
            double worstRC = INFTY;
            int item       = -1;
            for (int j = 0; j < inp.ri[i]; j++)
            {
                if (F0[i][j] == 1) continue;

                if (rc[i][j] < worstRC)
                {
                    worstRC = rc[i][j];
                    item    = j;
                }
            }
            if (item != -1)
                F0[i][item] = 1;
            else
            {
                break;
                cout << "something wrong here ... " << endl;
                exit(123);
            }
        }
    }
}



// Second version
void add_cover_inequality2(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        int * x, double Omega, double * sigma2)
{
    IloEnv env = model.getEnv();
#ifdef M_DEBUG
    cout << "Feasible solution is :: ";
    for (int i = 0; i < inp.nC; i++) 
        cout << " " << x[i];
    cout << endl;
#endif
    double * maxW  = new double[inp.nR];
    double * W     = new double[inp.nR];
    double * slack = new double[inp.nR];
    for (int k = 0; k < inp.nR; k++) 
    {
        W[k] = 0.0;
        maxW[k] = 0.0;
        for (int i = 0; i < inp.nC; i++)
        {
            W[k] += inp.w[i][x[i]][k];
            if (inp.w[i][x[i]][k] > maxW[k])
                maxW[k] = inp.w[i][x[i]][k];
        }
        // cout << "Checking constraint " << k << " LHS = " << W[k] + Omega*sqrt(Q[k])
        // << " vs RHS = " << inp.R[k] << endl;
        slack[k] = inp.R[k] - W[k];
        assert(slack[k] >= 0.0);
    }

    int nCovers = 0;
    int isCover = -1;
    for (int i = 0; i < inp.nC; i++) 
    {
        // cout << "COVERS for CLASS " << i << endl;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            isCover = -1;
            if (x[i] == j) continue;
            // cout << " ... trying with item " << j << endl;
            for (int k = 0; k < inp.nR; k++) 
            {
                double delta = -inp.w[i][x[i]][k] + inp.w[i][j][k];
                if (delta > slack[k])
                {
                    // is this cover minimal?
                    if (is_minimal_cover(inp, x, delta-slack[k], i, j, k) == 1)
                    {
                        // cout << "-----> constraint " << k << " is violated " << endl;
                        isCover = k; // constraint for which the cover is built
                        break;

                    }
                }       
            }
            // if infeasible, it is a cover
            if (isCover != -1)
            {
                // add cover
                nCovers++;
                IloExpr lhs(env);
                for (int l = 0; l < inp.nC; l++) 
                {
                    if (l == i) continue;
                    lhs += x_ilo[l][x[l]];
                    // cout << " " << x[l];
                }
                // cout << " + x[" << i << "," << j << " ] " << endl;
                lhs += x_ilo[i][j];

                // uplift the cover
                IloExpr uplifted(env);
                uplift_cover(inp, uplifted, x_ilo, isCover, x, i, j);
                model.add(lhs + uplifted <= inp.nC - 1);
            }
        }
    }
    cout << " ## ## ## ADDED " << nCovers << " covers. " << endl;
}

int is_minimal_cover(INSTANCE & inp, int * x, double excess, int group, int el,
        int constr)
{

    // cout << "IS IT A MINIMAL COVER ? " << endl;
    // cout << " .. excess is = " << excess << endl;
    int * xAux = new int[inp.nC];
    for (int i = 0; i < inp.nC; i++) 
        xAux[i] = x[i];
    xAux[group] = el;

    double maxInCover = 0.0;
    int maxClass   = -1;
    // find max value in cover for constraint "constr"
    for (int i = 0; i < inp.nC; i++) 
    {
        if (inp.w[i][xAux[i]][constr] > maxInCover)
        {
            maxInCover = inp.w[i][xAux[i]][constr]; 
            maxClass   = i;
        }
    }
    // cout << ".. max consumption is in class " << maxClass << " with value = "
    // << maxInCover << endl;

    // find max elements in maxClass
    double maxNotInCover = 0.0;
    for (int j = 0; j < inp.ri[maxClass]; j++) 
    {
        if (j == xAux[maxClass]) continue;
        if (inp.w[maxClass][j][constr] > maxNotInCover)
            maxNotInCover = inp.w[maxClass][j][constr];
    }
    // cout << ".. swap with max not in cover = " << maxNotInCover << endl;

    // if we swap, is it feasible?
    if (excess - maxInCover + maxNotInCover <= 0.0)
        return 1;
    else
        return -1;

}




void uplift_cover(INSTANCE inp, IloExpr & uplifted, TwoD & x_ilo, int isCover, 
        int * x, int group, int el)
{
    // cout << " CONSTRAINT " << isCover << endl;
    // sort coefficients from max to min
    int tempSol = x[group];
    x[group] = el; 
    /*
       cout << "The COVER is ";
       for (int i = 0; i < inp.nC; i++) 
       cout << " " << x[i];
       cout << endl;
       */ 
    double * w = new double[inp.nC];
    int * index = new int[inp.nC];
    for (int i = 0; i < inp.nC; i++) 
    {
        index[i] = i; 
        w[i] = inp.w[i][x[i]][isCover];
    }   
    quickSort(w, index, 0, inp.nC-1); // sort from min to max


    double * lower = new double[inp.nC-1];
    double * upper = new double[inp.nC-1];
    double progr = 0.0;
    for (int i = 0; i < inp.nC-1; i++) 
    {
        lower[i] = progr;
        upper[i] = progr + w[index[inp.nC-1-i]];
        progr += w[index[inp.nC-1-i]];
    }
#ifdef M_DEBUG
    for (int i = 0; i < inp.nC-1; i++) 
    {
        cout << "lower = " << lower[i] << " and upper " << upper[i] << endl;
    }
#endif

    std::vector < std::vector <int> > cols;
    std::vector < std::vector <int> > coeffs;

    for (int i = 0; i < inp.nC; i++) 
    {
        // cout << " CLASS " << i << " with x[" << i << "] = " << x[i] << endl;
        std::vector <int> auxi;
        std::vector <int> auxc;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            if (j == x[i]) continue;
            // cout << "... trying item " << j << endl;
            int l = 0;
            while (inp.w[i][j][isCover] >= upper[l])
                l++;
            if (l > 0)
            {
                auxi.push_back(j);
                auxc.push_back(l);
            }
        }

        // saving variables for current class
        cols.push_back(auxi);
        coeffs.push_back(auxc);
#ifdef M_DEBUG
        cout << "For class " << i << " we have ;; " << endl;
        int j = 0;
        for (std::vector<int>::iterator it = auxc.begin(); it != auxc.end(); ++it) 
        {
            cout << "item[" << i << "," << auxi[j] << "] = " << inp.w[i][auxi[j]][isCover] 
                << " in interval " << *it << endl;
            j++;
        }
#endif  

    }
    // cout << "SUMMARY for the VALID INEQUALITY " << endl;
    for (int i = 0; i < inp.nC; i++) 
        for (int j = 0; j < cols[i].size(); j++) 
        {
            // cout << " x[ " << i << "," << cols[i][j] << " ] with coeff " << coeffs[i][j] << endl;
            uplifted += x_ilo[i][cols[i][j]]*coeffs[i][j];

        }
    // is facet defining?
    double totS = 0.0;
    for (int i = 0; i < inp.nC; i++) 
        totS += w[i];

    int isFacet = 1;
    for (int h = 0; h < inp.nC; h++) 
    {

        totS -= w[index[inp.nC-1-h]];
        double rhs = inp.R[isCover] - w[index[inp.nC-1-h]];
        if (rhs < totS)
        {
            isFacet = 0;
            break;
        }
    }
    if (isFacet == 1)
    {
        cout << " FACET DEFINING " << endl;
        int aka;
        cin >> aka;
    }
    // restore solution
    x[group] = tempSol;
}

// Generate a cover inequality starting from the current feasible solution
void add_cover_inequality(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        int * x, double Omega, double * sigma2)
{
    IloEnv env = model.getEnv();
#ifdef M_DEBUG
    cout << "Feasible solution is :: ";
    for (int i = 0; i < inp.nC; i++) 
        cout << " " << x[i];
    cout << endl;
#endif
    double * W     = new double[inp.nR];
    double * Q     = new double[inp.nR];
    double * slack = new double[inp.nR];
    for (int k = 0; k < inp.nR; k++) 
    {
        W[k] = 0.0;
        Q[k] = 0.0;
        for (int i = 0; i < inp.nC; i++)
        {
            W[k] += inp.w[i][x[i]][k];
            Q[k] += sigma2[i]; 
        }
        // cout << "Checking constraint " << k << " LHS = " << W[k] + Omega*sqrt(Q[k])
        // << " vs RHS = " << inp.R[k] << endl;
        slack[k] = inp.R[k] - W[k] + Omega*sqrt(Q[k]);
        assert(slack[k] >= 0.0);
    }

    int nCovers = 0;
    int isCover = 0;
    for (int i = 0; i < inp.nC; i++) 
    {
        // cout << "COVERS for CLASS " << i << endl;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            isCover = 0;
            if (x[i] == j) continue;
            // cout << " ... trying with item " << j << endl;
            for (int k = 0; k < inp.nR; k++) 
            {
                // note: Q does not change because sigma_ij are all same
                double delta = -inp.w[i][x[i]][k] + inp.w[i][j][k] -sqrt(Q[k]) +
                    sqrt(Q[k] - sigma2[i] + sigma2[i]);
                if (delta > slack[k])
                {
                    // cout << "-----> constraint " << k << " is violated " << endl;
                    isCover = 1;
                    break;
                }

            }
            // if infeasible, it is a cover
            if (isCover)
            {
                // add cover
                nCovers++;
                IloExpr lhs(env);
                for (int l = 0; l < inp.nC; l++) 
                {
                    if (l == i) continue;
                    lhs += x_ilo[l][x[l]];
                    // cout << " " << x[l];
                }
                // cout << " + x[" << i << "," << j << " ] " << endl;
                lhs += x_ilo[i][j];
                model.add(lhs <= inp.nC - 1);
            }
        }
    }
    cout << " ## ## ## ADDED " << nCovers << " covers. " << endl;


}



/// Add cut to exclude current solution from feasible space
void cut_out_solution(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xIlo, IloRangeArray & cutSol)
{
    IloEnv env = model.getEnv();

    IloExpr lhs(env);
    //IloRangeArray cutSol1(env);
    for (int i = 0; i < inp.nC; i++)
        lhs += x_ilo[i][xIlo[i]];	 

    // for (int j = 0; j < inp.ri[i]; j++)
    //     if (cplex.getValue(x_ilo[i][j]) >= 1.0 - EPSI)
    //  	lhs += x_ilo[i][j];



    cutSol.add(lhs <= inp.nC-1);
    model.add(cutSol);

}

/// Add cut with obj function value of best found solution
void add_z_cut(IloModel & model, IloCplex & cplex, TwoD & x_ilo, double zBest)
{
    IloEnv env = model.getEnv();

    IloExpr lhs(env);
    //IloRangeArray cutSol1(env);
    for (int i = 0; i < inp.nC; i++)
        for (int j = 0; j < inp.ri[i]; j++)
            lhs += (double)inp.c[i][j]*x_ilo[i][j];

    model.add(lhs >= zBest);    
}
// sort vector x with respect to rc
void quickSort(double* rc, int* x, int low, int high)
{
    int  temp, left, right;
    long double median;

    if (high > low)
    {
        left   = low;
        right  = high;
        median = rc[x[low]];

        while (right > left)
        {
            while ( rc[x[left]] < median )
                left++;

            while (rc[x[right]] > median)
                right--;

            if (left > right) break;

            temp		= x[left];
            x[left]  	= x[right];
            x[right] 	= temp;

            left++;
            right--;
        }

        quickSort(rc, x, low, right);
        quickSort(rc, x, left,  high);
    }
}
// sort vector x with respect to rc
void quickSort(int* rc, int* x, int low, int high)
{
    int  temp, left, right, median;

    if (high > low)
    {
        left   = low;
        right  = high;
        median = rc[x[low]];

        while (right > left)
        {
            while ( rc[x[left]] < median )
                left++;

            while (rc[x[right]] > median)
                right--;

            if (left > right) break;

            temp		= x[left];
            x[left]  	= x[right];
            x[right] 	= temp;

            left++;
            right--;
        }

        quickSort(rc, x, low, right);
        quickSort(rc, x, left,  high);
    }
}

// Simulate the creation of nRuns random vectors w[i][j][k] and compute how 
// often the robust and the nominal solutions lead to infeasible situations.
/** It is worth noting that, in reality, we do not need to compute the number
 * of times the nominal solution leads to infeasibility, since the nominal 
 * solution never changes. However, since for every value of Omega we generate
 * different w vectors, we recompute the number of infeasibility for the 
 * nominal solution as well. This number should be pretty much constant.
 */
void robust_simulation(INSTANCE inp, int * xBest, int * xNominal, double zNominal,
        double zBest, int nRuns, double SSsigma)
{

    cout << "This is the current solution ;; " << endl;
    cout << "z^N vs z^R = " << zNominal << " " << zBest << endl;
    // cout << "CM solution :: ";
    // for (int i = 0; i < inp.nC; i++)
    // {
    // cout << " " << xBest[i];
    // }
    // cout << endl;
    // cout << "NOMINAL solution :: ";
    // for (int i = 0; i < inp.nC; i++)
    // {
    // cout << " " << xNominal[i];
    // }

    // now evaluate solution with random generated 
    int ** wN = new int*[inp.nC];
    for (int i = 0; i < inp.nC; i++)
        wN[i] = new int[inp.nR];

    int ** wR = new int*[inp.nC];
    for (int i = 0; i < inp.nC; i++)
        wR[i] = new int[inp.nR];

    // =============== EVAL ======================
    int totInfeasible        = 0;
    int totInfeasibleNominal = 0;
    for (int cc = 0; cc < nRuns; cc++)
    {
        // cout << "Iter " << cc << endl;

        for (int i = 0; i < inp.nC; i++)
        {
            for (int k = 0; k < inp.nR; k++)
            {
                double rr = (double)(rand()+1)/((double)(RAND_MAX)+1.0);
                if (rr <= 0.5)
                {
                    wR[i][k] = inp.w[i][xBest[i]][k] - sigma2[i];
                    wN[i][k] = inp.w[i][xNominal[i]][k] - sigma2[i];
                }
                else
                {
                    wN[i][k] = inp.w[i][xNominal[i]][k] + sigma2[i];
                    wR[i][k] = inp.w[i][xBest[i]][k] + sigma2[i];
                }
            }
        }

        // this is what we have done
        /* cout << "Comparison of nominal vs real :: " << endl;
         * for (int i = 0; i < inp.nC; i++)
         * {
         *     cout << "Item " << xBest[i] << " ... " << endl;
         *     for (int k = 0; k < inp.nR; k++)
         *     {
         *         cout << inp.w[i][xBest[i]][k] << " vs " << wR[i][k] << endl;
         *     }
         * } */

        // is robust solution feasible?
        for (int k = 0; k < inp.nR; k++)
        {
            double tot = 0.0;
            for (int i = 0; i < inp.nC; i++) 
                tot += wR[i][k];
            // cout << "lhs = " << tot << " vs rhs = " << inp.R[k] << endl;
            if (tot > inp.R[k])
            {
                // cout << "Infeasible in constraint " << k << endl;
                totInfeasible++;

                break;
            }
        }

        // is nominal solution feasible?
        for (int k = 0; k < inp.nR; k++)
        {
            double tot = 0.0;
            for (int i = 0; i < inp.nC; i++) 
                tot += wN[i][k];
            // cout << "lhs = " << tot << " vs rhs = " << inp.R[k] << endl;
            if (tot > inp.R[k])
            {
                // cout << "Infeasible in constraint " << k << endl;
                totInfeasibleNominal++;
                break;
            }
        }
    }

    cout << " [" << totInfeasible << "/" << nRuns << "]" 
        << endl;

    cout << " [" << totInfeasibleNominal << "/" << nRuns << "]" 
        << endl;

    string sRatio;
    if (zBest != -1)
    {
        double ratio = (zNominal - zBest)/zNominal;
        stringstream ss;
        ss << ratio;
        sRatio = ss.str();
    }
    else
        sRatio = "-";



    // write info to diskfile
    ofstream fRandom("robust.txt", ios::out);
    fRandom << zNominal << "\t" << totInfeasibleNominal << "\t" 
        << zBest << "\t" << sRatio << "\t" << Omega << "\t" << 
        SSsigma/(double)inp.nC << "\t" 
        << totInfeasible << "\t" << nRuns << endl; 
    fRandom.close();
}

// Find second best lagrangean cost for each class
void compute_stability(INSTANCE & inp, double ** rc, double * delta, int * xL)
{

    for (int i = 0; i < inp.nC; i++) 
    {
        double second = -INFTY;
        for (int j = 0; j < inp.ri[i]; j++) 
        {
            if (xL[i] == j) continue;
            if (rc[i][j] > second)
                second = rc[i][j];
        }
        delta[i] = rc[i][xL[i]] - second;
    }
    // cout << "Stability : " << endl;
    // for (int i = 0; i < inp.nC; i++)
    // {
    // cout << "s[" << i << " ] = " << delta[i] << endl;
    // }
}
