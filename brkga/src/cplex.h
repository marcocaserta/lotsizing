#ifndef CPLEX_H
#define CPLEX_H

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

typedef IloArray < IloNumVarArray > TwoD;

double defineModel(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj);
double defineRobustModel(IloModel & model, IloCplex & cplex, TwoD & x_ilo, 
        IloObjective & obj, IloNumVar & Q_ilo, double Omega, double * sigma2);
double defineRobustDet(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, double Omega, double sigma, double u);
double defineRobustBertsimas(IloModel & model, IloCplex & cplex, TwoD & x_ilo, IloObjective & obj, double Omega, double sigma, double u);
double solve_KNAP(IloModel model, IloCplex cplex, int solLimit, int displayLimit, int timeLim);
double get_cplex_sol(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xIlo);
void add_cut_corridor(IloModel & model, IloCplex & cplex, TwoD & x_ilo, int * xCorridor, int rhs, IloRangeArray & corridor, IloExpr lhs);

#endif
