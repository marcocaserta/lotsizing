#!/usr/bin/python
"""
===========================================================================

:Filename: clspBenders.py
:Author: marco caserta
:Date: 09.03.2017
:Last Update: |date|

.. |date| date:: %d.%m.%y
.. |time| date:: %H:%M

Copyright (C) 2017 by Marco Caserta  (marco dot caserta at ie dot edu)

(This document was generated on |date| at |time|.)

.. this is just a comment

===========================================================================

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, write to the Free Software Foundation, Inc., 59 Temple
   Place - Suite 330, Boston, MA  02111-1307, USA.

===========================================================================

Introduction
------------

**NOTE**: A completely modified version of this code is presented below.

This code implements a simple benders decomposition scheme for the
multi-item multi-period capacitated lot sizing problem. Benders cuts are added
via LazyConstraintCallback(). The case of an infeasibility cut is not
considered, since the subproblem is always feasible (due to the possibility
to choose an arbitrarily large value for the initial inventory level.) This
should be improved adding a cut based on extreme rays.

The following is a formulation for the MIMPLS:

.. math ::
    :nowrap:

    \\begin{eqnarray}
       \max     & & \displaystyle \sum_{j=1}^n \sum_{t=1}^T (f_{jt}y_{jt} +
        c_{jt}x_{jt} + h_{jt}s_{jt}) \\\\
                & s.t. & \\nonumber \\\\
                && \displaystyle \sum_{j=1}^n a_{jt}x_{jt}+m_{jt}y_{jt} \leq
                b_t, \quad \\forall t \\\\
                && x_{jt} + s_{jt-1} - s_{jt} = d_{jt}, \quad \\forall j,t \\\\
                && x_{jt} \leq M y_{jt}, \quad \\forall j,t \\\\
                && y_{jt} \in \left\{0,1\\right\}, x_{jt}, s_{jt} \geq 0
    \end{eqnarray}

To address the problem using Benders decomposition, we define a master, which
includes the *difficult* binary variables :math:`y_{jt}`, and a subproblem,
which deals with the continuous variables :math:`x_{jt}, s_{jt}`.

Master Problem
~~~~~~~~~~~~~~

The master is:

.. math ::
    :nowrap:

    \\begin{eqnarray}
       \max     && \displaystyle \sum_{j=1}^n \sum_{t=1}^T f_{jt}y_{jt} + z\\\\
                & s.t. & \\nonumber \\\\
                && z \geq \phi(\mathbf{y}) \\\\
                && y_{jt} \in \left\{0,1\\right\}
    \end{eqnarray}

where :math:`\phi(\mathbf{y})` is the best possible cost obtained for a given
:math:`\mathbf{y}`.

Subproblem
~~~~~~~~~~

Given a solution to the master problem :math:`\mathbf{y}^*`, we define a primal
subproblem as:

.. math ::
    :nowrap:

    \\begin{eqnarray}
       \max     & & \displaystyle \sum_{j=1}^n \sum_{t=1}^T c_{jt}x_{jt} + h_{jt}s_{jt} \\\\
                & s.t. & \\nonumber \\\\
                && \displaystyle \sum_{j=1}^n a_{jt}x_{jt}\leq
                b_t - m_{jt}y_{jt}^* , \quad \\forall t \\\\
                && x_{jt} + s_{jt-1} - s_{jt} = d_{jt}, \quad \\forall j,t \\\\
                && x_{jt} \leq M y_{jt}^*, \quad \\forall j,t \\\\
                && x_{jt}, s_{jt} \geq 0
    \end{eqnarray}

Note that the subproblem is an LP and, therefore, can easily be solved using
cplex. Once the subproblem is solved, we obtain the optimal dual values.

Assume we have:

* :math:`\lambda_t`: The dual values of the capacity constraints (2)
* :math:`\omega_{jt}`: The dual values of the demand constraints (3)
* :math:`\\nu{jt}`: The dual values of the logical constraints (4)

Then, we write the benders cut as follows:

.. math ::

    \displaystyle z \geq \sum_{t=1}^T \left( b_t - \sum_{j=1}^n m_{jt}y_{jt} \\right)
    + \sum_{j=1}^n \sum_{t=1}^T d_{jt}\omega_{jt} + \sum_{j=1}^n \sum_{t=1}^T
    M\\nu_{jt}y_{jt}

    \sum_{j=1}^n \sum_{t=1}^T y_{jt} \left(M\\nu_{jt} - \lambda_t m_{jt} \\right) -z \leq
    - \sum_{t=1}^t \lambda_t b_t - \sum_{j=1}^n \sum_{t=1}^T \omega_{jt} d_{jt}


Alternative Formulation Based on the Support Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(Based on Floudas, Non-linear and Mixed-integer Programming, Ch.6)

The general master program is:

.. math ::
    :nowrap:

    \\begin{eqnarray}
       \max     && \displaystyle \sum_{j=1}^n \sum_{t=1}^T f_{jt}y_{jt} + z\\\\
                & s.t. & \\nonumber \\\\
                && z \geq L(x,y,s,\lambda, \\nu) \quad \\forall \lambda, \\nu \\\\
                && 0 \geq L(x,y,s,\\bar{\lambda}, \\bar{\\nu}) \quad \\forall
                \\bar{\lambda}, \\bar{\\nu} \\\\
                && y_{jt} \in \left\{0,1\\right\}
    \end{eqnarray}

where constraints (2) are *optimality cuts*, while constraints (3) are
*feasibility cuts*. In this work, we ignore the feasility cuts, due to the use
of initial inventory variables ``sI``. Now, since the above formulation
requires the definition of all the possible value of :math:`\lambda` and
:math:`\mu`, we solve a **relaxation** of the master problem. Therefore, at
every iteration, the master provides a **lower bound** of the optimal value.

The master contains, as constraints, two inner optimization problems. We
express the inner minimization problem in terms of **support function**, i.e


.. math ::

    \\xi (y,\lambda, \\nu) = \min_{x \in X} L(x,y,\lambda, \\nu)

Using ideas from Geoffrion (1972), we consider the case of linearly separable
functions in :math:`x` and :math:`y`, i.e., the objective function
:math:`f(x,y)` and the constraints :math:`g_k(x,y)` can be separated. In this
case, the optimalit cut is:

.. math ::

    z \geq \displaystyle \sum_{j=1}^n \sum_{t=1}^T \left(-m_{jt}\lambda_t +
    M\\nu_{jt} \\right)\\times \left( y_{jt} - y_{jt}^*\\right) + z_{s}

where :math:`z_{s}` is the optimal objective function value of the last
subproblem, and :math:`y^*` indicates the current optimal solution of the
master problem.

.. note::

    The dual values have positive sign here, since they are treated as
    Lagrangean multipliers (this is the reason why the sign is reversed.)

How to Run This Code
--------------------


See :func:`parseCommandLine()`.


History
-------

15.03.17:

    The problem seems to be that the lower bound provided by the master is not
    tight. I attempted to tighten the bound, by adding constraints linking
    (somehow) the variables :math:`y_{jt}` with :math:`x_{jt}`, but no
    improvement is observed. I tried a constraint that ensures that the number
    of :math:`y` variables set to 1 is enough to covere the total demand. The
    constraint seems to be correct, but no improvement is observed with respect
    to the lower bound of the master.

23.03.17:

    Added user cuts, to get cuts for fractional values of y* as well.

Redefinition based on SPL
-------------------------

A modified version of this code is presented here. Basically, two important
changes have been introduced:

1. We use the Simple Plant Location (SPL) reformulation for the Lot Sizing
problem 
2. Rather than relying on LazyConstraintCallback(), we now implement a cycle
to define Benders' scheme, with iterative calls to master and subproblems.

The following is the SPL reformulation:

.. math ::
    :nowrap:

    \\begin{eqnarray}
      & \min z = &   \sum_{j=1}^n \sum_{t=1}^T f_{jt}y_{jt} +
      \sum_{j=1}^n \sum_{r=1}^T \sum_{t=1}^{r-1} h_{jtr}z_{jtr}
      \label{eq:SPL-obj}\\\\
      &\mbox{s.t} & \sum_{j=1}^n \left(a_{jt}\sum_{r=t}^T z_{jtr} +
      m_{jt}y_{jt}\\right) \leq b_t,
      \quad t = 1, \ldots, T \label{eq:SPL-capacity-constr} \\\\
      &&  \sum_{t=1}^r z_{jtr} = d_{jr},
      \quad j = 1, \ldots, n , \quad r = 1, \dots, T \label{eq:SPL-demand-constr}
      \\\\
      && z_{jtr} \leq d_{jr}y_{jt},
      \quad j = 1, \ldots, n , \quad t = 1, \dots, T, \quad r=t,\dots,T
      \label{eq:SPL-logic-constr} \\\\
      && \sum_{r=t}^T z_{jtr} \leq M y_{jt},
      \quad j = 1, \ldots, n , \quad t = 1, \dots, T \\\\
      && y_{jt} \in \left\{0,1\\right\}, \quad j = 1, \ldots, n, \quad
      t=1, \ldots, T \label{eq:SPL-y-binary}\\\\
      && z_{jtr} \geq 0, \quad j = 1, \ldots, n, \quad
      t=1, \ldots, T, \quad r = t,\dots,T \label{eq:SPL-z-cont}
    \end{eqnarray}

where we define :math:`h_{jtr} = \sum_{t'=t}^{r-1} h_{jt'}` as the
cumulative cost of keeping in inventory a unit of item *j* from period *t* to
period *r-1*. 

The separation scheme and, consequently, the construction of master and
subproblems is as before. However, the main change here concerns the use of a
cycle to iteratively solve the master and the subproblems. 

We call :func:`benderAlgorithm()`, which is a function that implements the cycle.
Before doing that, we might apply some fixing schemes (both to zero and to one)
based on the LP relaxation. See functions :meth:`MIP.solveLPOne()` and
:meth:`MIP.solveLPZero()`.

"""

from __future__ import print_function

import sys, getopt
import math
import time

import cplex
from cplex.callbacks import UserCutCallback, LazyConstraintCallback
from cplex.callbacks import SolveCallback, SimplexCallback

from cplex.exceptions import CplexError

from lagrange import *
from dw2 import *
#  from dw import *

lbSummary = "lowerBounds.txt"

#  lbFile = open(lbSummary, "a")

_INFTY    = sys.float_info.max
_EPSI     = sys.float_info.epsilon
inputfile = ""
userCuts  = "0"
cPercent  = 0.0
cZero     = 1.0 #  soft-fixing to zero is inactive
cOne      = 0.0 #  soft-fixing to one is inactive
algo      = -1
fixToZero = []
fixToOne  = []

#  we set the master variables as global, while the subproblem vars are local
z_ilo     = -1
y_ilo     = []
lCapacity = []
lLogic    = []
lDemand   = []

inout     = []
yRef      = []
yPool     = []
nPool     = 0
ubBest    = 0.0
startTime = time()

class SolveNodeCallback(SimplexCallback):
    def __call__(self):

        print("I AM here ")
        #  self.solve()

        print("OBJ ", self.get_objective_value())
        print("RC", self.get_reduced_costs())
        input(" ... ")

class BendersLazyConsCallback(LazyConstraintCallback):
    """
    This is the LazyConstraintCallback of cplex. We implement Benders algorithm
    via callback. That it, the master is solved within a branch and bound
    framework and, every time a new master solution is obtained, the callback
    is used to:

    * get the master solution :math:`y`
    * pass it to the subproblem to obtain :math:`\phi(y)`
    * get the dual values and define Benders cut
    * add the cut to the master
    * give control back to cplex

    """

    def __call__(self):
        """
        Define the actions to be carried out at every callback.

        Note that the ``separate`` function of the subproblem is called here.

        """

        # get data structure (self is the master)
        cpx    = self.cpx
        worker = self.worker
        y_ilo  = self.y_ilo
        z_ilo  = self.z_ilo
        inp    = self.inp
        yFixed = self.yFixed

        #  get current master solution
        zHat = self.get_values(z_ilo)
        ySol = []
        for j in range(inp.nI):
            ySol.append([])
            ySol[j] = self.get_values(y_ilo[j])

        #  flatten =  [item for sublist in ySol for item in sublist]
        #  benders cut separation
        cutType = worker.separate(inp, ySol, zHat, y_ilo, z_ilo)
        if cutType > 0:
            #  a = [float(worker.cutLhs.val[i]) for i in range(inp.nI*inp.nP)]
            #  lhsSum = sum([a[i]*flatten[i] for i in range(inp.nI*inp.nP)])
            #  print("LhsSum = ", lhsSum , " vs ", worker.cutRhs)
            #  print(lhsSum <= worker.cutRhs)
            #  violated = (lhsSum - worker.cutRhs) > 0.1
            #  print(" violated ? ", violated)
            violated = 1
            if violated:
            # add Benders cut to the master
                self.add(constraint = worker.cutLhs,
                         sense     = "L",
                         rhs        = worker.cutRhs,
                         use        = 0)


        zLP = self.get_best_objective_value()
        #  print("here ", zLP)
        #  input("...")
        if self.solved == 0:
            cpxCloneLP = cplex.Cplex(cpx)
            cpxCloneLP.set_problem_type(cpxCloneLP.problem_type.LP)
            cpxCloneLP.solve()
            self.solved = 1

            for j in range(inp.nI):
                self.rc.append(cpxCloneLP.solution.get_reduced_costs(y_ilo[j]))
            #  add cut here ??

            #  print("Before adding cut to master : ", cpx.linear_constraints.get_num())
            #  cutType = worker.separate(inp, yRef, 0.0, y_ilo, z_ilo)
            #  print(worker.cutLhs, " <= ", worker.cutRhs)
            #  if cutType > 0:
            #
            #      self.add(constraint = worker.cutLhs,
            #               sense     = "L",
            #               rhs        = worker.cutRhs,
            #               use        = 0)
            #
            #  print("Cut added to master : ", cpx.linear_constraints.get_num())
            #  input("....")


        #  nRowsMaster = cpx.linear_constraints.get_num()
        #  nRows = cpxClone.linear_constraints.get_num()
        #  print(" entering with ", nRowsMaster, " rows in master and ", nRows, "\
        #  rows in clone ... ")
        #  if nRowsMaster <= nRows:
        #      cpxClone.linear_constraints.add(lin_expr=[worker.cutLhs],
        #                                              senses  =["L"],
        #                                              rhs     =[worker.cutRhs])
        #      return
        #
        #  index = [i for i in range(nRows, nRowsMaster)]
        #
        #  #  print("ROWS ARE = ", cpx.linear_constraints.get_rows())
        #  #  print("rhs are  = ", cpx.linear_constraints.get_rhs())
        #  allConstr = cpx.linear_constraints.get_rows(index)
        #  allRhs    = cpx.linear_constraints.get_rhs(index)
        #  for i,j in enumerate(allRhs):
        #      #  print(i,j, allConstr[i])
        #
        #      cpxClone.linear_constraints.add(lin_expr = [allConstr[i]],
        #                                      senses   = ["L"],
        #                                      rhs      = [j])
        #
        #  #  cpx.set_problem_type(cpx.problem_type.LP)
        #  #  cpx.solve()
        #  #  for j in range(inp.nI):
        #  #      rc = cpx.solution.get_reduced_costs(y_ilo[j])
        #  #      print("REAL RC = ", rc)
        #
        #  #  solve Master LP
        #  cpxClone.solve()
        #  #  print("LP sol Master is ", cpx.solution.get_objective_value())
        #  zClone = cpxClone.solution.get_objective_value()
        #  slack = cpxClone.solution.get_linear_slacks()
        #  remove = [i for i in range(nRows) if slack[i] > _EPSI]
        #  print(" ... due to SLACK, removing ", len(remove), " constraints.")
        ub = self.get_objective_value()
        #  zClone = cpxCloneLP.solution.get_objective_value()
        #  print("CLONE z = ", zClone, " vs UB = ", ub)
        #  print("ubBes is ", ubBest, " vs ub = ", ub)
        #  from here
        fixInClone = []
        for j in range(inp.nI):
            #  rc = cpxCloneLP.solution.get_reduced_costs(y_ilo[j])
            for t in range(inp.nP):
                #  if yFixed[j][t] == 0 and (zLP + rc[t]) > ub:
                if yFixed[j][t] == 0 and (zLP + self.rc[j][t]) > ubBest:
                    yFixed[j][t] = 1
                    print(" [", self.nIter,"] ** ** ** ** fixing to zero ", y_ilo[j][t])
                    fixInClone.append(y_ilo[j][t])
                    self.add(constraint=cplex.SparsePair(ind=[y_ilo[j][t]],val=[1.0]),
                             sense = "E",
                             rhs   = 0.0)
                    #  cpxClone.variables.set_upper_bounds(y_ilo[j][t], 0.0)

        self.nIter += 1
        self.yFixed = yFixed


class BendersUserCutCallback(UserCutCallback):

    def __call__(self):
        """
        Define the actions to be carried out at every callback.

        Note that the ``separate`` function of the subproblem is called here.

        """
        # Skip the separation if not at the end of the cut loop
        print("SOL USER CUT ", self.get_best_objective_value())
        cpxClone = cplex.Cplex(self)
        cpxClone.set_problem_type(cpxClone.problem_type.LP)
        cpxClone.solve()
        print("SOL ", cpxClone.solution.get_objective_value())

        input("....")
        return

        if not self.is_after_cut_loop():
            return

        # get data structure (self is the master)
        worker = self.worker
        y_ilo  = self.y_ilo
        z_ilo  = self.z_ilo
        inp    = self.inp

        #  get current master solution
        zHat = self.get_values(z_ilo)
        ySol = []
        for j in range(inp.nI):
            ySol.append([])
            ySol[j] = self.get_values(y_ilo[j])

        flatten =  [item for sublist in ySol for item in sublist]
        #  benders cut separation
        cutType = worker.separate(inp, ySol, zHat, y_ilo, z_ilo)
        if cutType > 0:
            a = [float(worker.cutLhs.val[i]) for i in range(inp.nI*inp.nP)]
            lhsSum = sum([a[i]*flatten[i] for i in range(inp.nI*inp.nP)])
            print("LhsSum = ", lhsSum , " vs ", worker.cutRhs)
            violated = (lhsSum - worker.cutRhs) > 0.1
            print("VIOLATED = ", violated)
            #  if (lhsSum > worker.cutRhs):
            if violated:
                # add Benders cut to the master
                self.add(cut       = worker.cutLhs,
                         sense     = "L",
                         rhs       = worker.cutRhs,
                         use       = 0)




def parseCommandLine(argv):
    """
    .. func:parseCommandLine()

    Parse command line. Options are:

    -h help         usage help

    -i inputfile    instance file name

    -u userCuts    activate user cuts (fractional values)

    -c cpercent     corridor width

    -z zeros        soft fixing to zero
    
    -o ones        soft fixing to one

    -a algorithm   type of algorithm used:

    With respect to the type of algorithms that can be used, we have:

        1.  Benders Decomposition
        2.  Lagrangean Relaxation
        3.  Dantzig-Wolfe
        4.  Cplex MIP solver

    """
    global inputfile
    global userCuts
    global cPercent
    global cZero
    global cOne
    global algo

    try:
        opts, args = getopt.getopt(argv, "hi:u:c:z:o:a:",
        ["help","ifile=","ucuts","cpercent","zeros","ones","algorithm"])
    except getopt.GetoptError:
        print("Command Line Error. Usage : python cflp.py -i <inputfile> -u\
        <usercuts> -c <corridor width> -z <fix to zero> -o <fix to one> \
        -a <algorithm BD, LR, DW, Cplex>")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("Usage : python cflp.py -i <inputfile> -u <usercuts> -c \
            <corridor width> -z <fix to zero> -o <fix to one> \
            -a <algorithm - BD, LR, DW, Cplex>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-u", "--ucuts"):
            userCuts = arg
        elif opt in ("-c", "--cpercent"):
            cPercent = float(arg)
        elif opt in ("-z", "--zeros"):
            cZero = float(arg)
        elif opt in ("-o", "--ones"):
            cOne  = float(arg)
        elif opt in ("-a", "--algorithm"):
            algo  = float(arg)

class Instance:
    """
    Class used to read the instance from a disk file.
    Instances are obtained from Trigeiro. We also compute a tight value for the
    big M constant, as well as a cumulative demand value for t to T.

    """
    def __init__(self, inputfile):
       self.d        = []
       self.c        = []
       self.f        = []
       self.h        = []
       self.a        = []
       self.m        = []
       self.max_prod = []
       self.cap      = []
       self.dcum     = []
       with open(inputfile) as ff:
            data = ff.readline()
            self.nI, self.nP = [int(v) for v in data.split()]

            print("Nr. items and nr. periods = ", self.nI, " ", self.nP)

            data = ff.readline()
            self.cap = [float(v) for v in data.split()]*self.nP

            for i in range(self.nI):
                data = ff.readline()
                temp = [float(v) for v in data.split()]
                #  a = temp[0]*self.nP;
                a = [1.0]*self.nP
                h = [temp[1]]*self.nP
                m = [temp[2]]*self.nP
                f = [temp[3]]*self.nP
                c = [0.0]*self.nP

                self.a.append(a)
                self.h.append(h)
                self.m.append(m)
                self.f.append(f)
                self.c.append(c)
            self.d = [ [] for i in range(self.nI)]
            for t in range(self.nP):
                data = ff.readline()
                temp = [float(v) for v in data.split()]
                [self.d[i].append(temp[i]) for i in range(self.nI)]


       #  compute cumulative demand
       for j in range(self.nI):
           aux = []
           aux.append(self.d[j][self.nP-1])
           progr = 0
           for t in range(self.nP-2,-1,-1):
               aux.append(aux[progr]+self.d[j][t])
               progr += 1
           self.dcum.append(list(reversed(aux)))


       # max production of item j in period t is the minimum between
       # the limit set by capacity and the cumulative demand
       for j in range(self.nI):
           aa = [( (self.cap[t] - self.m[j][t])/self.a[j][t]) for t in range(self.nP)]
           self.max_prod.append([min(aa[t],self.dcum[j][t]) for t in \
           range(self.nP)])


class MIP:
    """
    .. class:MIP()

    Define the full model and solve it using cplex. We use this class to
    compute optimal values of the full MIP models and compare them, along with
    the achieve performance, with the Benders approach.

    This is the Standard Lot Sizing implementation, using variables
    :math:`y_{jt}, x_{jt}, s_{jt}`. This formulation is no longer used in the
    current version of the code, since the SPL reformulation provides tighter
    relaxations.

    """
    def __init__(self, inp):
        y_ilo = []
        x_ilo = []
        s_ilo = []
        sI    = []

        cpx = cplex.Cplex()
        cpx.objective.set_sense(cpx.objective.sense.minimize)

        #  create variables y_jt
        for j in range(inp.nI):
            y_ilo.append([])
            for t in range(inp.nP):
                varName = "x." + str(j) + "." + str(t)
                y_ilo[j].append(cpx.variables.get_num())
                cpx.variables.add(obj   = [inp.f[j][t]],
                                  lb    = [0],
                                  ub    = [1],
                                  types = ["B"],
                                  names = [varName])

        #  create variables x_jt
        for j in range(inp.nI):
            x_ilo.append([])
            for t in range(inp.nP):
                varName = "x." + str(j) + "." + str(t)
                x_ilo[j].append(cpx.variables.get_num())
                cpx.variables.add(obj   = [inp.c[j][t]],
                                  lb    = [0.0],
                                  ub    = [cplex.infinity],
                                  types = ["C"],
                                  names = [varName])
        #  create variables s_jt
        for j in range(inp.nI):
            s_ilo.append([])
            for t in range(inp.nP):
                varName = "s." + str(j) + "." + str(t)
                s_ilo[j].append(cpx.variables.get_num())
                cpx.variables.add(obj   = [inp.h[j][t]],
                                  lb    = [0.0],
                                  ub    = [cplex.infinity],
                                  types = ["C"],
                                  names = [varName])

        #  initial inventory level (avoid infeasibility)
        for j in range(inp.nI):
            varName = "sI." + str(j)
            sI.append(cpx.variables.get_num())
            cpx.variables.add(obj   = [100000],
                              lb    = [0.0],
                              ub    = [cplex.infinity],
                              types = ["C"],
                              names = [varName])

        #  demand constraints
        for j in range(inp.nI):
            #  first period
            index = [x_ilo[j][0], sI[j], s_ilo[j][0]]
            value = [1.0, 1.0, -1.0]
            demand_constraint = cplex.SparsePair(ind=index, val=value)
            cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                       senses   = ["E"],
                                       rhs      = [inp.d[j][0]])
            #  periods 2 to T-1
            for t in range(1,inp.nP-1):
                index = [x_ilo[j][t], s_ilo[j][t-1], s_ilo[j][t]]
                value = [1.0, 1.0, -1.0]
                demand_constraint = cplex.SparsePair(ind=index, val=value)
                cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                           senses   = ["E"],
                                           rhs      = [inp.d[j][t]])

            #  last period
            index = [x_ilo[j][inp.nP-1], s_ilo[j][inp.nP-2]]
            value = [1.0, 1.0]
            demand_constraint = cplex.SparsePair(ind=index, val=value)
            cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                       senses   = ["E"],
                                       rhs      = [inp.d[j][inp.nP-1]])

        #  capacity constraints
        for t in range(inp.nP):
            index = [x_ilo[j][t] for j in range(inp.nI)]
            value = [inp.a[j][t] for j in range(inp.nI)]
            index = index + [y_ilo[j][t] for j in range(inp.nI)]
            value = value + [inp.m[j][t] for j in range(inp.nI)]
            capacity_constraint = cplex.SparsePair(ind=index, val=value)
            cpx.linear_constraints.add(lin_expr = [capacity_constraint],
                                       senses   = ["L"],
                                       rhs      = [inp.cap[t]])

        #  logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                index = [x_ilo[j][t], y_ilo[j][t]]
                value = [1.0, -inp.max_prod[j][t]]
                logic_constraint = cplex.SparsePair(ind =index, val=value)
                cpx.linear_constraints.add(lin_expr = [logic_constraint],
                                           senses   = ["L"],
                                           rhs      = [0.0])

        self.cpx   = cpx
        self.y_ilo = y_ilo
        self.x_ilo = x_ilo
        self.s_ilo = s_ilo
        self.sI    = sI

    def solveLPZero(self, inp):
        """
        .. method:solveLPZero()

        Solve LP relaxation of original MIP twice:
        
        - the first time, we solve the LP relaxation of the whole problem, and
          we store the variables whose value is zero in the LP solution
          (indexLP1)
        - the second time, we add a "corridor" type of constraint to the LP,
          enforcing that at least a given number of variables in indexLP1 will
          change value, i.e., will take a value above 0. We store in indexLP2
          the variables that take value zero in the LP-constrained model.
        - the intersection between indexLP1 and indexLP2 gives the set of
          variables we want to keep fixed to zero.
        """
        cpx = self.cpx
        y_ilo = self.y_ilo

        #  transform MIP into LP and solve it
        cpx.set_problem_type(cpx.problem_type.LP)
        self.solve(inp)
        yLP = []
        for j in range(inp.nI):
            yLP.append(cpx.solution.get_values(y_ilo[j]))

        #  add "corridor" constraint (and solve LP again)
        indexLP1 = [y_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP) if
        yLP[j][t] <= _EPSI]
        value = [1.0]*len(indexLP1)
        rhsVal = 0.25*len(indexLP1) #  change at least rhsVal variables
        zero_constraint = cplex.SparsePair(ind=indexLP1,val=value)
        cpx.linear_constraints.add(lin_expr  = [zero_constraint],
                                   senses    = ["G"],
                                   rhs       = [rhsVal],
                                   names     = ["fixLP"])

        self.solve(inp) #  solve LP-constrained version
        yLP = []
        for j in range(inp.nI):
            yLP.append(cpx.solution.get_values(y_ilo[j]))
        indexLP2 = [y_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP) if
        yLP[j][t] <= _EPSI]
        #  fix to zero vars that are at zero in both LPs
        fixToZero = list(set(indexLP1).intersection(indexLP2))
        #  print("INTERSECTION = ", fixToZero)

        #  restore MIP
        cpx.set_problem_type(cpx.problem_type.MILP)
        for j in range(inp.nI):
            for t in range(inp.nP):
                cpx.variables.set_types(y_ilo[j][t], cpx.variables.type.binary)

        #  remove the corridor constraint from the main model
        cpx.linear_constraints.delete("fixLP")

        return fixToZero


    def solveLPOne(self, inp):
        """
        .. method:solveLPOne()

        Same idea presented in :func:`solveLPZero()`. See comment above.
        """
        cpx = self.cpx
        y_ilo = self.y_ilo

        #  transform into LP and solve it
        cpx.set_problem_type(cpx.problem_type.LP)
        self.solve(inp)
        yLP = []
        for j in range(inp.nI):
            yLP.append(cpx.solution.get_values(y_ilo[j]))

        #  add "corridor" contraint
        indexLP1 = [y_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP) if
        yLP[j][t] >= (1.0 -_EPSI)]
        value = [1.0]*len(indexLP1)
        rhsVal = (1.0-0.25)*len(indexLP1)
        zero_constraint = cplex.SparsePair(ind=indexLP1,val=value)
        cpx.linear_constraints.add(lin_expr  = [zero_constraint],
                                   senses    = ["L"],
                                   rhs       = [rhsVal],
                                   names     = ["fixLP"])

        #  solve constrained version of LP
        self.solve(inp)
        yLP = []
        for j in range(inp.nI):
            yLP.append(cpx.solution.get_values(y_ilo[j]))

        indexLP2 = [y_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP) if
        yLP[j][t] >= (1.0-_EPSI)]

        #  get interception between indexLP1 and indexLP2
        fixToOne = list(set(indexLP1).intersection(indexLP2))
        #  print("INTERSECTION = ", fixToOne)

        #  restore original MIP
        cpx.set_problem_type(cpx.problem_type.MILP)
        for j in range(inp.nI):
            for t in range(inp.nP):
                cpx.variables.set_types(y_ilo[j][t], cpx.variables.type.binary)

        #  eliminate corridor constraint
        cpx.linear_constraints.delete("fixLP")

        return fixToOne

    def solve(self, inp, nSol=99999, withPool=0, withPrinting=0, display=0,
            timeLimit = 10000):
        """
        .. method:solve()

        Solve the original MIMPLS using cplex branch and bound. A number of
        flags can be activate, to control the behavior of the solver:

        * nSol : maximum number of solutions to be visited
        * withPool : collect a pool of solutions during the optimization phase
        * withPrinting: control the output
        * display : control cplex output
        * timeLimit : set a maximum time limit
        """
        global yRef
        global ubBest
        global yPool
        global nPool
        global startTime

        y_ilo = self.y_ilo
        #  z_ilo = self.z_ilo
        cpx   = self.cpx

        #  cpx.set_results_stream(None)
        #  cpx.set_log_stream(None)
        cpx.parameters.timelimit.set(timeLimit)
        cpx.parameters.mip.limits.solutions.set(nSol)
        cpx.parameters.mip.display.set(display)
        #  cpx.parameters.mip.interval.set(500) # how often to print info
        cpx.parameters.mip.tolerances.mipgap.set(0.000000001)
        cpx.parameters.mip.tolerances.absmipgap.set(0.000000001)


        cpx.solve()

        if withPrinting == 1:
            print("STATUS = ", cpx.solution.status[cpx.solution.get_status()])
            print("OPT SOL found = ", cpx.solution.get_objective_value())
            print("Time          = ", time() - startTime)
        #  if cpx.solution.get_status() == cpx.solution.status.optimal_tolerance\
            #  or cpx.solution.get_status() == cpx.solution.status.optimal:
        ubBest = cpx.solution.get_objective_value()

        if withPrinting == 2:
            for j in range(inp.nI):
                yRef.append(cpx.solution.get_values(y_ilo[j]))
                print("y(",j,") = ", yRef[j])
                for t in range(inp.nP):
                    print("   z(",t,") = ",
                    [cpx.solution.get_values(z_ilo[j][t][r]) for r in
                    range(inp.nP)])

        if withPool==1:
            names = cpx.solution.pool.get_names()
            nPool = len(names)
            for n in names:
                print("z(",n,") = ", cpx.solution.pool.get_objective_value(n))

                yAux = []
                for j in range(inp.nI):
                    yAux.append(cpx.solution.pool.get_values(n, y_ilo[j]))
                yPool.append(yAux)

        return ubBest



class MIPReformulation(MIP):
    """
    .. class:MIPReformulation()

    This class implements the SPL reformulation, which is the one currently
    used in the code. The reformulation has been presented in the introduction
    of this code and makes use of two sets of variables, i.e.:

    * :math:`y_{jt}` : setup variables
    * :math:`z_{jtr}`: production variables, indicating the amount of production of item *j* produced in period *t* to satisfy the demand of period *r*. Obviously, :math:`z_{jtr} = 0` for all *t>r*.


    """

    def __init__(self, inp):
        y_ilo = []
        z_ilo = []

        cpx = cplex.Cplex()
        cpx.objective.set_sense(cpx.objective.sense.minimize)
        #  cpx.set_results_stream(None)
        #  cpx.set_log_stream(None)

        #  create variables y_jt
        for j in range(inp.nI):
            y_ilo.append([])
            for t in range(inp.nP):
                varName = "y." + str(j) + "." + str(t)
                y_ilo[j].append(cpx.variables.get_num())
                cpx.variables.add(obj   = [inp.f[j][t]],
                                  lb    = [0],
                                  ub    = [1],
                                  types = ["B"],
                                  names = [varName])
        #  create variables z_jts
        for j in range(inp.nI):
            z_ilo.append([])
            for t in range(inp.nP):
                z_ilo[j].append([])
                for r in range(inp.nP):
                    varName = "z." + str(j) + "." + str(t) + "." + str(r)
                    z_ilo[j][t].append(cpx.variables.get_num())
                    cpx.variables.add(obj   = [(r-t)*inp.h[j][t]],
                                      lb    = [0.0],
                                      ub    = [cplex.infinity],
                                      types = ["C"],
                                      names = [varName])

        #  demand constraints
        for j in range(inp.nI):
            for r in range(inp.nP):
                index = [z_ilo[j][t][r] for t in range(r+1)]
                value = [1.0]*(r+1)
                demand_constraint = cplex.SparsePair(ind=index, val=value)
                cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                           senses   = ["E"],
                                           rhs      = [inp.d[j][r]])

        #  capacity constraint
        for t in range(inp.nP):
            index = []
            value = []
            for j in range(inp.nI):
                index += [z_ilo[j][t][r] for r in range(t,inp.nP)]
                value += [inp.a[j][t]]*(inp.nP-t)
                index += [y_ilo[j][t]]
                value += [inp.m[j][t]]
            capacity_constraint = cplex.SparsePair(ind=index,val=value)
            cpx.linear_constraints.add(lin_expr  = [capacity_constraint],
                                       senses    = ["L"],
                                       rhs       = [inp.cap[t]])

        #  logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    index = [z_ilo[j][t][r], y_ilo[j][t]]
                    value = [1.0, -inp.d[j][r]]
                    logic_constraint = cplex.SparsePair(ind =index, val=value)
                    cpx.linear_constraints.add(lin_expr = [logic_constraint],
                                               senses   = ["L"],
                                               rhs      = [0.0])
        #  cumulative logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                index = [z_ilo[j][t][r] for r in range(t,inp.nP)]
                value = [1.0]*(inp.nP-t)
                index += [y_ilo[j][t]]
                value += [-inp.max_prod[j][t]]
                logic_constraint_cum = cplex.SparsePair(ind =index, val=value)
                cpx.linear_constraints.add(lin_expr = [logic_constraint_cum],
                                           senses   = ["L"],
                                           rhs      = [0.0])


        #  set to zero unused variables
        for j in range(inp.nI):
            for r in range(inp.nP-1):
                for t in range(r+1, inp.nP):
                    cpx.variables.set_upper_bounds(z_ilo[j][t][r], 0.0)

        self.cpx   = cpx
        self.y_ilo = y_ilo
        self.z_ilo = z_ilo


class WorkerLP:
    """
    .. class:WorkerLP()

    Define and solve the subproblem. We initilize the subproblem with the right
    hand side values of the constraints to zero, since we assume the initial
    values of :math:`y_{jt}` to be equal to zero. Next, within the
    :meth:`WorkerLP.separate()` function, we define the rhs values to the 
    correct values, depending on the solution obtained from the master.

    Cplex requires the presolve reductions to be turned off, using::

        cpx.parameters.preprocessing.reduce.set(0)

    In addition, we need to ensure that the LP is *not* solved using the
    interior point method, otherwise dual values won't be available. We can
    either use primal simplex or dual simplex.

    .. note ::

        The subproblem constraints should be defined with a name, in
        order to be able to recall the precise name of each constraint when we 
        want to obtain the dual values. Briefly, we need to:

        * define a unique name for each constraint, e.g., ``capacity.t`` for each t
        * store such names in a vector of names, e.g,::

            lCapacity = ["capacity." + str(t) for t in range(inp.nP)]

        * get the dual values using::

           dCapacity = cpx.solution.get_dual_values(lCapacity)

    """
    def __init__(self, inp):

        x_ilo = []
        s_ilo = []
        sI    = []

        cpx = cplex.Cplex()
        cpx.set_results_stream(None)
        cpx.set_log_stream(None)


        # Turn off the presolve reductions and set the CPLEX optimizer
        # to solve the worker LP with primal simplex method.
        cpx.parameters.preprocessing.presolve.set(
                cpx.parameters.preprocessing.presolve.values.off)
        cpx.parameters.preprocessing.reduce.set(0)
        #  cpx.parameters.lpmethod.set(cpx.parameters.lpmethod.values.primal)
        cpx.parameters.lpmethod.set(cpx.parameters.lpmethod.values.dual)

        cpx.objective.set_sense(cpx.objective.sense.minimize)

        for j in range(inp.nI):
            x_ilo.append([])
            for t in range(inp.nP):
                varName = "x." + str(j) + "." + str(t)
                x_ilo[j].append(cpx.variables.get_num())
                cpx.variables.add(obj   = [inp.c[j][t]],
                                  lb    = [0.0],
                                  ub    = [cplex.infinity],
                                  #  types = ["C"],
                                  names = [varName])

        for j in range(inp.nI):
            s_ilo.append([])
            for t in range(inp.nP):
                varName = "s." + str(j) + "." + str(t)
                s_ilo[j].append(cpx.variables.get_num())
                cpx.variables.add(obj   = [inp.h[j][t]],
                                  lb    = [0.0],
                                  ub    = [cplex.infinity],
                                  #  types = ["C"],
                                  names = [varName])


        #  this variables are used to allow for unlimited initial inventory
        #  for j in range(inp.nI):
        #      varName = "sI." + str(j)
        #      sI.append(cpx.variables.get_num())
        #      cpx.variables.add(obj   = [1000],
        #                        lb    = [0.0],
        #                        ub    = [cplex.infinity],
        #                        #  ub    = [0.0],
        #                        #  types = ["C"],
        #                        names = [varName])

        #  capacity constraints
        for t in range(inp.nP):
            index = [x_ilo[j][t] for j in range(inp.nI)]
            value = [inp.a[j][t] for j in range(inp.nI)]

            capacity_constraint = cplex.SparsePair(ind=index, val=value)
            constrName = "capacity." + str(t)
            cpx.linear_constraints.add(lin_expr = [capacity_constraint],
                                       senses   = ["L"],
                                       rhs      = [0.0],
                                       names    = [constrName])
        #  demand constraints
        #  note: If we want to include unlimited initial inventory, change this
        #  constraint, including sI
        for j in range(inp.nI):
            #  first period
            #  index = [x_ilo[j][0], sI[j], s_ilo[j][0]]
            #  value = [1.0, 1.0, -1.0]
            index = [x_ilo[j][0], s_ilo[j][0]]
            value = [1.0, -1.0]
#
            demand_constraint = cplex.SparsePair(ind=index, val=value)
            constrName = "demand." + str(j) + ".0"
            cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                       senses   = ["G"],
                                       rhs      = [inp.d[j][0]],
                                       names    = [constrName])
            #  periods 2 to T-1
            for t in range(1,inp.nP-1):
                index = [x_ilo[j][t], s_ilo[j][t-1], s_ilo[j][t]]
                value = [1.0, 1.0, -1.0]
                demand_constraint = cplex.SparsePair(ind=index, val=value)
                constrName = "demand." + str(j) + "." + str(t)
                cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                           senses   = ["G"],
                                           rhs      = [inp.d[j][t]],
                                           names    = [constrName])

            #  last period
            index = [x_ilo[j][inp.nP-1], s_ilo[j][inp.nP-2]]
            value = [1.0, 1.0]
            demand_constraint = cplex.SparsePair(ind=index, val=value)
            constrName = "demand." + str(j) + "." + str(inp.nP-1)
            cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                       senses   = ["G"],
                                       rhs      = [inp.d[j][inp.nP-1]],
                                       names    = [constrName])


        #  logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                constrName = "logic." + str(j) + "." + str(t)
                logic_constraint = cplex.SparsePair(ind = [x_ilo[j][t]],
                                                    val=[1.0])
                cpx.linear_constraints.add(lin_expr = [logic_constraint],
                                           senses   = ["L"],
                                           rhs      = [0.0],
                                           names    = [constrName])


        # define labels for constraints
        lCapacity = ["capacity." + str(t) for t in range(inp.nP)]
        lLogic    = ["logic." + str(j) + "." + str(t) for j in range(inp.nI)
        for t in range(inp.nP)]
        lDemand   = ["demand." + str(j) + "." + str(t) for j in range(inp.nI)
        for t in range(inp.nP)]

        self.cpx       = cpx
        self.x_ilo     = x_ilo
        self.s_ilo     = s_ilo
        self.sI        = sI
        self.lCapacity = lCapacity
        self.lLogic    = lLogic
        self.lDemand   = lDemand

    def separate(self, inp, ySol, zHat, y_ilo, z_ilo):
        """
        .. method:separate()

        This is the *old* implementation, based on the standard Lot Sizing
        formulation. Here is were Benders cut is obtained and passed to the master.

        The following steps describe the algorithm:

        * update rhs values of the subproblem, i.e., using the current optimal solution of the master problem :math:`y^*`
        * solve the subproblem
        * get the dual values :math:`\lambda, \\nu`
        * generate cut (lhs and rhs) and store them in a constraint structure
        * pass the cut to the master

        """

        cutType   = 0
        cpx       = self.cpx
        x_ilo     = self.x_ilo
        s_ilo     = self.s_ilo
        sI        = self.sI
        lCapacity = self.lCapacity
        lLogic    = self.lLogic
        lDemand   = self.lDemand

        #  update rhs values : logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                constrName = "logic." + str(j) + "." + str(t)
                #  if ySol[j][t] >= (1.0-_EPSI):
                rhsValue = inp.max_prod[j][t]*ySol[j][t]
                #  else:
                    #  rhsValue = 0.0
                cpx.linear_constraints.set_rhs(constrName, rhsValue)

        #  update rhs values : capacity constraints
        for t in range(inp.nP):
            sumY = sum([inp.m[j][t]*ySol[j][t] for j in range(inp.nI)])
            constrName = "capacity." + str(t)
            cpx.linear_constraints.set_rhs(constrName, inp.cap[t]-sumY)

        cpx.solve()


        #  print("STATUS = ", cpx.solution.status[cpx.solution.get_status()])
        #  print("z_sub  = ", cpx.solution.get_objective_value())
        #  xSol = cpx.solution.get_values(x_ilo)
        #  sISol = cpx.solution.get_values(sI)
        #  print("X : ", xSol)
        #  print("SI : ", sISol)
        #  vars and values for the cut lhs
        cutVars = []
        cutVals = []

        zSub = cpx.solution.get_objective_value()
        if cpx.solution.get_status() == cpx.solution.status.infeasible:
            cutType = 1
            #  add extreme ray
            #  cpx.solution.advanced.get_ray()
            #  note: we cannot get the rays without explicitly formulating the
            #  dual. We use Farkas certificate
            farkas, pp = cpx.solution.advanced.dual_farkas()
            #  print("Farkas are : ", farkas, " ** ", pp)
            progr     = inp.nP
            index     = range(progr)
            fCapacity = [farkas[i] for i in index]
            index     = range(progr,progr+inp.nP*inp.nI)
            fDemand   = [farkas[i] for i in index]
            progr    += inp.nP*inp.nI
            index     = range(progr, progr+inp.nP*inp.nI)
            fLogic    = [farkas[i] for i in index]

            #  rhsConstr = cpx.linear_constraints.get_rhs()
            #  progr = 0
            #  for i in enumerate(rhsConstr):
            #      print("i is ", i, " with farkas = ", farkas[progr])
            #      progr += 1

            #  feasibility cut
            cutRhs  = 0.0
            for j in range(inp.nI):
                for t in range(inp.nP):
                    progr = j*inp.nP + t
                    coeff = inp.max_prod[j][t]*fLogic[progr] - \
                    inp.m[j][t]*fCapacity[t]

                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)
                    cutRhs -= inp.d[j][t]*fDemand[progr]

            for t in range(inp.nP):
                cutRhs -= inp.cap[t]*fCapacity[t]
        elif cpx.solution.get_status() == cpx.solution.status.optimal:
            cutType = 2
            zSub = cpx.solution.get_objective_value()
            #  the solution is feasible (it will always be, due to vars sI)
            dCapacity = cpx.solution.get_dual_values(lCapacity)
            dLogic    = cpx.solution.get_dual_values(lLogic)
            dDemand   = cpx.solution.get_dual_values(lDemand)


            #  method 1: generalized benders
            cutRhs = -zSub
            progr  = 0
            for j in range(inp.nI):
                for t in range(inp.nP):
                    coeff = -inp.m[j][t]*dCapacity[t] +\
                    dLogic[progr]*inp.max_prod[j][t]

                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)

                    cutRhs += coeff*ySol[j][t]
                    progr += 1

            cutVars.append(z_ilo)
            cutVals.append(-1.0)


            # method 2: the "standard" benders
            #  cutRhs = 0.0
            #  progr  = 0
            #  for j in range(inp.nI):
            #      for t in range(inp.nP):
            #          coeff = inp.max_prod[j][t]*dLogic[progr] - \
            #          inp.m[j][t]*dCapacity[t]
            #
            #          cutVars.append(y_ilo[j][t])
            #          cutVals.append(coeff)
            #
            #          cutRhs -= inp.d[j][t]*dDemand[progr]
            #
            #          progr += 1
            #
            #  cutVars.append(z_ilo)
            #  cutVals.append(-1.0)
            #  for t in range(inp.nP):
            #      cutRhs -= inp.cap[t]*dCapacity[t]


        #  return cut and type
        cutLhs = cplex.SparsePair(ind=cutVars, val=cutVals)

        self.cutLhs = cutLhs
        self.cutRhs = cutRhs
        self.zSub   = zSub

        return cutType


class WorkerLPReformulation:
    """
    .. class:WorkerLPReformulation()

        Formulation of the subproblem (the worker) using the SPL reformulation.
        This is what is currently used in Benders' algorithm.
    """
    def __init__(self, inp):
        """

        **Note**: The rhs value of both the logic constraints and the capacity
        constraints, as well as the cumulative capacity constraints, will be
        update during the separation phase, implemented in
        :meth:`WorkerLPReformulation.separate()`, since the value of variables
        :math:`y_{jt}` will be provided by the master. Here we thus initialize
        these rhs values to either 0 or a convenient value.

        Note the use of labels for constraints. These labels are needed to
        access the dual values of those constraints in the separation phase.

        """

        z_ilo = []

        cpx = cplex.Cplex()
        cpx.set_results_stream(None)
        cpx.set_log_stream(None)

        # Turn off the presolve reductions and set the CPLEX optimizer
        # to solve the worker LP with primal simplex method.
        cpx.parameters.preprocessing.presolve.set(
                cpx.parameters.preprocessing.presolve.values.off)
        cpx.parameters.preprocessing.reduce.set(0)
        #  cpx.parameters.lpmethod.set(cpx.parameters.lpmethod.values.primal)
        cpx.parameters.lpmethod.set(cpx.parameters.lpmethod.values.dual)
        cpx.objective.set_sense(cpx.objective.sense.minimize)

        #  create variables z_jts
        for j in range(inp.nI):
            z_ilo.append([])
            for t in range(inp.nP):
                z_ilo[j].append([])
                for r in range(inp.nP):
                    varName = "z." + str(j) + "." + str(t) + "." + str(r)
                    z_ilo[j][t].append(cpx.variables.get_num())
                    cpx.variables.add(obj   = [(r-t)*inp.h[j][t]],
                                      lb    = [0.0],
                                      #  ub    = [cplex.infinity],
                                      #  ub    = [inp.d[j][r]],
                                      names = [varName])


        #  capacity constraint
        for t in range(inp.nP):
            index = []
            value = []
            for j in range(inp.nI):
                index += [z_ilo[j][t][r] for r in range(t,inp.nP)]
                value += [inp.a[j][t]]*(inp.nP-t)
            constrName = "capacity." + str(t)
            capacity_constraint = cplex.SparsePair(ind=index,val=value)
            cpx.linear_constraints.add(lin_expr  = [capacity_constraint],
                                       senses    = ["L"],
                                       rhs       = [inp.cap[t]],
                                       names     = [constrName])

        #  demand constraints
        for j in range(inp.nI):
            for r in range(inp.nP):
                index = [z_ilo[j][t][r] for t in range(r+1)]
                value = [1.0]*(r+1)
                constrName = "demand." + str(j) + "." + str(r)
                demand_constraint = cplex.SparsePair(ind=index, val=value)
                cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                           senses   = ["E"],
                                           rhs      = [inp.d[j][r]],
                                           names    = [constrName])
        #  logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    index = [z_ilo[j][t][r]]
                    value = [1.0]
                    logic_constraint = cplex.SparsePair(ind =index, val=value)
                    constrName = "logic."+str(j)+"."+str(t)+"."+str(r)
                    cpx.linear_constraints.add(lin_expr = [logic_constraint],
                                               senses   = ["L"],
                                               rhs      = [0.0],
                                               names    = [constrName])

        #  cumulative logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                index = [z_ilo[j][t][r] for r in range(t,inp.nP)]
                value = [1.0]*(inp.nP-t)
                constrName = "cumLogic."+str(j)+"."+str(t)
                logic_constraint_cum = cplex.SparsePair(ind =index, val=value)
                cpx.linear_constraints.add(lin_expr = [logic_constraint_cum],
                                           senses   = ["L"],
                                           rhs      = [0.0],
                                           names    = [constrName])


        #  set to zero unused variables
        for j in range(inp.nI):
            for r in range(inp.nP-1):
                for t in range(r+1, inp.nP):
                    cpx.variables.set_upper_bounds(z_ilo[j][t][r], 0.0)


        #  for j in range(inp.nI):
            #  for r in range(inp.nP):
                #  for t in range(r):
                    #  cpx.variables.set_upper_bounds(z_ilo[j][t][r], inp.d[j][r])

        # define labels for constraints (used in separation to get dual values)
        lCapacity = ["capacity." + str(t) for t in range(inp.nP)]
        lDemand   = ["demand." + str(j) + "." + str(t) for j in range(inp.nI)
        for t in range(inp.nP)]
        lLogic    = ["logic." + str(j) + "." + str(t) + "." + str(r) \
        for j in range(inp.nI) for t in range(inp.nP) for r in range(t,inp.nP)]
        lcumLogic = ["cumLogic." + str(j) + "." + str(t) for j in range(inp.nI)
        for t in range(inp.nP)]

        self.cpx       = cpx
        self.z_ilo     = z_ilo
        self.lCapacity = lCapacity
        self.lDemand   = lDemand
        self.lLogic    = lLogic
        self.lcumLogic = lcumLogic

    def formulateDual(self, cpx, inp, ySol, w_ilo, l_ilo, v_ilo, e_ilo):
        """
        Formulation of the subproblem dual. This is no longer used, but I leave
        it here for the sake of completeness. I checked that the formulation is
        correct, since the objective function value of the optimal dual is
        identical to that of the optimal primal.

        """

        for j in range(inp.nI):
            w_ilo.append([])
            for t in range(inp.nP):
                w_ilo[j].append(cpx.variables.get_num())
                varName = "w." + str(j) + "." + str(t)
                cpx.variables.add(obj   = [inp.d[j][t]],
                                  lb    = [-cplex.infinity],
                                  ub    = [cplex.infinity],
                                  names = [varName])

        for t in range(inp.nP):
            l_ilo.append(cpx.variables.get_num())
            varName = "l." + str(t)
            coeff = inp.cap[t]
            for j in range(inp.nI):
                coeff -= inp.m[j][t]*ySol[j][t]
            cpx.variables.add(obj   = [coeff],
                              lb    = [-cplex.infinity],
                              ub    = [0.0],
                              names = [varName])

        for j in range(inp.nI):
            v_ilo.append([])
            for t in range(inp.nP):
                v_ilo[j].append(cpx.variables.get_num())
                varName = "v." + str(j) + "." + str(t)
                cpx.variables.add(obj   = [inp.max_prod[j][t]*ySol[j][t]],
                                  lb    = [-cplex.infinity],
                                  ub    = [0.0],
                                  names = [varName])

        for j in range(inp.nI):
            e_ilo.append([])
            for t in range(inp.nP):
                e_ilo[j].append([])
                for r in range(inp.nP):
                    e_ilo[j][t].append(cpx.variables.get_num())
                    varName = "e." + str(j) + "." + str(t) + "." + str(r)
                    cpx.variables.add(obj    = [inp.d[j][r]*ySol[j][t]],
                                      lb     = [-cplex.infinity],
                                      ub     = [0.0],
                                      names  = [varName])

        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t+1):
                    cpx.variables.set_upper_bounds(e_ilo[j][t][r],0.0)
                    cpx.variables.set_lower_bounds(e_ilo[j][t][r],0.0)

        cpx.objective.set_sense(cpx.objective.sense.maximize)
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    constrName = "dual." + str(j) + "." + str(t) + "." + str(r)
                    index = [w_ilo[j][r], l_ilo[t], v_ilo[j][t], e_ilo[j][t][r]]
                    value = [1.0, inp.a[j][t], 1.0, 1.0]
                    dual_constraint = cplex.SparsePair(ind=index,val=value)
                    cpx.linear_constraints.add(lin_expr = [dual_constraint],
                                               senses   = ["L"],
                                               rhs      = [(r-t)*inp.h[j][t]],
                                               names    = [constrName])

    def paretoOptimal(self, inp, ySol, zDual, zHat):
        """
        Get pareto optimal cuts. See paper for an explanation. In principle,
        this should be helpful, since the dual problem is degenerate and,
        consequently, multiple optimal solutions should exist. This means that
        the "right" selection of dual values should make a difference in the
        strength of the cut. However, the extra effort required to get pareto
        optimal cuts does not seem to be compensated by the minor improvement
        obtained here.

        """
        cpx = cplex.Cplex()

        cpx.set_results_stream(None)
        cpx.set_log_stream(None)
        y0 = [ [0.5]*inp.nP for i in range(inp.nI)]
        #  y0 = yRef
        #  y0 = ySol
        w_ilo = []
        l_ilo = []
        v_ilo = []
        e_ilo = []
        self.formulateDual(cpx, inp, y0, w_ilo, l_ilo, v_ilo, e_ilo)

        index = [w_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP)]
        value = [inp.d[j][t] for j in range(inp.nI) for t in range(inp.nP)]
        index += [l_ilo[t] for t in range(inp.nP)]
        coeffs = []
        for t in range(inp.nP):
            aux = inp.cap[t]
            for j in range(inp.nI):
                aux -= ySol[j][t]*inp.m[j][t]
            coeffs.append(aux)
        value += coeffs
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    index += [e_ilo[j][t][r]]
                    value += [ySol[j][t]*inp.d[j][r]]

        index += [v_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP)]
        value += [inp.max_prod[j][t]*ySol[j][t] for j in range(inp.nI) for t in range(inp.nP)]

        obj_constraint = cplex.SparsePair(ind=index, val=value)
        cpx.linear_constraints.add(lin_expr  = [obj_constraint],
                                   senses    = ["E"],
                                   rhs       = [zDual],
                                   names     = ["obj_constraint"])


        #  cpx.write("pareto.lp")
        cpx.solve()

        dDemand = [cpx.solution.get_values(w_ilo[j][t]) for j in range(inp.nI)
        for t in range(inp.nP)]
        dCapacity = [cpx.solution.get_values(l_ilo[t]) for t in range(inp.nP)]
        dcumLogic = [cpx.solution.get_values(v_ilo[j][t]) for j in range(inp.nI)
        for t in range(inp.nP)]
        dLogic = [cpx.solution.get_values(e_ilo[j][t][r]) for j in
        range(inp.nI) for t in range(inp.nP) for r in range(t,inp.nP)]
        
        return dDemand, dCapacity, dcumLogic, dLogic

    def solveDual(self, inp, ySol, zHat):
        """
        Solve the dual problem (i.e., an explicit way of obtained the dual
        values needed to create the optimality cut. It is no longer used.

        """

        cpx = cplex.Cplex()
        cpx.set_results_stream(None)
        cpx.set_log_stream(None)
        w_ilo = []
        l_ilo = []
        v_ilo = []
        e_ilo = []
        self.formulateDual(cpx, inp, ySol, w_ilo, l_ilo, v_ilo, e_ilo)

        cpx.solve()
        zDual = cpx.solution.get_objective_value()

        dDemand, dCapacity, dcumLogic, dLogic = self.paretoOptimal(inp, ySol,
        zDual, zHat)

        return dDemand, dCapacity, dcumLogic, dLogic



    def separate(self, inp, ySol, zHat, y_ilo, z_ilo_m):
        """
        .. method:separate()

        This is the separation scheme for the WorkerLPReformulation, i.e., the
        SPL reformulation. Here is were Benders cut is obtained and passed to 
        the master.

        The following steps describe the algorithm:

        * update rhs values of the subproblem, i.e., using the current optimal solution of the master problem :math:`y^*`
        * solve the subproblem
        * get the dual values :math:`\lambda, \\nu`
        * generate cut (lhs and rhs) and store them in a constraint structure
        * pass the cut to the master

        Advanced cplex functions are used here, e.g.::
        
            farkas, pp = cpx.solution.advanced.dual_farkas()

        to get Farkas certificates used to get the extreme rays, and::

            dCapacity = cpx.solution.get_dual_values(lCapacity)

        to get the dual values, which are used in the creation of optimality
        cuts. 

        """

        cutType   = 0
        cpx       = self.cpx
        z_ilo     = self.z_ilo
        lCapacity = self.lCapacity
        lDemand   = self.lDemand
        lLogic    = self.lLogic
        lcumLogic = self.lcumLogic

        cpx.set_results_stream(None)
        cpx.set_log_stream(None)

        #  update rhs values : logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    constrName = "logic." + str(j) + "." + str(t) + "." + str(r)
                    rhsValue = inp.d[j][r]*ySol[j][t]
                    cpx.linear_constraints.set_rhs(constrName, rhsValue)

        #  update rhs values : capacity constraints
        for t in range(inp.nP):
            sumY = sum([inp.m[j][t]*ySol[j][t] for j in range(inp.nI)])
            constrName = "capacity." + str(t)
            cpx.linear_constraints.set_rhs(constrName, inp.cap[t]-sumY)

        #  update cumulative capacity
        for j in range(inp.nI):
            for t in range(inp.nP):
                constrName = "cumLogic." + str(j) + "." + str(t)
                rhsValue = inp.max_prod[j][t]*ySol[j][t]
                cpx.linear_constraints.set_rhs(constrName, rhsValue)


        cpx.solve()
        zSub = cpx.solution.get_objective_value()


        #  print("STATUS = ", cpx.solution.status[cpx.solution.get_status()])
        #  print("z_sub  = ", cpx.solution.get_objective_value())
        #  vars and values for the cut lhs
        cutVars = []
        cutVals = []

        if cpx.solution.get_status() == cpx.solution.status.infeasible:
            cutType = 1
            #  add extreme ray using Farkas certificate
            farkas, pp = cpx.solution.advanced.dual_farkas()
            #  print("[",len(farkas),"] Farkas are : ", farkas, " ** ", pp)
            progr     = inp.nP
            index     = range(progr)
            fCapacity = [farkas[i] for i in index]
            index     = range(progr,progr+inp.nP*inp.nI)
            fDemand   = [farkas[i] for i in index]
            progr    += inp.nP*inp.nI
            nr        = int( inp.nI*(inp.nP*(inp.nP+1))/2)
            index     = range(progr, progr+nr)
            fLogic    = [farkas[i] for i in index]
            progr    += nr
            index     = range(progr, progr+inp.nI*inp.nP)
            fcumLogic = [farkas[i] for i in index]


            #  feasibility cut
            cutRhs  = 0.0
            counter = 0
            for j in range(inp.nI):
                for t in range(inp.nP):
                    progr = j*inp.nP + t
                    coeff = inp.max_prod[j][t]*fcumLogic[progr] - \
                    inp.m[j][t]*fCapacity[t]

                    for r in range(t, inp.nP):
                        coeff += fLogic[counter]*inp.d[j][r]
                        counter += 1

                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)
                    cutRhs -= inp.d[j][t]*fDemand[progr]

            for t in range(inp.nP):
                cutRhs -= inp.cap[t]*fCapacity[t]


        elif cpx.solution.get_status() == cpx.solution.status.optimal:
            cutType = 2
            zSub = cpx.solution.get_objective_value()
            #  print("... ... zSUB = ", zSub)
            dCapacity = cpx.solution.get_dual_values(lCapacity)
            dDemand   = cpx.solution.get_dual_values(lDemand)
            dLogic    = cpx.solution.get_dual_values(lLogic)
            dcumLogic = cpx.solution.get_dual_values(lcumLogic)

            # method 2: the "standard" benders
            cutRhs = 0.0
            progr  = 0
            counter = 0
            for j in range(inp.nI):
                for t in range(inp.nP):
                    coeff = inp.max_prod[j][t]*dcumLogic[progr] - \
                    inp.m[j][t]*dCapacity[t]

                    for r in range(t, inp.nP):
                        coeff += dLogic[counter]*inp.d[j][r]
                        counter += 1

                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)

                    cutRhs -= inp.d[j][t]*dDemand[progr]

                    progr += 1

            cutVars.append(z_ilo_m)
            cutVals.append(-1.0)
            for t in range(inp.nP):
                cutRhs -= inp.cap[t]*dCapacity[t]


        #  return cut and type
        cutLhs = cplex.SparsePair(ind=cutVars, val=cutVals)

        self.cutLhs = cutLhs
        self.cutRhs = cutRhs
        self.zSub   = zSub

        return cutType

def findNext(j, t, inp):
    #  print("Item ", j, " from period ", t)
    tp = t
    maxProd = inp.max_prod[j][t]
    #  print("MAX = ", maxProd)
    dCum = inp.d[j][t]
    while dCum <= maxProd:
        #  print("dCum vs maxProd ", dCum, " ", maxProd, " up to period ",tp)
        tp += 1
        if tp == inp.nP:
            return tp
        dCum += inp.d[j][tp]
    return tp

def createMaster(inp, cpx):
    """
    Create benders master problem. In reality, this is a *relaxation* of the
    master, to which we progressively add cuts.

    .. note ::

        The master problem provides a **lower bound** to the optimal
        solution.

    Here we try different methods to tighen the lower bound.

    """

    global z_ilo
    global y_ilo
    cpx.objective.set_sense(cpx.objective.sense.minimize)

    #  create variables y_jt
    for j in range(inp.nI):
        y_ilo.append([])
        for t in range(inp.nP):
            varName = "y." + str(j) + "." + str(t)
            y_ilo[j].append(cpx.variables.get_num())
            cpx.variables.add(obj   = [inp.f[j][t]],
                              lb    = [0],
                              ub    = [1],
                              types = ["B"],
                              names = [varName])

    #  z_ilo.append(cpx.variables.get_num())
    z_ilo = cpx.variables.get_num()
    cpx.variables.add(obj   = [1.0],
                      lb    = [0.0],
                      ub    = [cplex.infinity],
                      types = ["C"],
                      names = ["zHat"])

    #  for j in range(inp.nI):
    #      if inp.d[j][0] > 0.0:
    #          cpx.variables.set_lower_bounds(y_ilo[j][0], 1)
    #          #  for t in range(inp.nP-1):
    #          for t in range(1):
    #              #  hop constraint
    #              tp = findNext(j,t, inp)
    #              print("... For item ",j," we go from ",t," to ",tp)
    #              #  input("...aka")
    #              if tp < inp.nP:
    #                  index = [y_ilo[j][t] for t in range(t+1,tp+1)]
    #                  value = [1.0]*len(range(t,tp))
    #                  hop_constraint = cplex.SparsePair(ind=index,val=value)
    #                  cpx.linear_constraints.add(lin_expr = [hop_constraint],
    #                                             senses   = ["G"],
    #                                             rhs      = [1])

    for j in range(inp.nI):
        index = [y_ilo[j][t] for t in range(inp.nP)]
        value = [inp.max_prod[j][t] for t in range(inp.nP)]
        c_constraint = cplex.SparsePair(ind=index,val=value)
        cpx.linear_constraints.add(lin_expr = [c_constraint],
                                   senses   = ["G"],
                                   rhs      = [inp.dcum[j][0]])





def barrierInit(inp):

    cpx = cplex.Cplex()
    cpx.objective.set_sense(cpx.objective.sense.maximize)
    cpx.parameters.lpmethod.set(cpx.parameters.lpmethod.values.barrier)
    #  cpx.parameters.solutiontype.set(2)
    cpx.parameters.barrier.crossover.set(-1) # no crossover
    yB = []
    for j in range(inp.nI):
        yB.append([])
        for t in range(inp.nP):
            yB[j].append(cpx.variables.get_num())
            cpx.variables.add(obj = [1.0],
                              lb  = [0.0],
                              ub  = [1.0])
    cpx.solve()

    #  print("BARRIER NO CROSSOVER SOLVED : ", cpx.solution.get_objective_value())
    ySol = []
    for j in range(inp.nI):
        ySol.append(cpx.solution.get_values(yB[j]))
    #
    #  print(ySol)
    #  input(" ... barrier ... ")
    return ySol

def inOutCycle(cpx, worker, y_ilo, z_ilo, inp, globalProgr):

    global inout
    _lambda = 0.1
    _alpha  = 0.9
    #  progr   = cpx.linear_constraints.get_num()


    cpx.set_problem_type(cpx.problem_type.LP)
    print("Problem type is ", cpx.problem_type[cpx.get_problem_type()])

    #  now solve simple problem to get interior point
    yt = barrierInit(inp)
    iter           = 0
    stopping       = 0
    noImprovement  = 0
    activateKelley = 0
    bestLP         = 0.0
    while not stopping:
        cpx.solve()
        zLP = cpx.solution.get_objective_value()
        yLP = []
        for j in range(inp.nI):
            yLP.append(cpx.solution.get_values(y_ilo[j]))
        zHat = cpx.solution.get_values(z_ilo)

        #  print("iter ", iter, " with z = ", zLP, " vs best ", bestLP)
        #  print("   .. current status : No Improv = ", noImprovement, \
                     #  " Kelley = ", activateKelley, " lambda = ", _lambda)
        #  if zLP > bestLP:
        if (zLP - bestLP) > 0.5:
            bestLP = zLP
            with open(lbSummary,"a") as ff:
                ff.write("{0:20.5f} {1:20.5} \n".format(bestLP,
                time()-startTime))
            noImprovement = 0
        else:
            noImprovement += 1

        if noImprovement >= 5:
            if activateKelley == 0:
                activateKelley = 1
                noImprovement = 0
                _lambda = 1.0
            else:
                stopping = 1
            continue

        yt = [ [_alpha*yt[j][t] + (1-_alpha)*yLP[j][t] for t in
        range(inp.nP)] for j in range(inp.nI)]
        ySep = [ [_lambda*yLP[j][t] + (1-_lambda)*yt[j][t] for t in
        range(inp.nP)] for j in range(inp.nI)]
        #  ySep = [[ ySep[j][t] + 2*_EPSI for t in
        #  range(inp.nP)] for j in range(inp.nI)]


        cutType = worker.separate(inp, ySep, zHat, y_ilo, z_ilo)
        iter += 1
        if iter % 5 == 0:
            check = [i for i in inout]
            #  print("checking these constraints ", check)

            slacks = cpx.solution.get_linear_slacks(check)
            #  nrConstr = cpx.linear_constraints.get_num()
            remove = [check[i] for i in range(len(slacks)) if slacks[i] > _EPSI]
            #  print("removing ", remove)
            cpx.linear_constraints.delete(remove)
            for i in remove:
                inout.remove(i)
            #  print("Removed ", remove)
            #  print("This is inout = ", inout)

        if cutType > 0:
            nrConstr = cpx.linear_constraints.get_num()
            constrName = "inout." + str(globalProgr)
            #  print("adding ", constrName)
            inout.append(constrName)
            globalProgr += 1
            cpx.linear_constraints.add(lin_expr = [worker.cutLhs],
                                       senses   = "L",
                                       rhs      = [worker.cutRhs],
                                       names    = [constrName])

    cpx.solve()

    check = [i for i in inout]
    slacks = cpx.solution.get_linear_slacks(check)
    #  print("[", nrConstr,"] SLACKS  == ", slacks)
    remove = [check[i] for i in range(len(slacks)) if slacks[i] > _EPSI]
    #  print("removing ", remove)
    cpx.linear_constraints.delete(remove)
    for i in remove:
        inout.remove(i)

    cpx.solve()
    #  print("OBJ AFTER REMOVING SLACKS = ", cpx.solution.get_objective_value())
    nrConstr = cpx.linear_constraints.get_num()
    print(" ... now constr is ", nrConstr)

    return globalProgr

def inOutCycle2(cpx, worker, y_ilo, z_ilo, inp, globalProgr):

    global inout
    _lambda = 0.1
    _alpha  = 0.9
    #  progr   = cpx.linear_constraints.get_num()

    cpx.set_problem_type(cpx.problem_type.LP)
    #  print("Problem type is ", cpx.problem_type[cpx.get_problem_type()])

    #  now solve simple problem to get interior point
    yIn = [[0.9]*inp.nP for i in range(inp.nI)]
    #  print("yIn = ", yIn)
    iter           = 0
    stopping       = 0
    noImprovement  = 0
    activateKelley = 0
    bestLP         = 0.0
    while not stopping:
        cpx.solve()
        zLP = cpx.solution.get_objective_value()
        print("zLP(",iter,") = ", zLP)
        yOut = []
        for j in range(inp.nI):
            yOut.append(cpx.solution.get_values(y_ilo[j]))
        #  print("yOut (LP) = ", yOut)
        zHat = cpx.solution.get_values(z_ilo)

        #  print("iter ", iter, " with z = ", zLP, " vs best ", bestLP)
        #  print("   .. current status : No Improv = ", noImprovement, \
                     #  " Kelley = ", activateKelley, " lambda = ", _lambda)
        #  if zLP > bestLP:
        if (zLP - bestLP) > 0.01:
            bestLP = zLP
            noImprovement = 0
        else:
            noImprovement += 1

        if noImprovement >= 5:
            if activateKelley == 0:
                activateKelley = 1
                noImprovement = 0
                _lambda = 1.0
            else:
                stopping = 1
            continue

        ySep = [ [_lambda*yOut[j][t] + (1-_lambda)*yIn[j][t] for t in
        range(inp.nP)] for j in range(inp.nI)]
        ySep = [[ ySep[j][t] + 2*_EPSI for t in
        range(inp.nP)] for j in range(inp.nI)]


        cutType = worker.separate(inp, ySep, zHat, y_ilo, z_ilo)
        flatten =  [item for sublist in ySep for item in sublist]
        a = [float(worker.cutLhs.val[i]) for i in range(inp.nI*inp.nP)]
        lhsSum = sum([a[i]*flatten[i] for i in range(inp.nI*inp.nP)])
        #  print("LhsSum = ", lhsSum , " vs ", worker.cutRhs)
        #  print(lhsSum <= worker.cutRhs)
        violated = (lhsSum - worker.cutRhs) > 0.1
        #  print(" violated ? ", violated)

        iter += 1
        if iter % 5 == 0:
            check = [i for i in inout]
            #  print("checking these constraints ", check)

            slacks = cpx.solution.get_linear_slacks(check)
            #  nrConstr = cpx.linear_constraints.get_num()
            remove = [check[i] for i in range(len(slacks)) if slacks[i] > _EPSI]
            #  print("removing ", remove)
            cpx.linear_constraints.delete(remove)
            for i in remove:
                inout.remove(i)
            #  print("Removed ", remove)
            #  print("This is inout = ", inout)

        if violated:
            #  progr = cpx.linear_constraints.get_num()
            constrName = "inout." + str(globalProgr)
            #  print("adding ", constrName)
            inout.append(constrName)
            globalProgr += 1
            cpx.linear_constraints.add(lin_expr = [worker.cutLhs],
                                       senses   = "L",
                                       rhs      = [worker.cutRhs],
                                       names    = [constrName])
        else:
            yIn = ySep
            print("... not violated ... moving 'in' point")

    cpx.solve()

    check = [i for i in inout]
    slacks = cpx.solution.get_linear_slacks(check)
    #  print("[", nrConstr,"] SLACKS  == ", slacks)
    remove = [check[i] for i in range(len(slacks)) if slacks[i] > _EPSI]
    cpx.linear_constraints.delete(remove)
    for i in remove:
        inout.remove(i)

    cpx.solve()
    #  print("OBJ AFTER REMOVING SLACKS = ", cpx.solution.get_objective_value())
    nrConstr = cpx.linear_constraints.get_num()
    print(" ... now constr is ", nrConstr)
    return globalProgr

def fixingToZero(inp, cpx, bestLB, bestUB, y_ilo, yFixed):
    cpxClone = cplex.Cplex(cpx)
    cpxClone.set_problem_type(cpxClone.problem_type.LP)
    cpxClone.solve()
    for j in range(inp.nI):
        rc = cpxClone.solution.get_reduced_costs(y_ilo[j])
        for t in range(inp.nP):
            if cpxClone.solution.get_values(y_ilo[j][t]) <= _EPSI and \
                yFixed[j][t] == 0 and (bestLB + rc[t]) > bestUB:
                yFixed[j][t] = 1

                print(" +** ** ** ** fixing to zero ", y_ilo[j][t])
                cpx.variables.set_upper_bounds(y_ilo[j][t], 0.0)

def addCorridor(inp, cpx, ySol, y_ilo, cWidth):
    """
    .. func:addCorridor()

    This function implements the corridor method for Benders' decomposition.
    The idea can be described as follows: Every time a new incumbent solution
    is found by the subproblem (i.e., a feasible solution which improves the
    upper bound), a new best solution is obtained. Let us indicate with :math:`y^I` 
    such solution. After adding the corresponding optimality cut to the master,
    we also impose the following corridor constraint:

    .. math ::
        :nowrap:
        
        \\begin{equation}
        \sum_{t=1}^T y^I_{jt}\left(1-y_{jt}\\right) + \left(1-y^I_{jt}\\right)y_{jt} \leq
        \gamma nT
        \end{equation}

    which ensure that no more than :math:`\gamma` percent of the setup
    variables will change value. This approach imposes a maximum hamming
    distance with respect to the incumbent solution, thus
    restricting the solution space of the master problem. This allows to speed up
    the convergence of Benders' algorithm. 
 
    """
    index = [y_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP)]
    value = [1-2*ySol[j][t] for j in range(inp.nI) for t in range(inp.nP)]
    rhsVal = cWidth - sum([ySol[j][t] for j in range(inp.nI) for t in
    range(inp.nP)])
    corridor_constr = cplex.SparsePair(ind=index,val=value)
    cpx.linear_constraints.add(lin_expr = [corridor_constr],
                               senses   = "L",
                               rhs      = [rhsVal],
                               names    = ["corridor"])
    print("... [",cWidth,"] After adding corridor constr : ", cpx.linear_constraints.get_num())

def setCpxParameters(cpx):
    cpx.parameters.preprocessing.presolve.set(cpx.parameters.preprocessing.presolve.values.off)
    cpx.parameters.threads.set(1)
    cpx.parameters.mip.strategy.search.set(cpx.parameters.mip.strategy.search.values.traditional)
    cpx.parameters.timelimit.set(30)

def getSolution(inp, cpx, y_ilo, z_ilo):
    ySol = []
    zHat = cpx.solution.get_values(z_ilo)
    for j in range(inp.nI):
        ySol.append(cpx.solution.get_values(y_ilo[j]))
    return zHat, ySol

def getUB(inp, zSub, ySol):
    z = zSub
    for j in range(inp.nI):
        for t in range(inp.nP):
            z += inp.f[j][t] * ySol[j][t]
    return z


def benderAlgorithm(inp, fixToZero, fixToOne, cPercent, cZero, cOne):
    """
    .. func:benderAlgorithm()

    This function implements Benders' scheme without the use of callbacks. We
    thus define a cycle that iteratively solves the master and the subproblems,
    until a termination criterion is reached.
    
    The current implementation is heuristic in nature, since fixing schemes as
    well as a corridor are used. Therefore, there is no guarantee that the
    solution returned by this algorithm is optimal. To get an exact approach,
    deactivate the fixing schemes (to zero and to one) and the corridor scheme.

    We make use of an in-out cycle, in a fashion similar to what has been done
    by Fischetti for the uncapacitated facility location problem (see their MS
    paper.) Probably, this part could be improved, to make it faster. However,
    it seems important to reach a good lower bound quickly. If the in-out cycle
    is not used, the lower bound is much worse and the convergence is extremely
    slow.

    Another feature of the Benders' implementation is connected with the use of
    a pool of solutions. Before calling Benders' algorithm, we heuristically
    solve the original CLSP with a maximum number of solutions to be reached.
    During this call to cplex, we collect a *pool* of solutions, which can be
    used to tighten the master. The idea is that we then pass these solutions
    to the separation scheme. Alternative approaches have been explored, e.g.,
    producing feasible solutions using a Cross Entropy scheme and, then, apply
    the separation mechanism to each of these solutions. To get the pool of
    solution, it suffices to call :meth:`MIP.solve()` with the flag ``withPool=1``.

    Here is a high-level description of the algorithm:

    1. create master and worker
    2. in-out cycle to quickly get a good lower bound
    3. fixing schemes (optionals) to zero and one
    4. add  corridor w.r.t. incumbent solution
    5. solve current master and get master solution yHat
    6. solve worker subproblem induced by yHat
    7. add cut induced by subproblem solution to current master
    8. if new incumbent is obtained, update bounds and corridor
    9. if gap < :math:`\epsilon`, STOP. Else, go to s5.

    """

    globalProgr = 0

    #  these are parameters for the corridor and the fixing schemes
    nWidth = max(cPercent*inp.nI*inp.nP, 1)
    nZero  = cZero*len(fixToZero)
    nOne   = cOne*len(fixToOne)

    cpx = cplex.Cplex()
    #  cpx.set_results_stream(None)
    #  cpx.set_log_stream(None)
    createMaster(inp, cpx)
    worker = WorkerLPReformulation(inp)
    setCpxParameters(cpx)

    globalProgr = inOutCycle2(cpx, worker, y_ilo, z_ilo, inp, globalProgr)
    #  cpx.write("inout-12-15.lp")
    #  cpx.read("inout-12-15.lp")
    #  cpx.read("inout-6-15.lp")
    cpx.set_problem_type(cpx.problem_type.MILP)
    for j in range(inp.nI):
        for t in range(inp.nP):
            cpx.variables.set_types(y_ilo[j][t], cpx.variables.type.binary)

    if len(fixToZero) > 0:
        #  add some fixing scheme here
        print("FIX TO ZERO HERE = ", fixToZero)
        zero_cut = cplex.SparsePair(ind=fixToZero, val=[1.0]*len(fixToZero))
        cpx.linear_constraints.add(lin_expr = [zero_cut],
                                   senses   = ["L"],
                                   rhs      = [nZero])
    if len(fixToOne) > 0:
        print("FIX TO ONE HERE = ", fixToOne)
        one_cut = cplex.SparsePair(ind=fixToOne, val=[1.0]*len(fixToOne))
        cpx.linear_constraints.add(lin_expr = [one_cut],
                                   senses   = ["G"],
                                   rhs      = [nOne])

    if nPool > 0:
        print("Before adding cut to master : ", cpx.linear_constraints.get_num())
        print("Adding solutions from pool of size ", nPool)
        for i in range(nPool):
        #  for i in range(0):
            cutType = worker.separate(inp, yPool[i], 0.0, y_ilo, z_ilo)
            if cutType > 0:
                constrName = "heur." + str(i)
                cpx.linear_constraints.add(lin_expr = [worker.cutLhs],
                                           senses   = "L",
                                           rhs      = [worker.cutRhs],
                                           names    = [constrName])
        print("After adding Pool cut to master : ", cpx.linear_constraints.get_num())

    #  #  add corridor constraint
    if nPool > 0:
        addCorridor(inp, cpx, yPool[0], y_ilo, nWidth)

    #  initialize data structure
    yFixed = [ [0 for i in range(inp.nP)] for t in range(inp.nI)]
    ySol   = []
    zHat   = -1
    for j in range(inp.nI):
        ySol.append([])

    nIter        = 0
    stopping     = 0
    nrCuts       = 0
    bestLB       = 0.1
    maxTolerance = 0.01
    bestUB       = cplex.infinity
    gap          = cplex.infinity
    # Benders' main cycle
    while not stopping:
        if nIter % 5 == 0:
            globalProgr = inOutCycle(cpx, worker, y_ilo, z_ilo, inp, globalProgr)
            cpx.set_problem_type(cpx.problem_type.MILP)
            for j in range(inp.nI):
                for t in range(inp.nP):
                    cpx.variables.set_types(y_ilo[j][t], cpx.variables.type.binary)

        cpx.solve()
        bestLB = cpx.solution.get_objective_value()
        with open(lbSummary,"a") as ff:
            ff.write("{0:20.5f} {1:20.5f} \n".format(bestLB, time()-startTime))
        
        if (bestUB - bestLB) <= maxTolerance:
            stopping = 1
            continue

        #  get master solution
        zHat, ySol = getSolution(inp, cpx, y_ilo, z_ilo)

        #  solve subproblem and add separation cut
        cutType = worker.separate(inp, ySol, zHat, y_ilo, z_ilo)
        if cutType > 0:
            cutName = "cut." + str(nrCuts)
            nrCuts += 1
            cpx.linear_constraints.add(lin_expr = [worker.cutLhs],
                                       senses   = ["L"],
                                       rhs      = [worker.cutRhs],
                                       names    = [cutName])
        if cutType == 2: # optimality cut
            ub = getUB(inp, worker.zSub, ySol)
            if ub < bestUB or ((nIter % 5 ) == 0 and gap < 0.1):
                if ub < bestUB:
                    bestUB = ub
                    #  add corridor constraint
                    #  addCorridor(inp, cpx, yPool[0], y_ilo, nWidth)
                    if cPercent > 0.0:
                        addCorridor(inp, cpx, ySol, y_ilo, nWidth)

                fixingToZero(inp, cpx, bestLB, bestUB, y_ilo, yFixed)

            if (bestUB - bestLB) <= maxTolerance:
                stopping = 1
                continue

        nIter +=1
        if bestLB > 0.0:
            gap = (bestUB-bestLB)/bestLB
        print("Iter {0:5d} :: GAP = {1:8.3f} ... LB = {2:8.3f} vs UB =\
        {3:8.3f}. ### Added {4:5d} cuts.".format(nIter, gap, bestLB, bestUB, \
        nrCuts))

    print("Final Result ====== ")
    print(" ** LB = ", bestLB)
    print(" ** UB = ", bestUB)


def printParameters():

    if algo == 1:
        algoType = "** Benders"
    elif algo == 2:
        algoType = "** Lagrange"
    elif algo == 3:
        algoType = "** Dantzig-Wolfe"
    elif algo == 4:
        algoType = "** CPLEX "
    else:
        algoType = "--> ERROR in algorithm selection !"

    print("==================================================================")
    print(" Marco Caserta 2017 - Benders, Lagrange, Dantzig-Wolfe ")
    print("==================================================================")
    print("\n")
    print("  Instance file \t :: {0:30s}".format(inputfile))
    print("  Algorithm     \t :: {0:30s}".format(algoType))
    print("  Parameters    \t :: ")
    print("   - Corridor Width \t = {0:20.5f}".format(cPercent))
    print("   - Soft Fix Zero  \t = {0:20.5f}".format(cZero))
    print("   - Soft Fix One   \t = {0:20.5f}".format(cOne))
    print("==================================================================")
    print("\n")



def main(argv):
    """
    Entry point.

    We first parse the command line, then reading the instance from a disk
    file.

    The algorithm to be used is selected via command line using flag -a. See
    :func:`parseCommandLine()` for more details on flags and options for this
    code.

    With respect to the type of algorithms that can be used, we have:

        1.  Benders Decomposition
        2.  Lagrangean Relaxation
        3.  Dantzig-Wolfe
        4.  Cplex MIP solver

    """
    global fixToZero
    global fixToOne
    global startTime
        #  outfile.write("{0:20s} {1:20.5f} {2:25s} {3:20.5f} {4:20.7f} {5:20.7f}\n".
        #                format(inputfile, zOpt, stat, lb, gap, zTime))

    parseCommandLine(argv)
    inp = Instance(inputfile)
    startTime = time()
    printParameters()
    if algo == 1:
        mip       = MIPReformulation(inp)
        #  fixToZero = mip.solveLPZero(inp)
        #  fixToOne  = mip.solveLPOne(inp)
        #  zHeur     = mip.solve(inp,nSol    = 100, display = 0, withPool = 1)
        #  print("zHeur = ", zHeur)
        benderAlgorithm(inp, fixToZero, fixToOne, cPercent, cZero, cOne)
        exit(101)
    if algo == 2:
        mip       = MIPReformulation(inp)
        fixToZero = mip.solveLPZero(inp)
        fixToOne  = mip.solveLPOne(inp)
        lagrange = Lagrange(inp, 99999999)
        lagrange.lagrangeanPhase(inp, mip, fixToZero, fixToOne, cPercent, cZero,
        cOne)
        exit(102)
    if algo == 3:
        dw = DantzigWolfe(inp)
        dw.dw_cycle(inp, lbSummary, startTime)
        exit(103)
    if algo == 4: #  Cplex MIP solver
        #  mip       = MIPReformulation(inp)
        mip       = MIP(inp)
        mip.solve(inp, withPrinting=1)
        exit(104)

    print("Algorithm Type not defined. Choose between 1 and 4.")
    exit(105)

    #  ======================================================================
    #  All this stuff below is used to define cplex callback. This was the first
    #  version of Benders implementation. It is no longer needed.
    #  ======================================================================
    cpx = cplex.Cplex()
    cpxClone = cplex.Cplex()
    # create master and worker (subproblem)
    createMaster(inp, cpx)
    createMaster(inp, cpxClone)
    cpxClone.set_problem_type(cpx.problem_type.LP)
    #  worker = WorkerLP(inp)
    worker = WorkerLPReformulation(inp)

    # Set up cplex parameters to use the cut callback for separating
    # Benders' cuts
    cpx.parameters.preprocessing.presolve.set(
        cpx.parameters.preprocessing.presolve.values.off)
    # Set the maximum number of threads to 1.
    cpx.parameters.threads.set(1)
    # Turn on traditional search for use with control callbacks
    cpx.parameters.mip.strategy.search.set(
        cpx.parameters.mip.strategy.search.values.traditional)

    #  inOutCycle(cpx, worker, y_ilo, z_ilo, inp)
    #  cpx.write("inout-6-15.lp")

    cpx.read("inout-6-15.lp")
    cpx.set_problem_type(cpx.problem_type.MILP)
    print("Type here = ", cpx.problem_type[cpx.get_problem_type()])

    #  binary variables must be re-specified
    for j in range(inp.nI):
        for t in range(inp.nP):
            cpx.variables.set_types(y_ilo[j][t], cpx.variables.type.binary)

    print("Before adding cut to master : ", cpx.linear_constraints.get_num())
    cutType = worker.separate(inp, yRef, 0.0, y_ilo, z_ilo)
    if cutType > 0:
        constrName = "heur." + str(0)
        cpx.linear_constraints.add(lin_expr = [worker.cutLhs],
                                   senses   = "L",
                                   rhs      = [worker.cutRhs],
                                   inames    = [constrName])
    print("Cut added to master : ", cpx.linear_constraints.get_num())
    #  solveCall = cpx.register_callback(SolveNodeCallback)

    # register LAZY callback
    lazyBenders        = cpx.register_callback(BendersLazyConsCallback)
    lazyBenders.cpx    = cpx
    lazyBenders.inp    = inp
    lazyBenders.z_ilo  = z_ilo
    lazyBenders.y_ilo  = y_ilo
    lazyBenders.yFixed = [ [0 for i in range(inp.nP)] for t in range(inp.nI)]
    lazyBenders.worker = worker
    lazyBenders.solved = 0
    lazyBenders.rc     = []
    lazyBenders.nIter  = 0
    if userCuts == "1":
        # register USER callback
        userBenders        = cpx.register_callback(BendersUserCutCallback)
        userBenders.inp    = inp
        userBenders.z_ilo  = z_ilo
        userBenders.y_ilo  = y_ilo
        userBenders.worker = worker

    startTime = time.time()
    # Solve the model
    cpx.solve()

    solution = cpx.solution
    print()
    print("Solution status: ", solution.status[solution.get_status()])
    print("Objective value: ", solution.get_objective_value())

    print("Thus time is ", time.time() - startTime)


if __name__ == "__main__":
    main(sys.argv[1:])

