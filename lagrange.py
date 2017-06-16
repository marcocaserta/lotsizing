
#!/usr/bin/python
"""
===========================================================================

:Filename: lagrange.py
:Author: marco caserta
:Date: 09.03.2017
:Last Update: |date|

.. |date| date:: %d.%m.%y
.. |time| date:: %H:%M

Copyright (C) 2017 by Marco Caserta  (marco dot caserta at ie dot edu)

(This document was generated on |date| at |time|.)

.. this is just a comment (it won't appear in the documentation)

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

This code implements a Lagrangean relaxation scheme for the multi-item
multi-period capacitated lot sizing problem.

"""

import sys, getopt
import math
import cplex
from time import time
import numpy as np
import random
random.seed(27)

_INFTY = sys.float_info.max
_EPSI  = sys.float_info.epsilon
maxIter = 500000

from clspBenders import lbSummary, startTime

class Lagrange:
    """
    This class implements the Lagrangean scheme for the CLSP. We relax in a
    Lagrangean fashion the capacity constraints. Therefore, the relaxed problem
    can be separated over the items.

    The same formulation used for the other techniques, i.e., Benders and
    Dantzig-Wolfe, is used here: The Single Plant Location formulation.
    
    Initialization of the Lagrangean class. We define a cplex object, which
    contains the **relaxed** problem, along with the decision variables
    :math:`y_{jt}` and :math:`z_{jtr}`.

    The formulation of the relaxed problem is given here (**note**: The problem
    can now be separated over the items):

    .. math::
        :nowrap:

        \\begin{eqnarray}
          & \min L(\\textbf{u}) = &   \sum_{j=1}^n \sum_{t=1}^T f'_{jt}y_{jt} +
          \sum_{j=1}^n \sum_{r=1}^T \sum_{t=1}^{r-1} h'_{jtr}z_{jtr} -\sum_{t=1}^T u_t
          \label{eq:LR-obj}\\\\
          &\mbox{s.t}&  \sum_{t=1}^r z_{jtr} = d_{jr},
          \quad j = 1, \ldots, n , \quad r = 1, \dots, T
          \label{eq:LR-demand-constr} \\\\
          && z_{jtr} \leq d_{jr}y_{jt},
          \quad j = 1, \ldots, n , \quad t = 1, \dots, T, \quad r=t,\dots,T
          \label{eq:LR-logic-constr} \\\\
          && \sum_{r=t}^T z_{jtr} \leq M y_{jt},
          \quad j = 1, \ldots, n , \quad t = 1, \dots, T
          \label{eq:LR-cumLogic-constr} \\\\
          && y_{jt} \in \left\{0,1\\right\}, \quad j = 1, \ldots, n, \quad
          t=1, \ldots, T \label{eq:LR-y-binary}\\\\
          && z_{jtr} \geq 0, \quad j = 1, \ldots, n, \quad
          t=1, \ldots, T, \quad r = t,\dots,T \label{eq:LR-z-cont}
        \end{eqnarray}
    """
    def __init__(self, inp, ub):
        z_ilo  = []
        y_ilo  = []

        cpx = cplex.Cplex()
        cpx.objective.set_sense(cpx.objective.sense.minimize)

        cpx.set_results_stream(None)
        cpx.set_log_stream(None)

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

        #  initialize lambda multipliers
        self.cpx    = cpx
        self.zL     = 0.0
        self.y_ilo  = y_ilo
        self.z_ilo  = z_ilo
        self.lmult  = [0.0]*inp.nP
        self.hVal0  = [0.0]*inp.nP
        self.hVal1  = [0.0]*inp.nP
        self.hVal2  = [0.0]*inp.nP
        self.iterL  = 0
        self.bestLB = 0.0
        self.normH  = _INFTY
        self.lbStar = _EPSI
        self.ubStar = ub
        self.delta  = 0.0


    def stoppingCriteria(self, maxIter, converged):
        """
        Stop when:

        * a maximum number of iterations is reached, or
        * the lagrangean bound is no longer improving (converged was reached)
        """
        iterL = self.iterL
        if iterL > maxIter or converged:
            return True
        else:
            return False
            

    def multRandomPerturbation(self, inp):
        """
        We run three Lagrangean cycles. At the end of each cycle, we randomly
        perturb the multipliers and restart. The idea is to introduce some
        "shaking" into the multipliers, without starting completely from
        scratch.
        """
        lmult = self.lmult

        for t in range(inp.nP):
            rr = random.uniform(-0.1,0.1)
            lmult[t] *= rr

        self.lmult = lmult


    def lagrangeanStep(self, inp):
        """
        Each Lagrangean step receives as input the current multipliers. Based
        on the value of :math:`\\mathbf{u}`, i.e., ``lmult``, we recompute the
        objective function coefficients as:

        .. math ::
            :nowrap:

            \\begin{eqnarray}
            f'_{jt} &=& f_{jt} + u_t m_{jt} \\\\
            h'_{jrt} &=& h_{jtr} + u_t a_{jt}
            \end{eqnarray}

        Next, we solve the current Lagrangean problem. Note that the Lagrangean
        problem obtained after relaxing the *capacity* constraints, could be
        separated over the item and solved using efficient methods, e.g.,
        Wagner-Within. 

        We indicate with ``zL``, ``yL``, and ``zL`` the objective function of the
        Lagrangean solution, the setup schedule, and the production plan,
        respectively. 
        
        Note that each Lagrangean solution provides a lower bound. Therefore,
        if the newly found value ``zL`` is better (higher) than the best known
        lower bound so far, we update the lower bound.
        """
        global startTime
        global lbSummary

        lmult  = self.lmult
        cpx    = self.cpx
        iterL  = self.iterL
        normH  = self.normH
        lbStar = self.lbStar
        ubStar = self.ubStar
        y_ilo  = self.y_ilo
        z_ilo  = self.z_ilo

        #  recompute obj function coefficients with the current multipliers
        for j in range(inp.nI):
            for t in range(inp.nP):
                coeff = inp.f[j][t] + inp.m[j][t]*lmult[t]
                cpx.objective.set_linear(y_ilo[j][t], coeff)
                for r in range(t,inp.nP):
                    coeff = (r-t)*inp.h[j][t] + lmult[t]*inp.a[j][t]
                    cpx.objective.set_linear(z_ilo[j][t][r], coeff)

        #  solve current Lagrangean problem
        cpx.solve()

        #  get Lagrangean solution
        zL = cpx.solution.get_objective_value()
        zL -= sum([inp.cap[t]*lmult[t] for t in range(inp.nP)])
        if zL > lbStar:
            improvement = (zL - lbStar)/lbStar
            lbStar = zL
            with open(lbSummary,"a") as ff:
                ff.write("{0:20.5f} {1:20.5} \n".format(lbStar,
                time()-startTime))
        else:
            improvement = 0.0

        if zL > self.lb50:
            self.lb50 = zL

        yL = [[cpx.solution.get_values(y_ilo[j][t]) for t in range(inp.nP)] \
              for j in range(inp.nI)]
        zVar = []
        for j in range(inp.nI):
            zVar.append([])
            for t in range(inp.nP):
                zVar[j].append([])
                for r in range(inp.nP):
                    zVar[j][t].append(cpx.solution.get_values(z_ilo[j][t][r]))

        if (iterL % 10 ) == 0:
            print("lb ({0:3d}) = {1:10.3f} (vs lb* = {2:10.3f} ub* =\
            {3:10.3f}) \t Status = {4:>10s} \t norm = {5:10.2f}".\
            format(iterL, zL, lbStar, ubStar, \
            cpx.solution.status[cpx.solution.get_status()], normH))


        self.zL          = zL
        self.yL          = yL
        self.zVar        = zVar
        self.lbStar      = lbStar
        self.improvement = improvement

    def lmultUpdate(self, inp):
        """
        Update Lagrangean multipliers. We use a subgradient deflection
        technique to update the value of the multipliers (see paper for
        details). A subgradient is:

        .. math::
            :nowrap:

              \\begin{equation}
              s_t (\\textbf{u}) = \sum_{j=1}^n \left( a_{jt}\sum_{r=t}^T z_{jtr}^L +
              m_{jt}y_{jt}^L\\right) - b_t,  \quad t = 1, \ldots, T
              \end{equation}

        where the subgradient :math:`s_t (\\textbf{u})` is a function of
        :math:`\\textbf{u}` via both :math:`\mathbf{y}^L` and
        :math:`\mathbf{z}^L`.

        """
        delta    = self.delta
        zL       = self.zL
        yL       = self.yL
        zVar     = self.zVar
        hVal0    = self.hVal0
        hVal1    = self.hVal1
        hVal2    = self.hVal2
        lmult    = self.lmult
        ub       = self.ubStar

        #  use subgradient deflection technique
        normSub = 0.0
        sub     = [0.0]*inp.nP
        for t in range(inp.nP):
            #  sub[t] = -inp.cap[t]
            #  sub[t] += sum([inp.m[j][t]*yL[j][t] for j in range(inp.nI)])
            #  sub[t] += sum([inp.a[j][t]*zVar[j][t][r] for j in range(inp.nI) \
            #  for t in range(inp.nP) for r in range(t,inp.nP)])
            sub[t] = -inp.cap[t]
            for j in range(inp.nI):
                sub[t] += inp.m[j][t]*yL[j][t]
                for r in range(t,inp.nP):
                    sub[t] += inp.a[j][t]*zVar[j][t][r]

            hVal0[t] = (sub[t] + 0.3*hVal1[t] + 0.1*hVal2[t])/1.4
            normSub += math.pow(hVal0[t],2)

        normSub = math.sqrt(normSub)

        if normSub < _EPSI:
            print("Norm of Subgradient is zero")
            exit(123)

        step = delta*(ub - zL)/normSub
        for t in range(inp.nP):
            lmult[t] = max(lmult[t] + step*hVal0[t], 0.0)

        hVal2 = hVal1[:]
        hVal1 = hVal0[:]
        delta = delta/1.01

        self.lmult    = lmult
        self.delta    = delta
        self.hVal0    = hVal0
        self.hVal1    = hVal1
        self.hVal2    = hVal2
        self.normH    = normSub


    def refineMIPSolution(self, inp, mip, cPercent):
        """
        Primal scheme based on the corridor method and the Lagrangean
        relaxation. The original MIP formulation is provided as input here
        (object ``mip``). We then build a corridor around the incumbent
        Lagrangean solution ``yL`` obtained within
        :meth:`lagrangeanStep()`. Once the corridor constraint is
        added to the original formulation, we solve the constrained CLSP.

        A solution found by this primal scheme provides a valid *upper bound*
        to the optimal solution. Therefore, if a new best value is found, we
        update the bound. 

        **Note**: The corridor constraint is removed at the end of the primal
        scheme, to avoid *overconstraining* the original formulation.
        """
        ubStar = self.ubStar
        yL    = self.yL
        y_ilo = mip.y_ilo
        cpx   = mip.cpx

        nrOnes = np.sum([yL[j][t] > 1.0-_EPSI for j in range(inp.nI) for t in
        range(inp.nP)])
        print("CORRIDOR = ", cPercent*inp.nP*inp.nI)
        rhs = cPercent*inp.nP*inp.nI - nrOnes
        index = [y_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP)]
        value = [(1.0 - 2.0*yL[j][t]) for j in range(inp.nI) for t in range(inp.nP)]
        corridor = cplex.SparsePair(ind=index, val=value)
        cpx.linear_constraints.add(lin_expr = [corridor],
                                   senses   = ["L"],
                                   rhs      = [rhs],
                                   names    = ["corridor"])
        mip.solve(inp, nSol=3, timeLimit=10)
        zCM = mip.cpx.solution.get_objective_value()
        if zCM < ubStar:
            ubStar = zCM
            self.bestTime = time() - self.initTime
            print(" *** CM({0:3d}) = {1:10.2f} after {2:7.3f} seconds.".\
                    format(self.iterL, ubStar, self.bestTime))


        cpx.linear_constraints.delete("corridor")

        self.ubStar = ubStar

    def FixVarsLP(self, mip, fixToZero, fixToOne, cZero, cOne):
        """
        Fixing scheme. Some variables are soft-fixed to 0 or 1, based on the LP
        scheme described in :meth:`clspBenders.MIP.solveLPZero()`. The two vectors
        ``fixToOne`` and ``fixToZero`` are obtained solving twice the LP
        relaxation of the original problem, with the goal of finding a subset
        of variables whose values can "safely" be set to 0 or 1.
        """
        cpx = mip.cpx
        y_ilo = mip.y_ilo

        rhsVal = cZero*len(fixToZero)
        print("Rhs ZERO [{0}/{1}] = {2}".format(rhsVal, len(fixToZero),
        fixToZero))
        zero_cut = cplex.SparsePair(ind=fixToZero, val=[1.0]*len(fixToZero))
        cpx.linear_constraints.add(lin_expr = [zero_cut],
                                   senses   = ["L"],
                                   rhs      = [rhsVal])

        #  rhsVal = (1-0.1)*len(fixToOne)
        rhsVal = cOne*len(fixToOne)
        print("Rhs ONE  [{0}/{1}] = {2}".format(rhsVal, len(fixToOne),
        fixToOne))
        one_cut = cplex.SparsePair(ind=fixToOne, val=[1.0]*len(fixToOne))
        cpx.linear_constraints.add(lin_expr = [one_cut],
                                   senses   = ["G"],
                                   rhs      = [rhsVal])
    def checkConvergence(self):
        """
        One of the two stopping criteria. If the improvement in the lower bound
        in the last 50 iterations is below a certain value, we assume the
        Lagrangean phase has converged and, therefore, we stop (or restart).
        """
        if (self.lb50 - self.lbInit)/self.lbInit < 0.0001:
        #  if (self.lb50 - self.lbInit)/self.lbInit < _EPSI:
            return True
        else:
            self.lb50   = self.zL
            self.lbInit = self.zL
            return False


    def lagrangeanPhase(self, inp, mip, fixToZero, fixToOne, cPercent, cZero,
    cOne):
        """
        We employ a primal-dual scheme, in which the Lagrangean relaxation
        provides a sequence of (non-monotonically) increasing lower bound,
        while a corridor-based primal scheme gives upper bounds. The basic
        steps of the algorithm are as folllows:

        1. soft-fix some of the variables to 0 or 1 (LP-based scheme)
        2. initialize Lagrangean multipliers
        3. perform one Lagrangean step, i.e., solve the relaxed problems for the given set of multipliers :math:`\\rightarrow` we obtain a Lagrangean lower bound ``zL`` and a solution ``yL``
        4. update Lagrangean multipliers (subgradient optimization)
        5. periodically, apply primal scheme (add corridor around ``yL`` and solve constrained CLSP)
        6. if stoppingCriteria are reached, either restart or stop (the outer cycle is repeated three times, after which the algorithm stops.)

        """

        #  soft-fix some variables to 0 or 1
        self.FixVarsLP(mip, fixToZero, fixToOne, cZero, cOne)
        
        print("*** LAGRANGEAN PHASE STARTS ***")
        self.initTime = time()
        mainCycle     = 0
        lbStar        = 0.0

        while mainCycle < 3:
            self.iterL  = 1
            self.lbStar = _EPSI
            normH       = 0.0
            hVal1       = [0.0]*inp.nP
            hVal2       = [0.0]*inp.nP
            self.delta  = 0.00000001
            #  self.delta  = 0.0001
            converged   = False
            self.lb50   = _EPSI
            self.lbInit = _EPSI

            print(" .. Lagrangean Cycle {0:4d} ..".format(mainCycle+1))
            while self.stoppingCriteria(maxIter, converged) is not True:
                self.lagrangeanStep(inp)
                self.lmultUpdate(inp)

                #  if self.iterL > 10 and self.improvement > 0.001:
                if (self.iterL % 10) == 0 :
                    if cPercent > 0.0:
                        self.refineMIPSolution(inp, mip, cPercent)
                self.iterL += 1

                if (self.iterL % 50) == 0:
                    converged = self.checkConvergence()

            mainCycle += 1
            if self.lbStar > self.bestLB:
                self.bestLB = self.lbStar
            if mainCycle < 3:
                self.multRandomPerturbation(inp)

            print(" .. End Lagrangean Cycle {0:4d} .. Best LB = {1:10.2f}".\
                    format(mainCycle, self.bestLB))

        print("Best UB Lagrange = {0:20.5f} in {1:20.5f}".format(self.ubStar, self.bestTime))
        #  write solution to disk
        with open("solution.txt", "w") as solFile:
            solFile.write("L   :: {0:20.5f} {1:20.5f}\
            {2:20.5f} {3:10.5f} {4:10.5f} {5:10.5f}\n".format(self.bestLB,
            self.ubStar, self.bestTime, cPercent, cZero, cOne))
        

    

    

    
    
