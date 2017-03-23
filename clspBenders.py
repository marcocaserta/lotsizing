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

.. note:: The dual values have positive sign here, since they are treated as Lagrangean multipliers (this is the reason why the sign is reversed.)

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


"""

from __future__ import print_function

import sys, getopt 
import math 
import time 
import cplex 
from cplex.callbacks import UserCutCallback, LazyConstraintCallback 
from cplex.exceptions import CplexError

_INFTY = sys.float_info.max 
_EPSI  = sys.float_info.epsilon

inputfile = ""
usercuts  = ""

#  we set the master variables as global, while the subproblem vars are local
z_ilo     = -1
y_ilo     = []
#  xc_ilo    = []
lCapacity = []
lLogic    = []
lDemand   = []

startTime = -1;

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
        worker = self.worker
        y_ilo  = self.y_ilo
        z_ilo  = self.z_ilo
        #  xc_ilo = self.xc_ilo
        inp    = self.inp


        #  get current master solution
        zHat = self.get_values(z_ilo)
        ySol = []
        for j in range(inp.nI):
            ySol.append([])
            ySol[j] = self.get_values(y_ilo[j])
        #  xcSol = self.get_values(xc_ilo)


        #  flatten =  [item for sublist in ySol for item in sublist]
    
        #  benders cut separation
        #  cutType = worker.separate(inp, ySol, zHat, y_ilo, z_ilo, xcSol, xc_ilo)
        cutType = worker.separate(inp, ySol, zHat, y_ilo, z_ilo)
        if cutType > 0:
            #  a = [float(worker.cutLhs.val[i]) for i in range(inp.nI*inp.nP)]
            #  lhsSum = sum([a[i]*flatten[i] for i in range(inp.nI*inp.nP)])
            # add Benders cut to the master
            self.add(constraint = worker.cutLhs,
                     sense     = "L",
                     rhs        = worker.cutRhs)


class BendersUserCutCallback(UserCutCallback):

    def __call__(self):
        """
        Define the actions to be carried out at every callback.

        Note that the ``separate`` function of the subproblem is called here.

        """
        # Skip the separation if not at the end of the cut loop
        if not self.is_after_cut_loop():
            return

        # get data structure (self is the master)
        worker = self.worker
        y_ilo  = self.y_ilo
        z_ilo  = self.z_ilo
        #  xc_ilo = self.xc_ilo
        inp    = self.inp


        #  get current master solution
        zHat = self.get_values(z_ilo)
        ySol = []
        for j in range(inp.nI):
            ySol.append([])
            ySol[j] = self.get_values(y_ilo[j])
        #  xcSol = self.get_values(xc_ilo)


        #  flatten =  [item for sublist in ySol for item in sublist]
    
        #  benders cut separation
        #  cutType = worker.separate(inp, ySol, zHat, y_ilo, z_ilo, xcSol, xc_ilo)
        cutType = worker.separate(inp, ySol, zHat, y_ilo, z_ilo)
        if cutType > 0:
            #  a = [float(worker.cutLhs.val[i]) for i in range(inp.nI*inp.nP)]
            #  lhsSum = sum([a[i]*flatten[i] for i in range(inp.nI*inp.nP)])
            # add Benders cut to the master
            self.add(cut       = worker.cutLhs,
                     sense     = "L",
                     rhs       = worker.cutRhs)




def parseCommandLine(argv):
    """
    .. func:parseCommandLine()

    Parse command line. Options are:

    -h help         usage help

    -i inputfile    instance file name

    -u userCuts    activate user cuts (fractional values)

    """
    global inputfile
    global userCuts
    
    try:
        opts, args = getopt.getopt(argv, "hi:u:", ["help","ifile=","ucuts"])
    except getopt.GetoptError:
        print("Command Line Error. Usage : python cflp.py -i <inputfile> -u\
        <usercuts>")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("Usage : python cflp.py -i <inputfile> -u <usercuts>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-u", "--ucuts"):
            userCuts = arg

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
       self.dj       = []

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
           
           #  aa = [math.floor( (self.cap[t] - self.m[j][t])/self.a[j][t]) for t in range(self.nP)]
           aa = [( (self.cap[t] - self.m[j][t])/self.a[j][t]) for t in range(self.nP)]
           print("AA vs CAP ", aa, " ", self.dcum[j])
           #  self.max_prod.append([min(aa[t],self.dcum[j][t]) for t in \
           #  range(self.nP)])
           self.max_prod.append([self.dcum[j][t] for t in \
           range(self.nP)])
           #  self.max_prod.append(aa)
           self.dj.append(sum([self.d[j][t] for t in range(self.nP)]))

           print(self.max_prod[j])


class MIP:
    """
    Define the full model and solve it using cplex. We use this class to
    compute optimal values of the full MIP models and compare them, along with
    the achieve performance, with the Benders approach.
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
        


    def solve(self, inp):
        """
        Solve MIMPLS using cplex branch and bound. 
        """
        
        cpx = self.cpx
        
        cpx.solve()

        print("STATUS = ", cpx.solution.status[cpx.solution.get_status()])
        if cpx.solution.get_status() == cpx.solution.status.optimal_tolerance:
            print("OPT SOL found ")
        

            
class WorkerLP:
    """
    Define and solve the subproblem. We initilize the subproblem with the right
    hand side values of the constraints to zero, since we assume the initial 
    values of :math:`y_{jt}` to be equal to zero. Next, within the
    ``separate()`` function, we define the rhs values to the correct values,
    depending on the solution obtained from the master.
    
    Cplex requires the presolve reductions to be turned off, using::

        cpx.parameters.preprocessing.reduce.set(0)

    In addition, we need to ensure that the LP is *not* solved using the
    interior point method, otherwise dual values won't be available. We can
    either use primal simplex or dual simplex.

    .. note :: 

        The subproblem constraints should be defined with a name, in
        order to be able to recall the precise name of each constraint when we want
        to obtain the dual values. Briefly, we need to:
    
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

        #  for i in range(20):
        #      constr = cpx.linear_constraints.get_rows(i)
        #      print(cpx.variables.get_names(constr.ind), constr.val, " =" ,
        #      cpx.linear_constraints.get_rhs(i))

            

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
        Here is were Benders cut is obtained and passed to the master.

        The following steps describe the algorithm:

        * update rhs values of the subproblem, i.e., using the current optimal
            solution of the master problem :math:`y^*`
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

        #  print("MASTER SOL :: ")
        #  for j in range(inp.nI):
        #      aux = [t for t in range(inp.nP) if ySol[j][t] >= (1.0-_EPSI)]
        #      print("item ", j, " :: ", aux)

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
            sumY = sum([inp.m[j][t]*ySol[j][t] for j in range(inp.nI) if ySol[j][t] >=
            (1.0-_EPSI)])
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

        return cutType
        


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
    #  global xc_ilo
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
    #      xc_ilo.append(cpx.variables.get_num())
    #      varName = "xc." + str(j)
    #      cpx.variables.add(obj   = [0.0],
    #                        lb    = [inp.dcum[j][0]],
    #                        ub    = [cplex.infinity],
    #                        types = ["C"],
    #                        names = [varName])

    # add cumulative constraints
    #  for j in range(inp.nI):
    #      index = [y_ilo[j][t] for t in range(inp.nP)]
    #      value = [inp.max_prod[j][t] for t in range(inp.nP)]
    #      cumulative_constraint = cplex.SparsePair(ind=index, val=value)
    #      constrName = "cumulative." + str(j)
    #      cpx.linear_constraints.add(lin_expr = [cumulative_constraint],
    #                                 senses   = ["G"],
    #                                 rhs      = [inp.dj[j]],
    #                                 names    = [constrName])
    #  for j in range(inp.nI):
    #      index = [y_ilo[j][t] for t in range(inp.nP)]
    #      value = [1.0]*inp.nP
    #      ones_constraint = cplex.SparsePair(ind=index, val=value)
    #      constrName = "sum." + str(j)
    #      cpx.linear_constraints.add(lin_expr = [ones_constraint],
    #                                 senses   = ["G"],
    #                                 rhs      = [2.0],
    #                                 names    = [constrName])



    
def main(argv):
    """
    Entry point. 

    We first parse the command line, then reading the instance from a disk
    file.

    The relevant part is the registration of the callback, obtained via::

         lazyBenders = cpx.register_callback(BendersLazyConsCallback)

    """

    parseCommandLine(argv)
    inp = Instance(inputfile)

    cpx = cplex.Cplex()

    # activate this part if we want to solve the original MIP via cplex
    #  mip = MIP(inp)
    #  mip.solve(inp)

    # create master and worker (subproblem)
    createMaster(inp, cpx)
    worker = WorkerLP(inp)

    # Set up cplex parameters to use the cut callback for separating
    # Benders' cuts
    cpx.parameters.preprocessing.presolve.set(
        cpx.parameters.preprocessing.presolve.values.off)

    # Set the maximum number of threads to 1.
    cpx.parameters.threads.set(1)

    # Turn on traditional search for use with control callbacks
    cpx.parameters.mip.strategy.search.set(
        cpx.parameters.mip.strategy.search.values.traditional)

    # register LAZY callback
    lazyBenders = cpx.register_callback(BendersLazyConsCallback)
    lazyBenders.inp   = inp
    lazyBenders.z_ilo = z_ilo
    lazyBenders.y_ilo = y_ilo
    #  lazyBenders.xc_ilo = xc_ilo
    lazyBenders.worker = worker
    if userCuts == "1":
        # register USER callback
        userBenders = cpx.register_callback(BendersUserCutCallback)
        userBenders.inp   = inp
        userBenders.z_ilo = z_ilo
        userBenders.y_ilo = y_ilo
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

