
#!/usr/bin/python
"""
===========================================================================

:Filename: bendersCallbacks.py
:Author: marco caserta
:Date: 05.09.2018
:Last Update: |date|

.. |date| date:: %d.%m.%y
.. |time| date:: %H:%M

Copyright (C) 2018 by Marco Caserta  (marco dot caserta at ie dot edu)

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
"""
import sys, getopt
import math
import time
import csv

import cplex
from cplex.exceptions import CplexError
from cplex.callbacks import IncumbentCallback
from cplex.callbacks import HeuristicCallback
from cplex.callbacks import UserCutCallback, LazyConstraintCallback

largeInstances = False

_INFTY    = sys.float_info.max
_EPSI     = sys.float_info.epsilon
inputfile = ""
lbSummary = "lowerBounds.txt"
ubSummary = "upperBounds.txt"

#  we set the master variables as global, while the subproblem vars are local
z_ilo     = -1
y_ilo     = []
u_ilo     = []
lCapacity = []
lLogic    = []
lDemand   = []
inout     = []

startTime = time.time()


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

        mip    = self.mip
        inp    = self.inp

        #  print("NR NODES = ", self.get_num_nodes())
        #  print("GAP = ", self.get_MIP_relative_gap())

        cpx    = self.cpx
        worker = self.worker
        y_ilo  = self.y_ilo
        z_ilo  = self.z_ilo
        bestLB = self.bestLB
        bestUB = self.bestUB
        bestCM = self.bestCM
        bestTime = self.bestTime


        self.nIter += 1
        #  get current master solution
        objMaster = self.get_objective_value()
        lb        = self.get_best_objective_value()
        if lb > bestLB:
            bestLB = lb
        zHat = self.get_values(z_ilo)
        ySol = [self.get_values(y_ilo[j]) for j in range(inp.nI)]
        #  for j in range(inp.nI):
        #      ySol.append([])
        #      ySol[j] = self.get_values(y_ilo[j])
        #  print("[{0:5d}] lb = {1:9.2f}; ub = {2:9.2f}".format(self.nIter,
        #  bestLB, bestCM))

        # add cut to the master
        # NOTE: cplex makes a working copy of the master (to which I am not able to
        # gain access.) The "cpx" object remains empty, without the extra cuts
        #  cutType, zSub = worker.solveSubDual(inp, ySol, zHat, y_ilo, z_ilo)
        #  cutType, zSub = worker.solveSubPrimal(inp, ySol, zHat, y_ilo, z_ilo)
        cutType, zSub = worker.solveSubQ(inp, ySol, zHat, y_ilo, z_ilo)

        # NOTE: Remove this part when activating CM
        zUB = (objMaster - zHat) + zSub
        if zUB < bestCM:
            bestCM = zUB
            print("** New UB found by Benders = ", bestCM)
            bestTime = time.time() - self.startTime



        self.add(constraint = worker.cutLhs,
                 sense      = "L",
                 rhs        = worker.cutRhs,
                 use        = self.use_constraint.purge)


        self.gap = self.get_MIP_relative_gap()
        if self.gap < 0.03:
            # call corridor method
            delta = int(inp.nI*inp.nP*0.05)
            delta = 20
            corridorLhs, corridorRhs = addCorridor(inp, delta, ySol, y_ilo)
            constrName = "corridor." + str(self.nCM)
            #  print("** Adding corridor Nr. ", constrName)
            mip.addConstraint(corridorLhs, corridorRhs, "L", constrName)
            self.nCM += 1
            #  mip.solve(inp, nSol=10)

            self.bestUB = self.get_best_objective_value()
            self.gap = self.get_MIP_relative_gap()

            MIPstatus = mip.solve(inp)
            print("MIP Status inside callback ", MIPstatus)
            #  if self.nIter >5:
            #      MIPstatus = -1
            if MIPstatus == -1:
                print("getting out of callbacks")
                self.abort()
                return

            mip.getSolution(inp, withPrinting=1)
            if mip.z < bestCM:
                bestCM = mip.z
                print("** New UB found by CM = ", bestCM)
                bestTime = time.time() - self.startTime
            #  mip.removeConstraint(constrName)
            mip.reverseConstraint(corridorLhs, corridorRhs, "G", constrName)

            # now add cut found by CM
            cutType, zSub = worker.solveSubQ(inp, mip.ySol, zHat, y_ilo, z_ilo)

            myUB = (objMaster - zHat) + zSub
            if myUB < bestUB:
                bestUB = myUB

            self.add(constraint = worker.cutLhs,
                     sense      = "L",
                     rhs        = worker.cutRhs,
                     use        = self.use_constraint.purge)

            #  print("Adding corridor", corridorLhs, corridorRhs)
            #  .add_local(constraint = corridorLhs,
            #           sense      = "L",
            #           rhs        = corridorRhs)
            #
            #  self.Lhs.append(worker.cutLhs)
            #  self.Rhs.append(worker.cutRhs)
            





        self.nIter += 1

        self.bestLB = bestLB
        self.bestUB = bestUB
        self.bestCM = bestCM
        self.bestTime = bestTime
        self.cpx    = cpx

def addCorridor(inp, delta, ySol, y_ilo):
    
    index = [y_ilo[j][t] for j in range(inp.nI) for t in range(inp.nP)]
    value = [1-2*ySol[j][t] for j in range(inp.nI) for t in range(inp.nP)]

    ySolSum = sum([1 for j in range(inp.nI) for t in range(inp.nP)
    if ySol[j][t] >= 1.0-_EPSI])

    corridorRhs = delta - ySolSum

    corridorLhs = cplex.SparsePair(ind=index,val=value)

    return corridorLhs, corridorRhs
    


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

    -p pool         number of solutions in pool (used in DW)

    -a algorithm   type of algorithm used:

    With respect to the type of algorithms that can be used, we have:

        1.  Benders Decomposition
        2.  Lagrangean Relaxation
        3.  Dantzig-Wolfe
        4.  Cplex MIP solver

    This is how the corridor and fixing flags work:
    -c : the constraint is of type <=. Thus, the smaller the value, the
    tighter the corridor. The value indicates the percentage of variables 
    that can change value with respect to the incumbent. To deactivate
    the corridor, we use -c 1.0

    -z : fixing to zero. The constraint is of type <=. Thus, this
    the percentage of variables currently to zero than can change
    value (to one). The smaller the value, the strongest the fixing scheme.
    If we set -z 0.0, this implies "hard fixing," i.e., all the variables
    currently at zero are going to be fixed to zero. If we want to 
    deactivate the scheme, we use -z 1.0 (all the variables currently
    at zero can take value 1.)

    -o : fixing to one. The constraint is of type >=. Thus, the right
    hand side value indicates the percentage of variables currently at one
    that shall be maintained at one. If we fix -o 1.0, this implies hard
    fixing, since it means that all the variables currently at one should
    be kept at the same value. If we want to deactivate the scheme, we
    use -o 0.0.

    """
    global inputfile
    global userCuts
    global cPercent
    global cZero
    global cOne
    global nSolInPool
    global algo

    try:
        opts, args = getopt.getopt(argv, "hi:u:c:z:o:p:a:",
        ["help","ifile=","ucuts","cpercent","zeros","ones","pool", "algorithm"])
    except getopt.GetoptError:
        print("Command Line Error. Usage : python cflp.py -i <inputfile> -u\
        <usercuts> -c <corridor width> -z <fix to zero> -o <fix to one> \
        -p <pool> -a <algorithm BD, LR, DW, Cplex>")
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
        elif opt in ("-p", "--pool"):
            nSolInPool  = int(arg)
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

       #  largeInstances = True
       if largeInstances == True:
           self.readLargeInstances(inputfile)
       else:
           self.readSmallInstances(inputfile)

       #  compute cumulative demand
       for j in range(self.nI):
           aux = []
           aux.append(self.d[j][self.nP-1])
           progr = 0
           for t in range(self.nP-2,-1,-1):
               aux.append(aux[progr]+self.d[j][t])
               progr += 1
           self.dcum.append(list(reversed(aux)))
           #  print("cum dem item ", j , " = ", self.dcum[j])


       # max production of item j in period t is the minimum between
       # the limit set by capacity and the cumulative demand
       for j in range(self.nI):
           aa = [( (self.cap[t] - self.m[j][t])/self.a[j][t]) for t in range(self.nP)]
           self.max_prod.append([min(aa[t],self.dcum[j][t]) for t in \
           range(self.nP)])
           
    
    def readSmallInstances(self, inputfile):
        """
        read Trigeiro instances 
        """

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
            # note : cum[v] gives the cumulative demand from v to T
            # thus, cum[v] - cum[t] gives the cumulative demand between 
            # v (included) and t (excluded)
            for j in range(self.nI):
                aux = []
                aux.append(self.d[j][self.nP-1])
                progr = 0
                for t in range(self.nP-2,-1,-1):
                    aux.append(aux[progr]+self.d[j][t])
                    progr += 1
                self.dcum.append(list(reversed(aux)))

        getInstanceStats(self)

    def readLargeInstances(self, inputfile):
        with open(inputfile) as ff:
            data = ff.readline()
            self.nI, self.nP = [int(v) for v in data.split()]

            self.d = [ [] for i in range(self.nI)]

            print("We have ", self.nI, self.nP)
            #  setup costs
            for i in range(self.nI):
                data = ff.readline()
                temp = [float(v) for v in data.split()]
                self.f.append([temp[1]]*self.nP)
            #  inventory holding costs
            for i in range(self.nI):
                data = ff.readline()
                temp = [float(v) for v in data.split()]
                self.h.append([temp[1]]*self.nP)
            #  demand for each item in each period
            for i in range(self.nI):
                for t in range(self.nP):
                    data = ff.readline()
                    temp = [float(v) for v in data.split()]
                    self.d[i].append(temp[2])

            #  skip cumulative demand
            for i in range(self.nI):
                for t in range(self.nP):
                    for tp in range(t, self.nP):
                        data = ff.readline()
            #  skip cumulative holding costs
            for i in range(self.nI):
                for t in range(self.nP):
                    for tp in range(t, self.nP):
                        data = ff.readline()

            #  resource usage a[j][t]
            for i in range(self.nI):
                data = ff.readline()
                temp = [float(v) for v in data.split()]
                a = [temp[1]]*self.nP
                self.a.append(a)
                self.m.append([0.0]*self.nP)
                self.c.append([0.0]*self.nP)


            for t in range(self.nP):
                data = ff.readline()
                temp = [float(v) for v in data.split()]
                self.cap.append(temp[1])

       #  compute cumulative demand
        #  for j in range(self.nI):
        #      aux = []
        #      aux.append(self.d[j][self.nP-1])
        #      progr = 0
        #      for t in range(self.nP-2,-1,-1):
        #          aux.append(aux[progr]+self.d[j][t])
        #          progr += 1
        #      self.dcum.append(list(reversed(aux)))
        #
        #  print("Cumulative demands are ", self.dcum)
        #  exit(23)

def separateTrigeiro(inputfile):

    print("Open file ", inputfile)
    namefile = "g"
    with open(inputfile) as ff:
        for i in range(5): # 5 instances per type

            print("** CYCLE = ", i," **")
            print("*"*80)
            data = ff.readline()
            nI, nP = [int(v) for v in data.split()]

            print("Nr. items and nr. periods = ", nI, " ", nP)
            namefile = "dataTrigeiroG/g" + str(nI) + "-" + str(nP) + "."  + str(i)
            fdata = open(namefile, "w")
            print("Writing in file", namefile)
            fdata.write(data)
            

            ff.readline() # skip this line (I do not know what this is)

            # capacity
            data = ff.readline()
            fdata.write(data)

            # a, h, m, f, c (c is set to 0)
            for i in range(nI):
                data = ff.readline()
                fdata.write(data)
                print(data)
                
            # demand (for each period, for each item)
            for t in range(nP):
                data = ff.readline()
                fdata.write(data)

            ff.readline()
        # done with this file




def getInstanceStats(inp):
    totUsed = 0
    totAvail = 0
    for t in range(inp.nP):
        usedCap = 0
        for j in range(inp.nI):
            usedCap += inp.m[j][t] + inp.a[j][t]*inp.d[j][t]
        totUsed += usedCap
        totAvail += inp.cap[t]
        print("Period {0} :: {1} vs {2}".format(t,usedCap, inp.cap[t]))

    print("Tot USED = {0} vs Tot AVAIL = {1}".format(totUsed, totAvail))
    #  exit(123)

                

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
        self.zBest = cplex.infinity #  best solution so far
        self.z     = cplex.infinity # current solution
        self.ySol = []
        y_ilo = []
        x_ilo = []
        s_ilo = []
        sI    = []

        cpx = cplex.Cplex()
        cpx.objective.set_sense(cpx.objective.sense.minimize)
        benders = cpx.long_annotations.add("cpxBendersPartition",1)
        objtype = cpx.long_annotations.object_type.variable
        
        cpx.parameters.benders.strategy.set(1)
        cpx.write_benders_annotation("benders.ann")

        #  create variables y_jt
        for j in range(inp.nI):
            y_ilo.append([])
            for t in range(inp.nP):
                varName = "y." + str(j) + "." + str(t)
                y_ilo[j].append(cpx.variables.get_num())
                y_var = cpx.variables.add(obj   = [inp.f[j][t]],
                                  lb    = [0],
                                  ub    = [1],
                                  types = ["B"],
                                  names = [varName])
                cpx.long_annotations.set_values(benders, objtype, [(i,0) for i in y_var])



        #  create variables x_jt
        for j in range(inp.nI):
            x_ilo.append([])
            for t in range(inp.nP):
                varName = "x." + str(j) + "." + str(t)
                x_ilo[j].append(cpx.variables.get_num())
                x_var= cpx.variables.add(obj   = [inp.c[j][t]],
                                  lb    = [0.0],
                                  ub    = [cplex.infinity],
                                  types = ["C"],
                                  names = [varName])

        print("ANNOTATIONS :: ", cpx.long_annotations.get_values(benders,
        objtype))

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

    def addConstraint(self, Lhs, Rhs, direction="L", constrName=""):
        cpx = self.cpx
        cpx.linear_constraints.add(lin_expr = [Lhs],
                                   senses   = [direction],
                                   rhs      = [Rhs],
                                   names    = [constrName])

        self.cpx = cpx

    def removeConstraint(self, constrName):
        cpx = self.cpx
        print("Nr. Constraints Before= ", cpx.linear_constraints.get_num())
        cpx.linear_constraints.delete(constrName)
        print("Nr. Constraints After= ", cpx.linear_constraints.get_num())
        self.cpx = cpx

    def reverseConstraint(self, Lhs, Rhs, direction = "G", constrName=""):
        cpx = self.cpx
        cpx.linear_constraints.set_rhs(constrName, Rhs+1)
        cpx.linear_constraints.set_senses(constrName, direction)
        self.cpx = cpx



        
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

        try:

            y_ilo = self.y_ilo
            #  z_ilo = self.z_ilo
            cpx   = self.cpx

            #  cpx.set_results_stream(None)
            #  cpx.set_log_stream(None)
            cpx.parameters.timelimit.set(timeLimit)
            cpx.parameters.mip.limits.solutions.set(nSol)
            cpx.parameters.mip.display.set(display)
            cpx.parameters.mip.interval.set(5000) # how often to print info
            cpx.parameters.mip.tolerances.mipgap.set(0.000000001)
            cpx.parameters.mip.tolerances.absmipgap.set(0.000000001)


            cpx.solve()
            if cpx.solution.get_status() == cpx.solution.status.MIP_infeasible:
                return -1
            else:
                return 1
        
        except CplexError as exc:
            print("CPLEX ERROR =", exc)
            exit(999)

    def getSolution(self, inp, withPrinting=0):

        cpx = self.cpx
        y_ilo = self.y_ilo

        if withPrinting >= 1:
            print("STATUS \t= ", cpx.solution.status[cpx.solution.get_status()])
            print("OPT SOL found \t= ", cpx.solution.get_objective_value())
            print("BEST BOUND found \t= ", cpx.solution.MIP.get_best_objective())
            print("Time          \t= ", time.time() - startTime)
        #  if cpx.solution.get_status() == cpx.solution.status.optimal_tolerance\
            #  or cpx.solution.get_status() == cpx.solution.status.optimal:
        ubBest = cpx.solution.get_objective_value()
        self.z = cpx.solution.get_objective_value()

        fcsv = csv.writer(open("ySol.csv", "w"))
        yRef = []
        for j in range(inp.nI):
            yRef.append(cpx.solution.get_values(y_ilo[j]))

        if withPrinting == 2:
            for j in range(inp.nI):
                print("y(",j,") = ", yRef[j])
                fcsv.writerow(yRef[j])
                #  for t in range(np.nP):
                #      print("   z(",t,") = ",
                #      [cpx.solution.get_values(z_ilo[j][t][r-t]) for r in
                #      range(t,inp.nP)])

        # if self.z < self.zBest:
        #     self.zBest = self.z
        self.ySol = yRef

        return ubBest

class MIPReformulation(MIP):
    """
    .. class:MIPReformulation()

    This class implements the SPL reformulation, which is the one currently
    used in the code. The reformulation has been presented in the introduction
    of this code and makes use of two sets of variables, i.e.:

    * :math:`y_{jt}` : setup variables
    * :math:`z_{jtr}`: production variables, indicating the amount of production of item *j* produced in period *t* to satisfy the demand of period *r*. Obviously, :math:`z_{jtr} = 0` for all *t>r*.


    NOTE: A modification w.r.t. the :math `z_{jtr}` variables has been
    introduced. Rather than creating all the variables, including those for
    which :math `r < t` (which make no sense, since demand cannot be satisfied
    retroactively), we now define the :math `z_{jtr}` variables only for :math
    `r \geq t`. The difference is now that z_ilo is a vector of vectors of
    unequal size. For example:

    * for t = 0, we define variables :math `z_{j00}` to `z_{j0T}`
    * for t = 1, we only define from :math `z_{j11}` to `z_{j1T}`
    * for t = 2, we only define from :math `z_{j22}` to `z_{j2T}`

    Thus, to access variable :math `z_{jtr}` in cplex, we need to scale back
    the value of r, i.e., access *x_ilo[j][t][r-t]*
    """

    def __init__(self, inp):

        self.zBest = cplex.infinity
        self.ySol = []
        y_ilo = []
        z_ilo = []

        cpx = cplex.Cplex()
        cpx.objective.set_sense(cpx.objective.sense.minimize)
        #  cpx.set_results_stream(None)
        #  cpx.set_log_stream(None)
        cpx.parameters.benders.strategy.set(1) #  use -1 to deactivate BENDERS
        cpx.write_benders_annotation("benders.ann")
        benders = cpx.long_annotations.add("cpxBendersPartition",1)
        objtype = cpx.long_annotations.object_type.variable

        #  create variables y_jt
        for j in range(inp.nI):
            y_ilo.append([])
            for t in range(inp.nP):
                varName = "y." + str(j) + "." + str(t)
                y_ilo[j].append(cpx.variables.get_num())
                y_var = cpx.variables.add(obj   = [inp.f[j][t]],
                                  lb    = [0],
                                  ub    = [1],
                                  types = ["B"],
                                  names = [varName])
                cpx.long_annotations.set_values(benders, objtype, [(i,0) for i in y_var])

        #  create variables z_jts
        for j in range(inp.nI):
            z_ilo.append([])
            for t in range(inp.nP):
                z_ilo[j].append([])
                for r in range(t,inp.nP):
                    varName = "z." + str(j) + "." + str(t) + "." + str(r)
                    z_ilo[j][t].append(cpx.variables.get_num())
                    z_var = cpx.variables.add(obj   = [(r-t)*inp.h[j][t]],
                                      lb    = [0.0],
                                      #  ub    = [cplex.infinity],
                                      #  types = ["C"],
                                      names = [varName])

        #  demand constraints
        for j in range(inp.nI):
            for r in range(inp.nP):
                index = [z_ilo[j][t][r-t] for t in range(r+1)]
                value = [1.0]*(r+1)
                demand_constraint = cplex.SparsePair(ind=index, val=value)
                cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                           senses   = ["G"],
                                           rhs      = [inp.d[j][r]])

        #  capacity constraint
        for t in range(inp.nP):
            index = []
            value = []
            for j in range(inp.nI):
                index += [z_ilo[j][t][r-t] for r in range(t,inp.nP)]
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
                    index = [z_ilo[j][t][r-t], y_ilo[j][t]]
                    value = [1.0, -inp.d[j][r]]
                    logic_constraint = cplex.SparsePair(ind =index, val=value)
                    cpx.linear_constraints.add(lin_expr = [logic_constraint],
                                               senses   = ["L"],
                                               rhs      = [0.0])
        #  #  cumulative logic constraints
        #  for j in range(inp.nI):
        #      for t in range(inp.nP):
        #          index = [z_ilo[j][t][r] for r in range(t,inp.nP)]
        #          value = [1.0]*(inp.nP-t)
        #          index += [y_ilo[j][t]]
        #          value += [-inp.max_prod[j][t]]
        #          logic_constraint_cum = cplex.SparsePair(ind =index, val=value)
        #          cpx.linear_constraints.add(lin_expr = [logic_constraint_cum],
        #                                     senses   = ["L"],
        #                                     rhs      = [0.0])


        #  for t in range(inp.nP):
        #      index = [y_ilo[j][t] for j in range(inp.nI)]
        #      value = [inp.m[j][t]  + inp.a[j][t]*inp.d[j][t] for j in range(inp.nI)]
        #
        #      hop_constraint = cplex.SparsePair(ind=index, val=value)
        #      cpx.linear_constraints.add(lin_expr = [hop_constraint],
        #                                 senses   = ["L"],
        #                                 rhs      = [inp.cap[t]])

        #  self.addSetupCuts(inp)

        self.cpx   = cpx
        self.y_ilo = y_ilo
        self.z_ilo = z_ilo
        self.benders = benders

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
    print("==================================================================")
    print("\n")



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


class WorkerLPPrimal:
    def __init__(self, inp):

        cpx = cplex.Cplex()
        z_ilo = []

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

        #  create variables z_jtr
        for j in range(inp.nI):
            z_ilo.append([])
            for t in range(inp.nP):
                z_ilo[j].append([])
                for r in range(t,inp.nP):
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
                index += [z_ilo[j][t][r-t] for r in range(t,inp.nP)]
                value += [inp.a[j][t]]*(inp.nP-t)

            constrName = "capacity." + str(t)
            capacity_constraint = cplex.SparsePair(ind=index,val=value)
            cpx.linear_constraints.add(lin_expr  = [capacity_constraint],
                                       senses    = ["L"],
                                       rhs       = [0.0],
                                       names     = [constrName])

        #  demand constraints
        for j in range(inp.nI):
            for r in range(inp.nP):
                index = [z_ilo[j][t][r-t] for t in range(r+1)]
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
                    index = [z_ilo[j][t][r-t]]
                    value = [1.0]
                    logic_constraint = cplex.SparsePair(ind =index, val=value)
                    constrName = "logic."+str(j)+"."+str(t)+"."+str(r)
                    cpx.linear_constraints.add(lin_expr = [logic_constraint],
                                               senses   = ["L"],
                                               rhs      = [0.0],
                                               names    = [constrName])

        # define labels for constraints (used in separation to get dual values)
        lCapacity = ["capacity." + str(t) for t in range(inp.nP)]
        #  lDemand   = ["demand." + str(j) + "." + str(r) for j in range(inp.nI)
        #  for r in range(inp.nP)]
        ll = [["demand." + str(j) + "." + str(r) for r in range(inp.nP)] for j in range(inp.nI)]

        lLogic    = ["logic." + str(j) + "." + str(t) + "." + str(r) \
        for j in range(inp.nI) for t in range(inp.nP) for r in range(t,inp.nP)]

        self.cpx       = cpx
        self.z_ilo     = z_ilo
        self.lCapacity = lCapacity
        #  self.lDemand   = lDemand
        self.lLogic    = lLogic
        self.ll = ll

    def solveSubPrimal(self, inp, ySol, zHat, y_ilo, z_ilo_master):
        
        cutType   = 0
        cpx       = self.cpx
        z_ilo     = self.z_ilo
        lCapacity = self.lCapacity
        #  lDemand   = self.lDemand
        lLogic    = self.lLogic
        ll = self.ll

        cpx.set_results_stream(None)
        cpx.set_log_stream(None)

        #  update rhs values : capacity constraints
        for t in range(inp.nP):
            sumY = sum([inp.m[j][t]*ySol[j][t] for j in range(inp.nI)])
            constrName = "capacity." + str(t)
            cpx.linear_constraints.set_rhs(constrName, (inp.cap[t]-sumY))

        #  update rhs values : logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    constrName = "logic." + str(j) + "." + str(t) + "." + str(r)
                    rhsValue = inp.d[j][r]*ySol[j][t]
                    cpx.linear_constraints.set_rhs(constrName, rhsValue)

        cpx.solve()

        cutVars = []
        cutVals = []
        zSub    = -1
        cutRhs = 0.0

        if cpx.solution.get_status() == cpx.solution.status.infeasible:
            cutType = 1
            #  add extreme ray using Farkas certificate
            farkas, pp = cpx.solution.advanced.dual_farkas()

            progr     = inp.nP
            #  index     = range(progr)
            #  fCapacity = [farkas[i] for i in index]
            lambda_sol = [farkas[t] for t in range(inp.nP)]

            index     = range(progr,progr+inp.nP*inp.nI)
            fDemand   = [farkas[i] for i in index]
            omega_sol = []
            for j in range(inp.nI):
                omega_sol.append([fDemand[k] for k in range(j*inp.nP, (j+1)*inp.nP)])

            progr    += inp.nP*inp.nI
            nr        = int( inp.nI*(inp.nP*(inp.nP+1))/2)
            index     = range(progr, progr+nr)
            fLogic    = [farkas[i] for i in index]
            epsilon_sol = []
            pp = 0
            for j in range(inp.nI):
                aux2D = []
                for t in range(inp.nP):
                    aux = []
                    for r in range(t,inp.nP):
                        aux.append(fLogic[pp])
                        pp += 1
                    aux2D.append(aux)
                epsilon_sol.append(aux2D)

            #  feasibility cut
            cutRhs  = 0.0
            for j in range(inp.nI):
                for t in range(inp.nP):
                    coeff  = -inp.m[j][t]*lambda_sol[t]
                    cutRhs -= inp.d[j][t]*omega_sol[j][t]

                    for r in range(t, inp.nP):
                        coeff += inp.d[j][r]*epsilon_sol[j][t][r-t]

                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)

            for t in range(inp.nP):
                cutRhs -= inp.cap[t]*lambda_sol[t]

        elif cpx.solution.get_status() == cpx.solution.status.optimal:
            cutType = 2
            zSub = cpx.solution.get_objective_value()

            lambda_sol = cpx.solution.get_dual_values(lCapacity)
            omega_sol = [ cpx.solution.get_dual_values(ll[j]) for j in range(inp.nI)]


            #  fDemand = cpx.solution.get_dual_values(lDemand)
            
            #  omega_sol = []
            #  for j in range(inp.nI):
            #      omega_sol.append([fDemand[k] for k in range(j*inp.nP, (j+1)*inp.nP)])

            fLogic  = cpx.solution.get_dual_values(lLogic)
            epsilon_sol = []
            pp = 0
            for j in range(inp.nI):
                aux2D = []
                for t in range(inp.nP):
                    aux = []
                    for r in range(t,inp.nP):
                        aux.append(fLogic[pp])
                        pp += 1
                    aux2D.append(aux)
                epsilon_sol.append(aux2D)

            #  epsilon_sol = [[ cpx.solution.get_dual_values(lL[j][t]) for t in range(inp.nP)] for j in range(inp.nI)]
            #
            cutRhs = -zSub
            progr = 0
            for j in range(inp.nI):
                for t in range(inp.nP):
                    coeff = inp.m[j][t]*lambda_sol[t]
                    for r in range(t, inp.nP):
                        coeff -= -epsilon_sol[j][t][r-t]*inp.d[j][r]

                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)

                    cutRhs += coeff*ySol[j][t]
                    progr += 1

            cutVars.append(z_ilo_master)
            cutVals.append(-1.0)

            #  # define cut
            #  for j in range(inp.nI):
            #      for t in range(inp.nP):
            #          coeff = inp.m[j][t]*lambda_sol[t]
            #
            #          for r in range(t, inp.nP):
            #              coeff -= inp.d[j][r]*epsilon_sol[j][t][r-t]
            #
            #          cutVars.append(y_ilo[j][t])
            #          cutVals.append(coeff)
            #
            #  cutVars.append(z_ilo_master)
            #  cutVals.append(-1.0)
            #
            #  cutRhs = 0.0
            #  for j in range(inp.nI):
            #      cutRhs = inp.cap[t]*lambda_sol[t]
            #      for t in range(inp.nP):
            #          cutRhs -= inp.d[j][t]*omega_sol[j][t]

        #  return cut and type
        cutLhs = cplex.SparsePair(ind=cutVars, val=cutVals)

        self.cutLhs = cutLhs
        self.cutRhs = cutRhs

        return cutType, zSub

class WorkerLPQ:
    def __init__(self, inp):
        cpxQ = cplex.Cplex()
        z_ilo = []
        q_ilo = []
        sI    = []
        cpxQ.set_results_stream(None)
        cpxQ.set_log_stream(None)

        #  cpxQ.parameters.preprocessing.presolve.set(
        #          cpxQ.parameters.preprocessing.presolve.values.off)
        #  cpxQ.parameters.preprocessing.reduce.set(0)
        cpxQ.parameters.lpmethod.set(cpxQ.parameters.lpmethod.values.dual)
        cpxQ.objective.set_sense(cpxQ.objective.sense.minimize)

        #  create variables q_jt
        for j in range(inp.nI):
            q_ilo.append([])
            for t in range(inp.nP):
                varName = "q." + str(j) + "." + str(t)
                q_ilo[j].append(cpxQ.variables.get_num())
                cpxQ.variables.add(obj   = [0.0],
                                  #  lb    = [0],
                                  #  ub    = [1],
                                  names = [varName])

        #  create variables z_jtr
        for j in range(inp.nI):
            z_ilo.append([])
            for t in range(inp.nP):
                z_ilo[j].append([])
                for r in range(t,inp.nP):
                    varName = "z." + str(j) + "." + str(t) + "." + str(r)
                    z_ilo[j][t].append(cpxQ.variables.get_num())
                    cpxQ.variables.add(obj   = [(r-t)*inp.h[j][t]],
                                      lb    = [0.0],
                                      #  ub    = [cplex.infinity],
                                      #  ub    = [inp.d[j][r]],
                                      names = [varName])

        # create initial inventory for each period
        for j in range(inp.nI):
            sI.append([])
            for t in range(inp.nP):
                varName = "sI." + str(j) + "." + str(t)
                sI[j].append(cpxQ.variables.get_num())
                cpxQ.variables.add(obj  = [50],
                                   lb   = [0.0],
                                   names = [varName])


        
        #  capacity constraint
        for t in range(inp.nP):
            index = []
            value = []
            for j in range(inp.nI):
                index += [z_ilo[j][t][r-t] for r in range(t,inp.nP)]
                value += [inp.a[j][t]]*(inp.nP-t)
                index += [q_ilo[j][t]]
                value += [inp.m[j][t]]

            constrName = "capacity." + str(t)
            capacity_constraint = cplex.SparsePair(ind=index,val=value)
            cpxQ.linear_constraints.add(lin_expr  = [capacity_constraint],
                                       senses    = ["L"],
                                       rhs       = [inp.cap[t]],
                                       names     = [constrName])

        #  demand constraints
        for j in range(inp.nI):
            for r in range(inp.nP):
                index = [z_ilo[j][t][r-t] for t in range(r+1)]
                value = [1.0]*(r+1)
                index += [sI[j][r]]
                value += [1.0]
                constrName = "demand." + str(j) + "." + str(r)
                demand_constraint = cplex.SparsePair(ind=index, val=value)
                cpxQ.linear_constraints.add(lin_expr = [demand_constraint],
                                           senses   = ["E"],
                                           rhs      = [inp.d[j][r]],
                                           names    = [constrName])
        #  logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    index = [z_ilo[j][t][r-t], q_ilo[j][t]]
                    value = [1.0, -inp.d[j][r]]
                    logic_constraint = cplex.SparsePair(ind =index, val=value)
                    constrName = "logic."+str(j)+"."+str(t)+"."+str(r)
                    cpxQ.linear_constraints.add(lin_expr = [logic_constraint],
                                               senses   = ["L"],
                                               rhs      = [0.0],
                                               names    = [constrName])
        # y-q equivalence
        for j in range(inp.nI):
            for t in range(inp.nP):
                index = [q_ilo[j][t]]
                value = [1.0]
                equiv_constraint = cplex.SparsePair(ind=index, val=value)
                constrName = "equiv." + str(j) + "." + str(t)
                cpxQ.linear_constraints.add(lin_expr = [equiv_constraint],
                                           senses   = ["E"],
                                           rhs      = [0.0],
                                           names    = [constrName])

        lEquiv   = [["equiv." + str(j) + "." + str(t) for t in range(inp.nP)] for j in range(inp.nI)]

        self.cpxQ = cpxQ
        self.z_ilo = z_ilo
        self.q_ilo = q_ilo
        self.sI = sI
        self.lEquiv = lEquiv


    def solveSubQ(self, inp, ySol, zHat, y_ilo, z_ilo_master):

        cpxQ = self.cpxQ
        z_ilo = self.z_ilo
        q_ilo = self.q_ilo
        lEquiv = self.lEquiv
        cutType = 2
        
        #  update rhs values : equivalence constraints
        #  update rhs values : logic constraints
        for j in range(inp.nI):
            for t in range(inp.nP):
                    constrName = "equiv." + str(j) + "." + str(t) 
                    rhsValue = ySol[j][t]
                    cpxQ.linear_constraints.set_rhs(constrName, rhsValue)
        #  for j in range(inp.nI):
            #  for t in range(inp.nP):
                #  cpxQ.variables.set_upper_bounds(q_ilo[j][t], ySol[j][t])
                #  cpxQ.variables.set_lower_bounds(q_ilo[j][t], ySol[j][t])

        cpxQ.solve()
        zSub = cpxQ.solution.get_objective_value()
        cutVars = []
        cutVals = []
        cutRhs = 0.0
        if cpxQ.solution.get_status() == cpxQ.solution.status.optimal:
            rr = [cpxQ.solution.get_dual_values(lEquiv[j]) for j in
            range(inp.nI)]

            cutRhs = -zSub
            for j in range(inp.nI):
                for t in range(inp.nP):

                    coeff = rr[j][t]
                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)

                    cutRhs += coeff*ySol[j][t]

            cutVars.append(z_ilo_master)
            cutVals.append(-1.0)

        #  return cut and type
        cutLhs = cplex.SparsePair(ind=cutVars, val=cutVals)

        self.cutLhs = cutLhs
        self.cutRhs = cutRhs

        return cutType, zSub

class WorkerLPDual:

    def __init__(self, inp):
        """
        Formulation of the subproblem dual. This is no longer used, but I leave
        it here for the sake of completeness. I checked that the formulation is
        correct, since the objective function value of the optimal dual is
        identical to that of the optimal primal.

        """
        cpx = cplex.Cplex()
        l_ilo = []
        w_ilo = []
        e_ilo = []
        #  v_ilo = []

        #  cpx.set_results_stream(None)
        #  cpx.set_log_stream(None)

        # Turn off the presolve reductions and set the CPLEX optimizer
        # to solve the worker LP with primal simplex method.
        cpx.parameters.preprocessing.presolve.set(
                cpx.parameters.preprocessing.presolve.values.off)
        cpx.parameters.preprocessing.reduce.set(0)
        cpx.parameters.lpmethod.set(cpx.parameters.lpmethod.values.primal)
        #  cpx.parameters.lpmethod.set(cpx.parameters.lpmethod.values.dual)
        
        cpx.objective.set_sense(cpx.objective.sense.maximize)

        # lambda variables (capacity constraints)
        for t in range(inp.nP):
            l_ilo.append(cpx.variables.get_num())
            varName = "l." + str(t)
            cpx.variables.add(obj   = [0.0], # to be changed
                              lb    = [-cplex.infinity],
                              ub    = [0.0],
                              names = [varName])
            
        # omega variables (demand constraints)
        for j in range(inp.nI):
            w_ilo.append([])
            for t in range(inp.nP):
                w_ilo[j].append(cpx.variables.get_num())
                varName = "w." + str(j) + "." + str(t)
                cpx.variables.add(obj   = [inp.d[j][t]],
                                  lb    = [0.0],
                                  names = [varName])


        # epsilon variables (logical constraints)
        for j in range(inp.nI):
            e_ilo.append([])
            for t in range(inp.nP):
                e_ilo[j].append([])
                for r in range(t,inp.nP): #  NOTE: from t
                    e_ilo[j][t].append(cpx.variables.get_num())
                    varName = "e." + str(j) + "." + str(t) + "." + str(r)
                    cpx.variables.add(obj    = [0.0], # to be  changed
                              lb    = [-cplex.infinity],
                                      ub     = [0.0],
                                      names  = [varName])

        #  # nu variables (cumulative logical constraints)
        #  for j in range(inp.nI):
        #      v_ilo.append([])
        #      for t in range(inp.nP):
        #          v_ilo[j].append(cpx.variables.get_num())
        #          varName = "v." + str(j) + "." + str(t)
        #          cpx.variables.add(obj   = [0.0], # to be changed
        #                        lb    = [-cplex.infinity],
        #                            ub    = [0.0],
        #                            names = [varName])


        # NOTE: Here w_ilo is with "r", not "j" !!!
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    constrName = "dual." + str(j) + "." + str(t) + "." + str(r)
                    #  index = [w_ilo[j][r], l_ilo[t], v_ilo[j][t], e_ilo[j][t][r-t]]
                    #  value = [1.0, inp.a[j][t], 1.0, 1.0]
                    index = [w_ilo[j][r], l_ilo[t], e_ilo[j][t][r-t]]
                    value = [1.0, inp.a[j][t], 1.0]
                    dual_constraint = cplex.SparsePair(ind=index,val=value)
                    cpx.linear_constraints.add(lin_expr = [dual_constraint],
                                               senses   = ["L"],
                                               rhs      = [(r-t)*inp.h[j][t]],
                                               names    = [constrName])

        self.cpx   = cpx
        self.l_ilo = l_ilo
        self.w_ilo = w_ilo
        self.e_ilo = e_ilo
        #  self.v_ilo = v_ilo


    def solveSubDual(self, inp, ySol, zHat, y_ilo, z_ilo):
        """
        Solve the dual problem (i.e., an explicit way of obtained the dual
        values needed to create the optimality cut. It is no longer used.

        """

        cpx = self.cpx
        l_ilo = self.l_ilo
        w_ilo = self.w_ilo
        e_ilo = self.e_ilo
        #  v_ilo = self.v_ilo

        # update objective function coefficients of dual problem
        # lambda vars
        for t in range(inp.nP):
            coeff = inp.cap[t] 
            for j in range(inp.nI):
                coeff -= inp.m[j][t]*ySol[j][t]

            cpx.objective.set_linear(l_ilo[t], coeff)

        # epsilon vars
        for j in range(inp.nI):
            for t in range(inp.nP):
                for r in range(t, inp.nP):
                    coeff = inp.d[j][r]*ySol[j][t]#NOTE: here it is "r", not "t"
                    #  if ySol[j][t] >= 1.0-_EPSI:
                    #      coeff = inp.d[j][r] #  NOTE: here it is "r", not "t"
                    #  else:
                    #      coeff = 0.0
                    cpx.objective.set_linear(e_ilo[j][t][r-t], coeff)
        #  # nu vars
        #  for j in range(inp.nI):
        #      for t in range(inp.nP):
        #          coeff = inp.max_prod[j][t]*ySol[j][t]
        #          #  if ySol[j][t] >= 1.0-_EPSI:
        #          #      coeff = inp.max_prod[j][t]
        #          #  else:
        #          #      coeff = 0.0
        #          cpx.objective.set_linear(v_ilo[j][t], coeff)

        #  cpx.write("dual.lp")
        cpx.solve()

        cutVars = []
        cutVals = []
        cutType = -1
        zSub = _INFTY
        cutRhs = 0.0

        if cpx.solution.get_status() == cpx.solution.status.unbounded:
            cutType = 1
            print("Dual unbounded, getting FEASIBILITY cut")
            rays = cpx.solution.advanced.get_ray()
            print("Rays = ", rays)

            progr     = inp.nP
            index     = range(progr)
            fCapacity = [rays[i] for i in index]
            lambda_sol = [fCapacity[t] for t in range(inp.nP)]

            
            index     = range(progr,progr+inp.nP*inp.nI)
            fDemand   = [rays[i] for i in index]
            omega_sol = []
            for j in range(inp.nI):
                omega_sol.append([fDemand[k] for k in range(j*inp.nP, (j+1)*inp.nP)])

            progr    += inp.nP*inp.nI
            nr        = int( inp.nI*(inp.nP*(inp.nP+1))/2)
            #  index     = range(progr, progr+inp.nI*inp.nP*inp.nP)
            index = range(progr, progr + nr)
            fLogic    = [rays[i] for i in index]
            epsilon_sol = []
            pp = 0
            for j in range(inp.nI):
                aux2D = []
                for t in range(inp.nP):
                    aux = []
                    for r in range(t,inp.nP):
                        aux.append(fLogic[pp])
                        pp += 1
                    aux2D.append(aux)
                epsilon_sol.append(aux2D)

            #  #  progr    += inp.nI*inp.nP*inp.nP
            #  progr += nr
            #  index     = range(progr, progr+inp.nI*inp.nP)
            #  fcumLogic = [rays[i] for i in index]
            #  nu_sol = []
            #  for j in range(inp.nI):
            #      nu_sol.append([fcumLogic[k] for k in range(j*inp.nP, (j+1)*inp.nP)])

            # define cut
            for j in range(inp.nI):
                for t in range(inp.nP):
                    #  coeff = inp.max_prod[j][t]*nu_sol[j][t] -\
                    coeff = -inp.m[j][t]*lambda_sol[t]

                    for r in range(t, inp.nP):
                        coeff += inp.d[j][r]*epsilon_sol[j][t][r-t]

                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)

            cutRhs = 0.0
            for t in range(inp.nP):
                cutRhs -= inp.cap[t]*lambda_sol[t]
                for j in range(inp.nI):
                    cutRhs -= inp.d[j][t]*omega_sol[j][t]

        elif cpx.solution.get_status() == cpx.solution.status.optimal:
            print("Getting OPTIMALITY cut")
            cutType = 2
            zSub = cpx.solution.get_objective_value()
            print("zDual* = ", zSub)
            lambda_sol = [cpx.solution.get_values(l_ilo[t]) for t in range(inp.nP)]
            omega_sol = [ [cpx.solution.get_values(w_ilo[j][t]) for t in range(inp.nP)] for j in range(inp.nI)]
            epsilon_sol = []
            for j in range(inp.nI):
                aux2D = []
                for t in range(inp.nP):
                    aux = []
                    for r in range(t, inp.nP):
                        aux.append(cpx.solution.get_values(e_ilo[j][t][r-t]))
                    aux2D.append(aux)
                epsilon_sol.append(aux2D)

            #  nu_sol = [ [cpx.solution.get_values(v_ilo[j][t]) for t in range(inp.nP)] for j in range(inp.nI)]
        
            # define cut
            for j in range(inp.nI):
                for t in range(inp.nP):
                    progr = j*inp.nP + t
                    #  coeff = inp.max_prod[j][t]*nu_sol[j][t] -\
                    coeff = -inp.m[j][t]*lambda_sol[t]

                    for r in range(t, inp.nP):
                        coeff += inp.d[j][r]*epsilon_sol[j][t][r-t]

                    cutVars.append(y_ilo[j][t])
                    cutVals.append(coeff)

            cutVars.append(z_ilo)
            cutVals.append(-1.0)

            cutRhs = 0.0
            for j in range(inp.nI):
                cutRhs -= inp.cap[t]*lambda_sol[t]
                for t in range(inp.nP):
                    cutRhs -= inp.d[j][t]*omega_sol[j][t]


        #  return cut and type
        cutLhs = cplex.SparsePair(ind=cutVars, val=cutVals)

        self.cutLhs = cutLhs
        self.cutRhs = cutRhs
        
        return cutType, zSub


def mipLPInterior(inp, mip):
    y_ilo = mip.y_ilo
    cpxLP = cplex.Cplex(mip.cpx)
    print("Solving here for interior ... ")
    cpxLP.set_problem_type(cpxLP.problem_type.LP)
    print("Problem type is ", cpxLP.problem_type[cpxLP.get_problem_type()])
    cpxLP.parameters.lpmethod.set(cpxLP.parameters.lpmethod.values.barrier)
    #  cpx.parameters.solutiontype.set(2)
    cpxLP.parameters.barrier.crossover.set(-1) # no crossover
    cpxLP.parameters.barrier.convergetol.set(1e-3) # no crossover

    cpxLP.solve()
    ySol = []
    for j in range(inp.nI):
        ySol.append(cpxLP.solution.get_values(y_ilo[j]))
    #
    #  print(ySol)
    #  input(" ... barrier ... ")
    return ySol

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

def newIOCycle(cpx, worker, y_ilo, z_ilo, inp):

    _alpha = 0.9
    _lambda = 0.1
    _delta = 0.00002
    _delta = 0.000000004

    globalProgr = 1

    # solve current master LP
    cpx.set_problem_type(cpx.problem_type.LP)
    cpx.parameters.simplex.display.set(0)
    yTilde = []
    for j in range(inp.nI):
        aux = [1]*inp.nP
        yTilde.append(aux)

    stopping = False
    noImprov = 0
    bestLB   = 0.0
    while not stopping:

        cpx.solve()
        zLP = cpx.solution.get_objective_value()
        if globalProgr % 25 == 0:
            print("[{0:5d}] Best LB = {1:9.2f}; zLP = {2:9.2f}" .format(globalProgr, bestLB, zLP))

        if zLP > bestLB:
            bestLB = zLP
        elif noImprov <= 5:
            noImprov += 1
        elif _delta == 0.0:
            stopping = True
        elif _lambda == 1.0:
            _delta = 0.0
            noImprov = 0
        else:
            _lambda = 1.0
            noImprov = 0

        zHat = cpx.solution.get_values(z_ilo)
        yStar = []
        for j in range(inp.nI):
            yStar.append(cpx.solution.get_values(y_ilo[j]))

        #  print("1 - Y* = ", yStar)

        # get interior point

        for j in range(inp.nI):
            for t in range(inp.nP):
                yTilde[j][t] = _alpha*yTilde[j][t] + (1.0-_alpha)*yStar[j][t]
        #  print("2 - yTilde = ", yTilde)

        yy = []
        for j in range(inp.nI):
            aux = []
            for t in range(inp.nP):
                aux.append(min(1.0, _lambda*yStar[j][t] + (1.0-_lambda)*yTilde[j][t] + _delta)) 
            yy.append(aux)

        #  print("3 - Separator 3 - Separator 3 - Separator :: ", yy)
        #  cutType, zSub = worker.solveSubPrimal(inp, yy, zHat, y_ilo, z_ilo)
        #  cutType, zSub = worker.solveSubDual(inp, yy, zHat, y_ilo, z_ilo)
        cutType, zSub = worker.solveSubQ(inp, yy, zHat, y_ilo, z_ilo)

        #  if globalProgr% 5 == 0:
        #  NOTE: removing these constraints does not seem useful
        #  if stopping:
        #      check = [i for i in inout]
        #      #  print("checking these constraints ", check)
        #
        #      slacks = cpx.solution.get_linear_slacks(check)
        #      #  nrConstr = cpx.linear_constraints.get_num()
        #      remove = [check[i] for i in range(len(slacks)) if slacks[i] > _EPSI]
        #      #  print("removing ", remove)
        #      cpx.linear_constraints.delete(remove)
        #      for i in remove:
        #          inout.remove(i)
        #      print("Removed ", len(remove), "constraints.")
        #      #  print("This is inout = ", inout)

        if cutType > 0:
            nrConstr = cpx.linear_constraints.get_num()
            constrName = "inout." + str(globalProgr)
            globalProgr += 1
            #  print("adding ", constrName)
            inout.append(constrName)
            cpx.linear_constraints.add(lin_expr = [worker.cutLhs],
                                       senses   = "L",
                                       rhs      = [worker.cutRhs],
                                       names    = [constrName])

        else:
            print("something wrong with cuttype !!!")
            input("aka")


def setCpxParameters(cpx):
    
    #  cpx.parameters.mip.limits.solutions.set(3)
    cpx.parameters.mip.display.set(1)
    cpx.parameters.preprocessing.presolve.set(cpx.parameters.preprocessing.presolve.values.off)
    cpx.parameters.threads.set(1)
    cpx.parameters.mip.strategy.search.set(cpx.parameters.mip.strategy.search.values.traditional)
    #  cpx.parameters.timelimit.set(5)
    # call RINS at every node
    #  cpx.parameters.mip.strategy.rinsheur.set(1)

def getUB(inp, zSub, ySol, zSol):
    z = zSub
    for j in range(inp.nI):
        for t in range(inp.nP):
            if ySol[j][t] > 1.0-_EPSI:
                z += inp.f[j][t]

    # is this feasible?
    for j in range(inp.nI):
        for t in range(inp.nP):
            for r in range(t, inp.nP):
                if zSol[j][t][r] > _EPSI:
                    assert(ySol[j][t] > 1.0-_EPSI)
    return z

def getSolution(inp, cpx, y_ilo, z_ilo):
    ySol = []
    zHat = cpx.solution.get_values(z_ilo)
    lb   = cpx.solution.get_objective_value()
    for j in range(inp.nI):
        aux = []
        for t in range(inp.nP):
            aux.append(cpx.solution.get_values(y_ilo[j][t]))
        ySol.append(aux)

    return lb, zHat, ySol

def bendersCallbackScheme(inp,mip):

    cpx = cplex.Cplex()
    createMaster(inp, cpx)
    #  worker = WorkerLPPrimal(inp)
    #  worker = WorkerLPDual(inp)
    workerPrimal = WorkerLPPrimal(inp)
    worker = WorkerLPQ(inp)
    # Set up cplex parameters to use the cut callback for separating
    # Benders' cuts
    #  cpx.parameters.preprocessing.presolve.set(
    #      cpx.parameters.preprocessing.presolve.values.off)
    # Set the maximum number of threads to 1.
    cpx.parameters.threads.set(4)
    cpx.parameters.parallel.set(-1)
    cpx.parameters.mip.cuts.cliques.set(3)
    cpx.parameters.mip.cuts.covers.set(3)
    cpx.parameters.mip.cuts.disjunctive.set(3)
    cpx.parameters.mip.cuts.flowcovers.set(2)
    cpx.parameters.mip.cuts.pathcut.set(2)
    cpx.parameters.mip.cuts.gomory.set(2)
    cpx.parameters.mip.cuts.gubcovers.set(2)
    cpx.parameters.mip.cuts.implied.set(2)
    cpx.parameters.mip.cuts.mircut.set(2)
    cpx.parameters.mip.cuts.zerohalfcut.set(2)
    # Turn on traditional search for use with control callbacks
    cpx.parameters.mip.strategy.search.set(
        cpx.parameters.mip.strategy.search.values.traditional)
    # call RINS at every node
    cpx.parameters.mip.strategy.rinsheur.set(1)

    cpx.parameters.mip.strategy.variableselect.set(3)
    cpx.parameters.mip.strategy.nodeselect.set(1)
    cpx.parameters.mip.strategy.bbinterval.set(1)
    cpx.parameters.mip.strategy.backtrack.set(0)

    print("Problem type is ", cpx.problem_type[cpx.get_problem_type()])
    #  newIOCycle(cpx, workerPrimal, y_ilo, z_ilo, inp)
    newIOCycle(cpx, worker, y_ilo, z_ilo, inp)
    #  cpx.write("inout-6-15.lp")
    #  cpx.read("inout-6-15.lp")
    #
    cpx.set_problem_type(cpx.problem_type.MILP)
    for j in range(inp.nI):
        for t in range(inp.nP):
            cpx.variables.set_types(y_ilo[j][t], cpx.variables.type.binary)
    print("Problem type is ", cpx.problem_type[cpx.get_problem_type()])

    startTime = time.time()

    # register LAZY callback
    lazyBenders           = cpx.register_callback(BendersLazyConsCallback)
    lazyBenders.cpx       = cpx
    lazyBenders.inp       = inp
    lazyBenders.mip       = mip
    lazyBenders.z_ilo     = z_ilo
    lazyBenders.y_ilo     = y_ilo
    lazyBenders.worker    = worker
    lazyBenders.solved    = 0
    lazyBenders.rc        = []
    lazyBenders.nIter     = 0
    lazyBenders.bestLB    = -_INFTY
    lazyBenders.bestUB    = _INFTY
    lazyBenders.Lhs       = []
    lazyBenders.Rhs       = []
    lazyBenders.nCM       = 0
    lazyBenders.bestCM    = _INFTY
    lazyBenders.bestTime  = -1
    lazyBenders.startTime = startTime



    # Solve the model
    cpx.parameters.mip.limits.nodes.set(1)
    cpx.solve()

    # get cuts here
    print("Cuts got here ", len(lazyBenders.Lhs))
    progr = 0
    for i in range(len(lazyBenders.Lhs)):
        cpx.linear_constraints.add(lin_expr = [lazyBenders.Lhs[i]],
                                   senses   = ["L"],
                                   rhs      = [lazyBenders.Rhs[i]],
                                   names    = ["staticCut."+str(progr)])
        progr += 1

    lazyBenders.Lhs    = []
    lazyBenders.Rhs    = []
    cpx.parameters.mip.limits.nodes.set(1)
    cpx.solve()

    # get cuts here
    print("Cuts got here ", len(lazyBenders.Lhs))
    for i in range(len(lazyBenders.Lhs)):
        cpx.linear_constraints.add(lin_expr = [lazyBenders.Lhs[i]],
                                   senses   = ["L"],
                                   rhs      = [lazyBenders.Rhs[i]],
                                   names    = ["staticCut."+str(progr)])
        progr += 1

    cpx.parameters.mip.limits.nodes.set(9223372036800000000)
    cpx.parameters.timelimit.set(10)
    cpx.solve()
                                   
    print("MAIN status = ", cpx.solution.get_status())
    print("\n")
    print("*"*80)
    if cpx.solution.get_status() == cpx.solution.status.optimal:
        solution = cpx.solution
        print()
        print("Solution status: ", solution.status[solution.get_status()])
        print("Objective value: ", solution.get_objective_value())
        print("Time to best   : ", lazyBenders.bestTime)

    else:
        print("Cplex lb \t=", lazyBenders.bestUB)
        print("Cplex gap \t=", lazyBenders.gap)
        print("Best CM \t=", lazyBenders.bestCM)
        print("Time to best \t=", lazyBenders.bestTime)

        print("Total Time\t=", time.time() - startTime)
    
    writeBendersSol2Disk("result.txt", inputfile, cpx, lazyBenders)

    print("*"*80)
    


def bendersDual(inp):

    cpx = cplex.Cplex()
    cpx.set_results_stream(None)
    cpx.set_log_stream(None)
    createMaster(inp, cpx)
    #  worker = WorkerLPDual(inp)
    worker = WorkerLPPrimal(inp)
    setCpxParameters(cpx)

    #  cpx.read("inout-6-15.lp")

    stopping = False
    nIters = 0
    nrCuts = 0
    ubBest = _INFTY
    bestLB = -_INFTY
    print("*"*80)
    print("Staring BENDERS Cycle")
    print("*"*80)

    while not stopping:
        nIters += 1
        cpx.write("master.lp")
        cpx.solve() # solve current Master
        bestLB, zHat, ySol = getSolution(inp, cpx, y_ilo, z_ilo)
        print("[{0}] Current MASTER solution (lb) :: {1}".format(nIters, bestLB))

        #  cutType, zSub = worker.solveSubDual(inp, ySol, zHat, y_ilo, z_ilo)
        cutType, zSub = worker.solveSubPrimal(inp, ySol, zHat, y_ilo, z_ilo)
        cutName = "cut." + str(nrCuts)

        if cutType == 1:
            print("\t Adding feasibiity cut = ", cutName)
        elif cutType == 2:
            print("\t Adding optimality cut = ", cutName, "[zSub = ", zSub,"]")
            fixedCost = bestLB - zHat
            ub = fixedCost + zSub
            if ub < ubBest:
                ubBest = ub

        nrCuts += 1

        cpx.linear_constraints.add(lin_expr = [worker.cutLhs],
                                  senses   = ["L"],
                                  rhs      = [worker.cutRhs],
                                  names    = [cutName])


        if nIters % 10 == 0:

            print("-"*30)
            print("Summary BENDERS status :: ")
            print("Tot Nr. Cuts so far = ", nrCuts)
            print("Best ub = \t", ubBest)
            print("Best lb = \t", bestLB)
            print("-"*30)
            print("\n")

        if ubBest - bestLB < _EPSI:
            stopping = True

def writeToDisk(namefile, instance, cpx):
    with open(namefile, "w") as ftxt:
        status = cpx.solution.status[cpx.solution.get_status()]
        zStar = "{0:9.2f}".format(cpx.solution.get_objective_value())
        lb    = "{0:9.2f}".format(cpx.solution.MIP.get_best_objective())
        totTime = "{0:5.2f}".format(time.time() - startTime)
        gap = "{0:7.5f}".format(cpx.solution.MIP.get_mip_relative_gap())
        row = instance + "\t" + str(status) + "\t"  +  str(zStar) + "\t" + str(totTime) + "\t" + str(lb) + "\t" + str(gap) + "\n"
        ftxt.write(row)

def writeBendersSol2Disk(namefile, instance, cpx, lazyBenders):

    with open(namefile, "w") as ff:
        status = cpx.solution.status[cpx.solution.get_status()]
        if cpx.solution.get_status() == cpx.solution.status.optimal:
            zStar = "{0:9.2f}".format(cpx.solution.get_objective_value())
            lb    = "{0:9.2f}".format(cpx.solution.MIP.get_best_objective())
            totTime = "{0:5.2f}".format(time.time() - startTime)
            gap = "{0:7.5f}".format(cpx.solution.MIP.get_mip_relative_gap())
            row = instance + "\t" + str(status) + "\t"  +  str(zStar) + "\t" + str(totTime) + "\t" + str(lb) + "\t" + str(gap) + "\n"
            ff.write(row)
        else:
            zStar = "{0:9.2f}".format(lazyBenders.bestCM)
            lb    = "{0:9.2f}".format(lazyBenders.bestLB)
            totTime = "{0:5.2f}".format(lazyBenders.bestTime)
            gap     = "{0:7.5f}".format(lazyBenders.gap)

            row = instance + "\t" + str(status) + "\t"  +  str(zStar) + "\t" + str(totTime) + "\t" + str(lb) + "\t" + str(gap) + "\n"
            ff.write(row)




    


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
    global yPool #  contains all the solutions found by CM added as cuts
    yPool = []
    global ubBest
    global yBest
    yBest = []
    ubBest = _INFTY

    global startTime
        #  outfile.write("{0:20s} {1:20.5f} {2:25s} {3:20.5f} {4:20.7f} {5:20.7f}\n".
        #                format(inputfile, zOpt, stat, lb, gap, zTime))

    #  separateTrigeiro("../dantzig-wolfe/data/trigeiro/G2430.dat")
    #  exit(123)

    parseCommandLine(argv)
    inp = Instance(inputfile)
    startTime = time.time()
    printParameters()

    #  mip       = MIPReformulation(inp)
    #  mip.solve(inp, withPrinting = 1, display = 4)
    #  exit(104)

    if algo == 4: #  Cplex MIP solver
        #  mip       = MIP(inp)
        mip       = MIPReformulation(inp)
        status = mip.solve(inp, withPrinting = 1, display = 4, timeLimit=1800)
        if status == 1:
            mip.getSolution(inp, withPrinting=1)
            writeToDisk("result.txt", inputfile, mip.cpx)
        exit(104)
    if algo == 3: # Cplex with callbacks
        mip       = MIPReformulation(inp)
        #  bendersAlgo(inp, mip)
        bendersCallbackScheme(inp, mip)
    if algo == 1: # Cplex with callbacks
        mip       = MIPReformulation(inp)
        #  bendersAlgo(inp, mip)
        bendersDual(inp)


if __name__ == "__main__":
    main(sys.argv[1:])
