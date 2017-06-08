#!/usr/bin/python
"""
===========================================================================

:Filename: dw2.py
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

This code implements a Dantzig-Wolfe reformulation scheme for the multi-item
multi-period capacitated lot sizing problem.

"""

import sys, getopt
import math
import cplex
from time import time
import numpy as np
import random
from decimal import *
random.seed(27)

_INFTY   = sys.float_info.max
#  _EPSI = sys.float_info.epsilon
_EPSI    = 0.00001
bigM     = 1.0e8
nSolInPool =1000000
maxIter  = 40


class DantzigWolfe:
    """
    This class implements a Dantzig-Wolfe reformulation scheme for the CLSP.
    """
    def __init__(self, inp):
        self.l_ilo = []
        self.n_ilo = []
        self.y_ilo = []
        self.x_ilo = []
        self.s_ilo = []

        self.yWW   = []
        self.xWW   = []
        self.sWW   = []

        for j in range(inp.nI):
            self.l_ilo.append([])
            self.n_ilo.append([])

        self.cpxMaster = cplex.Cplex()
        self.cpxMaster.objective.set_sense(self.cpxMaster.objective.sense.minimize)
        self.cpxMaster.set_results_stream(None)
        self.cpxMaster.set_log_stream(None)
        self.define_master(inp)

        u, alpha, theta = self.initialize_multipliers(inp)

        self.cpxSub = cplex.Cplex()
        self.cpxSub.objective.set_sense(self.cpxSub.objective.sense.minimize)
        self.cpxSub.set_results_stream(None)
        self.cpxSub.set_log_stream(None)
        self.cpxSub.set_warning_stream(None)
        self.define_ulsp(inp)


    def initialize_multipliers(self, inp):
        u     = [0.0]*inp.nP
        alpha = [0.0]*inp.nI
        theta = 0.0

        return u, alpha, theta


    def define_master(self, inp):
        """
        Define the master problem and include the column corresponding to an
        arbitrarily high initial inventory (s_{j0}), with high cost. This ensures
        that a feasible solution is always available.
        """
        cpxMaster = self.cpxMaster
        l_ilo     = self.l_ilo
        n_ilo     = self.n_ilo
        k         = 0  #  column 0 corresponds to initial inventory

        for j in range(inp.nI):
            varName = "lambda." + str(j) + "." + str(k)
            l_ilo[j].append(cpxMaster.variables.get_num())
            cpxMaster.variables.add(obj   = [0.0],
                              lb    = [0],
                              ub    = [1],
                              # types = ["B"],
                              names = [varName])

            varName = "mu." + str(j) + "." + str(k) + "." + str(0)
            n_ilo[j].append([])
            n_ilo[j][k].append(cpxMaster.variables.get_num())

            cpxMaster.variables.add(obj = [bigM],
                lb                = [0.0],
                ub                = [1.0],
                # types             = ["C"],
                names             = [varName])
        
        #  add constraints and define names
        for t in range(inp.nP):
            constrName = "capacity." + str(t)
            index = []
            value = []
            capacity_constraint = cplex.SparsePair(ind=index, val=value)
            cpxMaster.linear_constraints.add(lin_expr = [capacity_constraint],
                                       senses   = ["L"],
                                       rhs      = [inp.cap[t]],
                                       names    = [constrName])

        for j in range(inp.nI):
            constrName = "convexity." + str(j)
            index = [l_ilo[j][k]]
            value = [1.0]
            convexity_constraint = cplex.SparsePair(ind=index, val=value)
            cpxMaster.linear_constraints.add(lin_expr = [convexity_constraint],
                                       senses   = ["E"],
                                       rhs      = [1.0],
                                       names    = [constrName])

        for j in range(inp.nI):
            constrName = "linking." + str(j) + "." + str(k)
            index = [n_ilo[j][k][0], l_ilo[j][k]]
            value = [1.0, -1.0]
            linking_constraint = cplex.SparsePair(ind=index, val=value)
            cpxMaster.linear_constraints.add(lin_expr = [linking_constraint],
                                       senses   = ["E"],
                                       rhs      = [0.0],
                                       names    = [constrName])

        #  define labels to access constraints
        lCapacity = ["capacity." + str(t) for t in range(inp.nP)]
        lConvexity = ["convexity." + str(j) for j in range(inp.nI)]
        lLinking = ["linking." + str(j) + "." + str(k) for j in
        range(inp.nI)] 

        self.cpxMaster  = cpxMaster
        self.lCapacity  = lCapacity
        self.lConvexity = lConvexity
        self.lLinking   = lLinking

    def get_master_solution(self, inp):
        cpxMaster = self.cpxMaster
        u         = cpxMaster.solution.get_dual_values(self.lCapacity)
        alpha     = cpxMaster.solution.get_dual_values(self.lConvexity)

        l_ilo = self.l_ilo
        theta = []
        for j in range(inp.nI):
            aux = []
            for k in range(len(l_ilo[j])):
                constrName = "linking." + str(j) + "." + str(k)
                #  print("alpha_",j,".",k,"=",
                #  cpxMaster.solution.get_dual_values(constrName))
                aux.append(cpxMaster.solution.get_dual_values(constrName))
            theta.append(aux)

        return u, alpha, theta

    def define_ulsp(self, inp):
        cpxSub = self.cpxSub
        y_ilo  = self.y_ilo
        x_ilo  = self.x_ilo
        s_ilo  = self.s_ilo

        for t in range(inp.nP):
            varName = "y." + str(t)
            y_ilo.append(cpxSub.variables.get_num())
            cpxSub.variables.add(obj   = [0.0],
                                 lb    = [0],
                                 ub    = [1],
                                 types = ["B"],
                                 names = [varName])

            varName = "x." + str(t)
            x_ilo.append(cpxSub.variables.get_num())
            cpxSub.variables.add(obj   = [0.0],
                                 lb    = [0.0],
                                 ub    = [cplex.infinity],
                                 types = ["C"],
                                 names = [varName])

            varName = "s." + str(t)
            s_ilo.append(cpxSub.variables.get_num())
            cpxSub.variables.add(obj   = [0.0],
                                 lb    = [0.0],
                                 ub    = [cplex.infinity],
                                 types = ["C"],
                                 names = [varName])

        # demand constraints
        #  first period
        constrName = "demand." + str(0)
        index = [x_ilo[0], s_ilo[0]]
        value = [1.0, -1.0]
        #  value = [1.0, 1.0, -1.0]
        demand_constraint = cplex.SparsePair(ind=index, val=value)
        cpxSub.linear_constraints.add(lin_expr = [demand_constraint],
                                   senses   = ["E"],
                                   rhs      = [0.0],
                                   names    = [constrName])
        #  periods 2 to T-1
        for t in range(1,inp.nP-1):
            constrName = "demand." + str(t)
            index = [x_ilo[t], s_ilo[t-1], s_ilo[t]]
            value = [1.0, 1.0, -1.0]
            demand_constraint = cplex.SparsePair(ind=index, val=value)
            cpxSub.linear_constraints.add(lin_expr = [demand_constraint],
                                       senses   = ["E"],
                                       rhs      = [0.0],
                                       names    = [constrName])

        #  last period
        constrName = "demand." + str(inp.nP-1)
        index = [x_ilo[inp.nP-1], s_ilo[inp.nP-2]]
        value = [1.0, 1.0]
        demand_constraint = cplex.SparsePair(ind=index, val=value)
        cpxSub.linear_constraints.add(lin_expr = [demand_constraint],
                                   senses   = ["E"],
                                   rhs      = [0.0],
                                   names    = [constrName])

        #  logic constraints
        for t in range(inp.nP):
            constrName = "logic." + str(t)
            index = [x_ilo[t], y_ilo[t]]
            value = [1.0, 0.0]
            logic_constraint = cplex.SparsePair(ind =index, val=value)
            cpxSub.linear_constraints.add(lin_expr = [logic_constraint],
                                       senses   = ["L"],
                                       rhs      = [0.0],
                                       names    = [constrName])

        self.cpxSub = cpxSub
        self.y_ilo = y_ilo
        self.x_ilo = x_ilo
        self.s_ilo = s_ilo

    def solve_ulsp(self, inp, u, alpha, theta, item):
        """
        Solve ULS for a given set of multipliers u and alpha.
        The item for which this problem is solved is given in "item".
        """
        cpxSub = self.cpxSub
        y_ilo  = self.y_ilo
        x_ilo  = self.x_ilo
        s_ilo  = self.s_ilo

        # change objective function coefficients and bounds
        for t in range(inp.nP):
            coeff = inp.f[item][t] - u[t]*inp.m[item][t]
            cpxSub.objective.set_linear(y_ilo[t], coeff)
            coeff = inp.c[item][t] - u[t]*inp.a[item][t]
            cpxSub.objective.set_linear(x_ilo[t], coeff)
            cpxSub.objective.set_linear(s_ilo[t], inp.h[item][t])
            cpxSub.variables.set_upper_bounds(x_ilo[t], inp.max_prod[item][t])
            cpxSub.variables.set_upper_bounds(s_ilo[t], inp.max_prod[item][t])

        # change rhs of demand  and logic constraints
        for t in range(inp.nP):
            constrName = "demand." + str(t)
            cpxSub.linear_constraints.set_rhs(constrName, inp.d[item][t])
            constrName = "logic." + str(t)
            cpxSub.linear_constraints.set_coefficients(constrName, y_ilo[t],
            -inp.max_prod[item][t])

        cpxSub.write("sub" + str(item) + ".lp")
        cpxSub.solve()
        rc  = cpxSub.solution.get_objective_value()
        rc -= alpha[item]
        #  print ("z_WW(",item,") = ", rc+alpha[item], " and rc = ", rc)

        yWW  = [cpxSub.solution.get_values(y_ilo[t]) for t in range(inp.nP)]
        xWW  = [cpxSub.solution.get_values(x_ilo[t]) for t in range(inp.nP)]
        sWW  = [cpxSub.solution.get_values(s_ilo[t]) for t in range(inp.nP)]

        #  recompute
        #  z = 0.0
        #  for t in range(inp.nP):
        #      if yWW[t] >= 1.0 - _EPSI:
        #          z += inp.f[item][t] - u[t]*inp.m[item][t]
        #      z += (inp.c[item][t] - u[t]*inp.a[item][t])*xWW[t]
        #      z += sWW[t]*inp.h[item][t]
        #  z -= alpha[item]
        #  print(" ... ... ... recomputed rc = ", z)
        #  if abs(z - rc) > 0.00001:
        #      print("problems here ... ")
        #      print("y = ", yWW, " with costs = ",
        #      [inp.f[item][t]-u[t]*inp.m[item][t] for t in range(inp.nP)])
        #      print("x = ", xWW, " with cots  = ",
        #      [inp.c[item][t]-u[t]*inp.a[item][t] for t in range(inp.nP)])
        #      print("s = ", sWW, " with costs = ", inp.h[item][t])
        #
        #      input (" aka ")




        # if item == 2:
        #     print("Sol item 2 = ")
        #     print(yWW)
        #     print("  ", xWW)
        #     print("  ", sWW)
        #     print("sI = ", sIWW)
        #     z = 0
        #     for t in range(inp.nP):
        #         z += (inp.f[item][t] - u[t]*inp.m[item][t])*yWW[t]
        #         z += (inp.c[item][t] - u[t]*inp.a[item][t])*xWW[t]
        #         z += sWW[t]*inp.h[item][t]
        #
        #     print(" --- > recomputed z = ", z)
        #
        #     zOpt = []
        #     yOpt =  [0.0, 1.0, 0.0, 1.0, 0.0]
        #     xOpt =  [0.0, 201.0, 0.0, 185.0, 0.0]
        #     sOpt =  [0.0, 124.0, 0.0, 79.0, 0.0]
        #
        #     z = 0.0
        #     for t in range(inp.nP):
        #         z += (inp.f[item][t] - u[t]*inp.m[item][t])*yOpt[t]
        #         z += (inp.c[item][t] - u[t]*inp.a[item][t])*xOpt[t]
        #         z += sOpt[t]*inp.h[item][t]
        #
        #     print(" --- comparing sol vs optimal : ", z)
        #     if z < rc:
        #         print("problem here ... ")
        #         exit(111)


        self.yWW  = yWW
        self.xWW  = xWW
        self.sWW  = sWW

        return rc


    def add_column_to_master(self, inp, item, q):
        """
        Note: "q" is the index of the convex set of dominated solutions.
        q = 0 --> undominated WW solution
        """
        
        lLinking  = self.lLinking
        l_ilo     = self.l_ilo
        n_ilo     = self.n_ilo
        cpxMaster = self.cpxMaster
        yWW       = self.yWW
        xWW       = self.xWW
        sWW       = self.sWW

        nrVars     = len(l_ilo[item])
        varName   = "lambda." + str(item) + "." + str(nrVars)

        l_ilo[item].append(cpxMaster.variables.get_num())
        cpxMaster.variables.add(obj   = [0.0],
                                lb    = [0],
                                ub    = [1],
                                names = [varName])

        varName = "mu." + str(item) + "." + str(nrVars) + "." + str(q)
        n_ilo[item].append([])
        n_ilo[item][nrVars].append(cpxMaster.variables.get_num())
        coeff = sum([inp.f[item][t]*yWW[t] + inp.c[item][t]*xWW[t] +
        inp.h[item][t]*sWW[t] for t in range(inp.nP)]) 
        cpxMaster.variables.add(obj   = [coeff],
                                lb    = [0],
                                ub    = [1],
                                names = [varName])

        # add variables to the corresponding constraints
        for t in range(inp.nP):
            constrName = "capacity." + str(t)
            if yWW[t] > 1.0 - _EPSI:
                cap = inp.m[item][t] + inp.a[item][t]*xWW[t]
            else:
                cap = 0.0
            cpxMaster.linear_constraints.set_coefficients(constrName,
            n_ilo[item][nrVars][q], cap)
        #     cpxMaster.linear_constraints.set_coefficients(zip(constrName, index,
        #     value))

        constrName = "convexity." + str(item)
        cpxMaster.linear_constraints.set_coefficients(constrName,
        l_ilo[item][nrVars],1.0)

        constrName = "linking." + str(item) + "." + str(nrVars)
        index = [n_ilo[item][nrVars][q], l_ilo[item][nrVars]]
        value = [1.0, -1.0]
        linking_constraint = cplex.SparsePair(ind=index, val=value)
        cpxMaster.linear_constraints.add(lin_expr = [linking_constraint],
                                   senses   = ["E"],
                                   rhs      = [0.0],
                                   names    = [constrName])
        lLinking.append(constrName)

        self.lLinking  = lLinking
        self.cpxMaster = cpxMaster
        self.l_ilo     = l_ilo
        self.n_ilo     = n_ilo

    def compute_rc(self, inp, item, y, x, s, u, alpha):
        rc = 0.0
        for t in range(inp.nP):
            if y[t] >= 1.0-_EPSI:
                rc += inp.f[item][t] - u[t]*inp.m[item][t]
                rc += (inp.c[item][t] - u[t]*inp.a[item][t])*x[t]
            if s[t] > _EPSI:
                rc += inp.h[item][t]*s[t]
        rc -= alpha[item]

        return rc
            
    def column_generation(self, inp, u, alpha, theta):
        cpxSub = self.cpxSub
        stopping = True
        for j in range(inp.nI):
            #  print("====================")
            #  print("* ITEM ", j)
            #  print("====================")
            rc = self.solve_ulsp(inp, u, alpha, theta, j)
            #  print("WW solution (q=0)")
            #  print("yWW :: ", self.yWW)
            #  print("xWW :: ", self.xWW)
            #  print("sWW :: ", self.sWW)
            #  if rc < -_EPSI or self.iterNr <= 1:
            if rc < 10000*_EPSI or self.iterNr <= 1:
                stopping = False
                self.add_column_to_master(inp, j, 0)
                self.get_columns(inp, j, u, alpha)
            #  input(".. next item ..")

        return stopping
        
    def add_dominated_column_to_master(self, inp, item, q):
        cpxMaster = self.cpxMaster
        lLinking  = self.lLinking
        l_ilo     = self.l_ilo
        n_ilo     = self.n_ilo
        yWW       = self.yWW
        xWW       = self.xWW
        sWW       = self.sWW
        #  print("Solution q = ", q)
        #  print("y = ", yWW)
        #  print("x = ", xWW)
        #  print("s = ", sWW)

        nrVars    = len(l_ilo[item])-1
        nrSubVar = len(n_ilo[item][nrVars])
        #  print("How many mu vars do we have for this lambda?")
        #  print("nr = ", nrSubVar)
        #  print("vs = ", n_ilo[item][nrVars])
        #  input("aaa")
        #  varName = "mu." + str(item) + "." + str(nrVars) + "." + str(nrSubVar)
        varName = "mu." + str(item) + "." + str(nrVars) + "." + str(q)
        n_ilo[item][nrVars].append(cpxMaster.variables.get_num())
        coeff = sum([inp.f[item][t]*yWW[t] + inp.c[item][t]*xWW[t] +
        inp.h[item][t]*sWW[t] for t in range(inp.nP)]) 

        cpxMaster.variables.add(obj   = [coeff],
                                lb    = [0],
                                ub    = [1],
                                names = [varName])

        for t in range(inp.nP):
            constrName = "capacity." + str(t)
            if yWW[t] > 1.0-_EPSI:
                cap = inp.m[item][t] + inp.a[item][t]*xWW[t]
            else:
                cap = 0.0
            cpxMaster.linear_constraints.set_coefficients(constrName,
            n_ilo[item][nrVars][nrSubVar], cap)
            #  n_ilo[item][nrVars][q], cap)

        constrName = "linking." + str(item) + "." + str(nrVars)
        cpxMaster.linear_constraints.set_coefficients(constrName,
        n_ilo[item][nrVars][nrSubVar], 1.0)
        #  n_ilo[item][nrVars][q], 1.0)

        self.cpxMaster = cpxMaster
        self.n_ilo    = n_ilo
        
    def get_columns(self, inp, item, u, alpha):
        cpxSub = self.cpxSub
        y_ilo  = self.y_ilo
        x_ilo  = self.x_ilo
        s_ilo  = self.s_ilo
        yWW    = self.yWW #  get the current WW solution

        #  add WW constraint to find "dominated" setup plans
        lDominated = []
        for t in range(inp.nP):
            constrName = "dominated." + str(t)
            lDominated.append(constrName)
            sub_contraint = cplex.SparsePair(ind=[y_ilo[t]], val=[1.0])
            cpxSub.linear_constraints.add(lin_expr = [sub_contraint],
                                       senses   = ["L"],
                                       rhs      = [yWW[t]],
                                       names    = [constrName])

        rhs = 1 - sum([yWW[t] for t in range(inp.nP)])
        index = [y_ilo[t] for t in range(inp.nP)]
        value = [(1-2*yWW[t]) for t in range(inp.nP)]
        diff_constraint = cplex.SparsePair(ind=index, val=value)
        cpxSub.linear_constraints.add(lin_expr = [diff_constraint],
                                   senses   = ["G"],
                                   rhs      = [rhs],
                                   names    = ["different"])

        #  first, we solve the ULSP, and then we call the populate solution
        #  self.cpxSub.parameters.mip.tolerances.absmipgap.set()
        #  self.cpxSub.parameters.mip.tolerances.mipgap.set(0.05)
        #  cpxSub.solve()
        cpxSub.parameters.mip.pool.intensity.set(4)
        cpxSub.parameters.mip.pool.capacity.set(nSolInPool)
        #  cpxSub.parameters.mip.pool.relgap.set(0.9)
        cpxSub.parameters.mip.limits.populate.set(nSolInPool)
        cpxSub.populate_solution_pool()

        # collecting pool solutions (one of them should be WW)
        names = cpxSub.solution.pool.get_names()
        nPool = len(names)
        #  print("POOL size = ", nPool)
        #  for n in names:
        #      print("z(",n,") = ", cpxSub.solution.pool.get_objective_value(n))
            
        q     = 1
        for n in names:
            #  rc = cpxSub.solution.pool.get_objective_value(n) - alpha[item]
            #  if rc <= -_EPSI:
                #  yAux      = [cpxSub.solution.get_values(y_ilo[t]) for t in range(inp.nP)]
            #  self.yWW  = [cpxSub.solution.get_values(y_ilo[t]) for t in range(inp.nP)]
            self.xWW  = [cpxSub.solution.pool.get_values(n, x_ilo[t]) for t in range(inp.nP)]
            self.sWW  = [cpxSub.solution.pool.get_values(n, s_ilo[t]) for t in range(inp.nP)]
                #  rc = self.compute_rc(inp, item, self.yWW, self.xWW, self.sWW, u, alpha)
                #  input("...aka")
                #  if rc <= -_EPSI:
            self.add_dominated_column_to_master(inp, item, q)
            q += 1

        cpxSub.linear_constraints.delete(lDominated)
        cpxSub.linear_constraints.delete("different")



    def dw_cycle(self, inp):
        
        u, alpha, theta = self.initialize_multipliers(inp)
        self.iterNr     = 0
        stopping        = False
        while not stopping and self.iterNr < maxIter:
            self.iterNr += 1
            stopping = self.column_generation(inp, u, alpha, theta)
            self.cpxMaster.solve()
            #  self.cpxMaster.write("master.lp")
            u, alpha, theta = self.get_master_solution(inp)
            nLambda = sum([len(self.l_ilo[j]) for j in range(inp.nI)])
            print("MASTER INFO :: z({0:3d}) = {1:15.2f} -- size :: lambda =\
            {2:5d}".format(self.iterNr, self.cpxMaster.solution.get_objective_value(), nLambda))
            
        
        # redefine master as MIP and solve it
        self.cpxMaster.set_problem_type(self.cpxMaster.problem_type.MILP)
        for j in range(inp.nI):
            for k in range(len(self.l_ilo[j])):
                self.cpxMaster.variables.set_types(self.l_ilo[j][k],
                self.cpxMaster.variables.type.binary)

        self.cpxMaster.write("MIP.lp")
        self.cpxMaster.parameters.mip.tolerances.absmipgap.set(_EPSI)
        self.cpxMaster.parameters.mip.tolerances.mipgap.set(_EPSI)
        self.cpxMaster.solve()
        print("MIP = ", self.cpxMaster.solution.get_objective_value())
        print("STATUS = ",
        self.cpxMaster.solution.status[self.cpxMaster.solution.get_status()])
        for j in range(inp.nI):
            lSol = self.cpxMaster.solution.get_values(self.l_ilo[j])
            print("lambda = ", lSol)
        
        

        


