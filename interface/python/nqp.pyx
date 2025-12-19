# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 16:15:49 2025

@author: Kirill Mitenkov
"""
import numpy as np
cimport nqp
"""
user-defined parameters
"""
LIN_SOLVER_TYPE = "linSolverType"
LOG_LEVEL = "logLevel"
LARGE_BONDS_PENALTY = "largeBoundsPenalty"
POS_DEF_CORRECTION = "posDefCorrection"
GAMMA_UPDATE = "gammaUpdate"

VALID_PARAMS = \
{
   LIN_SOLVER_TYPE : ("ldlt", "qr"),
   LOG_LEVEL  : (0, 1, 2, 3, 4, 5),
   LARGE_BONDS_PENALTY : (True, False),
   POS_DEF_CORRECTION : (True, False),
   GAMMA_UPDATE : (True, False)
}

def Solve(double[:,:]H, double[:]c, double[:,:]A = None, double[:]b = None, \
          double[:]lw = None, double[:]up = None, unsigned int nEq = 0, **params):
    actualParams = {}
    if params:
        for k, v in params.items():
            if k not in VALID_PARAMS:
                print("invalid key {k}")
                continue
            elif v not in VALID_PARAMS[k]:
                print("invalid value {v} for key {k}")
                continue
            elif k == LIN_SOLVER_TYPE:
                if v == "ldlt":
                    actualParams[k] = 0 
                elif v == "qr":
                    actualParams[k] = 1
                
    if H.shape[0] != H.shape[1]:
       print("ERROR: H must be symmetric {H.shape}")
       return None
    nVariables = H.shape[0]   
    if c.shape[0] != nVariables:
       print("ERROR: c size {s.shape[0]} must be equal to H size {nVariables}")
       return None
    nLwBnds = 0 if lw is None else lw.shape[0]
    nUpBnds = 0 if up is None else up.shape[0]
    if (nLwBnds > 0 and nUpBnds > 0 and (nLwBnds != nUpBnds or nUpBnds != nVariables)) or \
       (nLwBnds == 0 and nUpBnds > 0 and nUpBnds != nVariables) or \
       (nUpBnds == 0 and nLwBnds > 0 and nLwBnds != nVariables):    
        print("ERROR: Inconsistent lw up bounds sizes {nLwBnds} {nUpBnds}")
        return None
    nConstraints = 0 if A is None else A.shape[0] 
    bSize = 0 if b is None else b.shape[0]
    if bSize != nConstraints:
        print("ERROR: Inconsistent A and b sizes {nConstraints} {bSize}")
        return None
    if nEq > nConstraints:
        print("ERROR: number of equlaity constraints {nEq} must be <= total number of constraints {nConstraints}")
        return None 
    cdef int n_lw = nLwBnds
    cdef int n_up = nUpBnds
    cdef int n_x = nVariables
    cdef int n_c = nConstraints
    cdef Input problem
    cdef unsigned int ii, jj
    cdef Py_ssize_t i, j
    problem.A.resize(n_c)
    problem.b.resize(n_c)
    for i in range(nConstraints):
        problem.b[<unsigned int>i ] = b[i]
        problem.A[<unsigned int>i].resize(n_x)
        for j in range(nVariables):
            problem.A[<unsigned int>i][<unsigned int>j] = A[i][j]            
    problem.H.resize(n_x)
    problem.c.resize(n_x)
    problem.lw.resize(n_lw)
    problem.up.resize(n_up)
    for i in range(nVariables):
        problem.H[<unsigned int>i].resize(n_x)
        problem.c[<unsigned int>i] = c[i]
        if nLwBnds > 0:
            problem.lw[<unsigned int>i] = lw[i] 
        if nUpBnds > 0:
            problem.up[<unsigned int>i] = up[i] 
        for j in range(nVariables):
            problem.H[<unsigned int>i][<unsigned int>j] = H[i][j] 
    problem.nEqConstraints = nEq
    print("problem:", "nVariables", n_x, "nConstarints", n_c, "nLowerBounds", n_lw, "nUpperBounds", n_up)
    cdef Configuration configDefault
    for k, v in actualParams.items():
        if k == LIN_SOLVER_TYPE:
            configDefault.linSolverType = <unsigned int>v
        elif k == LOG_LEVEL:
            configDefault.logLevel = <unsigned int>v
        elif k == LARGE_BONDS_PENALTY:
            configDefault.largeBoundsPenalty = <bool>v
        elif k == POS_DEF_CORRECTION:
            configDefault.posDefCorrection = <bool>v
        elif k == GAMMA_UPDATE:
            configDefault.gammaUpdate = <bool>v
    cdef QPNNLS solver
    solver.Init(configDefault)
    if not solver.SetProblem(problem):
        return (100, None, None, None, None, None)
    solver.Solve()
    cdef Output out = solver.GetOutput()
    return (out.exitStatus, out.cost, np.array(out.x), np.array(out.lambdaC),
            np.array(out.lambdaLw), np.array(out.lambdaUp))
 