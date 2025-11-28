# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 16:15:49 2025

@author: m00829527
"""

import numpy as np
cimport nqp
def Solve(double[:,:]H, double[:]c, double[:,:]A = None, double[:]b = None, \
          double[:]lw = None, double[:]up = None, unsigned int nEq = 0):
    nn = H.shape
    if nn[0] != nn[1]:
       print("Invalid H")
       return None
    n = nn[0]   
    if c.shape[0] != n:
       print("Invalid c")
       return None
    if lw is None:
        nl = 0
    else:
        nl = lw.shape[0]
    if up is None:
        nup = 0
    else:
        nup = up.shape[0]   
    if nl > 0 and nup > 0 and nl != nup:
        print("Inconsistent lw up bounds")
        return None
    if  A is None:
        nc = 0
    else: 
        nc = A.shape[0]
    if b is None:
        nb = 0
    else:
        nb = b.shape[0]
    if nb != nc:
        print("Inconsistent A and b")
        return None
    if nEq > nc:
        print("Invalid nEq")
        return None 
    cdef int n_lw = nl
    cdef int n_up = nup
    cdef int n_x = n
    cdef int n_c = nc
    cdef Input problem
    problem.A.resize(n_c)
    cdef unsigned int ii, jj
    cdef Py_ssize_t i, j
    for i in range(nc):
        ii = <unsigned int>i
        problem.A[ii].resize(n_x)
        for j in range(n):
            jj = <unsigned int>j
            problem.A[ii][jj] = A[i][j]
            
    problem.H.resize(n_x)
    for i in range(n):
        ii = <unsigned int>i
        problem.H[ii].resize(n_x)
        for j in range(n):
            jj = <unsigned int>j
            problem.H[ii][jj] = H[i][j] 
            
    problem.c.resize(n_x)
    for i in range(n):
        ii = <unsigned int>i  
        problem.c[ii] = c[i]
        
    problem.b.resize(n_c)
    for i in range(nc):
        ii = <unsigned int>i  
        problem.b[ii] = b[i]  
        
    problem.lw.resize(n_lw)
    for i in range(nl):
        ii = <unsigned int>i  
        problem.lw[ii] = lw[i] 
    
    problem.up.resize(n_up)
    for i in range(nup):
        ii = <unsigned int>i  
        problem.up[ii] = up[i]       
    problem.nEqConstraints = nEq
    print("problem:", "nx",n,"nC",nc,"nLw",nl,"nUp",nup)
    cdef Configuration configDefault
    cdef QPNNLS solver
    solver.Init(configDefault)
    solver.SetProblem(problem)
    solver.Solve()
    cdef Output out = solver.GetOutput()
    return (out.dualExitStatus, out.cost, np.array(out.x), np.array(out.lambdaC),
            np.array(out.lambdaLw), np.array(out.lambdaUp))
 