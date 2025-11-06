# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 15:36:13 2025

@author: Kirill Mitenkov
"""
from libcpp.vector cimport vector
ctypedef unsigned char bool
cdef extern from "types.h":
    cdef struct Configuration:
        Configuratuion()
        unsigned char linSolverType
        unsigned char logLevel
        bool largeBoundsPenalty
        bool posDefCorrection
        bool gammaUpdate
        bool rejectSingular
        double nnlsResidNormFsb
        double origPrimalFsb 
        
    cdef struct Input:
        vector[vector[double]] A
        vector[double] b
        vector[double] c
        vector[double] up
        vector[double] lw
        unsigned int nEqConstraints
    
    cdef struct Output:
        bool isPositiveDefinite
        unsigned char dualExitStatus
        unsigned char primalExitStatus
        unsigned int nIterations
        double maxBViolation
        double maxCViolation
        double cost
        vector[double] x
        vector[double] lambdaC
        vector[double] lambdaLw
        vector[double] lambdaUp
        
cdef extern from "decorator.h":
    cdef cppclass Core:
        pass
    cdef cppclass Callback:
        pass
    cdef cppclass QPNNLS:
        QPNNLS()
        void Init(const Configuration& config)
        bool SetProblem(const Input& problem)
        void SetCallback(Callback* callback = nullptr)
        void Solve()
        unsigned char GetInitStatus()
        const Output& GetOutput()  
        #bool isInitialized
        #Core* core
        #Output output
        
    



