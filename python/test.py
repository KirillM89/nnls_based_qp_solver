# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 16:34:01 2025

@author: Kirill Mitenkov
"""

import unittest
from ctypes import c_double, c_int
import nqp
import numpy as np

class Testing(unittest.TestCase):
 
    def test_python_demo(self):
        c = np.array([2.17299e6, -13162.7, 1.36145e6, 557966.0, 1.46127e6, 6.23872e6,\
                      24368.4, -24368.4, -14567.7, -9040.55, 46731.1, 36536.2], dtype = np.float64)
        nc = c.shape[0]
        H = np.zeros((nc, nc), dtype = np.float64)
        for i in range(nc):
            H[i][i] = 1.0
        
        A = np.array([[1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], \
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],  \
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0]], dtype = np.float64)
        b = np.array([1.0e-7, 1.0e-7, 1.0e-7], dtype = np.float64) 
        lw = np.array([-10.00553, -5.9893, -0.0107047, 0.367332, 0.0107047, 0.0107047, \
                       -0.989295, -10.0, -1.0e-8, -1.0e-8, -1.0e-8, 1.0e-8], dtype = np.float64)
        up = np.array([9998.4, 0.0107047, 9998.69, 9998.73, 0.189295, 0.189295, \
                       0.0107047, 10.0, 50.0, 50.0, 50.0, 50.0], dtype = np.float64)
        output = nqp.Solve(H, c, A, b, lw, up, 0)
        #out.dualExitStatus, out.cost, out.x, out.lambdaC, out.lambdaLw, out.lambdaUp
        print('status ',output[0])
        print('objective ', output[1])
        print('x ', output[2])
        print('lamC ', output[3])
        print('lamLw ', output[4])
        print('lamUp ', output[5])
        
        
if __name__ == '__main__':
    unittest.main()