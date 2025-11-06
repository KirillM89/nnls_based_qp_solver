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
    """Testing class for daqp python interface."""

    def test_python_demo(self):
        nqp.Solve()
        
if __name__ == '__main__':
    unittest.main()