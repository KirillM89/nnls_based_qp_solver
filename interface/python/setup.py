# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 15:35:29 2025

@author: Kirill Mitenkov
"""

import sys
import os
from setuptools import setup, find_packages, Extension
from shutil import copyfile, copytree, rmtree
from Cython.Distutils import build_ext 
from Cython.Build import cythonize
cwd = os.getcwd()
src_copy = os.path.join(cwd, 'source')
if not os.path.exists(src_copy):
    root = os.path.join(cwd, '..', '..')
    src_orig = os.path.join(root, 'src')
    pub_orig = os.path.join(root, 'include')
    if not os.path.exists(src_orig) or not os.path.exists(pub_orig):
        print("ERROR: Directory with code didn't found")
        sys.exit(1)
    else:
        copytree(src_orig, os.path.join(src_copy, 'src'))
        copytree(pub_orig, os.path.join(src_copy, 'include'))
        copyfile(os.path.join(root, 'CMakeLists.txt'), os.path.join(src_copy, 'CMakeLists.txt'))

src_files = [os.path.join(src_copy, 'src', fl) for fl in os.listdir(os.path.join(src_copy, 'src')) \
             if os.path.splitext(os.path.join(src_copy, 'src', fl))[1] == '.cpp']  + ['nqp.pyx'] 
print(src_files)
try:         
    os.remove(os.path.join(src_copy, 'src', 'CMakeLists.txt'))
except:
    pass
cython_ext = Extension('nqp',
        sources = src_files,  
        language='c++',
        extra_compile_args=['-O3', '-DPROFILING', '-Wall'],
        include_dirs=['source/include', 'source/src'],
        library_dirs=['.'] )
    

setup(name='nqp',
      version='0.1',
      description='NQP solver',
      ext_modules=cythonize(cython_ext),
      cmdclass={'build_ext': build_ext},
      zip_safe=False,
      )

if os.path.exists(src_copy):
    rmtree(src_copy)