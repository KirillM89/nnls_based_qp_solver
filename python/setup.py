# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 15:35:29 2025

@author: Kirill Mitenkov
"""

import sys
import os
from setuptools import setup, find_packages, Extension
from shutil import copyfile, rmtree
from Cython.Distutils import build_ext 
from Cython.Build import cythonize

sources = ['callback.cpp', 'core.cpp', 'decorator.cpp', 'linSolvers.cpp', \
           'log.cpp', 'timers.cpp', 'utils.cpp']
priv_headers = ['core.h', 'linSolvers.h', 'log.h', 'timers.h', 'utils.h', 'core_types.h']
pub_headers = ['types.h', 'decorator.h', 'callback.h']
cwd = os.getcwd()
src_rel_path = os.path.join('..', 'src')
pub_rel_path = os.path.join('..', 'include')
try:
    src_path = os.path.join(cwd, src_rel_path)
    hdr_path = os.path.join(cwd, pub_rel_path)
except:
    src_path = []
    hdr_path = []
cp_src_path = os.path.join(cwd, 'src_copy')
if os.path.exists(cp_src_path):
    rmtree(cp_src_path)
if (src_path and os.path.exists(src_path)) and (hdr_path and os.path.exists(hdr_path)):
    os.mkdir(cp_src_path)
    for source in sources:
        copyfile(os.path.join(src_path, source), os.path.join(cp_src_path, source)) 
    for hdr in priv_headers:
        copyfile(os.path.join(src_path, hdr), os.path.join(cp_src_path, hdr))
    for hdr in pub_headers:
        copyfile(os.path.join(hdr_path, hdr), os.path.join(cp_src_path, hdr))    
else:   
    print("Could not find source directory")
    sys.exit(1)
sources_full_path = ['nqp.pyx']
for source in sources:
    sources_full_path.append(os.path.join(cp_src_path, source))
cython_ext = Extension('nqp',
        sources = sources_full_path,  
        language='c++',
        extra_compile_args=['-O3', '-DPROFILING', '-Wall'],
        include_dirs=[cp_src_path])

setup(name='nqp',
      version='0.1',
      description='NQP solver',
      ext_modules=cythonize(cython_ext),
      cmdclass={'build_ext': build_ext})

if os.path.exists(src_path):
    rmtree(cp_src_path)