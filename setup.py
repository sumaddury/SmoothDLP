# setup.py

import sys
import setuptools
from setuptools import Extension, setup
import pybind11
import os

# The build_ext command will know where to find pybind11 headers:
pybind_include = pybind11.get_include()

GMP_ROOT = "/opt/homebrew/opt/gmp"
BOOST_ROOT = "/opt/homebrew/opt/boost"


core_module = Extension(
    name="smooth._core",                 # the full Python package name
    sources=["smooth/core.cpp", "smooth/gauss_dream.cpp"],         # our single C++ source file
    include_dirs=[pybind_include, os.path.join(GMP_ROOT, "include")],       # so the compiler can find pybind11
    library_dirs=[os.path.join(GMP_ROOT, "lib")],
    libraries=["gmpxx", "gmp"],
    language="c++",
    extra_compile_args=["-std=c++17"],   # compile as C++17
)

setup(
    name="smooth",
    version="0.0.1",
    author="Your Name",
    description="A minimal smooth-numbers frontend with C++ backend",
    packages=["smooth"],                 # tells setuptools to include `smooth/`
    ext_modules=[core_module],
    zip_safe=False,
)
