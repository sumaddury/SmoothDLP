# setup.py

import sys
import setuptools
from setuptools import Extension, setup
import pybind11
import os

pybind_include = pybind11.get_include()

GMP_ROOT = "/opt/homebrew/opt/gmp"
BOOST_ROOT = "/opt/homebrew/opt/boost"


core_module = Extension(
    name="smooth._core",
    sources=["smooth/core.cpp", "smooth/gauss_dream.cpp", "smooth/smooth_algos.cpp"],
    include_dirs=[pybind_include, os.path.join(GMP_ROOT, "include")],
    library_dirs=[os.path.join(GMP_ROOT, "lib")],
    libraries=["gmpxx", "gmp"],
    language="c++",
    extra_compile_args=["-std=c++17"],
)

setup(
    name="smooth",
    version="0.0.1",
    author="Your Name",
    description="A minimal smooth-numbers frontend with C++ backend",
    packages=["smooth"],
    ext_modules=[core_module],
    zip_safe=False,
)
