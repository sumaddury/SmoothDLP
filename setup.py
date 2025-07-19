# setup.py

import sys
import setuptools
from setuptools import Extension, setup
import pybind11
import os
import platform

pybind_include = pybind11.get_include()


# GMP_ROOT = "/opt/homebrew/opt/gmp"
# BOOST_ROOT = "/opt/homebrew/opt/boost"
CONDA = os.environ.get("CONDA_PREFIX", sys.prefix)

host_arch = platform.machine()
arch_flags = ["-arch", host_arch]

core_module = Extension(
    name="smooth._core",
    sources=["smooth/core.cpp", "smooth/gauss_dream.cpp", "smooth/smooth_algos.cpp", "smooth/infra.cpp", 
             "smooth/bench_helpers.cpp"],
    include_dirs=[pybind_include, os.path.join(CONDA, "include"), 
                #   os.path.join(GMP_ROOT, "include"), os.path.join(BOOST_ROOT, "include")
                  ],
    library_dirs=[os.path.join(CONDA, "lib")
        # os.path.join(GMP_ROOT, "lib"), os.path.join(BOOST_ROOT, "lib")
        ],
    libraries=["gmpxx", "gmp", "ntl", "linbox", "givaro", "fflas", "ffpack"],
    language="c++",
    extra_compile_args=["-std=c++17", "-D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION"] + arch_flags,
    extra_link_args=arch_flags
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
