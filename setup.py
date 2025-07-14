# setup.py

import sys
import setuptools
from setuptools import Extension, setup
import pybind11
import os
import platform

pybind_include = pybind11.get_include()

GMP_ROOT = "/opt/homebrew/opt/gmp"
BOOST_ROOT = "/opt/homebrew/opt/boost"

host_arch = platform.machine()
arch_flags = ["-arch", host_arch]

core_module = Extension(
    name="smooth._core",
    sources=["smooth/core.cpp", "smooth/gauss_dream.cpp", "smooth/smooth_algos.cpp", "smooth/dlp_infra.cpp"],
    include_dirs=[pybind_include, os.path.join(GMP_ROOT, "include"), os.path.join(BOOST_ROOT, "include")],
    library_dirs=[os.path.join(GMP_ROOT, "lib"), os.path.join(BOOST_ROOT, "lib")],
    libraries=["gmpxx", "gmp", "boost_thread", "boost_system"],
    language="c++",
    extra_compile_args=["-std=c++17"] + arch_flags,
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
