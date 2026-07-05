import glob
import os
import platform
import sys

import pybind11
from setuptools import Extension, setup

pybind_include = pybind11.get_include()

CONDA = os.environ.get("CONDA_PREFIX", sys.prefix)

if "CC" not in os.environ and "CXX" not in os.environ:
    clangxx_candidates = glob.glob(os.path.join(CONDA, "bin", "*-apple-darwin*-clang++"))
    if clangxx_candidates:
        clangxx = clangxx_candidates[0]
        os.environ["CXX"] = clangxx
        os.environ["CC"] = clangxx[: -len("clang++")] + "clang"

host_arch = platform.machine()
arch_flags = ["-arch", host_arch]

core_module = Extension(
    name="smooth._core",
    sources=["src/core.cpp", "src/gauss_dream.cpp", "src/smooth_algos.cpp", "src/infra.cpp"],
    include_dirs=[pybind_include, os.path.join(CONDA, "include")],
    library_dirs=[os.path.join(CONDA, "lib")],
    libraries=["gmpxx", "gmp", "ntl", "linbox", "givaro", "fflas", "ffpack"],
    language="c++",
    extra_compile_args=["-std=c++17", "-D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION",
                         "-Wno-invalid-specialization", "-Wno-template-body"] + arch_flags,
    extra_link_args=arch_flags + [f"-Wl,-rpath,{os.path.join(CONDA, 'lib')}"]
)

setup(
    packages=["smooth"],
    package_dir={"smooth": "src"},
    ext_modules=[core_module],
)
