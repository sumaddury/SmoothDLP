# Standalone build for the internal montgomery test suite.
#
# montgomery.{h,cpp} has no dependencies on the rest of this project (no
# LinBox/Givaro/FFLAS-FFPACK, no gauss_dream/infra), so this only needs GMP
# (already required by the whole project -- see ../../environment.yml) plus a
# C++17 compiler. It never touches setup.py/pybind11, so it builds and runs
# independently of the Python extension.
#
# Usage:
#   conda env create -f ../../environment.yml   # once, if not already done
#   conda activate smoothnumbers-dev
#   make -C tests/internal test                 # build + run
#
# Add a new `test_<module>.cpp` + matching `../../src/<module>.cpp` to SRC in
# a new target as more internal suites are added.

ifndef CONDA_PREFIX
$(error CONDA_PREFIX is not set -- run "conda activate smoothnumbers-dev" first)
endif

CXX ?= c++
SRCDIR := ../../src

CXXFLAGS := -std=c++17 -O2 -g -Wall -Wextra \
            -I$(SRCDIR) -I$(CONDA_PREFIX)/include
LDFLAGS  := -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
LIBS     := -lgmpxx -lgmp

MONT_SRC := test_montgomery.cpp $(SRCDIR)/montgomery.cpp
MONT_BIN := test_montgomery

.PHONY: all test clean

all: $(MONT_BIN)

$(MONT_BIN): $(MONT_SRC) $(SRCDIR)/montgomery.h $(SRCDIR)/types.h
	$(CXX) $(CXXFLAGS) $(MONT_SRC) -o $(MONT_BIN) $(LDFLAGS) $(LIBS)

test: $(MONT_BIN)
	./$(MONT_BIN)

clean:
	rm -f $(MONT_BIN)
