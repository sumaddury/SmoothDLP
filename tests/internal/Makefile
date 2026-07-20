# Standalone builds for the internal (non-pybind) C++ test suites, one per
# header that has complete, headered functions: montgomery.h, gauss_dream.h,
# smooth_algos.h, infra.h.
#
# test_montgomery/test_gauss_dream/test_smooth_algos need only GMP (+ Boost
# for smooth_algos.h's Dickman-rho asymptotic tail) plus a C++17 compiler.
# test_infra and test_lin_alg both link LinBox/Givaro/NTL/FFLAS-FFPACK --
# infra.cpp's crtSolve calls into lin_alg (Hensel lift + CRT recombination
# over p-1's factors) -- so they share LIN_ALG_FLAGS/LIN_ALG_LIBS below and
# are both unity builds for the same reason (see LIN_ALG_SRC's comment).
# None of it touches setup.py/pybind11, so everything here builds and runs
# independently of the Python extension.
#
# Usage:
#   conda env create -f ../../environment.yml   # once, if not already done
#   conda activate smoothnumbers-dev
#   make -C tests/internal test                 # build + run all suites
#
# Add a new `test_<module>.cpp` + matching `../../src/<module>.cpp` source
# list as more internal suites are added.

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

# ProblemParams (in infra.cpp) computes its own smoothness bound and factor
# base, so it calls into smooth_algos.cpp (mp_ln/logDickman) and
# gauss_dream.cpp (sieveTo/factorize) -- both must link here too. test_infra.cpp
# itself unity-includes infra.cpp (which in turn unity-includes lin_alg.cpp --
# see infra.cpp's own top-of-file comment, and test_infra.cpp's), so it needs
# LIN_ALG_FLAGS/LIN_ALG_LIBS below, not the plain CXXFLAGS/LIBS the LinBox-free
# suites use; infra.cpp/lin_alg.cpp are listed as prerequisites only, same as
# test_lin_alg's rule, since they're pulled in via #include rather than
# compiled as their own translation units here.
INFRA_SRC := test_infra.cpp $(SRCDIR)/smooth_algos.cpp $(SRCDIR)/gauss_dream.cpp $(SRCDIR)/montgomery.cpp
INFRA_BIN := test_infra

GAUSS_SRC := test_gauss_dream.cpp $(SRCDIR)/gauss_dream.cpp $(SRCDIR)/montgomery.cpp
GAUSS_BIN := test_gauss_dream

SMOOTH_SRC := test_smooth_algos.cpp $(SRCDIR)/smooth_algos.cpp $(SRCDIR)/gauss_dream.cpp $(SRCDIR)/montgomery.cpp
SMOOTH_BIN := test_smooth_algos

# test_lin_alg is the only suite that links LinBox. Its extra flags are not
# optional:
#   -isystem for the conda includes, so LinBox/Givaro's own warnings don't
#     drown out -Wall -Wextra on our code
#   -Wno-invalid-specialization, because newer libc++ marks std::make_signed
#     as non-specializable and recint (Givaro's fixed-width big-int backend)
#     specializes it anyway
#   -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION, for LinBox's use of
#     the removed std::unary_function/binary_function
# setup.py passes the same two flags for the same reasons.
#
# This is a unity build -- test_lin_alg.cpp #includes lin_alg.cpp directly
# rather than compiling it as a second translation unit linked in. LinBox's
# commentator/args-parser machinery is pulled in transitively by even the
# base matrix/vector headers lin_alg.h needs, and defines several symbols
# out-of-line without `inline` in this LinBox build -- so two TUs that both
# include lin_alg.h and get linked together hit ~60 duplicate-symbol errors
# at link time. A single TU sidesteps that entirely. lin_alg.cpp is still
# listed as a prerequisite below so `make` rebuilds when it changes, even
# though it's not passed to the compiler as its own source argument.
LIN_ALG_SRC   := test_lin_alg.cpp
LIN_ALG_BIN   := test_lin_alg
LIN_ALG_FLAGS := -std=c++17 -O2 -g -Wall -Wextra \
                 -I$(SRCDIR) -isystem $(CONDA_PREFIX)/include \
                 -Wno-invalid-specialization \
                 -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION
LIN_ALG_LIBS  := -llinbox -lgivaro -lntl -lfflas -lffpack -lgmpxx -lgmp

# Unlike test_infra.cpp, dlp.cpp does NOT unity-include infra.cpp (it only
# includes infra.h -- dlp.cpp calls infra::* functions but never touches
# lin_alg:: directly, so it has no need for lin_alg.h's LinBox headers
# itself). That means infra.cpp needs to be compiled here as its own real
# source file (not just a Makefile prerequisite the way it is for
# INFRA_BIN above), so DLP::solve's calls into infra::addRelations/
# crtSolve/filterDetermined/modInv actually have something to link against.
# Still needs LIN_ALG_FLAGS/LIN_ALG_LIBS throughout, since infra.cpp itself
# unity-includes lin_alg.cpp internally and pulls in LinBox regardless of
# how it's reached.
DLP_SRC := test_dlp.cpp $(SRCDIR)/dlp.cpp $(SRCDIR)/infra.cpp $(SRCDIR)/smooth_algos.cpp $(SRCDIR)/gauss_dream.cpp $(SRCDIR)/montgomery.cpp
DLP_BIN := test_dlp

.PHONY: all test clean

all: $(MONT_BIN) $(INFRA_BIN) $(GAUSS_BIN) $(SMOOTH_BIN) $(LIN_ALG_BIN) $(DLP_BIN)

$(MONT_BIN): $(MONT_SRC) $(SRCDIR)/montgomery.h $(SRCDIR)/types.h
	$(CXX) $(CXXFLAGS) $(MONT_SRC) -o $(MONT_BIN) $(LDFLAGS) $(LIBS)

$(INFRA_BIN): $(INFRA_SRC) $(SRCDIR)/infra.cpp $(SRCDIR)/lin_alg.cpp $(SRCDIR)/lin_alg.h $(SRCDIR)/infra.h $(SRCDIR)/types.h $(SRCDIR)/montgomery.h $(SRCDIR)/smooth_algos.h $(SRCDIR)/gauss_dream.h
	$(CXX) $(LIN_ALG_FLAGS) $(INFRA_SRC) -o $(INFRA_BIN) $(LDFLAGS) $(LIN_ALG_LIBS)

$(GAUSS_BIN): $(GAUSS_SRC) $(SRCDIR)/gauss_dream.h $(SRCDIR)/types.h $(SRCDIR)/montgomery.h
	$(CXX) $(CXXFLAGS) $(GAUSS_SRC) -o $(GAUSS_BIN) $(LDFLAGS) $(LIBS)

$(SMOOTH_BIN): $(SMOOTH_SRC) $(SRCDIR)/smooth_algos.h $(SRCDIR)/gauss_dream.h $(SRCDIR)/types.h $(SRCDIR)/montgomery.h
	$(CXX) $(CXXFLAGS) $(SMOOTH_SRC) -o $(SMOOTH_BIN) $(LDFLAGS) $(LIBS)

$(LIN_ALG_BIN): $(LIN_ALG_SRC) $(SRCDIR)/lin_alg.cpp $(SRCDIR)/lin_alg.h $(SRCDIR)/types.h
	$(CXX) $(LIN_ALG_FLAGS) $(LIN_ALG_SRC) -o $(LIN_ALG_BIN) $(LDFLAGS) $(LIN_ALG_LIBS)

$(DLP_BIN): $(DLP_SRC) $(SRCDIR)/dlp.h $(SRCDIR)/infra.h $(SRCDIR)/lin_alg.cpp $(SRCDIR)/lin_alg.h $(SRCDIR)/types.h $(SRCDIR)/montgomery.h $(SRCDIR)/smooth_algos.h $(SRCDIR)/gauss_dream.h
	$(CXX) $(LIN_ALG_FLAGS) $(DLP_SRC) -o $(DLP_BIN) $(LDFLAGS) $(LIN_ALG_LIBS)

test: $(MONT_BIN) $(INFRA_BIN) $(GAUSS_BIN) $(SMOOTH_BIN) $(LIN_ALG_BIN) $(DLP_BIN)
	./$(MONT_BIN)
	./$(INFRA_BIN)
	./$(GAUSS_BIN)
	./$(SMOOTH_BIN)
	./$(LIN_ALG_BIN)
	./$(DLP_BIN)

clean:
	rm -f $(MONT_BIN) $(INFRA_BIN) $(GAUSS_BIN) $(SMOOTH_BIN) $(LIN_ALG_BIN) $(DLP_BIN)
	rm -rf $(MONT_BIN).dSYM $(INFRA_BIN).dSYM $(GAUSS_BIN).dSYM $(SMOOTH_BIN).dSYM $(LIN_ALG_BIN).dSYM $(DLP_BIN).dSYM
