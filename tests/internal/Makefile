# Standalone builds for the internal (non-pybind) C++ test suites, one per
# header that has complete, headered functions: montgomery.h, gauss_dream.h,
# smooth_algos.h, infra.h.
#
# Most of these suites need only GMP (+ Boost for smooth_algos.h's Dickman-rho
# asymptotic tail) plus a C++17 compiler: infra.cpp's old LinBox-based sparse
# solve (linSolve/crtSolve/rank_relation_gf2) is still disabled (#if 0'd out in
# infra.cpp). The one exception is test_lin_alg, which exercises src/lin_alg.cpp
# and therefore does link LinBox/Givaro/NTL/FFLAS-FFPACK -- it has its own flags
# and library list below. None of it touches setup.py/pybind11, so everything
# here builds and runs independently of the Python extension.
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

INFRA_SRC := test_infra.cpp $(SRCDIR)/infra.cpp $(SRCDIR)/montgomery.cpp
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
LIN_ALG_SRC   := test_lin_alg.cpp $(SRCDIR)/lin_alg.cpp
LIN_ALG_BIN   := test_lin_alg
LIN_ALG_FLAGS := -std=c++17 -O2 -g -Wall -Wextra \
                 -I$(SRCDIR) -isystem $(CONDA_PREFIX)/include \
                 -Wno-invalid-specialization \
                 -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION
LIN_ALG_LIBS  := -llinbox -lgivaro -lntl -lfflas -lffpack -lgmpxx -lgmp

.PHONY: all test clean

all: $(MONT_BIN) $(INFRA_BIN) $(GAUSS_BIN) $(SMOOTH_BIN) $(LIN_ALG_BIN)

$(MONT_BIN): $(MONT_SRC) $(SRCDIR)/montgomery.h $(SRCDIR)/types.h
	$(CXX) $(CXXFLAGS) $(MONT_SRC) -o $(MONT_BIN) $(LDFLAGS) $(LIBS)

$(INFRA_BIN): $(INFRA_SRC) $(SRCDIR)/infra.h $(SRCDIR)/types.h $(SRCDIR)/montgomery.h
	$(CXX) $(CXXFLAGS) $(INFRA_SRC) -o $(INFRA_BIN) $(LDFLAGS) $(LIBS)

$(GAUSS_BIN): $(GAUSS_SRC) $(SRCDIR)/gauss_dream.h $(SRCDIR)/types.h $(SRCDIR)/montgomery.h
	$(CXX) $(CXXFLAGS) $(GAUSS_SRC) -o $(GAUSS_BIN) $(LDFLAGS) $(LIBS)

$(SMOOTH_BIN): $(SMOOTH_SRC) $(SRCDIR)/smooth_algos.h $(SRCDIR)/gauss_dream.h $(SRCDIR)/types.h $(SRCDIR)/montgomery.h
	$(CXX) $(CXXFLAGS) $(SMOOTH_SRC) -o $(SMOOTH_BIN) $(LDFLAGS) $(LIBS)

$(LIN_ALG_BIN): $(LIN_ALG_SRC) $(SRCDIR)/lin_alg.h $(SRCDIR)/types.h
	$(CXX) $(LIN_ALG_FLAGS) $(LIN_ALG_SRC) -o $(LIN_ALG_BIN) $(LDFLAGS) $(LIN_ALG_LIBS)

test: $(MONT_BIN) $(INFRA_BIN) $(GAUSS_BIN) $(SMOOTH_BIN) $(LIN_ALG_BIN)
	./$(MONT_BIN)
	./$(INFRA_BIN)
	./$(GAUSS_BIN)
	./$(SMOOTH_BIN)
	./$(LIN_ALG_BIN)

clean:
	rm -f $(MONT_BIN) $(INFRA_BIN) $(GAUSS_BIN) $(SMOOTH_BIN) $(LIN_ALG_BIN)
	rm -rf $(MONT_BIN).dSYM $(INFRA_BIN).dSYM $(GAUSS_BIN).dSYM $(SMOOTH_BIN).dSYM $(LIN_ALG_BIN).dSYM
