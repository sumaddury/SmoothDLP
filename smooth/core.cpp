// smooth/core.cpp

#include <pybind11/pybind11.h>
#include "gauss_dream.h"
#include "smooth_algos.h"
#include <pybind11/stl.h>
#include "dlp_infra.h"
#include <vector>
#include <utility>

// python setup.py build_ext --inplace
// python -m pip install -e . --no-build-isolation

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "Very basic core module for smooth—powered by pybind11";

    m.def("sieve_to", [](py::object o) {
        std::string s = py::str(o);
        mpz_class n_mp(s);
        std::vector<mpz_class> primes = sieveTo(n_mp);
        py::list out(primes.size());
        for (size_t i = 0; i < primes.size(); ++i) {
            out[i] = py::int_(py::str(primes[i].get_str()));
        }
        return out;
    }, py::arg("n"), "Generate a list of all primes ≤ n (accepts and returns Python big-ints).");

    m.def("is_prime", [](py::object o){
        std::string s = py::str(o);
        mpz_class n_mp(s);
        return isPrime(n_mp);
    }, py::arg("n"), "Return true if n is prime (handles arbitrarily large integers).");

    m.def("factorize", [](py::object o) {
        mpz_class n_mp{ std::string(py::str(o)) };
        std::vector<std::pair<mpz_class, uint32_t>> vec = factorize(n_mp);

        py::list out(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            auto &pr = vec[i];
            py::int_ p_py( py::str(pr.first.get_str()) );
            py::int_ e_py( pr.second );
            out[i] = py::make_tuple(p_py, e_py);
        }
        return out;
    }, py::arg("n"), "Fully factor n into a list of (prime, exponent) pairs");

    m.def("is_smooth", [](py::object o, uint32_t y){
        std::string s = py::str(o);
        mpz_class x_mp(s);
        return isSmooth(x_mp, y);
    }, py::arg("x"), py::arg("y"), "Return true if integer n is y-smooth (all prime factors ≤ y).");

    m.def("log_dickman", [](double u) {
        return logDickman(u);
    }, py::arg("u"), "Compute the natural logarithm of the Dickman-ρ function at u.\n\n");

    m.def("psi_approx", [](py::object o, uint64_t y) {
        std::string xs = py::str(o);
        mpz_class x_mp(xs);
        mpz_class z = psiApprox(x_mp, y);
        return py::int_(py::str(z.get_str()));
    }, py::arg("x"), py::arg("y"), "Approximate ψ(x, y) ≈ floor(x * ρ(ln(x)/ln(y))) with ≤0.1% error.\n\n");
}
