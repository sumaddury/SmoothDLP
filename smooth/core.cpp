// smooth/core.cpp

#include <pybind11/pybind11.h>
#include "gauss_dream.h"
#include <pybind11/stl.h>

// python setup.py build_ext --inplace
// python -m pip install -e . --no-build-isolation

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "Very basic core module for smooth—powered by pybind11";

    m.def("sieve_to", [](py::object o) {
        std::string s = py::str(o);
        mpz_class n_mp(s);
        auto primes = sieveTo(n_mp);
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

    // m.def("factorize_naive", [](py::object o) {
    //     mpz_class n_mp{ std::string(py::str(o)) };
    //     auto result = factorize_naive(n_mp);
    //     auto &vec = result.first;
    //     auto  rem = result.second;

    //     py::list out(vec.size());
    //     for (size_t i = 0; i < vec.size(); ++i) {
    //         auto &pr = vec[i];
    //         py::int_  p_py( py::str(pr.first.get_str()) );
    //         py::int_  e_py(pr.second);
    //         out[i] = py::make_tuple(p_py, e_py);
    //     }
    //     py::int_ rem_py = py::int_( py::str(rem.get_str()) );
    //     return py::make_tuple(out, rem_py);
    // }, py::arg("n"), "Naïvely factor n into [(p,e)...], returns (list of (prime, exponent), remainder)");

    // m.def("squfof", [](py::object o) {
    //     mpz_class n_mp{ std::string(py::str(o)) };
    //     mpz_class factor = squfof(n_mp);
    //     return py::int_( py::str(factor.get_str()) );
    //   }, py::arg("n"), "Return a nontrivial factor of n using single‐instance SQUFOF");

    m.def("factorize", [](py::object o) {
            mpz_class n_mp{ std::string(py::str(o)) };
            auto vec = factorize(n_mp);

            py::list out(vec.size());
            for (size_t i = 0; i < vec.size(); ++i) {
                auto &pr = vec[i];
                py::int_ p_py( py::str(pr.first.get_str()) );
                py::int_ e_py( pr.second );
                out[i] = py::make_tuple(p_py, e_py);
            }
            return out;
        }, py::arg("n"), "Fully factor n into a list of (prime, exponent) pairs");
}
