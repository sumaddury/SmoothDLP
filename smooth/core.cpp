// smooth/core.cpp

#include <pybind11/pybind11.h>
#include "gauss_dream.h"
#include <pybind11/stl.h>


namespace py = pybind11;

// A trivial "add" function in C++:
int add(int a, int b) {
    return a + b;
}

// Optional: a function that prints something to stdout.
// We'll return a Python string so Python can print it too.
std::string say_hello(const std::string &name) {
    return "Hello from C++, " + name + "!";
}

// The PYBIND11_MODULE macro: name the extension "_core"
PYBIND11_MODULE(_core, m) {
    m.doc() = "Very basic core module for smooth—powered by pybind11";

    // Expose `add` to Python: call it smooth._core.add(...)
    m.def("add", &add, "Add two integers");

    // Expose `say_hello` to Python: call it smooth._core.say_hello(...)
    m.def("say_hello", &say_hello, "Greet someone from C++");

    m.def("sieve_to", &sieveTo, "Generate a vector of all primes ≤ n");

    m.def("is_prime", [](py::object o){
        std::string s = py::str(o);
        mpz_class n_mp(s);
        return isPrime(n_mp);
    }, py::arg("n"), "Return true if n is prime (handles arbitrarily large integers).");

}
