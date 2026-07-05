#include <pybind11/pybind11.h>
#include "gauss_dream.h"
#include "smooth_algos.h"
#include <pybind11/stl.h>
#include "infra.h"
#include "types.h"
#include <vector>
#include <utility>
#include <string>

namespace py = pybind11;

namespace {

u128 pyint_to_u128(const py::object& o) {
    py::bytes raw = o.attr("to_bytes")(16, "little", py::arg("signed") = false);
    std::string buf = py::cast<std::string>(raw);
    u128 result = 0;
    for (int i = 15; i >= 0; --i) {
        result = (result << 8) | static_cast<unsigned char>(buf[static_cast<size_t>(i)]);
    }
    return result;
}

py::int_ u128_to_pyint(u128 v) {
    char buf[16];
    for (int i = 0; i < 16; ++i) {
        buf[i] = static_cast<char>(static_cast<unsigned char>(v & 0xFF));
        v >>= 8;
    }
    py::bytes raw(buf, 16);
    static py::object int_type = py::module_::import("builtins").attr("int");
    return int_type.attr("from_bytes")(raw, "little").cast<py::int_>();
}

U128Vector pylist_to_u128vector(const py::iterable& o) {
    U128Vector v;
    for (auto item : o) v.push_back(pyint_to_u128(py::reinterpret_borrow<py::object>(item)));
    return v;
}

py::list u128vector_to_pylist(const U128Vector& v) {
    py::list out(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[i] = u128_to_pyint(v[i]);
    return out;
}

FactorList pylist_to_factorlist(const py::iterable& o) {
    FactorList v;
    for (auto item : o) {
        py::tuple t = py::reinterpret_borrow<py::tuple>(item);
        v.push_back({pyint_to_u128(t[0]), t[1].cast<uint32_t>()});
    }
    return v;
}

mpz_class pyint_to_mpz(const py::object& o) {
    return mpz_class(std::string(py::str(o)));
}

py::int_ mpz_to_pyint(const mpz_class& m) {
    return py::int_(py::str(m.get_str()));
}

MpzVector pylist_to_mpzvector(const py::iterable& o) {
    MpzVector v;
    for (auto item : o) v.push_back(pyint_to_mpz(py::reinterpret_borrow<py::object>(item)));
    return v;
}

py::list mpzvector_to_pylist(const MpzVector& v) {
    py::list out(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[i] = mpz_to_pyint(v[i]);
    return out;
}

std::vector<MpzVector> pylist_to_levels(const py::iterable& o) {
    std::vector<MpzVector> levels;
    for (auto item : o) levels.push_back(pylist_to_mpzvector(py::reinterpret_borrow<py::iterable>(item)));
    return levels;
}

py::list levels_to_pylist(const std::vector<MpzVector>& levels) {
    py::list out(levels.size());
    for (size_t i = 0; i < levels.size(); ++i) out[i] = mpzvector_to_pylist(levels[i]);
    return out;
}

} // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "";

    m.def("sieve_to", [](uint64_t n) {
        std::vector<uint32_t> primes = sieveTo(n);
        py::list out(primes.size());
        for (size_t i = 0; i < primes.size(); ++i) out[i] = py::int_(primes[i]);
        return out;
    }, py::arg("n"));

    m.def("is_prime", [](py::object o) {
        return isPrime(pyint_to_u128(o));
    }, py::arg("n"));

    m.def("factorize", [](py::object o) {
        auto vec = factorize(pyint_to_u128(o));
        py::list out(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            out[i] = py::make_tuple(u128_to_pyint(vec[i].first), py::int_(vec[i].second));
        }
        return out;
    }, py::arg("n"));

    m.def("factorize_naive", [](py::object o) {
        auto result = factorize_naive(pyint_to_u128(o));
        py::list out(result.first.size());
        for (size_t i = 0; i < result.first.size(); ++i) {
            out[i] = py::make_tuple(u128_to_pyint(result.first[i].first), py::int_(result.first[i].second));
        }
        return py::make_tuple(out, u128_to_pyint(result.second));
    }, py::arg("n"));

    m.def("squfof", [](py::object o) {
        return u128_to_pyint(squfof(pyint_to_u128(o)));
    }, py::arg("n"));

    m.def("is_smooth", [](py::object o, uint32_t y) {
        return isSmooth(pyint_to_u128(o), y);
    }, py::arg("x"), py::arg("y"));

    m.def("log_dickman", [](double u) {
        return logDickman(u);
    }, py::arg("u"));

    m.def("mp_ln", [](py::object o) {
        return mp_ln(pyint_to_u128(o));
    }, py::arg("x"));

    m.def("log_mul", [](py::object o, double log_rho) {
        return u128_to_pyint(log_mul(pyint_to_u128(o), log_rho));
    }, py::arg("x"), py::arg("log_rho"));

    m.def("psi_approx", [](py::object o, uint64_t y) {
        return u128_to_pyint(psiApprox(pyint_to_u128(o), y));
    }, py::arg("x"), py::arg("y"));

    m.def("build_product_tree", [](py::iterable level) {
        return levels_to_pylist(buildProductTree(pylist_to_mpzvector(level)));
    }, py::arg("level"));

    m.def("smooth_candidates", [](py::iterable p_levels, py::iterable X) {
        auto idx = smoothCandidates(pylist_to_levels(p_levels), pylist_to_u128vector(X));
        py::list out(idx.size());
        for (size_t i = 0; i < idx.size(); ++i) out[i] = py::int_(idx[i]);
        return out;
    }, py::arg("p_levels"), py::arg("X"));

    m.def("tree_factorize", [](py::iterable p_levels, py::object d) {
        SparseList result = treeFactorize(pylist_to_levels(p_levels), pyint_to_u128(d));
        py::list out(result.size());
        for (size_t i = 0; i < result.size(); ++i)
            out[i] = py::make_tuple(py::int_(result[i].first), py::int_(result[i].second));
        return out;
    }, py::arg("p_levels"), py::arg("d"));

}
