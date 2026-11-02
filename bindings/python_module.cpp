// bindings/python_bindings.cpp
#include <pybind11/pybind11.h>
#include "models/Heston.hpp"
#include "core/payoff.hpp"

namespace py = pybind11;
using namespace quant::pricing;

PYBIND11_MODULE(pricing_engine, m) {
    m.doc() = "High-Performance C++ Derivatives Pricing Engine via Pybind11";

    // Bind the Payoff Abstract Base Class
    py::class_<core::IPayoff>(m, "Payoff")
        .def("__call__", &core::IPayoff::operator());

    // Bind Concrete Payoff Classes
    py::class_<core::CallPayoff, core::IPayoff>(m, "CallPayoff")
        .def(py::init<double>(), py::arg("strike"));

    py::class_<core::PutPayoff, core::IPayoff>(m, "PutPayoff")
        .def(py::init<double>(), py::arg("strike"));

    // Bind Heston Parameters Struct
    py::class_<models::HestonParams>(m, "HestonParams")
        .def(py::init<double, double, double, double, double>(),
             py::arg("v0"), py::arg("kappa"), py::arg("theta"), 
             py::arg("sigma"), py::arg("rho"))
        .def_readwrite("v0", &models::HestonParams::v0)
        .def_readwrite("kappa", &models::HestonParams::kappa)
        .def_readwrite("theta", &models::HestonParams::theta)
        .def_readwrite("sigma", &models::HestonParams::sigma)
        .def_readwrite("rho", &models::HestonParams::rho);

    // Bind Heston Pricer Engine
    py::class_<models::HestonPricer>(m, "HestonPricer")
        .def(py::init<const models::HestonParams&, double>(), 
             py::arg("params"), py::arg("risk_free_rate"))
        .def("price_european_mc", &models::HestonPricer::price_european_mc,
             "Monte Carlo Pricing for European Options under Heston Dynamics",
             py::arg("spot"), py::arg("time_to_maturity"), 
             py::arg("num_paths"), py::arg("num_time_steps"), py::arg("payoff"));
}