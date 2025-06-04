#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "signal_lib.h"

namespace py = pybind11;

PYBIND11_MODULE(signal_lib, m) {
    m.def("generate_sin", &generate_sin, "Generate sine wave");
    m.def("generate_cos", &generate_cos, "Generate cosine wave");
    m.def("generate_square", &generate_square, "Generate square wave");
    m.def("generate_sawtooth", &generate_sawtooth, "Generate sawtooth wave");
    m.def("plot_signal", &plot_signal, "Plot signal using matplot++");
}
