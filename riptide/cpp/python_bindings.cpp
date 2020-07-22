#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cstring> // memset()
#include <chrono>

#include "kernels.hpp"
#include "block.hpp"
#include "transforms.hpp"


namespace py = pybind11;


py::array_t<float> rollback(py::array_t<float> arr_x, size_t shift)
{
    // Array proxies
    auto x = arr_x.unchecked<1>();
    const size_t size = x.size();

    auto arr_output = py::array_t<float, py::array::c_style>({size});
    riptide::rollback(x.data(0), size, shift, arr_output.mutable_data(0));
    return arr_output;
}


py::array_t<float> fused_rollback_add(py::array_t<float> arr_x, py::array_t<float> arr_y, size_t shift)
{
    if (arr_x.size() != arr_x.size())
        throw std::invalid_argument("Arrays must have the same number of elements");
    
    // Array proxies
    auto x = arr_x.unchecked<1>();
    auto y = arr_y.unchecked<1>();
    const size_t size = x.size();

    auto arr_output = py::array_t<float, py::array::c_style>({size});
    riptide::fused_rollback_add(x.data(0), y.data(0), size, shift, arr_output.mutable_data(0));
    return arr_output;
}


py::array_t<float> circular_prefix_sum(py::array_t<float> arr_x, size_t nsum)
{
    // Array proxies
    auto x = arr_x.unchecked<1>();
    const size_t size = x.size();
    auto arr_output = py::array_t<float, py::array::c_style>({nsum});
    riptide::circular_prefix_sum(x.data(0), size, nsum, arr_output.mutable_data(0));
    return arr_output;
}


py::array_t<float> ffa2(py::array_t<float> arr_input)
{
    auto input = arr_input.unchecked<2>();
    const size_t rows = input.shape(0);
    const size_t cols = input.shape(1);

    std::unique_ptr<float[]> temp(new float[rows * cols]);
    auto output = py::array_t<float, py::array::c_style>({rows, cols});

    riptide::transform(input.data(0, 0), rows, cols, temp.get(), output.mutable_data(0, 0));
    return output;
}

/* Benchmark the ffa2() function. Returns the time per loop in seconds */
double benchmark_ffa2(size_t rows, size_t cols, size_t loops)
{
    const size_t size = rows * cols;

    // NOTE: performance slightly increases when all buffers are contiguous
    // (better memory locality)
    std::unique_ptr<float[]> buffer(new float[3 * size]);
    float* input = buffer.get();
    float* temp = input + size;
    float* out = temp + size;
    memset(input, 0, size * sizeof(float));

    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < loops; ++i)
        riptide::transform(input, rows, cols, temp, out);

    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count() / loops;
}


PYBIND11_MODULE(libcpp, m)
{
    m.def(
        "rollback", &rollback, 
        "Rotate input array backwards by shift elements. shift must be positive. In numpy that would be equivalent to out = roll(x, -shift)"
    );

    m.def(
        "fused_rollback_add", &fused_rollback_add, 
        "Add x with y rolled backwards by shift elements, and store the result in z. shift must be positive. In numpy that would equivalent to: z = x + roll(y, -shift)"
    );

    m.def(
        "circular_prefix_sum", &circular_prefix_sum, 
        "Compute the circular prefix sum of the input array over nsum elements"
    );

    m.def(
        "ffa2", &ffa2, 
        "FFA transform a 2D input array"
    );

    m.def(
        "benchmark_ffa2", &benchmark_ffa2, 
        "Benchmark the ffa2() function. Returns the time per loop in seconds."
    );

} // PYBIND11_MODULE