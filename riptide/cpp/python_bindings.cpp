#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cstring> // memset()
#include <chrono>

#include "kernels.hpp"
#include "block.hpp"
#include "transforms.hpp"


namespace py = pybind11;


py::array_t<float> fused_rollback_add(py::array_t<float> arr_x, py::array_t<float> arr_y, ptrdiff_t shift)
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


py::array_t<float> ffa2(py::array_t<float> arr_input)
{
    auto input = arr_input.mutable_unchecked<2>();
    const size_t rows = input.shape(0);
    const size_t cols = input.shape(1);

    std::unique_ptr<float[]> temp(new float[rows * cols]);
    auto output = py::array_t<float, py::array::c_style>({rows, cols});

    riptide::transform(
        riptide::Block(input.mutable_data(0, 0), rows, cols), 
        temp.get(), 
        output.mutable_data(0, 0)
        );

    return output;
}

/* Benchmark the ffa2() function. Returns the time per loop in seconds */
double benchmark_ffa2(size_t rows, size_t cols, size_t loops)
{
    const size_t size = rows * cols;
    std::unique_ptr<float[]> input(new float[size]);
    std::unique_ptr<float[]> temp(new float[size]);
    std::unique_ptr<float[]> out(new float[size]);
    memset(input.get(), 0, size * sizeof(float));

    auto block = riptide::Block(input.get(), rows, cols);
    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < loops; ++i)
        riptide::transform(block, temp.get(), out.get());

    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count() / loops;
}


PYBIND11_MODULE(libcpp, m)
{
    
    m.def(
        "fused_rollback_add", &fused_rollback_add, 
        "Add x with y rolled by shift elements to the LEFT, and store the result in z. shift must be positive. In numpy that would equivalent to: z = x + roll(y, -shift)"
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