#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <algorithm>
#include <stdexcept>
#include <cstring> // memset()
#include <chrono>
#include <initializer_list>

#include "kernels.hpp"
#include "block.hpp"
#include "transforms.hpp"
#include "snr.hpp"
#include "downsample.hpp"
#include "periodogram.hpp"
#include "running_median.hpp"


namespace py = pybind11;


// Shorthand for allocating a new array to be returned as a numpy array to the
// python caller. Also a sneaky workaround the fact that the constructor of
// pybind11::array_t now expects an initializer list of ssize_t (signed),
// rather than size_t (unsigned). This avoids a narrowing conversion error.
// See: https://github.com/v-morello/riptide/issues/4
template <typename T>
py::array_t<T> new_cstyle_array(std::initializer_list<size_t> shape) {
    return py::array_t<T, py::array::c_style>(shape);
}


template<typename T>
void assert_c_contiguous(py::array_t<T> arr)
{
    const bool b = arr.flags() & py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_;
    if (!b)
    {
        throw std::invalid_argument("Input numpy array must be contiguous in memory");
    }
}


py::array_t<float> rollback(py::array_t<float> arr_x, size_t shift)
{
    assert_c_contiguous(arr_x);
    auto x = arr_x.unchecked<1>();
    const size_t size = x.size();
    auto arr_output = new_cstyle_array<float>({size});
    riptide::rollback(x.data(0), size, shift, arr_output.mutable_data(0));
    return arr_output;
}


py::array_t<float> fused_rollback_add(py::array_t<float> arr_x, py::array_t<float> arr_y, size_t shift)
{
    if (arr_x.size() != arr_y.size())
    {
        throw std::invalid_argument("Arrays must have the same number of elements");
    }
    assert_c_contiguous(arr_x);
    assert_c_contiguous(arr_y);
    
    auto x = arr_x.unchecked<1>();
    auto y = arr_y.unchecked<1>();
    const size_t size = x.size();
    auto arr_output = new_cstyle_array<float>({size});
    riptide::fused_rollback_add(x.data(0), y.data(0), size, shift, arr_output.mutable_data(0));
    return arr_output;
}


py::array_t<float> circular_prefix_sum(py::array_t<float> arr_x, size_t nsum)
{
    assert_c_contiguous(arr_x);
    auto x = arr_x.unchecked<1>();
    const size_t size = x.size();
    auto arr_output = new_cstyle_array<float>({size});
    riptide::circular_prefix_sum(x.data(0), size, nsum, arr_output.mutable_data(0));
    return arr_output;
}


py::array_t<float> ffa2(py::array_t<float> arr_input)
{
    assert_c_contiguous(arr_input);
    auto input = arr_input.unchecked<2>();
    const size_t rows = input.shape(0);
    const size_t cols = input.shape(1);

    std::unique_ptr<float[]> temp(new float[rows * cols]);
    auto output = new_cstyle_array<float>({rows, cols});

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


py::array_t<float> snr1(py::array_t<float> arr_x, py::array_t<size_t> arr_widths, float stdnoise)
{
    assert_c_contiguous(arr_x);
    assert_c_contiguous(arr_widths);
    auto x = arr_x.unchecked<1>();
    const size_t size = x.size();

    auto widths = arr_widths.unchecked<1>();
    const size_t num_widths = widths.size();

    riptide::check_stdnoise(stdnoise);
    riptide::check_trial_widths(widths.data(0), num_widths, size); 

    auto arr_output = new_cstyle_array<float>({num_widths});

    riptide::snr1(x.data(0), size, widths.data(0), num_widths, stdnoise, arr_output.mutable_data(0));
    return arr_output;
}


py::array_t<float> snr2(py::array_t<float> arr_x, py::array_t<size_t> arr_widths, float stdnoise)
{
    assert_c_contiguous(arr_x);
    assert_c_contiguous(arr_widths);
    auto x = arr_x.unchecked<2>();
    const size_t rows = x.shape(0);
    const size_t cols = x.shape(1);

    auto widths = arr_widths.unchecked<1>();
    const size_t num_widths = widths.size();

    riptide::check_stdnoise(stdnoise);
    riptide::check_trial_widths(widths.data(0), num_widths, cols); 

    auto arr_output = new_cstyle_array<float>({rows, num_widths});
    auto block = riptide::ConstBlock(x.data(0, 0), rows, cols);

    riptide::snr2(block, widths.data(0), num_widths, stdnoise, arr_output.mutable_data(0, 0));
    return arr_output;
}


py::array_t<float> downsample(py::array_t<float> arr_x, double f)
{
    assert_c_contiguous(arr_x);
    auto x = arr_x.unchecked<1>();
    const size_t size = x.size();

    riptide::check_downsampling_factor(size, f);

    // Allocate output array
    const size_t outsize = riptide::downsampled_size(size, f);
    auto output = new_cstyle_array<float>({outsize});

    riptide::downsample(x.data(0), size, f, output.mutable_data(0));
    return output;
}


std::tuple< py::array_t<double>, py::array_t<uint32_t>, py::array_t<float> > periodogram(
    py::array_t<float> arr_data,
    double tsamp,
    py::array_t<size_t> arr_widths,
    double period_min,
    double period_max,
    size_t bins_min,
    size_t bins_max)
{
    assert_c_contiguous(arr_data);
    assert_c_contiguous(arr_widths);
    auto data = arr_data.unchecked<1>();
    size_t size = data.size();
    auto widths = arr_widths.unchecked<1>();
    size_t num_widths = widths.size();

    size_t length = riptide::periodogram_length(size, tsamp, period_min, period_max, bins_min, bins_max);

    auto periods = new_cstyle_array<double>({length});
    auto foldbins = new_cstyle_array<uint32_t>({length});
    auto snrs = new_cstyle_array<float>({length, num_widths});

    riptide::periodogram(
        data.data(0), size, tsamp, widths.data(0), num_widths, 
        period_min, period_max, bins_min, bins_max, 
        periods.mutable_data(0), foldbins.mutable_data(0), snrs.mutable_data(0)
        );

    return std::make_tuple(periods, foldbins, snrs);
}


py::array_t<float> running_median(py::array_t<float> arr_x, size_t width)
{
    assert_c_contiguous(arr_x);
    auto x = arr_x.unchecked<1>();
    const size_t size = x.size();

    auto output = new_cstyle_array<float>({size});

    riptide::running_median<float>(x.data(0), size, width, output.mutable_data(0));
    return output;
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

    m.def(
        "snr1", &snr1, py::arg("data"), py::arg("widths"), py::arg("stdnoise") = 1.0,
        "S/N of a single pulse profile for multiple boxcar filter widths"
    );

    m.def(
        "snr2", &snr2, py::arg("data"), py::arg("widths"), py::arg("stdnoise") = 1.0,
        "S/N of multiple pulse profiles for multiple boxcar filter widths. 'data' must be a 2D array with shape (num_profiles, num_bins)."
    );

    m.def(
        "downsample", &downsample, py::arg("data"), py::arg("factor"),
        "Downsample data by a real-valued factor"
    );

    m.def(
        "periodogram", &periodogram,
        py::arg("data"), py::arg("tsamp"), py::arg("widths"), py::arg("period_min"), py::arg("period_max"), py::arg("bins_min"), py::arg("bins_max"),
        "Compute the periodogram of a time series. Returns a 3-tuple of arrays: trial periods, number of phase bins, S/N"
    );

    m.def(
        "running_median", &running_median, py::arg("data"), py::arg("width"),
        "Calculate the running median of a 1D array with a median window of 'width' elements.\n"
        "The data must be contiguous in memory, and width must be an odd number smaller than the input length.\n"
        "Throws std::invalid argument if any of the above conditions are not met."
    );

} // PYBIND11_MODULE