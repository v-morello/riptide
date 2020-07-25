#ifndef DOWNSAMPLE_HPP
#define DOWNSAMPLE_HPP

#include <cstddef> // size_t
#include <cmath>
#include <stdexcept>

namespace riptide {

// Check that 1 < f <= size; throw std::invalid_argument if that is not the case
void check_downsampling_factor(size_t size, double f)
    {
    bool valid = (f > 1.0) & (f <= size);
    if (!valid)
        throw std::invalid_argument("Downsampling factor must verify: 1 < f <= size");
    }

/*
Number of samples in a time series after it has been downsampled by a real-valued factor f
*/
size_t downsampled_size(size_t num_samples, double f)
    {
    return floor(num_samples / f);
    }

/*
Variance of background Gaussian noise after downsampling a time series with 'num_samples' elements by a real-valued factor f.
*/
double downsampled_variance(size_t num_samples, double f)
    {
    const double k = floor(f);
    const double r = f - k;
    const double x = downsampled_size(num_samples, f) * r;
    if (x > 1)
        return f - 1.0 / 3.0;
    else
        return pow(k-1, 2) + 2.0/3.0 * pow(x, 2) - x + 1.0;
    }


/* 
Downsample input array by a real-valued factor. Output must have capacity for floor((size - 1.0) / f) elements
*/
void downsample(const float* __restrict__ in, size_t size, double f, float* __restrict__ out)
    {
    check_downsampling_factor(size, f);

    // N = number of input samples
    // n = number of valid output samples
    const size_t N = size;
    const size_t n = downsampled_size(N, f);

    // k = output sample index
    for (size_t k = 0; k < n; ++k)
        {
        // minimum and maximum x-coordinate (real-valued)
        // of the input data range required to compute output sample k
        const double start = k * f;
        const double end = start + f;

        // minimum (inclusive) and maximum (inclusive)
        // input sample indices that must be read to compute output sample k
        const size_t imin = floor(start);
        // NOTE: floor() is OK in imax calculation, 
        // because the input sample index imax covers the input x-coordinate range [imax, imax+1]
        // NOTE 2: when f divides size, for the last output sample we have end = N,
        // hence the need to enforce imax < N
        const size_t imax = std::min(floor(end), N - 1.0);

        // Weights to be applied to samples imin and imax
        // Other input samples are weighted by 1
        const float wmin = (imin + 1) - start; // ceil(k*f) - k*f
        const float wmax = end - imax;         // (k+1)*f - floor((k+1)*f)

        float acc = wmin * in[imin];
        for (size_t i = imin + 1; i < imax; ++i)
            acc += in[i];
        acc += wmax * in[imax];

        out[k] = acc;
        }
    }

} // namespace riptide

#endif // DOWNSAMPLE_HPP
