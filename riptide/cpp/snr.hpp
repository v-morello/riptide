#ifndef SNR_HPP
#define SNR_HPP

#include <algorithm>
#include <cmath>

#include "block.hpp"
#include "kernels.hpp"


namespace riptide {

// Check that stdnoise > 0; throw std::invalid_argument if that is not the case
bool check_stdnoise(float stdnoise)
    {
    if (!(stdnoise > 0))
        throw std::invalid_argument("stdnoise must be > 0");
    }

// Check that all trial widths are > 0 and < bins; throw std::invalid_argument if that is not the case
bool check_trial_widths(const size_t* widths, size_t num_widths, size_t bins)
    {
    bool b = std::all_of(
        widths, 
        widths + num_widths, 
        [bins](size_t w){return (w > 0) & (w < bins);}
        );

    if (!b)
        throw std::invalid_argument("trial widths must be all > 0 and < columns");
    }

/*
Compute the S/N of a single pulse profile. Trial widths must be > 0 and < cols. stdnoise must be > 0.
output must have capacity for num_widths elements.
*/
void snr1(const float* __restrict__ arr, size_t size, const size_t* widths, size_t num_widths, float stdnoise, float* __restrict__ out)
    {
    const size_t wmax = *std::max_element(widths, widths + num_widths);
    float cpfsum[size + wmax];
    circular_prefix_sum(arr, size, size + wmax, cpfsum);
    const float sum = cpfsum[size - 1]; // sum of input array

    for (size_t iw = 0; iw < num_widths; ++iw)
        {
        const size_t w = widths[iw];
        // Height and baseline of a boxcar filter with width w bins
        // and zero mean and unit square sum
        const float h = sqrt((size - w) / float(size * w)); // boxcar height = +h
        const float b = w / float(size - w) * h; // boxcar baseline = -b
        const float dmax = diff_max(cpfsum + w, cpfsum, size);
        const float snr = ((h + b) * dmax - b * sum) / stdnoise;
        out[iw] = snr;
        }
    }


void snr2(ConstBlock block, const size_t* widths, size_t num_widths, float stdnoise, float* __restrict__ out)
    {
    for (size_t i = 0; i < block.rows; ++i)
        {
        snr1(block.rowptr(i), block.cols, widths, num_widths, stdnoise, out);
        out += num_widths;
        }
    }


} // namespace riptide


#endif