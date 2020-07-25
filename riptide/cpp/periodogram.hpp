#ifndef PERIODOGRAM_HPP
#define PERIODOGRAM_HPP

#include <cstddef> // size_t
#include <cstdint>
#include <cmath>
#include <memory>

#include "downsample.hpp"
#include "transforms.hpp"
#include "block.hpp"
#include "snr.hpp"


namespace riptide {


void periodogram_check_arg(bool condition, std::string const& errmsg)
    {
    if (!(condition))
        throw std::invalid_argument(errmsg);
    }


void periodogram_check_arguments(
    size_t size, 
    double tsamp,
    double period_min, 
    double period_max, 
    size_t bins_min, 
    size_t bins_max)
    {    
    periodogram_check_arg(tsamp > 0, "tsamp must be > 0");
    periodogram_check_arg(period_min > 0, "period_min must be > 0");
    periodogram_check_arg(period_max > period_min, "period_max must be > period_min");
    periodogram_check_arg(bins_min > 1, "bins_min must be > 1");
    periodogram_check_arg(bins_max >= bins_min, "bins_max must be >= bins_min");
    periodogram_check_arg(period_min >= tsamp * bins_min, "Must have: period_min >= tsamp * bins_min ");
    // NOTE: we don't check period_max, the search will automatically stop when the maximum allowable trial period is reached
    }


/*
Returns the first shift in an FFA transform that corresponds to a trial period
equal to, or greater than pmax. pmax must be expressed in units of the 
sampling interval.

This function is useful to calculate how many rows of an FFA transform should
be evaluated for S/N, since often we wish to only consider the rows that
correspond to period trials smaller than p + 1. The index of the last row
to evaluate is equal to ceilshift - 1, or equivalently, the total number
of rows to evaluate is equal to ceilshift.
*/
size_t ceilshift(size_t rows, size_t cols, double pmax)
    {
    return ceil(cols * (rows - 1.0) * (1.0 - cols / pmax));
    }


/*
Returns the total number of trial periods in a periodogram
*/
size_t periodogram_length(
    size_t size,
    double tsamp,
    double period_min,
    double period_max,
    size_t bins_min,
    size_t bins_max)
    {
    periodogram_check_arguments(size, tsamp, period_min, period_max, bins_min, bins_max);

    // Initial downsampling factor
    // We want: ds_ini * tsamp * bmin = period_min
    double ds_ini = period_min / (tsamp * bins_min);

    // Geometric growth factor for the downsampling factor
    double ds_geo = (bins_max + 1.0) / bins_min;

    // Number of required downsampling cycles
    size_t num_downsamplings = ceil(log(period_max / period_min) / log(ds_geo));
    size_t length = 0; // total number of period trials, to be calculated

    /* Downsampling loop */
    for (size_t ids = 0; ids < num_downsamplings; ++ids)
        {
        const double f = ds_ini * pow(ds_geo, ids); // current downsampling factor
        const double tau = f * tsamp; // current sampling time
        const double period_max_samples = period_max / tau;
        const size_t n = downsampled_size(size, f); // current number of input samples

        // Min and max number of bins with which to FFA transform in order to
        // cover all trial periods between period_min and period_max.
        // NOTE: bstop is INclusive
        // Also, we MUST enforce bstop <= n, to avoid doing an FFA transform with 0 rows
        const size_t bstart = bins_min;
        const size_t bstop = std::min({ bins_max, n, size_t(period_max_samples) });

        /* FFA transform loop */
        for (size_t bins = bstart; bins <= bstop; ++bins)
            {
            const size_t rows = n / bins;
            const double period_ceil = std::min(period_max_samples, bins + 1.0);
            const size_t rows_eval = std::min(rows, ceilshift(rows, bins, period_ceil));
            length += rows_eval;
            }
        }
    return length;
    }


/*
Compute the periodogram of a time series that has been normalised to zero mean and unit variance.
Outputs are: trial periods (num_periods elements), number of phase bins used in the fold (num_periods elements), 
and signal to noise ratio (num_periods * num_widths elements)
*/
void periodogram(
    const float* __restrict__ data,
    size_t size,
    double tsamp,
    const size_t* __restrict__ widths,
    size_t num_widths,
    double period_min,
    double period_max,
    size_t bins_min,
    size_t bins_max,
    double* __restrict__ periods,
    uint32_t* __restrict__ foldbins,
    float* __restrict__ snr)
    {
    periodogram_check_arguments(size, tsamp, period_min, period_max, bins_min, bins_max);

    // Initial downsampling factor
    // We want: ds_ini * tsamp * bmin = period_min
    double ds_ini = period_min / (tsamp * bins_min);

    // Geometric growth factor for the downsampling factor
    double ds_geo = (bins_max + 1.0) / bins_min;

    // Number of required downsampling cycles
    size_t num_downsamplings = ceil(log(period_max / period_min) / log(ds_geo));

    // Allocate buffers
    const size_t bufsize = downsampled_size(size, ds_ini);
    std::unique_ptr<float[]> input_mem(new float[bufsize]);
    std::unique_ptr<float[]> ffabuf_mem(new float[bufsize]);
    std::unique_ptr<float[]> ffaout_mem(new float[bufsize]);
    float* input = input_mem.get();
    float* ffabuf = ffabuf_mem.get();
    float* ffaout = ffaout_mem.get();

    /* Downsampling loop */
    for (size_t ids = 0; ids < num_downsamplings; ++ids)
        {
        const double f = ds_ini * pow(ds_geo, ids); // current downsampling factor
        const double tau = f * tsamp; // current sampling time
        const double period_max_samples = period_max / tau;
        const size_t n = downsampled_size(size, f); // current number of input samples

        downsample(data, size, f, input);

        // Min and max number of bins with which to FFA transform in order to
        // cover all trial periods between period_min and period_max.
        // NOTE: bstop is INclusive
        // Also, we MUST enforce bstop <= n, to avoid doing an FFA transform with 0 rows
        const size_t bstart = bins_min;
        const size_t bstop = std::min({ bins_max, n, size_t(period_max_samples) });

        /* FFA transform loop */
        for (size_t bins = bstart; bins <= bstop; ++bins)
            {
            const size_t rows = n / bins;
            const float stdnoise = sqrt(rows * downsampled_variance(size, f));
            const double period_ceil = std::min(period_max_samples, bins + 1.0);
            const size_t rows_eval = std::min(rows, ceilshift(rows, bins, period_ceil));

            transform(input, rows, bins, ffabuf, ffaout);
            
            auto block = ConstBlock(ffaout, rows_eval, bins);
            snr2(block, widths, num_widths, stdnoise, snr);

            for (size_t s = 0; s < rows_eval; ++s)
                {
                periods[s] = tau * bins * bins / (bins - s / (rows - 1.0));
                foldbins[s] = bins;
                }

            snr += rows_eval * num_widths;
            periods += rows_eval;
            foldbins += rows_eval;
            }
        }
    }


} // namespace riptide

#endif // PERIODOGRAM_HPP