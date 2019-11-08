#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // memcpy()

#include "arrayops.h"
#include "kernels.h"

// Downsample a time series by any real-valued factor
void downsample(
    const float* restrict in,
    size_t size,
    double factor,
    float* out
    )
    {
    size_t outsize = floor(size / factor);
    
    // ii = output sample index
    for (size_t ii=0; ii<outsize; ++ii)
        {
        const double start = factor * ii;
        const double end = start + factor;
        const size_t js = ceil(start);
        const size_t je = floor(end);
        
        float val = 0.0;

        // jj = input sample index
        for (size_t jj=js; jj<je; ++jj)
            val += in[jj];

        if (start < js)
            val += in[js-1] * (js - start);
        if (end > je)
            val += in[je] * (end - je);

        out[ii] = val;
        }
    }


// Compute the S/N ratio of a profile for a number of boxcar width trials
// This version of the function convolves profiles with a boxcar that has zero mean
// and unit power, as it should be.
void get_snr(
    const float* restrict in,        // input profile
    size_t size,                     // number of bins in profile
    const long int* restrict widths, // width trials, array of size nw
    size_t nw,                       // number of width trials
    float varnoise,                  // background noise variance
    float* restrict out              // output buffer of size nw
    )
    {
    const size_t wmax = widths[nw - 1];

    // Compute profile cumulative sum
    float csum[size + wmax];
    csum[0] = in[0];
    for (size_t ii = 1; ii < size; ++ii)
        csum[ii] = csum[ii-1] + in[ii];
    for (size_t ii = size; ii < size + wmax; ++ii)
        csum[ii] = csum[ii-1] + in[ii - size];

    const float stdnoise = sqrt(varnoise);
    const float profile_sum = csum[size - 1];
    const float nf = size;

    for (size_t iw=0; iw < nw; ++iw)
        {
        const long int w = widths[iw];

        // Height and baseline value of a boxcar filter with
        // width w bins, zero mean and unit square sum
        const float h = sqrt((nf - w) / (nf * w));
        const float b = -w / (nf - w) * h;
        const float H = h - b; // height measured from baseline

        // diff at start phase 0, width w
        // csum[k] = x0 + ... + xk
        float best_diff = csum[w - 1];
        const float* A = &csum[0];
        const float* B = &csum[w];

        for (size_t ii = 0; ii < size - 1; ++ii)
            {
            float diff = B[ii] - A[ii];
            if (diff > best_diff)
                best_diff = diff;
            }
        out[iw] = (H * best_diff + b * profile_sum) / stdnoise;
        }
    }


// Compute the S/N ratios of multiple profiles for a number of boxcar width trials 
void get_snr_2d(
    const float* restrict in,        // 2D input array with m lines, each representing a profile with b bins
    size_t m,                        // number of profiles
    size_t b,                        // number of bins in each profile
    const long int* restrict widths, // width trials, array of size nw
    size_t nw,                       // number of width trials
    double varnoise,                 // background noise variance. NOTE: has type double, for easier interfacing with python.
    float* restrict out              // 2D output buffer with m lines and nw cols
    )
    {
    for (size_t ii=0; ii<m; ++ii)
        get_snr(in + ii * b, b, widths, nw, varnoise, out + ii * nw);
    }


// FFA Transform function called from python
void py_transform(
    const float* restrict in, // 2D input array, size m x b
    size_t m,                 // number of profiles in input
    size_t b,                 // number of bins in each profile
    float* restrict out       // output of size m x b, pre-allocated in python
    )
    {
    float* buf = (float*)malloc(m * b * sizeof(float));
    transform(in, m, b, buf, out);
    free(buf); 
    }


// Periodogram function called from python
void py_periodogram(
    const float* restrict tseries,      // input time series, normalized to zero mean and unit variance
    size_t nsamp,                       // number of samples in input time series
    const double* dsfactor,             // sequence of downsampling factor, size = number of plan steps
    const long int* restrict bmin,      // sequence of min. number of bins, size = number of plan steps
    const long int* restrict bmax,      // sequence of max. number of bins, size = number of plan steps
    size_t nsteps,                      // number of plan steps
    const long int* restrict widths,    // sequence of pulse width trials, in number of bins
    size_t nw,                          // number of width trials
    float* restrict periods,            // period trials (output)
    float* restrict snr                 // S/N (output)
    )
    {
    // Allocate buffer for downsampled time series
    float* in = (float*)malloc(nsamp * sizeof(float)); 

    // Allocate FFA buffers
    float* out = (float*)malloc(nsamp * sizeof(float));
    float* buf = (float*)malloc(nsamp * sizeof(float));

    // MAIN LOOP
    for (size_t istep=0; istep<nsteps; ++istep)
        {
        // Downsample
        double ds = dsfactor[istep];         // current downsampling factor
        size_t ns = floor(nsamp / ds);       // number of samples after downsampling
        
        downsample(tseries, nsamp, ds, in);

        // Search the range of period numbers specified by current plan step
        for (size_t b=bmin[istep]; b<bmax[istep]; ++b)
            {
            // FFA transform
            size_t m = ns / b;
            transform(in, m, b, buf, out);

            // S/N evaluation
            float varnoise = m * ds;
            get_snr_2d(out, m, b, widths, nw, varnoise, snr);

            // Fill period trials
            for (size_t sh=0; sh<m; ++sh)
                periods[sh] = ds * b + ds * sh * b / (m * b - sh);

            // move output pointers forward
            snr = snr + m * nw;
            periods = periods + m;
            } // end current plan step
        
        } // END MAIN LOOP

    // Free all temporary arrays
    free(in);
    free(out);
    free(buf);
    }
