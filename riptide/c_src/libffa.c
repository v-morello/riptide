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
    const double epsilon = 1e-7;
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
        
        // NOTE: in some rare case where nsamp_in is an integer
        // multiple of nsamp_out, end might exceed je = floor(end) by
        // an infinitesimal amount. in[je] is out of bounds in that case,
        // causing bogus values to be read, or possibly a segfault.
        // Requiring end > (je + epsilon) protects against that case.
        // TODO: rewrite the whole downsample() function more cleanly
        if (end > (je + epsilon))
            val += in[je] * (end - je);

        out[ii] = val;
        }
    }

void cumsum(
    const float* restrict in, 
    size_t size,
    size_t nwrap,
    float* csum)
    {
    // csum[k] = x0 + ... + xk
    csum[0] = in[0];
    for (size_t ii = 1; ii < size; ++ii)
        csum[ii] = csum[ii-1] + in[ii];
    for (size_t ii = size; ii < size + nwrap; ++ii)
        csum[ii] = csum[ii-1] + in[ii - size];
    }


float get_snr_singlewidth(
    const float* restrict csum,
    size_t size,
    size_t width,
    float stdnoise)   
    {
    const float n = size;
    const float w = width;

    // Height and baseline value of a boxcar filter with
    // width w bins, zero mean and unit square sum
    const float h = sqrt((n - w) / (n * w));
    const float b = -w / (n - w) * h;
    const float H = h - b; // height measured from baseline
    const float profile_sum = csum[size - 1];

    // diff at start phase 0, width w
    // csum[k] = x0 + ... + xk
    float best_diff = csum[width - 1];

    for (size_t ii = 0; ii < size - 1; ++ii)
        {
        float diff = csum[ii+width] - csum[ii];
        if (diff > best_diff)
            best_diff = diff;
        }
    float best_snr = (H * best_diff + b * profile_sum) / stdnoise;
    return best_snr;
    }


// Compute the S/N ratio of a profile for a number of boxcar width trials
// This version of the function convolves profiles with a boxcar that has zero mean
// and unit power, as it should be.
void get_snr(
    const float* restrict in,        // input profile
    size_t size,                     // number of bins in profile
    const long int* restrict widths, // width trials, array of size nw
    size_t nw,                       // number of width trials
    float stdnoise,                  // background noise stddev
    float* restrict out              // output buffer of size nw
    )
    {
    const size_t wmax = widths[nw - 1];
    float csum[size + wmax];
    cumsum(in, size, wmax, &csum[0]);

    for (size_t iw=0; iw < nw; ++iw)
        out[iw] = get_snr_singlewidth(&csum[0], size, widths[iw], stdnoise);
    }


// Compute the S/N ratios of multiple profiles for a number of boxcar width trials 
void get_snr_2d(
    const float* restrict in,        // 2D input array with m lines, each representing a profile with b bins
    size_t m,                        // number of profiles
    size_t b,                        // number of bins in each profile
    const long int* restrict widths, // width trials, array of size nw
    size_t nw,                       // number of width trials
    double stdnoise,                 // background noise stddev. NOTE: has type double, for easier interfacing with python.
    float* restrict out              // 2D output buffer with m lines and nw cols
    )
    {
    for (size_t ii=0; ii<m; ++ii)
        get_snr(in + ii * b, b, widths, nw, stdnoise, out + ii * nw);
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

        // If ds is non-integer, then the noise variance in the downsampled data
        // is scaled by (ds - 1/3) rather than ds
        double ds_fracpart = ds - floor(ds);
        double varscale = ds;

        // The correction is applied only if ds "far enough" from integer number
        // i.e. if the downsampling window shifts by more than one new sample over
        // the input data
        if (ds_fracpart * (double)ns > 1.0)
            varscale = ds - 1.0/3.0;
        
        downsample(tseries, nsamp, ds, in);

        // Search the range of period numbers specified by current plan step
        for (size_t b=bmin[istep]; b<bmax[istep]; ++b)
            {
            // FFA transform
            size_t m = ns / b;
            transform(in, m, b, buf, out);

            // S/N evaluation
            float stdnoise = sqrt(m * varscale);
            get_snr_2d(out, m, b, widths, nw, stdnoise, snr);

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
