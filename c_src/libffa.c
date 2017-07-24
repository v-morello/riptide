#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // memcpy()
#include <omp.h>

typedef struct {
    const float* data;
    size_t m;
    size_t b;
    } CBlock;


typedef struct {
    size_t h;
    size_t t;
    size_t b;
    size_t g;
    } CShifts;




// Get the number of profiles in the 'head' part of an input data with m profiles
inline size_t get_nh(size_t m)
    {return pow(2.0, floor(log2(m - 1.0)));}

// Get the circular shifts to apply during the recombine operation of the FFA transform
inline CShifts get_shifts(
        size_t m,  // number of profiles in the input
        size_t nh, // number of profiles in the head part of the input
        size_t g   // global top to bottom shift
        )
    {
    size_t nt = m - nh; // number of profiles in the tail part of the input
    size_t h = (nh - 1.0) / (m - 1.0) * g + 0.5; // head shift
    size_t t = (nt - 1.0) / (m - 1.0) * g + 0.5; // tail shift
    CShifts sh = {
        .h = h,
        .t = t,
        .b = g - h - t,
        .g = g
        };
    return sh;
    }

// C = A + B
inline void float_array_add(const float* restrict aa, const float* restrict bb, size_t size, float* restrict cc)
    {
    for (size_t jj=0; jj<size; ++jj)
        cc[jj] = aa[jj] + bb[jj];
    }

// Compute the S/N ratio of a profile for a number of boxcar width trials 
void get_snr(
    const float* in,        // input profile, b bins
    size_t b,               // number of bins in profile
    const long int* widths, // width trials, array of size nw
    size_t nw,              // number of width trials
    float varnoise,         // background noise variance
    float* out              // output buffer of size nw
    )
    {
    const size_t wmax = widths[nw - 1];
    const size_t csize = b + wmax; // size of the cumsum buffer

    // Compute profile cumulative sum
    float csum[csize];
    csum[0] = in[0];
    for (size_t ii=1; ii<b; ++ii)
        csum[ii] = csum[ii-1] + in[ii];
    for (size_t ii=b; ii<csize; ++ii)
        csum[ii] = csum[ii-1] + in[ii - b];

    for (size_t iw=0; iw<nw; ++iw)
        {
        const long int w = widths[iw];
        const float norm = pow(w * varnoise, -0.5);

        // Initialize best S/N with first phase trial at current width
        float best_snr = norm * (csum[w] - csum[0]);

        for (size_t ii=0; ii<b; ++ii)
            {
            float snr = norm * (csum[ii+w] - csum[ii]);
            if (snr > best_snr)
                best_snr = snr;
            }
        out[iw] = best_snr;
        }
    }


// Compute the S/N ratios of multiple profiles for a number of boxcar width trials 
void get_snr_2d(
    const float* in,        // 2D input array with m lines, each representing a profile with b bins
    size_t m,               // number of profiles
    size_t b,               // number of bins in each profile
    const long int* widths, // width trials, array of size nw
    size_t nw,              // number of width trials
    double varnoise,        // background noise variance. NOTE: has type double, for easier interfacing with python.
    size_t threads,         // number of OpenMP threads to use
    float* out              // 2D output buffer with m lines and nw cols
    )
    {
    omp_set_num_threads(threads);
    #pragma omp parallel for
    for (size_t ii=0; ii<m; ++ii)
        get_snr(in + ii * b, b, widths, nw, varnoise, out + ii * nw);
    }


// Perform the FFA recombine operation
// head and tail are FFA transforms here
void recombine(const CBlock head, const CBlock tail, float* out)
    {
    size_t gmax = head.m + tail.m;
	size_t b = head.b;

    for (size_t g=0; g<gmax; ++g)
        {
        CShifts sh = get_shifts(gmax, head.m, g);
        
        // Pointers to relevant lines of the inputs and output
        const float* lh = head.data + sh.h * b;
        const float* lt = tail.data + sh.t * b;
        float* lo = out + sh.g * b;
        
        // number of bins by which lt has to be left-shifted
        size_t nleft = (sh.h + sh.b) % b;
        size_t nright = b - nleft;

        // Do the operation: lo = lh + circ_lshift(lt, nleft)
        float_array_add(lh, lt + nleft, nright, lo);
        float_array_add(lh + nright, lt, nleft, lo + nright);
        } 
    }


// FFA Transform C function
void transform(
    const CBlock in,  // Block of input data to FFA transform
    float* buf,       // Pre-allocated buffer used to store intermediate FFA transforms of head and tail part of input data block
    float* out        // Where the FFA transform output is stored
    )
    {
    if (in.m <= 1)
        {
        memcpy(out, in.data, in.m * in.b * sizeof(float));
        return;
        }

    size_t b = in.b;
    size_t m = in.m;
    size_t nh = get_nh(m);
    size_t nt = m - nh;

    CBlock head = {.data = in.data, .m = nh, .b = b};
    CBlock tail = {.data = in.data + nh * b, .m = nt, .b = b};

    transform(head, out, buf);
    transform(tail, out + nh * b, buf + nh * b);

    CBlock trhead = {.data = buf, .m = nh, .b = b};
    CBlock trtail = {.data = buf + nh * b, .m = nt, .b = b};
    recombine(trhead, trtail, out);
    }



// FFA Transform function called from python
void py_transform(
    const float* in, // 2D input array, size m x b
    size_t m,        // number of profiles in input
    size_t b,        // number of bins in each profile
    float* out       // output of size m x b, pre-allocated in python
    )
    {
    CBlock block = {.data = in, .m = m, .b = b};
    float* buf = (float*)malloc(m * b * sizeof(float));
    transform(block, buf, out);
    free(buf); 
    }
