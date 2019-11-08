#ifndef KERNELS_HPP
#define KERNELS_HPP

#include <string.h> // memcpy(), size_t
#include "arrayops.h"

// FFA Transform input with 2 lines and p columns
void transform_m2(
    const float* restrict in,
    const size_t p,
    float* restrict out)
{
    array_add(in, in + p, p, out);
    array_add_clshift(in, in + p, p, 1, out + p);
}


void merge(
    const float* restrict th,
    const float* restrict tt,
    const size_t hm, // lines in head
    const size_t tm, // lines in tail
    const size_t p,
    float* restrict out)
{
    const size_t m = hm + tm;
    const float kh = (hm - 1.0) / (m - 1.0);
    const float kt = (tm - 1.0) / (m - 1.0);

    for (size_t g = 0; g < m; ++g)
    {
        size_t h = kh * g + 0.5; // head shift
        size_t t = kt * g + 0.5; // tail shift
        size_t b = g - h - t;    // between shift

        const float* lh = th + h * p;
        const float* lt = tt + t * p;
        float* lo = out + g * p;

        array_add_clshift(lh, lt, p, h + b, lo);
    }
}


void transform(
    const float* restrict in,
    const size_t m,
    const size_t p,
    float* restrict buf,
    float* restrict out)
{
    if (m == 2)
    {
        transform_m2(in, p, out);
        return;
    }
    else if (m == 1)
    {
        memcpy(out, in, p * sizeof(float));
        return;
    }

    const size_t hm = m >> 1; // lines in head
    const size_t tm = m - hm; // lines in tail
    const size_t ofs = hm * p;

    transform(in, hm, p, out, buf);
    transform(in + ofs, tm, p, out + ofs, buf + ofs);
    merge(buf, buf + ofs, hm, tm, p, out);
}


#endif // KERNELS_HPP
