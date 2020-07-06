#ifndef KERNELS_HPP
#define KERNELS_HPP

#include <cstddef> // size_t

namespace riptide {

/* Add x with y, and store the result in z */
void add(const float* __restrict__ x, const float* __restrict__ y, size_t size, float* __restrict__ z)
    {
    for (size_t i = 0; i < size; ++i)
        z[i] = x[i] + y[i];
    }

/* 
Add x with y rolled by shift elements to the LEFT, and store the result in z. shift must be positive.
In numpy that would equivalent to: z = x + roll(y, -shift)
*/
void fused_rollback_add(const float* __restrict__ x, const float* __restrict__ y, size_t size, size_t shift, float* __restrict__ z)
    {
    const size_t p = shift % size;
    const size_t q = size - p;
    add(x, y + p, q, z);
    add(x + q, y, p, z + q);
    }


/*
Rotate x backwards by shift elements and stored the result in out. shift must be positive.
In numpy that would equivalent to: out = roll(x, -shift)
*/
void rollback(const float* __restrict__ x, size_t size, size_t shift, float* __restrict__ out)
    {
    const size_t p = shift % size;
    const size_t q = size - p;
    memcpy(out, x + p, q * sizeof(float));
    memcpy(out + q, x, p * sizeof(float));
    }


} // namespace riptide

#endif // KERNELS_HPP