#ifndef ARRAYOPS_H
#define ARRAYOPS_H

#include <string.h> // memcpy(), size_t
#include <stdint.h>

// C = A + B
void array_add(
    const float* restrict A,
    const float* restrict B,
    const size_t size,
    float* restrict C)
{
    // Using a 32-bit int as a loop index is a tiny bit faster
    const uint32_t imax = size;
    for (uint32_t i = 0; i < imax; ++i)
        C[i] = A[i] + B[i];
}

// C = A + circular_left_shift(B, shift)
// shift is a POSITIVE integer number
void array_add_clshift(
    const float* restrict A,
    const float* restrict B,
    const size_t size,
    const size_t shift,
    float* restrict C)
{
    const size_t p = shift % size;
    const size_t q = size - p;
    array_add(A, B + p, q, C);
    array_add(A + q, B, p, C + q);
}

#endif // ARRAYOPS_H