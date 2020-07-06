#ifndef KERNELS_HPP
#define KERNELS_HPP

#include <cstddef> // size_t

namespace riptide {

/* Add x with y, and store the result in z */
void add(const float *x, const float *y, size_t size, float *z)
    {
    for (size_t i = 0; i < size; ++i)
        z[i] = x[i] + y[i];
    }

/* 
Add x with y rolled by shift elements to the LEFT, and store the result in z. shift must be positive.
In numpy that would equivalent to: z = x + roll(y, -shift)
*/
void fused_rollback_add(const float *x, const float *y, size_t size, size_t shift, float *z)
    {
    const size_t p = shift % size;
    const size_t q = size - p;
    add(x, y + p, q, z);
    add(x + q, y, p, z + q);
    }

} // namespace riptide

#endif // KERNELS_HPP