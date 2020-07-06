#ifndef TRANSFORMS_HPP
#define TRANSFORMS_HPP

#include <cstddef> // size_t
#include <cstring> // memcpy()

#include "kernels.hpp"
#include "block.hpp"


namespace riptide {

void merge(Block thead, Block ttail, Block out)
    {
    const size_t m = out.rows;
    const size_t p = out.cols;
    const float kh = (thead.rows - 1.0f) / (m - 1.0f);
    const float kt = (ttail.rows - 1.0f) / (m - 1.0f);

    for (size_t s = 0; s < m; ++s)
        {
        const size_t h = kh * s + 0.5f;
        const size_t t = kt * s + 0.5f;
        const size_t b = s - (h + t);
        fused_rollback_add(thead.rowptr(h), ttail.rowptr(t), p, h+b, out.rowptr(s));
        }
    }


void transform(Block input, Block temp, Block out)
    {
    const size_t m = input.rows;
    const size_t p = input.cols;

    if (m == 2)
        {
        add(input.rowptr(0), input.rowptr(1), p, out.rowptr(0));
        fused_rollback_add(input.rowptr(0), input.rowptr(1), p, 1, out.rowptr(1));
        return;
        }
    else if (m == 1)
        {
        memcpy(out.data, input.data, p * sizeof(float));
        return;
        }

    transform(input.head(), out.head(), temp.head());
    transform(input.tail(), out.tail(), temp.tail());
    merge(temp.head(), temp.tail(), out);
    }


void transform(Block input, float* temp, float* out)
    {
    transform(
        input, 
        Block(temp, input.rows, input.cols), 
        Block(out, input.rows, input.cols)
        );
    }


} // namespace riptide

#endif // TRANSFORMS_HPP