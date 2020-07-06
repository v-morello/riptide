#ifndef BLOCK_HPP
#define BLOCK_HPP

#include <cstddef> // size_t
#include "kernels.hpp"


namespace riptide {

/* Represents a contiguous two-dimensional array */
class Block {

public:
    float* const data;
    const size_t rows;
    const size_t cols;

public:
    Block(float* ptr, size_t rows, size_t cols)
        : data(ptr), rows(rows), cols(cols)
        {}

    size_t head_size() const 
        {return rows >> 1;}

    float* rowptr(size_t irow) const 
        {return data + irow * cols;}

    Block head() const 
        {return Block(data, head_size(), cols);}

    Block tail() const 
        {
        const size_t h = head_size();
        return Block(rowptr(h), rows - h, cols);
        }

}; // class Block

} // namespace riptide

#endif // BLOCK_HPP