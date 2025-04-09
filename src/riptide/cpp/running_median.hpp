#ifndef RUNNING_MEDIAN_HPP
#define RUNNING_MEDIAN_HPP

#include <cstring> // size_t
#include <vector>
#include <stdexcept>
#include <cstdlib>


namespace riptide {


template <typename T>
size_t partition(T* __restrict__ data, size_t start, size_t end, size_t ipiv)
    {
    const T pivot = data[ipiv];

    // Temporarily place pivot at the end
    std::swap(data[ipiv], data[end]);

    // Place all elements strictly smaller than pivot value in front
    size_t j = start;
    for (size_t i = start; i < end; ++i)
        {
        if (data[i] < pivot)
            {
            std::swap(data[i], data[j]);
            ++j;
            }
        }
    
    // Place pivot in j-th position (where it belongs)
    std::swap(data[j], data[end]);
    return j;
    }


// Find the n-th element in data between indices start and end (inclusive)
template <typename T>
float nth_element(T* __restrict__ data, size_t start, size_t end, size_t n)
    {
    if (start == end)
        return data[start];
    
    // NOTE: Using end as pivot index is fastest
    const size_t j = partition(data, start, end, end);

    if (j == n)
        return data[j];

    else if (j < n)
        return nth_element(data, j + 1, end, n);
    
    else // j > n
        return nth_element(data, start, j - 1, n);
    }


template <typename T>
class RunningMedian {

private:
    const size_t width;

    // Ring buffer
    std::vector<T> ring;
    size_t writepos; 

    // Temporary buffer for locating median
    std::vector<T> buffer;

    void ringbuf_push(T val)
        {
        ring[writepos] = val;
        ++writepos;
        if (writepos >= width)
            writepos = 0;
        }

public:
    RunningMedian(size_t width)
        : width(width), ring(width), writepos(0), buffer(width)
        {}
    
    float push(float val)
        {
        // insert new value into ring buffer
        ringbuf_push(val);

        // copy ring buffer to quickselect buffer
        std::copy(ring.begin(), ring.end(), buffer.begin());

        // return new median
        return nth_element<T>(buffer.data(), 0, width - 1, width / 2);
        }

}; // RunningMedian


template <typename T>
void running_median(const T* __restrict__ data, size_t size, size_t width, T* __restrict__ out)
    {
    if (!(width % 2))
        throw std::invalid_argument("width must be an odd number");

    if (!(width < size))
        throw std::invalid_argument("width must be < size");

    const size_t half = width / 2;
    auto rmed = RunningMedian<T>(width);

    // Initialize by pushing:
    // * 'half' copies of data[0]
    // * data[0]
    // * data[1:half]  i.e. 1 <= i < half
    for (size_t i = 0; i < half + 1; ++i)
        rmed.push(data[0]);

    for (size_t i = 1; i < half; ++i)
        rmed.push(data[i]);

    // From this point we can start pushing elements from data[half]
    // When we push data[i] we collect the running median with the window
    // centered on data[i - half]
    for (size_t i = half; i < size; ++i)
        out[i - half] = rmed.push(data[i]);

    // And now we deal with the last 'half' elements by pushing copies
    // of the last element: data[size - 1]
    for (size_t i = 0; i < half; ++i)
        out[size - half + i] = rmed.push(data[size - 1]);
    }


} // namespace riptide


#endif // RUNNING_MEDIAN_HPP