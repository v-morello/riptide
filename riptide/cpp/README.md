# C++ kernels dev notes

#### Downsampling by a real-valued factor

Let `N` be the number of input samples, `n` the number of output samples, `f` the downsampling factor and `kmax = n - 1` the index of the last output sample. Indices start at 0. With respect to the input data, the x coordinate of the rising edge of output sample index `k` is `kf` and that of its falling edge is `(k+1)f`. Calculating output sample `k` requires reading input samples indices between `imin` and `imax` inclusive, given by:
```
imin = floor(k * f)
imax = floor((k+1) * f)
```
Using `floor()` in the `imax` formula is correct, because the input sample index `i` covers the input coordinate range `[i, i+1[`.  

We can calculate `n` as `n = floor(N / f)` but there is a **caveat**: when `f` divides `N`, for the last output sample `k + 1 = n` and we thus have `imax = N`. **In the code we need to explicitly enforce `imax < N`.**


#### Largest FFA shift corresponding to a period < p + 1 

From the formula in Morello 2020:
```
smax = floor(p * (m - 1) / (p + 1))
```
which is always stricly less than `m - 1`. It is thus not worth evaluating for S/N the folded profiles whose corresponding shift is strictly larger than `smax`.
Equivalently, the *ceiling shift* `ceilshift` is the lowest shift that corresponds to a trial period larger or equal to `p + 1`:
```
ceilshift = ceil(p * (m - 1) / (p + 1))
```


#### Largest FFA shift corresponding to a period < Pmax

```
smax = floor(p * (m - 1) * (1 - p * tau / Pmax))
```

Where `tau` is the sampling interval. For `Pmax = p + 1` and `tau = 1`, this yields the previous formula. Likewise:

```
ceilshift = ceil(p * (m - 1) * (1 - p * tau / Pmax))
```

### Height and baseline of a boxcar filter with zero mean and unit square sum

Let `n` be the total number of bins of the boxcar, and `w` its width expressed in number of bins. Its height `h` and baseline value `b` are given by:
```
h = sqrt((n - w) / (n * w))
b = - w / (n - w) * h
```