# C++ kernels dev notes

#### Downsampling by a real-valued factor

Let `N` be the number of input samples, `n` the number of output samples, `f` the downsampling factor and `kmax = n - 1` the index of the last output sample. Indices start at 0. With respect to the input data, the x coordinate of the rising edge of output sample index `k` is `kf` and that of its falling edge is `(k+1)f`. The falling edge of the last valid output sample must be located at a coordinate strictly lower than the falling edge of the last valid input sample, which is `N`. We must thus have:

```
        (kmax + 1) * f < N
i.e.             n * f < N
```

But we must be careful and note that `n = floor(N/f)` is **not** a sufficient condition ! If `f` divides `N`, then we have `floor(N/f) * f = N`. A sufficient condition is therefore:
```
n = floor((N - 1) / f)
```

Calculating output sample `k` requires reading input samples indices between `imin` and `imax` inclusive, given by:
```
imin = floor(k * f)
imax = floor((k+1) * f)
```
Using `floor()` in the `imax` formula is correct, because the input sample index `i` covers the input coordinate range `[i, i+1[`.

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