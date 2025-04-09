from setuptools import setup

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext

# WARNING: Enabling -freciprocal-math (which is also enabled by
# -funsafe-math-optimizations or -ffast-math) causes surprising inconsistencies
# when calculating the size of the output periodogram, and the size of the
# buffers to store data downsampled by a real-valued factor.
# In turn, this can cause segmentation faults (code attempts to write past the
# end of these arrays).
# The flags below provide the same speedups as -ffast-math, without the risks.
SAFE_FAST_MATH_FLAGS = [
    "-fassociative-math",
    "-fno-math-errno",
    "-ffinite-math-only",
    "-fno-rounding-math",
    "-fno-signed-zeros",
    "-fno-trapping-math",
]

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)
ext_modules = [
    Pybind11Extension(
        "riptide.libcpp",
        sorted(["src/riptide/cpp/python_bindings.cpp"]),
        extra_compile_args=["-O3", "-march=native"] + SAFE_FAST_MATH_FLAGS,
    ),
]


if __name__ == "__main__":
    setup(
        ext_modules=ext_modules,
        # Currently, build_ext only provides an optional "highest supported C++
        # level" feature, but in the future it may provide more features.
        cmdclass={"build_ext": build_ext},
    )
