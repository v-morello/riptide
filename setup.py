import os
import sys
import setuptools
import subprocess
import versioneer

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        return pybind11.get_include()


ext_modules = [
    Extension(
        # Must match the module name passed to the PYBIND11_MODULE macro in the C++ source
        'riptide.libcpp',

        # Sort input source files to ensure bit-for-bit reproducible builds
        # (https://github.com/pybind/python_example/pull/53)
        sorted(['riptide/cpp/python_bindings.cpp']),

        # Path to pybind11 headers
        include_dirs=[get_pybind_include()],
        language='c++'
    ),
]


# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    import os
    with tempfile.NamedTemporaryFile('w', suffix='.cpp', delete=False) as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        fname = f.name
    try:
        compiler.compile([fname], extra_postargs=[flagname])
    except setuptools.distutils.errors.CompileError:
        return False
    finally:
        try:
            os.remove(fname)
        except OSError:
            pass
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']
    for flag in flags:
        if has_flag(compiler, flag):
            return flag
    raise RuntimeError('Unsupported compiler -- at least C++11 support is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': ['-O3', '-march=native', '-ffast-math'],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')

        for ext in self.extensions:
            ext.define_macros = [('VERSION_INFO', '"{}"'.format(self.distribution.get_version()))]
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        build_ext.build_extensions(self)


with open("README.md", "r") as fh:
    long_description = fh.read()


install_requires = [
    # NOTE: order matters when pip resolves which numpy version
    # is required. astropy seems to have the strictest requirements right now
    # Apparently pip will get a real dependency resolver in the near future
    'astropy',
    'pandas',
    'schema',
    'pyyaml',
    'threadpoolctl', # in the pipeline, dynamically limit the number of threads used by numpy libs
    'pytest',
    'pytest-cov',
    'matplotlib',
]


setup_requires = [
    'pybind11>=2.5.0'
]


setup(
    name='riptide-ffa',
    version=versioneer.get_version(),
    url='https://github.com/v-morello/riptide',
    author='Vincent Morello',
    author_email='vmorello@gmail.com',
    description='Pulsar searching with the Fast Folding Algorithm (FFA)',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    setup_requires=setup_requires,
    license='MIT License',

    # NOTE (IMPORTANT): This means that everything mentioned in MANIFEST.in will be copied at install time 
    # to the package’s folder placed in 'site-packages'
    include_package_data=True,

    entry_points = {
        'console_scripts': [
            'rffa=riptide.pipeline.pipeline:main',
            'rseek=riptide.apps.rseek:main'
        ],
    },

    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: C",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering :: Astronomy"
        ],

    ext_modules=ext_modules,

    # NOTE TO DEVELOPERS ONLY: 
    # As of March 2020, the latest official release of versioneer (0.18) does
    # not support custom cmdclass in setup.py.
    # Setting up versioneer for riptide thus requires using the tip of the master branch on github.
    # pip install git+https://github.com/warner/python-versioneer@master
    cmdclass=versioneer.get_cmdclass({'build_ext': BuildExt}),
)