import os
import setuptools
import subprocess
from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.build_py import build_py


def build_libffa():
    """ A custom sequence of build actions to run """
    try:
        cwd = os.getcwd()
        thisdir, __ = os.path.split(__file__)
        srcdir = os.path.realpath(os.path.join(thisdir, 'riptide', 'c_src'))
        os.chdir(srcdir)
        subprocess.check_call('make clean', shell=True)
        subprocess.check_call('make all', shell=True)
    except:
        raise
    finally:
        os.chdir(cwd)


class CustomInstall(install):
    def run(self):
        self.announce("Building C library ...")
        build_libffa()
        install.run(self)


class CustomDevelop(develop):
    def run(self):
        self.announce("Building C library ...")
        build_libffa()
        develop.run(self)


# NOTE: It is necessary to also have a custom build_py command !
# During the first 'pip install' on a given machine, pip caches the packages
# with all its build files into a wheel (.whl). To do so, it runs the
# 'build_py' command and stores all resulting files in the .whl. If 'build_py'
# does not call build_actions(), then the .whl will NOT contain the compiled 
# binaries. And in this case, running pip install another time (in another 
# conda environment for example) actually installs that cached .whl which 
# contains no binaries, and the module does not work.
class CustomBuildPy(build_py):
    def run(self):
        self.announce("Building C library ...")
        build_libffa()
        build_py.run(self)


with open("README.md", "r") as fh:
    long_description = fh.read()


install_requires = [
    'numpy',
    'pandas',
    'astropy',
    'pyyaml',
    'h5py',
    'matplotlib'
]


# TODO: dynamic versioning
VERSION = '0.0.3'


setup(
    name='riptide-ffa',
    url='https://bitbucket.org/vmorello/riptide',
    author='Vincent Morello',
    author_email='vmorello@gmail.com',
    description='Pulsar searching with the Fast Folding Algorithm (FFA)',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=VERSION,
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    license='MIT License',

    # NOTE (IMPORTANT): This means that everything mentioned in MANIFEST.in will be copied at install time 
    # to the package’s folder placed in 'site-packages'
    include_package_data=True,

    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Fortran",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Astronomy"
        ],
    cmdclass={
        'install': CustomInstall,
        'develop': CustomDevelop,
        'build_py': CustomBuildPy
        }
)