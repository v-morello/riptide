import os
import setuptools
import subprocess
import versioneer

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
    'matplotlib',
    'threadpoolctl', # in the pipeline, dynamically limit the number of threads used by numpy libs
    'schema'
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
    license='MIT License',

    # NOTE (IMPORTANT): This means that everything mentioned in MANIFEST.in will be copied at install time 
    # to the package’s folder placed in 'site-packages'
    include_package_data=True,

    entry_points = {
        'console_scripts': ['rffa=riptide.pipeline.pipeline:main'],
    },

    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Astronomy"
        ],

    # NOTE TO DEVELOPERS ONLY: 
    # As of March 2020, the latest official release of versioneer (0.18) does
    # not support custom cmdclass in setup.py.
    # Setting up versioneer for riptide thus requires using the tip of the master branch on github.
    # pip install git+https://github.com/warner/python-versioneer@master
    cmdclass=versioneer.get_cmdclass({
        'install': CustomInstall,
        'develop': CustomDevelop,
        'build_py': CustomBuildPy
        })
)