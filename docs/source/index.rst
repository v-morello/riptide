riptide: Finding pulsars with the Fast Folding Algorithm (FFA)
==============================================================

.. image:: http://img.shields.io/badge/astro.ph-2004.03701-B31B1B.svg
   :target: https://arxiv.org/abs/2004.03701
   :alt: arXiv
  
.. image:: https://travis-ci.com/v-morello/riptide.svg?branch=master
   :target: https://travis-ci.com/v-morello/riptide
   :alt: Build status
  
.. image:: https://codecov.io/gh/v-morello/riptide/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/v-morello/riptide
   :alt: Coverage


``riptide`` ("sea\ **r**\ ch\ **i**\ ng for **p**\ ulsars in the **ti**\ me **d**\ omain") is a pulsar searching
package that implements the Fast Folding Algorithm (FFA), the theoretically optimal search method
for periodic signals. Its interface is entirely in python while the core algorithms are implemented
in C++. riptide provides:  

* A library of functions and classes to manipulate and search individual dedispersed time series  
* A pipeline executable to process a set of DM trials and output a list of candidate files, plots and other data products


Citation
--------

If using ``riptide`` contributes to a project that leads to a scientific publication, please cite the article:  
`Optimal periodicity searching: Revisiting the Fast Folding Algorithm for large scale pulsar surveys`__

__ https://arxiv.org/abs/2004.03701


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   quickstart
   kernfuncs
   pipeline
   docker
   reference

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


