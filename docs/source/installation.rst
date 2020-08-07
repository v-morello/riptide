Installation
============

The easiest method is to use pip install, which pulls the latest release from PyPI and installs all required dependencies:

.. code-block:: console

   pip install riptide-ffa


Alternatively you can clone the repository and run ``make install``

.. code-block:: console

   git clone https://github.com/v-morello/riptide
   cd riptide/
   make install

This simply runs ``pip install`` in `editable mode`_, which means you can freely edit the code. 
It also installs any required dependencies with ``pip`` that are not present already. You can check
that it all works by running the test suite in a Python or IPython console:

    >>> import riptide
    >>> riptide.test()

.. _`editable mode`: https://pip.pypa.io/en/latest/reference/pip_install/#editable-installs

There should also now be two command-line apps in your python environment:

* ``rseek``: A lightweight app to search a single time series and print significant candidates found, useful for quick data checks. See :ref:`The rseek command-line app`.
* ``rffa``: The full end-to-end pipeline to search multiple DM trials, see :ref:`Using the Pipeline`.
