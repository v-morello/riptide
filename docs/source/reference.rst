API Reference
=============


TimeSeries
----------

.. autoclass:: riptide.TimeSeries
   :members:


ffa_search
----------

.. autofunction:: riptide.ffa_search


find_peaks
----------

.. autofunction:: riptide.find_peaks


Periodogram
-----------

.. autoclass:: riptide.Periodogram
   :members:


Candidate
---------

.. autoclass:: riptide.Candidate
   :members:


Save / Load data
----------------

Most objects in ``riptide`` can be converted to/from JSON via their ``to_dict()`` and ``from_dict()`` methods. 
This includes in particular :class:`riptide.TimeSeries`, :class:`riptide.Periodogram`, and :class:`riptide.Candidate`. 
It is also possible to save / load a list of such objects, a dictionary containing such objects, and any other composition that is JSON-serializable.

.. autofunction:: riptide.load_json

.. autofunction:: riptide.save_json


Kernel Functions
----------------

.. automodule:: riptide.libffa
   :members: ffa1, ffa2, ffafreq, ffaprd, boxcar_snr
