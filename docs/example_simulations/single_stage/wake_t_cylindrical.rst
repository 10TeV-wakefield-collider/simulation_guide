Single-stage laser-driven acceleration in cylindrical geometry (using Wake-T)
=============================================================================

This page describes how to run a single-stage laser-driven acceleration simulation in cylindrical geometry,
using the quasi-static code `Wake-T <https://wake-t.readthedocs.io/en/latest/>`_. For more details on the
approximations made in Wake-T and their validity, see :doc:`../../validity_of_codes/index`.

Overview:
- Laser in a flat stage, with parabolic plasma channel for guiding
- Trapezoidal beam ; beam transverse size and length need to be carefully chosen (typically with optimization)
optimize matchign and beamloading

Installation
------------

Wake-T is open-source, and available on `Github <https://github.com/AngelFP/Wake-T>`_
and `PyPI <https://pypi.org/project/Wake-T/>`_.

Instructions for installing Wake-T can be found
`here <https://wake-t.readthedocs.io/en/latest/getting_started/installation.html>`_.

Simulation script
-----------------

You can download the script below by clicking on :download:`this link <wake_t_cylindrical.py>`, and then
run it in your terminal by typing::

   python wake_t_cylindrical.py

.. literalinclude:: wake_t_cylindrical.py
   :language: python
