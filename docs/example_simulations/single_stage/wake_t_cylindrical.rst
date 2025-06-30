Single-stage laser-driven acceleration in cylindrical geometry (using Wake-T)
=============================================================================

This page describes how to run a **single-stage laser-driven acceleration simulation in cylindrical geometry**,
using the quasi-static code `Wake-T <https://wake-t.readthedocs.io/en/latest/>`_.

For more details on the
approximations made in Wake-T, and how they compare to other codes, see :doc:`../../validity_of_codes/index`.

Overview
--------

This simulates the acceleration of a 200 pC electron beam in laser-driven plasma stage.
The plasma density profile is longitudinally uniform and transversally parabolic, so as to **guide the laser driver**.
The electron beam profile is transversally Gaussian and longitudinally trapezoidal.

- The transverse beam size has been tuned so as to be **matched to the transverse focusing** in the plasma wakefield.
  (This minimizes oscillations of the beam size during acceleration, thus ensuring better preservation of emittance.)

- The trapezoidal beam profile has been tuned so as to **optimize beamloading**, i.e. ensure that each longitudinal
  slice of the beam gains the same amount of energy over the length of the plasma stage, to avoid large energy spread.

.. warning::

   Is ion motion turned on by defauld in the latest version of Wake-T? It is not clear from the documentation.

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

The script takes about 12 minutes to run on a single CPU core, and is numerically converged.

.. literalinclude:: wake_t_cylindrical.py
   :language: python

Visualization and postprocessing
--------------------------------

.. warning::

   This section needs to be completed.