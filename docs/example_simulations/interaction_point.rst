Beam-beam simulations at the interaction point
==============================================

This section contains example scripts for simulations of 
beam-beam crossing at the interaction point using `WarpX <https://blast-warpx.github.io/>`_ and `Guinea-Pig <https://gitlab.cern.ch/clic-software/guinea-pig>`_.
As an arbitrary choice, we selected the :math:`C^3` collider parameters at 125 GeV.

Overview
--------

This simulates the collision between a beam of electrons and a beam of positrons traveling at 125 GeV. 
As usual, :math:`x, y, z` are the horizontal, vertical, and longitudinal coordinates, respectively.
The beams are initialized such that their centroids are :math:`4 \sigma_z` away from interaction point. 

Beamstrahlung and coherent pair production (though not actually important for these beam parameters) are
activated. Note that incoherent processes are not enabled in either WarpX and Guinea-Pig. 
We will add them in future versions of this example. 

WarpX
-----

Resources
^^^^^^^^^
Some useful WarpX's links if you get stuck or want to know more:

 * `GitHub repository <https://github.com/BLAST-WarpX/warpx>`_
 * `ReadTheDocs documentation <https://github.com/BLAST-WarpX/warpx>`_
 * `Beam-beam tutorial <https://blast-warpx.github.io/warpx-tutorials/a-beam-beam-collision.html>`_ 

Build
^^^^^

WarpX is open-source and available on `GitHub <https://github.com/BLAST-WarpX/warpx>`_.
You can install WarpX in several different ways (for example from source, using Conda, Spack),  
and on all the major operative systems (Linux, Windows, MacOS). 
You can compile it on CPU or on GPU, depending on your machine.
`Find the instructions here <https://warpx.readthedocs.io/en/latest/install/users.html>`_.

Run
^^^

Download or copy :download:`the input file below <interaction_point/warpx_inputs.txt>`.  

.. dropdown:: See WarpX inputs
   :color: primary
   :icon: eye
   :animate: fade-in-slide-down
   :chevron: down-up
   
   .. literalinclude:: interaction_point/warpx_inputs.txt
      :language: output

On your local machine, you can run it with:: 

   mpirun -np <n_ranks> ./<warpx_executable> <input_file>

The simulation takes about ~500 seconds on two cores on an ordinary not-too-old laptop. 

Guinea-Pig
----------

Build
^^^^^

First, install `FFTW3 <https://www.fftw.org/download.html>`_.
Then follow `the instructions here <https://gitlab.cern.ch/clic-software/guinea-pig/-/blob/master/README>`_.

Run
^^^

You can download or copy :download:`the input file below <interaction_point/acc.dat>`

.. dropdown:: See Guinea-Pig inputs
   :color: primary
   :icon: eye
   :animate: fade-in-slide-down
   :chevron: down-up
   
   .. literalinclude:: interaction_point/acc.dat  
   
To run, do::

   ./<guineapig_executable> C3_250_PS1 params acc.out  
   
Visualization and postprocessing
--------------------------------

You can download :download:`this Jupyter notebbok <interaction_point/tutorial.ipynb>`
that shows how to extract and make a video of the main particle data.

Soon we will add the luminosities.

.. dropdown:: Jupyter Notebook
   :color: primary
   :icon: eye
   :animate: fade-in-slide-down
   :chevron: down-up
   :open:
   
   .. include:: interaction_point/tutorial.rst
