.. image:: figures/pipnick-logo-with-text.png
   :width: 320px
   :align: center
   :alt: pipnick logo

The Photometric Imaging Pipeline for the Nickel Telescope
=========================================================

Welcome!  pipnick is open-source Python software that processes images taken by
the science cameras mounted on the 1-meter Nickel telescope at Lick Observatory
on Mt. Hamilton, CA.  The code should successfully reduce and analyze data taken
with the previous C2 camera and the new Andor camera.

Starting with the raw data obtained at the telescope, pipnick enables you to
perform:

- :ref:`reduce`
- :ref:`astrometry`
- :ref:`photometry`
- :ref:`calibrate`
- :ref:`mosaic`
- :ref:`rgb`

The processing is done by a sequence of command-line scripts.  Each script is a
simple wrapper of an underlying python function that can also be executed
directly from, say, a Jupyter notebook.  We also provide status and inspection
functions that can be executed at any time to determine the state and success of
each processing step.

To use pipnick with your data, start with :ref:`install` and the
:ref:`quickstart`.

Please `submit an issue <https://github.com/UCObservatories/pipnick/issues>`__
if you run into bugs or have specific feature requests.

.. toctree::
   :maxdepth: 2

   getting_started
   scripts
   using
   API <api/pipnick>
