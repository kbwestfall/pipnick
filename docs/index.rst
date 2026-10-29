pipnick: the data-reduction pipeline of the Nickel 1-m telescope
================================================================

.. image:: figures/pipnick-logo-with-text.png
   :width: 240px
   :align: center
   :alt: The pipnick logo: the Nickel dome on a pixel grid

pipnick reduces images taken with the science cameras of the 1-meter Nickel
telescope at Lick Observatory on Mt. Hamilton, CA.  It takes the raw frames
of a night to reduced images with astrometric solutions, calibrated source
catalogs, and image mosaics.

The processing is done by a sequence of scripts, each run on the table of raw
frames of one instrument setup:

.. code-block:: text

    pipnick_metadata     the table of raw frames: frame types, setups, groups
          |
    pipnick_reduce       bias, flat field, cosmic rays: the reduced images
          |
    pipnick_astrometry   the astrometric solution of each image
          |
    pipnick_photometry   PSF and aperture photometry: the source tables
          |
    pipnick_calibrate    the photometric calibration, from standard stars
          |
    pipnick_mosaic       mosaics of the images of each field and filter
          |
    pipnick_rgb          colour images from three mosaics

At any point, ``pipnick_status``, ``pipnick_report``, and ``pipnick_inspect``
show what has been done and what needs attention.  Start with
:doc:`install` and the :doc:`quickstart`.

Please `submit an issue <https://github.com/UCObservatories/pipnick/issues>`__
if you have any problems.

.. toctree::
   :maxdepth: 2

   getting_started
   scripts
   using
   reference
