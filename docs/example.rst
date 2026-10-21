Example
=======

Here is a fully worked example with plots showing the output.

.. _example_pipeline:

.. code:: python

    from pipnick.pipelines.reduction import reduce_all
    from pipnick.pipelines.astrometry import astrometry_all
    from pipnick.pipelines.photometry import photometry_all
    from pipnick.pipelines.final_calib import final_calib_all
    from pipnick.convenience.graphs import plot_sources

Initialize Logger
-----------------

Initialize logger to see output of functions, setting verbosity level to 
'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'. Logs are displayed 
wherever code is being run (terminal, console, etc.), and 'DEBUG' level 
logs are always saved to a .log file.

.. code:: python

    import logging
    from pipnick.convenience.log import adjust_global_logger
    adjust_global_logger('INFO', __name__)
    logger = logging.getLogger(__name__)

Astrometry Pipeline API Key
---------------------------

To run the astrometry pipeline, you need to create an account on 
https://nova.astrometry.net/ and grab an API key. This is a randomly 
generated string tied to your user account; everything you do using this 
pipeline will be through your account on the Astrometry.net web service 
(e.g., all image uploads will show up on your web profile).

Your key is available in the My Profile section of the dashboard.

.. code:: python

    api_key = 'exampleapikey'

Define Raw Images Directory
---------------------------

Define the directory where images are stored. All raw FITS files
to be calibrated must be in /maindir/raw/.

All results of the ``pipnick`` pipeline will be saved to this
directory: reduction products will be saved to a folder called
/reduced/ in ``maindir``. If ``save_inters`` is set to True,
intermediate products will be saved to /processing/ in ``maindir``

.. code:: python

    maindir = 'data_example'

Basic Reduction
---------------

.. code:: python

    red_files = reduce_all(maindir)

Output:

.. code-block:: text

    19:03:31 - INFO - ---- reduce_all() called on directory test_data_example\\raw
    19:03:31 - INFO - 19 raw files extracted from raw directory
    19:03:31 - INFO - Saving table of file data to test_data_example\\files_table.tbl
    19:03:31 - INFO - Manually excluding files with names []
    19:03:31 - INFO - Manually excluding files with [] in object name: []
    19:03:31 - INFO - Manually excluding files with [] filters: []
    19:03:31 - INFO - Automatically excluding files with 'Focus' in object name: []
    19:03:31 - INFO - Automatically excluding files already commented out in the table file: []
    19:03:31 - INFO - Initializing CCDData objects & removing cosmic rays
    19:03:38 - INFO - Combining bias files into master bias
    19:03:38 - INFO - Using 10 bias frames: ['d1001', 'd1002', 'd1003', 'd1004', 'd1005', 'd1006', 'd1007', 'd1008', 'd1009', 'd1010']
    19:03:39 - INFO - Saving master bias to test_data_example\\processing\\master_bias.fits
    19:03:39 - INFO - Combining flat files into master flat
    19:03:39 - INFO - Using 3 flat frames: ['d1035', 'd1036', 'd1037']
    19:03:39 - INFO - Saving V-band master flat to test_data_example\\processing\\master_flat_V.fits
    19:03:39 - INFO - Performing overscan subtraction & trimming on 6 science images
    19:03:39 - INFO - Saving 6 _over images overscan images to test_data_example\\processing\\overscan
    19:03:39 - INFO - Performing bias subtraction on 6 science images
    19:03:39 - INFO - Saving 6 _unbias images unbias images to test_data_example\\processing\\unbias
    19:03:40 - INFO - Performing flat division
    19:03:40 - INFO - Saving 1 _red images UGC_9837_V images to test_data_example\\reduced\\UGC_9837_V
    19:03:40 - INFO - Saving 1 _red images PG1323-085_V images to test_data_example\\reduced\\PG1323-085_V
    19:03:40 - INFO - Saving 1 _red images NGC_3982_V images to test_data_example\\reduced\\NGC_3982_V
    19:03:40 - INFO - Saving 1 _red images NGC_4147_V images to test_data_example\\reduced\\NGC_4147_V
    19:03:40 - INFO - Saving 1 _red images PG1530+057_V images to test_data_example\\reduced\\PG1530+057_V
    19:03:40 - INFO - Saving 1 _red images NGC_6205_V images to test_data_example\\reduced\\NGC_6205_V
    19:03:40 - INFO - Flat divided images saved to test_data_example\\reduced
    19:03:40 - INFO - ---- reduce_all() call ended

Astrometric Calibration
-----------------------

.. code:: python

    astro_calib_files = astrometry_all(maindir, api_key)

Output:

.. code-block:: text

    19:03:40 - INFO - ---- astrometry_all() called on directory test_data_example\\reduced
    19:03:40 - INFO - ---- astrometry_all() call ended
    19:03:40 - INFO - Zeroing out masked regions for faster astrometric solves
    19:03:40 - INFO - Connecting to astrometry.net
    19:03:41 - INFO - Returning local copy of d1052_red.fits's solution; astrometry.net not used
    19:03:41 - INFO - Returning local copy of d1048_red.fits's solution; astrometry.net not used
    19:03:41 - INFO - Returning local copy of d1059_red.fits's solution; astrometry.net not used
    19:03:41 - INFO - Returning local copy of d1047_red.fits's solution; astrometry.net not used
    19:03:41 - INFO - Returning local copy of d1083_red.fits's solution; astrometry.net not used
    19:03:41 - INFO - Returning local copy of d1074_red.fits's solution; astrometry.net not used

Photometric Calibration
-----------------------

.. code:: python

    src_catalog_paths = photometry_all(maindir, group=False, plot_final=False, plot_inters=False)

Output:

.. code-block:: text

    19:03:41 - INFO - ---- photometry_all() called on directory test_data_example\\reduced
    19:03:42 - INFO - Working on image d1052_red.fits (NGC_3982 - V) (11 stamps)
    19:03:43 - INFO - Stack FWHM = 6.499938121481804
    19:03:45 - WARNING - One or more fit(s) may not have converged. Please check the "flags" column in the output table.
    19:03:45 - INFO - Aperture photometry cannot handle masked pixels--sources with masked pixels have flux_aper = nan
    WARNING: One or more fit(s) may not have converged. Please check the "flags" column in the output table. [photutils.psf.photometry]
    19:03:48 - INFO - Working on image d1048_red.fits (NGC_4147 - V) (57 stamps)
    19:03:50 - INFO - Stack FWHM = 6.226742495381551
    19:03:59 - WARNING - One or more fit(s) may not have converged. Please check the "flags" column in the output table.
    19:03:59 - WARNING - Some groups have more than 25 sources. Fitting such groups may take a long time and be error-prone. You may want to consider using different `SourceGrouper` parameters or changing the "group_id" column in "init_params".
    19:03:59 - INFO - Aperture photometry cannot handle masked pixels--sources with masked pixels have flux_aper = nan
    WARNING: One or more fit(s) may not have converged. Please check the "flags" column in the output table. [photutils.psf.photometry]
    19:04:14 - INFO - Working on image d1059_red.fits (NGC_6205 - V) (186 stamps)
    19:04:22 - INFO - Clipped Avg FWHM = 6.880937236174258
    19:05:14 - WARNING - One or more fit(s) may not have converged. Please check the "flags" column in the output table.
    19:05:14 - INFO - Aperture photometry cannot handle masked pixels--sources with masked pixels have flux_aper = nan
    WARNING: One or more fit(s) may not have converged. Please check the "flags" column in the output table. [photutils.psf.photometry]
    19:05:19 - INFO - Working on image d1047_red.fits (PG1323-085 - V) (9 stamps)
    19:05:19 - INFO - Stack FWHM = 6.228083336601202
    19:05:21 - INFO - Aperture photometry cannot handle masked pixels--sources with masked pixels have flux_aper = nan
    19:05:22 - WARNING - One or more fit(s) may not have converged. Please check the "flags" column in the output table.
    19:05:22 - INFO - Aperture photometry cannot handle masked pixels--sources with masked pixels have flux_aper = nan
    WARNING: One or more fit(s) may not have converged. Please check the "flags" column in the output table. [photutils.psf.photometry]
    19:05:23 - INFO - Working on image d1083_red.fits (UGC_9837 - V) (12 stamps)
    19:05:23 - INFO - Stack FWHM = 7.189019259783787
    19:05:27 - WARNING - One or more fit(s) may not have converged. Please check the "flags" column in the output table.
    19:05:27 - INFO - Aperture photometry cannot handle masked pixels--sources with masked pixels have flux_aper = nan
    19:05:27 - INFO - ---- photometry_all() call ended

Final Calibration (Convert Pixel Coordinates -> RA/Dec)
-------------------------------------------------------

.. code:: python

    astrophot_data_tables = final_calib_all(maindir)

Output:

.. code-block:: text

    19:05:27 - INFO - Saving photometric source catalogs with sky coordinates (RA/Dec) to test_data_example\\final_calib\\astrophotsrcs\\NGC_3982_V
    19:05:27 - INFO - Saving photometric source catalogs with sky coordinates (RA/Dec) to test_data_example\\final_calib\\astrophotsrcs\\NGC_4147_V
    19:05:27 - INFO - Saving photometric source catalogs with sky coordinates (RA/Dec) to test_data_example\\final_calib\\astrophotsrcs\\NGC_6205_V
    19:05:27 - INFO - Saving photometric source catalogs with sky coordinates (RA/Dec) to test_data_example\\final_calib\\astrophotsrcs\\PG1323-085_V
    19:05:27 - INFO - Saving photometric source catalogs with sky coordinates (RA/Dec) to test_data_example\\final_calib\\astrophotsrcs\\PG1530+057_V
    19:05:27 - INFO - Saving photometric source catalogs with sky coordinates (RA/Dec) to test_data_example\\final_calib\\astrophotsrcs\\UGC_9837_V

Display Images & Annotate Sources
---------------------------------

.. code:: python

    for object, src_table_dict in astrophot_data_tables.items():
        for file_key, src_table in src_table_dict.items():
            plot_sources(src_table, given_fwhm=8.0, flux_name='flux_psf', scale=1.5)

Output:

.. code-block:: text

    19:12:40 - INFO - Plotting image d1052_red.fits (NGC_3982 - V)
    19:12:40 - INFO - Plotting image d1059_red.fits (NGC_6205 - V)
    19:12:40 - INFO - Plotting image d1047_red.fits (PG1323-085 - V)
    19:12:40 - INFO - Plotting image d1083_red.fits (UGC_9837 - V)
    19:12:40 - INFO - Plotting image d1048_red.fits (NGC_4147 - V)

Images:

.. image:: /figures/d1052_sources.png
    :alt: NGC 3982 - V
    :align: center
    :width: 600px

.. image:: /figures/d1059_sources.png
    :alt: NGC 6205 - V
    :align: center
    :width: 600px

.. image:: /figures/d1047_sources.png
    :alt: PG1323-085 - V
    :align: center
    :width: 600px

.. image:: /figures/d1083_sources.png
    :alt: UGC 9837 - V
    :align: center
    :width: 600px

.. image:: /figures/d1048_sources.png
    :alt: NGC 4147 - V
    :align: center
    :width: 600px
