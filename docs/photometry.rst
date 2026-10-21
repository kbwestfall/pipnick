Point-Source and Aperture Photometry
====================================

This tutorial demonstrates how to use the `photometry_all()` function
to perform photometric analysis on images that have undergone basic reduction.
Point-source photometry is performed using photutils.

Overview
--------

The `photometry_all()` function extracts and analyzes sources from
images and stores the results in `.csv` tables. These tables contain
data such as source positions and fluxes. Fluxes are calculated using PSF
photometry as well as aperture photometry. The function can also detect
multiple sources in one group, but groups can be consolidated.

The astrometry and photometry pipelines are independent and can usually
run simultaneously. However, conflicts may arise if both pipelines attempt
to access the same file concurrently.


Using the Command-Line Script
-----------------------------

For command-line usage, you can run the `pipnick_photometry`
script, which provides flexible options for executing photometric analysis.

**Basic Use**

To execute the script, use the following command:

.. code::

   pipnick_photometry <maindir> [options]

Replace `<maindir>` with the path to the directory containing all data.
The reduced directory containing images to be calibrated must be
at /maindir/reduced/.

Once you have run this script the first time, an ASCII Astropy table of
all files in /maindir/reduced/ will be saved for reference at 
/maindir/astrometry_files.tbl. You may 'comment out' files to be
ignored in future runs with a `#`.

.. code::

  pipnick_photometry <maindir> -t [options]

Toggle on use_table with `-t` to ignore all 'commented-out' files in the table.

**Using Optional Arguments**

- `-t` or `--use_table` (flag, optional)
  Whether to use the table file to automatically exclude files that have been commented-out.

- `--excl_files` (list, optional):
  List of file stem substrings to exclude (exact match not necessary).

- `--excl_objs` (list, optional):
  List of object substrings to exclude (exact match not necessary).

- `--excl_filts` (list, optional):
  List of filter substrings to exclude (exact match not necessary).

- `-th` or `--thresh` (float, optional)
  Threshold for source detection. The default is 8.0.

- `-g` or `--group` (flag, optional)
  Consolidates groups of sources detected together into one source.

- `-f` or `--fittype` (str, optional)
  Type of Moffat fit to use. Options are 'circ' (default) or 'ellip'.

- `-pf` or `--plot_final` (flag, optional)
  Displays images with sources and fluxes labeled.

- `-pi` or `--plot_inters` (flag, optional)
  Displays images with initial sources and source groups to determine
  by inspection which groups are valid.

- `-vv` or `--very_verbose` (flag, optional)
  Enables the most detailed logging. Overrides verbosity settings provided via `--verbosity`.

- `--verbosity` (int, optional)
  Sets the level of verbosity for logging. Acceptable values are 1 (CRITICAL),
  2 (ERROR), 3 (WARNING), 4 (INFO, default), and 5 (DEBUG). Overrides `--verbose`.

For example:

.. code::

   pipnick_photometry 'path/to/maindir' -t -th 10.0 -g -f ellip -pf -pi


Using the Photometry Function
-----------------------------

The same code can be run in a Python module with the same functionality.
To run the `photometry_all()` function, follow these steps:

**Basic Use**

1. **Import the Function**

   First, import the `photometry_all()` function from the
   `pipnick` package.

   .. code:: python

      from pipnick.pipelines.photometry import photometry_all

2. **Initialize Logging**

   Set up a logger to capture output from the function. You can
   configure the verbosity level to 'DEBUG', 'INFO', 'WARNING',
   'ERROR', or 'CRITICAL'. Logs are displayed in the terminal or
   console and are always saved to a `.log` file at the 'DEBUG' level.

   .. code:: python

      import logging
      from pipnick.convenience.log import adjust_global_logger

      adjust_global_logger('INFO', __name__)
      logger = logging.getLogger(__name__)

3. **Specify the Image Directory**

   Define the "main" directory containing all data. The reduced directory
   contains images to be calibrated, and must be at /maindir/reduced/.
   The .csv source tables will be saved to a directory called
   /maindir/photometric/[unconsolidated or consolidated] in the same
   parent as /maindir/reduced/. Tables will be organized by object name.

   .. code:: python

      maindir = 'path/to/data/'

4. **Run the Photometry Pipeline**

   Use the `photometry_all()` function to process the images.
   
   With default parameters, this call: 
    - Performs photometry on all images in /maindir/reduced/
    - Uses the default detection threshold (8.0 = detect only sources brighter than 8.0 x background STD)
    - Uses a circular Moffat fit
    - Preserves source groups
    - Uses the photutils setting mode = 'all' (recommended not to change--see https://photutils.readthedocs.io/en/stable/api/photutils.psf.IterativePSFPhotometry.html)

   It saves .csv source tables to /maindir/photometric/unconsolidated/,
   organized by object name.

   As with the reduction pipeline, this call also creates an ascii Astropy
   table, but of all files in /maindir/reduced/ at /maindir/astrometry_files.tbl,
   commenting out any files that were excluded.

   .. code:: python

      src_catalog_paths = photometry_all(maindir)

5. **Customizing Parameters**

   You can customize the function's behavior with various parameters.
   
   This example:
    - Uses the table produced by the first call to determine file exclusions
    - Additionally excludes any files using the ``'B'`` filter
    - Uses an elliptical Moffat fit
    - Consolidates groups of sources into one source
    - Generates matplotlib plots showing all detected sources & their fluxes
    - Generates matplotlib plots showing cutouts of all source groups to manually determine if a group should be consolidated

   .. code:: python

      src_catalog_paths = photometry_all(maindir, use_table=True, excl_filts=['B'],
                                         thresh=15.0, group=True, fittype='ellip',
                                         plot_final=True, plot_inters=True)

Viewing Results
---------------

If possible, you should run final calibrations using the
``pipnick.pipelines.final_calib`` module before exporting these
.csv files for further analysis.

The output `.csv` files contain tables of detected sources with their
positions and fluxes. These tables are organized by object name and
saved in the specified output directory. If plotting options were
enabled, Matplotlib plots will show detected sources and source groups
for further inspection.
