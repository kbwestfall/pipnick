.. _metadata:

The Metadata Table
==================

Most of the pipnick scripts require a table that provides the metadata for each
raw frame collected at the telescope that effectively defines how a night is
processed.  The table is intended to be directly edited by the user, which is
a primary reason for it being in the human-readable, `ECSV
<https://docs.astropy.org/en/stable/io/ascii/ecsv.html>`__ format.  The table
includes one row per raw frame and the ``ecsv`` file also contains metadata
about the reductions, such as the name of the directories with the raw data and
for the processed output.

This page provides additional detail on the content of the file and the script
used to construct it, ``pipnick_metadata``.

.. _comb_id:

Image combination groups
------------------------

A key component of the metadata table is that it indicates which frames should
be combined.  The table produced by ``pipnick_metadata`` includes a ``comb_id``
column that is a simple integer indicating the combination group to which each
frame belongs.  All frames with the same integer are combined into a single
processed image, and frames with the same ``comb_id`` **must** have the same
frame type.  By default, ``pipnick_metadata`` will group all biases, all dome
flats in a given filter, all sky flats in a given filter, and all dark frames
with a given exposure time.  Each science frame is, by default, processed
independently.  **An important aspect of editing the file is ensuring that the
``comb_id`` value for each frame is correct.**

Frames should only be combined if they are expected to have
*exactly* the same pointing.  That is, the images in a given combination group
are averaged pixel-by-pixel with no offsets.  To combine science images with
small on-sky offsets (dithers), use :ref:`mosaic`.

The algorithm
-------------

The algorithm used by ``pipnick_metadata`` is summarized as follows:

#. **Find the frames.**  Every file with the extension ``--ext`` (``.fits`` by
   default) in the raw directory is read, and its camera is identified from its
   header (or set with ``--camera``; frames from other cameras are then
   skipped).  Frames whose camera cannot be identified are commented out.  The
   supported cameras are the Nickel C2 and Andor cameras; see
   :mod:`pipnick.cameras`.

#. **Read the header quantities.**  The extracted data are (normalized across
   the two cameras): object, filter, exposure time, airmass, coordinates,
   observation time, binning, readout mode, image size, and readout window.

#. **Classify each frame.**  Based on the header data, each frame is classified
   as a ``bias``, ``dark``, ``domeflat``, ``skyflat``, ``science``, ``focus``,
   ``pointing``, or ``unknown``.  This classification is performed using a
   deterministic set of rules in the camera's configuration file (e.g.,
   closed-shutter frames shorter than 1 s are biases; an object name of "dome
   flat" marks a dome flat).  **This classification can be wrong and must be
   checked.**

#. **Check each frame.** Each frame is compared with the camera's ``[quality]``
   limits set in its configuration file.  From the header, frames taken when the
   detector is too warm or with an undefined readout setting are rejected.  If
   the ``--classify stats`` option is set, the image data are also read and an
   image is rejected if it is empty, has too many saturated pixels, indicates
   that the shutter is closed but the counts are significantly above the bias
   level, or was identified as a flat but has too few counts.  Frames labeled as
   calibrations are reclassified as biases if their level is that expected for a
   bias.  **Rejected frames are commented out of the table, with the reason in
   the** ``note`` **column.**

#. **Assign the setups.**  Frames with the same camera, binning, readout mode,
   and readout window share a *setup*, labeled ``A``, ``B``, etc. in order of
   first appearance.  Calibrations can only be applied to frames of the same
   setup.

#. **Write the file(s).** When the ``--split`` option is used, the script will
   write one table per setup, add a ``comb_id`` column that groups the frames to
   combine, and comment out the pointing and focus frames.  When the
   ``--assign_mosaics`` option is used, the table also assigns frames to a set
   of expected mosaics (see :ref:`mosaic`).

Inputs
------

The only input required for this step is a directory with a set of raw frames.

.. important::

    **Nothing should ever be written to the raw directory.**

Execution
---------

An example call looks like this:

.. code-block:: console

    pipnick_metadata raw n3 --split --assign_mosaics

The two required arguments are:

- ``rawdir``: The directory (either an absolute path or a path relative to the
  current working directory) with the raw frames of one night.

- ``raw_table``: The name of the table, with or without the ``.ecsv``
  extension: ``n3`` and ``n3.ecsv`` are both written to ``n3.ecsv``.  With
  ``--split``, it is the base name of the per-setup tables: e.g., ``n3`` (or
  ``n3.ecsv``) gives ``n3_2x2_fast_full.ecsv``.

The options are:

- ``--split``: Write one table per setup (see below).  The only reason to *not*
  set this option is if you first want to look at all frames collected into a
  single file.  All of the processing scripts expect the provided table to only
  contain data from a single instrument configuration, meaning you will commonly
  use this option.

- ``--assign_mosaics`` (requires ``--split``): Label the science frames with the
  mosaics they belong to, grouped by pointing and filter (see :ref:`mosaic`).
  Only frames without labels are assigned.

- ``--mosaic_radius``: The linking length (arcmin) of the pointing groups of
  ``--assign_mosaics``; by default, half the field of view.

- ``--camera``: The camera (``NickelC2`` or ``NickelAndor``).  By default, it
  is identified from the headers, and all frames must be from one camera.

- ``--ext``: The extension of the raw files (default ``.fits``).

- ``--classify``: ``header`` (the default; fast) or ``stats``, which also reads
  the image data to apply the checks that need it (see algorithm step 4).

- ``-o``, ``--overwrite``: Rebuild the tables from scratch, **discarding any
  edits**.

- ``-a``, ``--append``: Add the new frames in the raw directory to existing
  tables, keeping the existing rows, including edits and commented-out rows.
  Use this if you're processing data as it is collected and you want to update
  the table with newly acquired data.

A table is never written to a file that already exists without setting
``--overwrite`` or ``--append``.  The logging options (``-v``, ``--log_file``,
``--log_level``) are common to every script; see :ref:`common`.

Outputs
-------

The primary output is one or more metadata tables.  When ``--split`` is not set,
the name is as provided on the command line.  When using ``--split``, the file
name is ``<base>_<xbin>x<ybin>_<readmode>_<window>.ecsv``; e.g.,
``n3_2x2_fast_full.ecsv``.

- ``<base>`` is ``raw_table`` without the ``.ecsv`` extension

- ``<xbin>`` and ``<ybin>`` are the number of pixels binned during readout,

- ``<readmode>`` is ``fast`` or ``slow``,

- and ``<window>`` is ``full`` for the full readout or ``win00``, ``win01``,
  etc. for custom readout windows.

Each `ECSV <https://docs.astropy.org/en/stable/io/ascii/ecsv.html>`__ table has
one row per frame with the following columns:

- ``file``: the raw file name;

- ``frametype``: the frame type (algorithm step 3);

- ``setup``: the setup label (algorithm step 5);

- ``object``, ``filter``, ``exptime``, ``airmass``, ``ra``, ``dec`` ``date``,
  ``mjd``: details regarding the target and exposure all pulled (or constructed)
  from the header;

- ``xbin``, ``ybin``, ``readmode``, ``naxis1``, ``naxis2``, ``window``: the
  setup values (in a per-setup table, these are stored once, in the file header;
  see below);

- ``note``: why a frame was rejected or reclassified;

- For per-setup tables only: ``comb_id``, which groups the frames to combine
  (see :ref:`above<comb_id>`), and ``mosaic`` and ``wcs_group``, the labels of
  the mosaics (empty unless assigned; see :ref:`mosaic`).

The ECSV table also contains a commented YAML format heaader with setup
metadata.  This metadata can be edited via various scripts (e.g.,
``pipnick_reduce`` adds the data with the reduction output).  The base-level
information added by ``pypeit_metadata`` includes the camera, the raw directory,
the setups, the classification method, and the pipnick version.

As with all pipnick scripts, a log is also produced in the parent of the raw
directory (``log/``); see :ref:`common`.

Editing the table
-----------------

The ECSV table produced by the ``pypeit_metadata`` script is really just a
starting point.  **You should always check the file and correct it.**  The table
is used by nearly all of the remaining scripts and becomes the authoritative
version of the metadata.  That is, for any quantities in the table that differ
from header values, pipnick always takes the value in the table.  This allows
you to correct values (such as the object name) that are real errors in the raw
file headers.  You can edit the file using any simple text editor.

Some typical edits you will make are:

- **Exclude a frame:** Force pipnick to ignore any frame by commenting out its
  row: add a ``#`` at the beginning of the line.  Restore it by removing the
  ``#``.

- **Correct any value:**  The table, not the raw header, is the authority for
  all data in the file, as discussed above.

- **Change the combination groups:**  As discussed :ref:`above<comb_id>`, frames
  with the same ``comb_id`` are combined pixel-by-pixel.  By default, the biases
  are combined; the flats are combined by type and filter; darks are combined if
  their exposure times agree to within 0.1 s; and every science frame is
  processed on its own.  To combine multiple science frames, give them the same
  ``comb_id``.  The integers used for the ``comb_id`` are largely irrelevant
  (they don't need to take on any given value or be in sequence), but, in terms
  of the naming convention (see :ref:`reduce_outputs`), they are expected to be
  999 or lower.

- **Define the mosaics** by editing the ``mosaic`` and ``wcs_group`` labels;
  see :ref:`mosaic`.

.. note::
    
    Edits are preserved when using the ``--append`` optin, but they are ignored
    when using ``--overwrite``

If you edit the table after running any of the processing scripts, pipnick will
incorporate the edit into subsequent processing steps and reprocesses only the
outputs that rely on frames whose data was changed.

Troubleshooting
---------------

The following are a set of issues that may occur, their possible causes, and our
suggestion for how to deal with the issue.

- **"Files in ... are from more than one camera"**:  This message should be rare
  and indicates that frames from more than one camera were found.  Fix this by
  using the ``--camera`` option to set the camera by hand.  Frames not from this
  camera will be ignored and left out of the table.

- **"Unable to determine the camera used to obtain the data"**: The header of a
  frame matches no supported camera.  It will be commented out in the table.  We
  **highly** recommend that you only keep the raw observing files in the "raw
  directory."

- **A good frame is commented out**: pipnick will provide the reason the frame
  was commented out in the ``note`` column.  If the reason is incorrect, you can
  simply uncomment the line, but you may need to correct metadata that led to
  the frame being commented out.

- **Frames are misclassified:**  For example, a science frame could be labeled
  as a flat, or a bias as a dark.  We expect this kind of misclassification to
  happen often.  **It is critical that you catch and correct these errors.**
  When you change the frame type, you may also need to correct the ``comb_id``.
  Remember, frames with the same ``comb_id`` **must** have the same frame type;
  an error is raised otherwise.

- **Flats with little to no flux will not typically be caught:**  The default
  approach of classifying frames uses header data.  You may get better automatic
  classifications if you use ``--classify stats``.

- **The C2 fast readout cannot yet be classified:** The header value of the C2's
  fast readout is not yet confirmed, so C2 frames taken in this mode may be
  rejected as having an undefined readout setting; see :ref:`issues`.

Using the function in Python
----------------------------

.. code-block:: python

    from pipnick.metadata.build import build_metadata

    files, tbl = build_metadata(rawdir='raw', filename='n3.ecsv', split=True, assign=True)

See :func:`~pipnick.metadata.build.build_metadata`, and
:func:`~pipnick.metadata.table.read_metadata` to read a table (with its
commented-out rows, if needed).
