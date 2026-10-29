.. _metadata:

The Metadata Table
==================

Most of the pipnick scripts require a table that provides the metadata for each
raw frame collected at the telescope that effectively defines how a night is
processed.  The table is intended to be directly edited by the user, which is
primary reason for it being in the human-readable, `ECSV
<https://docs.astropy.org/en/stable/io/ascii/ecsv.html>`__ format.  The table
includes one row per raw frame and the ``ecsv`` file also contains metadata
about the reductions, such as the name of the directories with the raw data and
for the processed output.

This page provides additional detail on the content of the file and the script
used to construct it, ``pipnick_metadata``.

The algorithm
-------------

The algorithm used by ``pipnick_metadata`` is summarized as follows:

#. **Find the frames.**  Every file with the extension ``--ext`` (``.fits`` by
   default) in the raw directory is read, and its camera is identified from its
   header (or set with ``--camera``; frames from other cameras are then
   skipped).  The supported cameras are the Nickel C2 and Andor cameras; see
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
   that the closed-shutter frames but the counts are significantly above the
   bias level, or was identified as a flats but has too few counts.  Frames
   labeled as calibrations are reclassified as biases if their level is that
   expected for a bias.  **Rejected frames are commented out of the table, with
   the reason in the** ``note`` **column.**

#. **Assign the setups.**  Frames with the same camera, binning, readout mode,
   and readout window share a *setup*, labeled ``A``, ``B``, etc. in order of
   first appearance.  Calibrations can only be applied to frames of the same
   setup.

#. **Write the file(s).** When the ``--split`` option is used, the script will
   write one table per setup, add a ``comb_id`` column that groups the frames to
   combine, and comment out the pointing and focus frames.  When the
   ``--assign_mosaics`` option is used, the table also assigns frames to a set
   of expected mosaics (see :doc:`mosaic`).

Inputs
------

The only inputs required for this step is a directory with a set of raw frames.

.. important::

    **Nothing should ever be written to the raw directory.**

Execution
---------

An example call looks like this:

.. code-block:: console

    pipnick_metadata -d raw -t n3.ecsv --split --assign_mosaics

The options are:

- ``-d``, ``--rawdir``: The directory with the raw frames of one night
  (required to build a table).
- ``-t``, ``--raw_table``: The table, or, with ``--split``, the base name of
  the per-setup tables.  If not provided, a name is made from the camera and
  the time.
- ``--split``: Write one table per setup (see below).  The reduction requires
  per-setup tables.
- ``--assign_mosaics`` (requires ``--split``): Label the science frames with
  the mosaics they belong to, grouped by pointing and filter (see
  :doc:`mosaic`).  Only frames without labels are assigned.
- ``--mosaic_radius``: The linking length (arcmin) of the pointing groups of
  ``--assign_mosaics``; by default, half the field of view.
- ``--camera``: The camera (``NickelC2`` or ``NickelAndor``).  By default, it
  is identified from the headers, and all frames must be from one camera.
- ``--ext``: The extension of the raw files (default ``.fits``).
- ``--classify``: ``header`` (the default; fast) or ``stats``, which also
  reads the image data to apply the checks that need it (step 4).
- ``-o``, ``--overwrite``: Rebuild the tables from scratch, **discarding any
  edits**.
- ``-a``, ``--append``: Add the new frames in the raw directory to existing
  tables, keeping the existing rows, including edits and commented-out rows.
  Use this when new frames are taken during the night.

A table that exists is never overwritten without ``--overwrite`` or
``--append``.  The logging options (``-v``, ``--log_file``, ``--log_level``)
are common to every script; see :doc:`common`.

Outputs
-------

- **The table**, or, with ``--split``, **one table per setup**, named
  ``<base>_<xbin>x<ybin>_<readmode>_<window>.ecsv``: e.g.,
  ``n3_2x2_fast_full.ecsv``, where the window is ``full`` for the full readout
  or ``win00``, ``win01``, etc. for custom readout windows.  Each is an ECSV
  file (a text file with a YAML header), with one row per frame and the
  columns:

  - ``file``: the raw file name;
  - ``frametype``: the frame type (step 3);
  - ``setup``: the setup label;
  - ``object``, ``filter``, ``exptime`` (s), ``airmass``, ``ra``, ``dec``
    (sexagesimal), ``date`` (UT), ``mjd``: from the header;
  - ``xbin``, ``ybin``, ``readmode``, ``naxis1``, ``naxis2``, ``window``: the
    setup values (in a per-setup table, these are stored once, in the
    ``setup`` entry of its metadata);
  - ``note``: why a frame was rejected or reclassified;
  - per-setup tables only: ``comb_id``, which groups the frames to combine
    (see below), and ``mosaic`` and ``wcs_group``, the labels of the mosaics
    (empty unless assigned; see :doc:`mosaic`).

  The table metadata record the camera, the raw directory, the setups, the
  classification method, the pipnick version, and, once the data are reduced,
  the reduction directory (``rdxdir``).
- **The log**, in the parent of the raw directory (``log/``); see
  :doc:`common`.

Editing the table
-----------------

The reduction does exactly what the table says, and the table can be edited
with any text editor:

- **Exclude a frame** by commenting out its row: add a ``#`` at the
  beginning of the line.  Restore it by removing the ``#``.
- **Correct a value**, e.g., an object name or a filter.  The table, not the
  raw header, is the authority for the object, filter, exposure time,
  coordinates, and airmass, and the reduction records the header value it
  replaced.
- **Change the combination groups**, ``comb_id``: frames with the same
  ``comb_id`` are combined into one product.  By default, the biases are
  combined; the flats are combined by type and filter; darks are combined if
  their exposure times agree to within 0.1 s; and every science frame is
  processed on its own.  Give several science frames the same ``comb_id`` to
  combine them (without registration, so they must share a pointing).
- **Define the mosaics** by editing the ``mosaic`` and ``wcs_group`` labels;
  see :doc:`mosaic`.

Edits are preserved by ``--append``, and discarded by ``--overwrite``.  After
an edit, rerun the reduction: it reprocesses only the products whose frames
or values changed.

Troubleshooting
---------------

- **"Files in ... are from more than one camera"**: give ``--camera``, and
  the frames from the other cameras are skipped (or separate the files).
- **"Unable to determine the camera used to obtain the data"**: the header of
  a frame matches no supported camera (e.g., a frame from another
  instrument).  Remove it from the raw directory.
- **A good frame is commented out**: its reason is in the ``note`` column.
  If the check is wrong for your data, restore the row.
- **Frames are misclassified** (e.g., a science frame labeled with a flat's
  object name): correct its ``frametype`` in the table, and check its
  ``comb_id``.
- **Flats with too little light, or empty frames, are not caught**: use
  ``--classify stats``.
- **The C2 fast readout**: the header value of the C2's fast readout is not
  yet confirmed, so C2 frames taken with it may be rejected as having an
  undefined readout setting; see :doc:`issues`.

Using the function in Python
----------------------------

.. code-block:: python

    from pipnick.metadata.build import build_metadata

    files, tbl = build_metadata(rawdir='raw', filename='n3.ecsv', split=True, assign=True)

See :func:`~pipnick.metadata.build.build_metadata`, and
:func:`~pipnick.metadata.table.read_metadata` to read a table (with its
commented-out rows, if needed).
