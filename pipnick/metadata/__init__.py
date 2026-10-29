"""
Construct and manage the table of metadata for a set of observations.

The metadata table is the first product of pipnick (see the ``pipnick_metadata``
script) and defines how the data are processed.  It has one row per frame,
with the frame type, the instrument setup, and the header quantities needed
for the reduction, normalized across the Nickel cameras (see
:mod:`pipnick.cameras`).  The table is written as an ECSV file.

The package is organized as follows: the table input and output
(:mod:`~pipnick.metadata.table`), the setups and per-setup files
(:mod:`~pipnick.metadata.setups`), the combination groups and mosaics
(:mod:`~pipnick.metadata.groups`), the names derived from the metadata
(:mod:`~pipnick.metadata.names`), and the construction of the table
(:mod:`~pipnick.metadata.build`).

Setups
------

Frames obtained with the same camera, binning, readout mode, and readout
window share a *setup*, labeled ``A``, ``B``, etc. in order of first
appearance.  A typical night uses a single setup.  Calibrations (biases,
darks, and flats) can only be applied to frames with the same setup.  The
definition of each setup is stored in the ``setups`` entry of the table
``meta`` dictionary.

Per-setup files
---------------

The metadata can instead be written to one file per setup (see
``pipnick_metadata --split`` and
:func:`~pipnick.metadata.build.build_metadata`), named
``<base>_<xbin>x<ybin>_<readmode>_<window>.ecsv``, where the window is
``full`` for the full readout or ``win00``, ``win01``, etc. for custom readout
windows.  Because the setup values and the image size are the same for every
frame in the file (see :data:`~pipnick.metadata.setups.setup_file_keys`),
they are stored in the ``setup`` entry of the file ``meta`` dictionary instead
of as columns.  :func:`~pipnick.metadata.table.read_metadata` restores them as
columns, so a per-setup file can be used anywhere a table with all setups
can.

Per-setup files also differ from a table with all setups in two ways:

- They have a ``comb_id`` column with integers that identify the frames to
  combine; see :func:`~pipnick.metadata.groups.assign_comb_ids`.  By default,
  biases are combined, flats are combined by type and filter, darks are
  combined if their exposure times agree to within 0.1 s, and all other frames,
  including science frames, are processed independently.  Edit the column to
  change the groups.

- Pointing and focus frames are commented out by default.

- They have ``mosaic`` and ``wcs_group`` columns with string labels that
  identify the reduced science images to combine into mosaics, and the mosaics
  that share an output WCS; see :data:`~pipnick.metadata.groups.mosaic_columns`
  and :func:`~pipnick.metadata.groups.assign_mosaics`.  The labels are empty
  (no mosaics) unless they are assigned automatically (``pipnick_metadata
  --assign_mosaics``) or by editing the table.

Excluding frames
----------------

Frames are excluded from the reduction by *commenting out* their rows in the
ECSV file, i.e., by adding a ``#`` to the beginning of the line.  Frames that
fail the quality checks (see
:func:`~pipnick.cameras.cameras.NickelCamera.assess`) are commented out
automatically, with the reason given in the ``note`` column.  Users can
comment out (or restore) any row, and can edit the table, e.g., to correct a
frame type.  These edits are preserved when new frames are appended to the
table; see :func:`~pipnick.metadata.build.build_metadata`.
"""
