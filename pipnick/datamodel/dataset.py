"""
The data set of one setup: its metadata table, reduction directory, and
products.

A per-setup metadata table (see :func:`~pipnick.metadata.build.build_metadata`)
defines the frames of one setup and their combination groups (``comb_id``).
Each group to be reduced produces one product (:class:`Product`), e.g., a
combined bias or flat, or a reduced science image.  :class:`DataSet` reads
and checks the table (:meth:`DataSet.from_file`), and provides the products,
the metadata table values written to them, and the reduced science images
(:class:`~pipnick.datamodel.image.ReducedImage`) used by the later steps.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.time import Time

from pipnick import cameras, log
from pipnick.datamodel.image import ReducedImage
from pipnick.metadata.build import build_metadata
from pipnick.metadata.names import filter_name, safe_name
from pipnick.metadata.setups import setup_tag
from pipnick.metadata.table import frame_path, read_metadata, update_metadata_meta
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickPathError
from pipnick.pkg.paths import step_dir


@dataclass(frozen=True, eq=False)
class Product:
    """
    One combination group of a data set to be reduced, and its output file.

    Products are constructed by :meth:`DataSet.products`.

    Parameters
    ----------
    frametype : str
        The frame type of the group, e.g., ``'bias'`` or ``'science'``.
    comb_id : int
        The combination group.
    indx : :class:`numpy.ndarray`
        The rows of its frames in the metadata table, sorted by MJD; see
        :func:`_group_indices`.
    filter : str
        The filter of its first frame (an empty string if unknown).
    exptime : float, None
        The exposure time of its first frame (None if unknown).
    file : :class:`pathlib.Path`
        The output file.
    """
    frametype: str
    comb_id: int
    indx: np.ndarray
    filter: str
    exptime: float | None
    file: Path

    @property
    def name(self):
        """The name of the output file."""
        return self.file.name

    def matches(self, item):
        """
        Determine if the product is selected by an item: its file name, its
        file name without the extension, or its combination group
        (``comb_id``).

        Parameters
        ----------
        item : str, int
            The item.

        Returns
        -------
        bool
            True if the item selects the product.
        """
        return (
            str(item) in [self.file.name, self.file.stem]
            or (str(item).isdigit() and self.comb_id == int(item))
        )


class DataSet:
    """
    The data set of one setup: its metadata table, reduction directory, and
    camera.

    Use :meth:`from_file` to read a per-setup metadata file.

    Parameters
    ----------
    metafile : :class:`pathlib.Path`
        The metadata file.
    table : :class:`astropy.table.Table`
        The metadata table, excluding commented-out rows; see
        :func:`~pipnick.metadata.table.read_metadata`.
    rdxdir : :class:`pathlib.Path`
        The reduction directory, which may not yet exist.
    tag : str
        The setup tag; see :func:`~pipnick.metadata.setups.setup_tag`.
    camera : :class:`~pipnick.cameras.cameras.NickelCamera`
        The camera used to obtain the data.
    """
    def __init__(self, metafile, table, rdxdir, tag, camera):
        self.metafile = metafile
        self.table = table
        self.rdxdir = rdxdir
        self.tag = tag
        self.camera = camera
        self._products = None

    @classmethod
    def from_file(cls, table, rdxdir=None, check_mosaics=False, read_only=False):
        """
        Read and check a per-setup metadata file.

        Parameters
        ----------
        table : str, :class:`pathlib.Path`
            A per-setup metadata file, written by ``pipnick_metadata
            --split``.
        rdxdir : str, :class:`pathlib.Path`, optional
            Directory for the reduced data.  If None, the directory saved in
            the metadata file is used or, if there is none, the parent
            directory of the raw data directory, which is then saved to the
            metadata file.  A directory that differs from the saved one is
            also saved; see :func:`~pipnick.metadata.build.build_metadata`.
        check_mosaics : bool, optional
            Also check the mosaic labels; see :func:`check_mosaic_groups`.
        read_only : bool, optional
            Never change the metadata file, e.g., to save the reduction
            directory; see :func:`reduction_dir`.

        Returns
        -------
        :class:`DataSet`
            The data set.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickPathError`
            Raised if the metadata file does not exist.
        :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
            Raised if the table includes more than one setup, has no frames,
            or has invalid combination groups (or mosaic labels).
        """
        _table = Path(table).absolute()
        if not _table.is_file():
            raise PipnickPathError(f'{_table} does not exist.')
        if read_only:
            metafile, metadata = _table, read_metadata(_table)
        else:
            metafile, metadata = build_metadata(filename=_table, rdxdir=rdxdir)
        check_single_setup(metadata, metafile)
        if len(metadata) == 0:
            raise PipnickDataModelError(f'{metafile} has no frames to reduce.')
        check_comb_groups(metadata, metafile)
        if check_mosaics:
            check_mosaic_groups(metadata, metafile)
        tag = setup_tag(metadata.meta['setups'][metadata['setup'][0]])
        camera = cameras.get_camera(metadata.meta['camera'])

        if read_only:
            _rdxdir = reduction_dir(_table, rdxdir=rdxdir)
        elif 'rdxdir' in metadata.meta:
            _rdxdir = Path(metadata.meta['rdxdir'])
        else:
            _rdxdir = Path(metadata.meta['rawdir']).absolute().parent
            metadata.meta['rdxdir'] = str(_rdxdir)
            # Update the file to include the reduction directory
            update_metadata_meta(metafile, rdxdir=str(_rdxdir))
        return cls(metafile, metadata, _rdxdir, tag, camera)

    def products(self, frametype=None):
        """
        Return the products of the data set: one for each combination group
        to be reduced.

        Parameters
        ----------
        frametype : str, optional
            Only return the products of this frame type, e.g.,
            ``'science'``.

        Returns
        -------
        list
            The :class:`Product` objects, ordered by ``comb_id``.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
            Raised if two groups would have the same output file, e.g., for
            flats with filter names that only differ in characters that cannot
            be used in a file name.
        """
        if self._products is None:
            self._products = _plan_products(self.table, self.rdxdir)
        if frametype is None:
            return list(self._products)
        return [p for p in self._products if p.frametype == frametype]

    def select(self, selection, frametype='science'):
        """
        Select products.

        Parameters
        ----------
        selection : list
            The products to select, by file name, file name without the
            extension, or combination group (``comb_id``); see
            :meth:`Product.matches`.
        frametype : str, optional
            The frame type of the products to select from.

        Returns
        -------
        list
            The selected products, in the order of :meth:`products`.

        Raises
        ------
        ValueError
            Raised if an entry does not match any product.
        """
        products = self.products(frametype)
        unknown = [str(i) for i in selection if not any(p.matches(i) for p in products)]
        if len(unknown) > 0:
            raise ValueError(f'No reduced science image matches: {", ".join(unknown)}.')
        return [p for p in products if any(p.matches(i) for i in selection)]

    def science_images(self, selection=None, include_missing=False):
        """
        Return the reduced science images.

        Parameters
        ----------
        selection : list, optional
            Only return the selected images; see :meth:`select`.
        include_missing : bool, optional
            Include the images that have not been reduced.  Otherwise, they
            are skipped with a warning.

        Returns
        -------
        list
            The :class:`~pipnick.datamodel.image.ReducedImage` objects,
            ordered by ``comb_id``.

        Raises
        ------
        :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
            Raised if an image is from a different setup; see
            :meth:`~pipnick.datamodel.image.ReducedImage.check_setup`.
        """
        products = self.products('science') if selection is None else self.select(selection)
        images = []
        for product in products:
            image = ReducedImage(product, self)
            if image.exists:
                image.check_setup()
            elif not include_missing:
                log.warning(f'{image.name} does not exist; reduce the data first.  Skipping.')
                continue
            images.append(image)
        return images

    def table_values(self, product):
        """
        Return the metadata table values to write to the header of a product.

        The table is the authority for these values, so that edits to it
        (e.g., a corrected object name) propagate to the reduced images.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        dict
            The values keyed by the quantity (``'filter'`` and ``'exptime'``
            for all frame types, and ``'object'``, ``'ra'``, ``'dec'``, and
            ``'airmass'`` for science frames).  Unknown values are None.  The
            object is the most common in the group (see
            :func:`_group_object`), the exposure time is the median, the
            airmass is the mean (of the frames with an airmass), and the
            other values are from the earliest frame.
        """
        return _table_values(self.table, product.indx, product.frametype)

    def calibrations(self, product):
        """
        Return the calibrations used to reduce a product, as the reduction
        selects them: the bias for darks, flats, and science frames, and the
        sky flat of the same filter (or the dome flat, if there is no sky
        flat) for science frames.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        dict
            The calibration :class:`Product` objects, keyed by their type
            (``'bias'``, ``'flat'``), with None for any that the data set does
            not have.
        """
        if product.frametype == 'bias':
            return {}
        bias = self.products('bias')
        calibs = {'bias': None if len(bias) == 0 else bias[0]}
        if product.frametype != 'science':
            return calibs
        calibs['flat'] = None
        for frametype in ['skyflat', 'domeflat']:
            flats = [p for p in self.products(frametype) if p.filter == product.filter]
            if len(flats) > 0:
                calibs['flat'] = flats[0]
                break
        return calibs

    def raw_files(self, product):
        """
        Return the raw files of a product.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        list
            The :class:`pathlib.Path` of each raw frame, in the order of
            :attr:`Product.indx`.
        """
        return frame_path(self.table, select=product.indx, path_key='rawdir')

    def exptimes(self, product):
        """
        Return the exposure time of each frame of a product.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        list
            The exposure time (s) of each frame, None if unknown, in the order
            of :attr:`Product.indx`.
        """
        return [_exptime(self.table, i) for i in product.indx]

    def pixel_scale(self, product):
        """
        Return the nominal pixel scale of a product.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        float
            The plate scale of the camera times the binning of the first
            frame (arcsec per pixel).
        """
        return self.camera.plate_scale * int(self.table['xbin'][product.indx[0]])

    def pointing(self, product):
        """
        Return the telescope pointing of a product.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        :class:`astropy.coordinates.SkyCoord`
            The coordinates of the first frame, or None if they are unknown.
        """
        i = product.indx[0]
        ra = _masked_to_none(self.table['ra'][i])
        dec = _masked_to_none(self.table['dec'][i])
        if ra is None or dec is None:
            return None
        return SkyCoord(ra, dec, unit=(units.hourangle, units.deg))

    def mjd(self, product):
        """
        Return the MJD of a product.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        float
            The MJD of the first frame, or None if it is unknown.
        """
        mjd = _masked_to_none(self.table['mjd'][product.indx[0]])
        return None if mjd is None else float(mjd)

    def night(self, product):
        """
        Return the night of a product.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        str
            The UT date (``YYYY-MM-DD``) of the first frame, or None if it is
            unknown; see :func:`night`.
        """
        mjd = self.mjd(product)
        return None if mjd is None else night(mjd)

    def mosaic_labels(self, product):
        """
        Return the mosaic labels of a product.

        Parameters
        ----------
        product : :class:`Product`
            The product.

        Returns
        -------
        tuple
            The ``mosaic`` and ``wcs_group`` labels (see
            :data:`~pipnick.metadata.groups.mosaic_columns`); empty strings if the
            product is not in a mosaic, or the table has no mosaic columns.
            The WCS group defaults to the mosaic label.
        """
        if 'mosaic' not in self.table.colnames:
            return '', ''
        i = product.indx[0]
        label = str(self.table['mosaic'][i])
        group = str(self.table['wcs_group'][i])
        return label, label if len(group) == 0 else group

    def outdir(self, step):
        """
        Return the directory for the products of a step; see
        :func:`~pipnick.pkg.paths.step_dir`.
        """
        return step_dir(self.rdxdir, step)


def reduction_dir(table, rdxdir=None):
    """
    Return the reduction directory of a per-setup metadata file.

    Unlike :meth:`DataSet.from_file`, this only reads the metadata of the
    table, and never changes the file; it is used to find the log directory
    before a step starts.

    Parameters
    ----------
    table : str, :class:`pathlib.Path`
        A per-setup metadata file.
    rdxdir : str, :class:`pathlib.Path`, optional
        The reduction directory given by the user, which takes precedence.

    Returns
    -------
    :class:`pathlib.Path`
        The reduction directory: ``rdxdir``, the directory saved in the table,
        or the parent of its raw directory; see :meth:`DataSet.from_file`.
        None if it cannot be determined (e.g., the table does not exist).
    """
    if rdxdir is not None:
        return Path(rdxdir).absolute()
    _table = Path(table)
    if not _table.is_file():
        return None
    meta = read_metadata(_table).meta
    if 'rdxdir' in meta:
        return Path(meta['rdxdir'])
    if 'rawdir' in meta:
        return Path(meta['rawdir']).absolute().parent
    return None


def night(mjd):
    """
    Return the night (UT date, ``YYYY-MM-DD``) of an MJD.

    Parameters
    ----------
    mjd : float
        The MJD.

    Returns
    -------
    str
        The UT date.
    """
    return Time(mjd, format='mjd', scale='utc').to_value('iso', subfmt='date')


def check_single_setup(metadata, metafile):
    """
    Check that a metadata table includes frames from only one setup.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table, excluding commented-out rows; see
        :func:`~pipnick.metadata.table.read_metadata`.
    metafile : str, :class:`pathlib.Path`
        The file the table was read from, used in error messages.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the table has no ``setup`` column, or includes frames from
        more than one setup.
    """
    if 'setup' not in metadata.colnames:
        raise PipnickDataModelError(
            f'{metafile} has no setup column; it was likely built by an older version of '
            'pipnick.  Rebuild it using pipnick_metadata.'
        )
    setups, counts = np.unique(metadata['setup'], return_counts=True)
    if len(setups) > 1:
        definitions = metadata.meta.get('setups', {})
        summary = '\n'.join(
            f'    {s} ({n} frames): {definitions.get(s)}' for s, n in zip(setups, counts)
        )
        raise PipnickDataModelError(
            f'{metafile} includes frames from more than one setup:\n{summary}\nData from only '
            'one setup can be reduced at a time.  Comment out the rows for the other setups, '
            'or use pipnick_metadata --split to write one file per setup.'
        )


def check_comb_groups(metadata, metafile):
    """
    Check that the combination groups of a metadata table can be reduced.

    The groups are set by the ``comb_id`` column.  So that the calibrations
    used for each frame are unambiguous, the groups must follow these rules:

        - all frames in a group have the same frame type;
        - there is only one group of biases;
        - the frames in a group of darks have the same exposure time (to
          within 0.1 s), and there is only one group of darks for each
          exposure time;
        - the frames in a group of dome or sky flats have the same filter, and
          there is only one group of each flat type for each filter;
        - the frames in a group of science frames have the same filter.

    Science frames in a group with different exposure times are allowed, but a
    warning is issued.  Frames with other types (e.g., focus or pointing
    frames) are not reduced, and are only noted.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table, excluding commented-out rows; see
        :func:`~pipnick.metadata.table.read_metadata`.
    metafile : str, :class:`pathlib.Path`
        The file the table was read from, used in error messages.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the table has no ``comb_id`` column, or if any of the rules
        above are broken.  All broken rules are listed.
    """
    if 'comb_id' not in metadata.colnames:
        raise PipnickDataModelError(
            f'{metafile} has no comb_id column.  The reduction requires a per-setup metadata '
            'file, which defines the frames to combine; use pipnick_metadata --split to write '
            'one.'
        )
    violations = []
    bias_groups = []
    darks = {}
    flats = {}
    ignored = []
    for comb_id, indx in _group_indices(metadata).items():
        files = ', '.join(metadata['file'][indx])
        frametypes = sorted(set(metadata['frametype'][indx]))
        if len(frametypes) > 1:
            violations.append(
                f'Group {comb_id} mixes frame types ({", ".join(frametypes)}): {files}'
            )
            continue
        frametype = frametypes[0]
        filters = sorted(set(_filter(metadata, i) for i in indx))
        exptimes = set(_exptime(metadata, i) for i in indx)
        exptimes = set(None if e is None else round(e, 1) for e in exptimes)
        if frametype == 'bias':
            bias_groups.append(comb_id)
        elif frametype == 'dark':
            if len(exptimes) > 1:
                violations.append(
                    f'Dark group {comb_id} mixes exposure times '
                    f'({", ".join(sorted(str(e) for e in exptimes))}): {files}'
                )
            elif None in exptimes:
                violations.append(f'Dark group {comb_id} has no exposure time: {files}')
            else:
                darks.setdefault(exptimes.pop(), []).append(comb_id)
        elif frametype in ['domeflat', 'skyflat']:
            if len(filters) > 1:
                violations.append(
                    f'{frametype} group {comb_id} mixes filters ({", ".join(filters)}): {files}'
                )
            else:
                flats.setdefault((frametype, filters[0]), []).append(comb_id)
        elif frametype == 'science':
            if len(filters) > 1:
                violations.append(
                    f'Science group {comb_id} mixes filters ({", ".join(filters)}): {files}'
                )
            if len(exptimes) > 1:
                log.warning(
                    f'Science group {comb_id} mixes exposure times; the frames will be scaled '
                    f'to the median exposure time before they are combined: {files}'
                )
        else:
            ignored.extend(metadata['file'][indx])

    if len(bias_groups) > 1:
        violations.append(
            f'There is more than one bias group ({", ".join(str(c) for c in bias_groups)})'
        )
    for exptime, ids in darks.items():
        if len(ids) > 1:
            violations.append(
                f'There is more than one group of {exptime:g} s darks '
                f'({", ".join(str(c) for c in ids)})'
            )
    for (frametype, filt), ids in flats.items():
        if len(ids) > 1:
            violations.append(
                f'There is more than one group of {filt} {frametype} frames '
                f'({", ".join(str(c) for c in ids)})'
            )
    if len(ignored) > 0:
        log.info(f'These frames have types that are not reduced: {", ".join(ignored)}')
    if len(violations) > 0:
        raise PipnickDataModelError(
            f'{metafile} has invalid combination groups (comb_id):\n'
            + '\n'.join(f'    {v}' for v in violations)
        )


def check_mosaic_groups(metadata, metafile):
    """
    Check the mosaic labels of a per-setup metadata table.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table; see :meth:`DataSet.from_file`.
    metafile : :class:`pathlib.Path`
        The metadata file, used in messages.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if the frames of one combination group (``comb_id``) have
        different labels, or one mosaic has several WCS groups.
    """
    if 'mosaic' not in metadata.colnames:
        return
    errors = []
    for comb_id in np.unique(metadata['comb_id']):
        rows = metadata[metadata['comb_id'] == comb_id]
        for name in ['mosaic', 'wcs_group']:
            if len(set(str(v) for v in rows[name])) > 1:
                errors.append(
                    f'the frames of combination group {comb_id} have different {name} labels'
                )
    for label in set(str(v) for v in metadata['mosaic']) - {''}:
        groups = set(str(v) for v in metadata['wcs_group'][metadata['mosaic'] == label])
        if len(groups) > 1:
            errors.append(f'mosaic {label} has several WCS groups: {", ".join(sorted(groups))}')
    if len(errors) > 0:
        raise PipnickDataModelError(f'Invalid mosaics in {metafile.name}: ' + '; '.join(errors))


def _masked_to_none(value):
    """
    Return None if a table value is masked, or the value otherwise.
    """
    return None if value is np.ma.masked else value


def _filter(metadata, i):
    """
    Return the filter of a frame in the metadata table, or an empty string if
    it is unknown.
    """
    value = _masked_to_none(metadata['filter'][i])
    return '' if value is None else str(value)


def _exptime(metadata, i):
    """
    Return the exposure time of a frame in the metadata table, or None if it is
    unknown.
    """
    value = _masked_to_none(metadata['exptime'][i])
    return None if value is None else float(value)


def _group_indices(metadata):
    """
    Find the rows of the metadata table in each combination group.

    Returns
    -------
    dict
        The row indices of each group, keyed by the ``comb_id`` and ordered by
        it.  The indices are sorted by the MJD of the frame (frames without an
        MJD are last), and then by the row order.
    """
    mjd = np.ma.filled(metadata['mjd'], np.inf)
    srt = np.lexsort((np.arange(len(metadata)), mjd))
    comb_id = np.asarray(metadata['comb_id'])
    return {int(c): srt[comb_id[srt] == c] for c in np.unique(comb_id)}


def _dark_name(exptime):
    """
    Return the name of the combined dark for an exposure time, which is rounded
    to 0.1 s; e.g., ``Dark_120s.fits`` or ``Dark_4.3s.fits``.
    """
    return f'Dark_{round(exptime, 1):g}s.fits'


def _flat_name(frametype, filt):
    """
    Return the name of the combined flat of a given type (``'skyflat'`` or
    ``'domeflat'``) and filter; e.g., ``SkyFlat_V.fits`` or
    ``SkyFlat_ip.fits`` for the ``i'`` filter (see :func:`~pipnick.metadata.names.filter_name`).
    """
    title = 'SkyFlat' if frametype == 'skyflat' else 'DomeFlat'
    return f'{title}_{filter_name(filt)}.fits'


def _group_object(metadata, indx):
    """
    Return the most common object name in a combination group.

    Ties go to the earliest frame.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table.
    indx : :class:`numpy.ndarray`
        The rows in the group, sorted by MJD; see :func:`_group_indices`.

    Returns
    -------
    str, None
        The object name, as given in the table, or None if it is unknown.
    """
    objects = [_masked_to_none(metadata['object'][i]) for i in indx]
    # max returns the first of the most common objects, i.e., the earliest
    return max(objects, key=objects.count)


def _table_values(metadata, indx, frametype):
    """
    Return the metadata table values to write to the header of an output.

    The table is the authority for these values, so that edits to it (e.g., a
    corrected object name) propagate to the reduced images.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table.
    indx : :class:`numpy.ndarray`
        The rows in the group, sorted by MJD; see :func:`_group_indices`.
    frametype : str
        The frame type of the group.

    Returns
    -------
    dict
        The values keyed by the quantity (``'filter'`` and ``'exptime'`` for
        all frame types, and ``'object'``, ``'ra'``, ``'dec'``, and
        ``'airmass'`` for science frames).  Unknown values are None.  The
        object is the most common in the group (see :func:`_group_object`),
        the exposure time is the median, the airmass is the mean (of the
        frames with an airmass), and the other values are from the earliest
        frame.
    """
    exptimes = [_exptime(metadata, i) for i in indx]
    values = {
        'filter': _masked_to_none(metadata['filter'][indx[0]]),
        'exptime': None if None in exptimes else float(np.median(exptimes))
    }
    if frametype == 'science':
        values['object'] = _group_object(metadata, indx)
        airmass = [_masked_to_none(metadata['airmass'][i]) for i in indx]
        airmass = [float(a) for a in airmass if a is not None]
        values['airmass'] = None if len(airmass) == 0 else float(np.mean(airmass))
        values['ra'] = _masked_to_none(metadata['ra'][indx[0]])
        values['dec'] = _masked_to_none(metadata['dec'][indx[0]])
    return values


def _science_name(metadata, indx, comb_id):
    """
    Return the name of the output science image for a combination group.

    The name is ``<object>_<filter>_C<comb_id>_<date>.fits``, where
    ``<object>`` is the most common object name in the group (ties go to the
    earliest frame), ``<filter>`` is the filter of the earliest frame (all
    frames in a group have the same filter; see :func:`check_comb_groups` and
    :func:`~pipnick.metadata.names.filter_name`), ``<comb_id>`` is zero-padded to three digits, and
    ``<date>`` is the UT date of the earliest frame.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table.
    indx : :class:`numpy.ndarray`
        The rows in the group, sorted by MJD; see :func:`_group_indices`.
    comb_id : int
        The combination group.

    Returns
    -------
    str
        The file name.
    """
    obj = _group_object(metadata, indx)
    obj = '' if obj is None else safe_name(obj)
    date = _masked_to_none(metadata['date'][indx[0]])
    date = 'nodate' if date is None else str(date).split('T')[0]
    filt = filter_name(_filter(metadata, indx[0]))
    return f'{"unknown" if len(obj) == 0 else obj}_{filt}_C{comb_id:03d}_{date}.fits'


def _plan_products(metadata, rdxdir):
    """
    Set the output file for each combination group to be reduced.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        The metadata table, with valid combination groups; see
        :func:`check_comb_groups`.
    rdxdir : :class:`pathlib.Path`
        The reduction directory.

    Returns
    -------
    list
        A :class:`Product` for each group to be reduced, ordered by
        ``comb_id``.

    Raises
    ------
    :exc:`~pipnick.pkg.exceptions.PipnickDataModelError`
        Raised if two groups would have the same output file, e.g., for flats
        with filter names that only differ in characters that cannot be used
        in a file name.
    """
    products = []
    for comb_id, indx in _group_indices(metadata).items():
        frametype = str(metadata['frametype'][indx[0]])
        filt = _filter(metadata, indx[0])
        exptime = _exptime(metadata, indx[0])
        if frametype == 'bias':
            name = 'Bias.fits'
        elif frametype == 'dark':
            name = _dark_name(exptime)
        elif frametype in ['domeflat', 'skyflat']:
            name = _flat_name(frametype, filt)
        elif frametype == 'science':
            name = _science_name(metadata, indx, comb_id)
        else:
            continue
        products.append(
            Product(
                frametype=frametype, comb_id=comb_id, indx=indx, filter=filt, exptime=exptime,
                file=rdxdir / name
            )
        )
    names = [p.name for p in products]
    duplicates = sorted(set(n for n in names if names.count(n) > 1))
    if len(duplicates) > 0:
        raise PipnickDataModelError(
            f'More than one combination group would be written to: {", ".join(duplicates)}.'
        )
    return products
