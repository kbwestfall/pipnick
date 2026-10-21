
from pathlib import Path
import warnings

from IPython import embed

import numpy as np

from astropy.nddata import CCDData
from astropy.stats import sigma_clip
from astropy.wcs.wcs import FITSFixedWarning
from astropy import units
import ccdproc

from pipnick import cameras
from pipnick.utils import dir_nav
from pipnick import log


def reduce_all(rawdir=None, table=None, rdxdir=None, overwrite=False, append=False):
    """
    Perform reduction of raw astronomical data frames (overscan subtraction,
    bias subtraction, flat division, cosmic ray masking).

    The processing steps use :func:`ccdproc.subtract_overscan`,
    :func:`ccdproc.trim_image`, :func:`ccdproc.subtract_bias`,
    :func:`ccdproc.flat_correct`, and :func:`ccdproc.cosmicray_lacosmic`.

    Parameters
    ----------
    rawdir : str, :class:`pathlib.Path`, optional
        Path to the parent directory of the raw directory containing the raw
        FITS files to be reduced.  If None, ``table`` *must* be provided.
    table : str, :class:`pathlib.Path`, optional
        A file with the tabulated data to reduce.  If None, this file is
        automatically generated based on the data found and reduced in
        ``rawdir``; see :func:`~pipnick.utils.dir_nav.build_metadata`.
    rdxdir : str, :class:`pathlib.Path`, optional
        Top-level directory for the reduced data.  If None, this is the parent
        directory of ``rawdir``.  I.e., if the raw directory is
        ``/User/janedoe/Nickel/raw``, ``rdxdir`` will be set to
        ``/User/janedoe/Nickel/``.
    overwrite : bool, optional
        Overwrite any existing files.
    append : bool, optional
        If the data-reduction table already exists, add any new science frames
        found in the raw data directory to the provided table.  In effect, this
        simply overwrites the existing ``table`` with a new listing of files in
        ``rawdir``.  The difference with ``overwrite`` is that the *products* of
        the reduction are *not* overwritten.  This means that new raw
        calibration frames are ignored.  To incorporate new calibration files in
        the reduction, you should set ``overwrite=True``.

    Returns
    -------
    list
        Paths to the reduced images written to disk.
    """
    # Organize raw files based on input directory or table
    metafile, metadata = dir_nav.build_metadata(
        rawdir=rawdir, filename=table, overwrite=overwrite or append, rdxdir=rdxdir
    )

    # Get the output directory
    if 'rdxdir' in metadata.meta:
        _rdxdir = Path(metadata.meta['rdxdir'])
    else:
        _rdxdir = Path(metadata.meta['rawdir']).absolute().parent
        metadata.meta['rdxdir'] = str(_rdxdir)
        # Overwrite the file to include the reduction directory
        metadata.write(metafile, format='ascii.ecsv', overwrite=True)
    if not _rdxdir.is_dir():
        _rdxdir.mkdir(parents=True)

    log.info(f'Reducing data:')
    log.info(f'    Metadata file: {metafile}')
    log.info(f'    Reduction directory: {_rdxdir}')

    # Get the super bias
    super_bias = get_super_bias(metadata, rdxdir=metadata.meta['rdxdir'], overwrite=overwrite)
    if super_bias is None:
        log.warning('No bias frames available. Attempting to continue anyway...')

    # Get the super flat.  Try to get sky flats first:
    super_flat = get_super_flats(metadata, rdxdir=metadata.meta['rdxdir'],
                                 flattype='skyflat', bias=super_bias, overwrite=overwrite)
    # If there were no sky flats, try finding some dome flats
    if super_flat is None:
        super_flat = get_super_flats(metadata, rdxdir=metadata.meta['rdxdir'],
                                     flattype='domeflat', bias=super_bias, overwrite=overwrite)
    # Still no flats.  Warn the user and keep going
    if super_flat is None:
        log.warning('No flat frames available.  Attempting to continue anyway...')
    
    # Assume all unknown frame types are on-sky science observations (but this
    # will include pointing and focus frames!)
    is_science = metadata['frametype'] == 'None'
    if np.sum(is_science) == 0:
        log.warning('No science frames available.  Reduction will end.')
        return metadata
    
    # Perform the basic processing of all science frames
    rdx_file = np.empty(len(metadata), dtype=object)
    frame = dir_nav.frame_path(metadata, path_key='rawdir')
    nfiles = len(metadata)
    for i in range(nfiles):
        if not is_science[i]:
            continue

        if super_flat is None or metadata['filter'][i] not in super_flat.keys():
            flat = None
        else:
            flat = super_flat[metadata['filter'][i]]

        # Trim and overscan correct it
        rdx_file[i] = _rdxdir / dir_nav.default_reduced_file_name(frame[i])
        if rdx_file[i].is_file() and not overwrite:
            log.info(f'{rdx_file[i]} exists and overwrite is False.  Continuing...')
            continue

        rdx_frame = process_frame(frame[i], flag_cosmics=True, bias=super_bias, flat=flat)
        rdx_frame.write(rdx_file[i], overwrite=True)
        log.info(f'Saved processed frame to {rdx_file[i]}')

    # Return the metadata table
    return metadata


def init_ccddata(frame, bpm=None, saturation=None, gain=None, readnoise=None):
    """
    Read data from a provided file.
    
    This function also uses the header information to identify the camera used
    to acquire the data.  If the camera is successfully identified, the
    camera information is used to 

        - apply any known bad-pixel mask
        - mask saturated pixels
        - convert the units from DN to counts by multiplying by the camera gain

    The relevant keyword arguments can be used to supplement information not
    available for a given camera or when the camera cannot be identified.  Note
    that values provided via a keyword *always* take precedence.

    Parameters
    ----------
    frame : str, :class:`pathlib.Path`
        Path to the FITS file.

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Initialized and processed CCDData object.
    """
#    warnings.simplefilter('ignore', FITSFixedWarning)
    _frame = Path(frame).absolute()

    log.info(f'Initializing CCD imaging data in {_frame.name}.')
    ccd = CCDData.read(_frame, unit=units.adu)

    # Determine the camera used
    try:
        camera = cameras.identify_camera(ccd.header)
    except:
        warnings.warn(f'Unable to determine camera used for {_frame.name}.')
        camera = None
    else:
        ccd.meta['camera'] = camera
        camera = eval(f'cameras.{camera}')

    # Mask bad pixels
    _bpm = camera.bpm() if bpm is None and camera is not None else bpm
    if _bpm is not None:
        ccd.mask = _bpm

    # Mask saturated pixels; saturation is in ADU
    _sat = camera.saturation if saturation is None and camera is not None else saturation
    if _sat is not None:
        ccd.mask[ccd.data >= _sat] = True

    # Get and save the RN (in electrons) and gain (in electrons/ADU)
    ccd.meta['PROCGAIN'] = camera.gain if gain is None else gain
    ccd.meta['PROCRN'] = camera.readnoise if readnoise is None else readnoise

    # Apply gain if available
    if ccd.meta['PROCGAIN'] is not None:
        ccd.data = ccd.data * ccd.meta['PROCGAIN']
        ccd.unit *= (units.electron / units.adu)

    # Done
    return ccd


def process_frame(raw_frame, bias=None, flat=None, flag_cosmics=True):
    """
    Perform basic image processing

    This includes overscan subtraction, trimming, bias subtraction, and
    field-flattening.    

    Parameters
    ----------
    raw_frame : str, :class:`pathlib.Path`
        Data file to process.
    flag_cosmics : bool, optional
        Identify and mask cosmic rays
    bias : :class:`astropy.nddata.CCDData`, optional
        Image to use for the bias subtraction
    flat : :class:`astropy.nddata.CCDData`, optional
        Image to use for the field-flattening

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Processed CCDData object with overscan subtracted and image trimmed.
    """
#    warnings.simplefilter("error", RuntimeWarning)
    ccd = init_ccddata(raw_frame)
    try:
        oscansec = eval(f'cameras.{ccd.meta["camera"]}').oscansec(ccd.header)
    except:
        log.warning('Could not determine overscan region!')
        oscansec = None
    else:
        log.info('Subtracting overscan and trimming image.')
        ccd = ccdproc.subtract_overscan(ccd, fits_section=oscansec, overscan_axis=1)
        ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['DATASEC'])
    if bias is not None:
        log.info('Subtracting super bias')
        ccd = ccdproc.subtract_bias(ccd, bias)
    if flat is not None:
        log.info('Applying flat-field correction')
        ccd = ccdproc.flat_correct(ccd, flat)
    if flag_cosmics:
        log.info(f'Removing cosmic rays from {raw_frame}')
        ccd = ccdproc.cosmicray_lacosmic(ccd, gain_apply=False, gain=1.0,
                                         readnoise=ccd.meta['PROCRN'], verbose=True)
    return ccd


def stack_frames(raw_frames, scale=False, bias=None, flat=None, flag_cosmics=True):
    """
    Process and stack a list of frames.

    Parameters
    ----------
    raw_frames : list
        List of file names to process and combine.
    scale : bool, optional
        Scale the image data by their (masked) mean before combining them?

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        Combined CCDData object.
    """
    # TODO: For *many* images, this may be too memory intensive.
    rdx_frames = [process_frame(frame, bias=bias, flat=flat, flag_cosmics=flag_cosmics) 
                    for frame in raw_frames]

    # TODO: Need to figure out how to handle the scaling; i.e., how is this
    # passed to the combiner?
    if scale:
        factor = np.array([float(f.mean().data) for f in rdx_frames])
        factor = factor[len(factor)//2] / factor
        rdx_frames = [f.multiply(s) for s,f in zip(factor, rdx_frames)]

    log.info('Identifying outliers using sigma clipping')
    combiner = ccdproc.Combiner(rdx_frames)

    old_n_masked = 0
    new_n_masked = 1
    while new_n_masked > old_n_masked:
#        combiner.sigma_clipping(low_thresh=3, high_thresh=3, func=np.ma.mean)
        # TODO:
        #   - does ccdproc provide an iterative sigma clipper?
        #   - allow user to change the clipping parameters
        combiner.data_arr = sigma_clip(combiner.data_arr, sigma=5., maxiters=1,
                                       cenfunc=np.ma.mean, stdfunc=np.ma.std,
                                       axis=0, masked=True, copy=False)
        old_n_masked = new_n_masked
        new_n_masked = combiner.data_arr.mask.sum()

#    if frame_type == 'flat':
#        combiner.scaling = lambda arr: 1/np.ma.mean(arr)

    log.info('Creating sigma-clipped average image.')
    mean = combiner.average_combine(scale_func=np.ma.mean)

    # TODO: Fill the masked regions?
    return mean


def get_super_bias(metadata, rdxdir=None, overwrite=False):
    """
    Create a combined (super) bias frame from individual bias frames.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        Astropy Table containing file information.
    rdxdir : str, :class:`pathlib.Path`, optional
        Directory to save the super bias frame.
    overwrite : bool, optional
        Overwrite any existing files.

    Returns
    -------
    :class:`astropy.nddata.CCDData`
        SuperBias CCDData object.
    """
    if rdxdir is not None:
        ofile = Path(rdxdir).absolute() / 'Bias.fits'
        if ofile.is_file() and not overwrite:
            log.info(f'{ofile} exists.  Loading existing data and continuing.')
            super_bias = CCDData.read(ofile)
            return super_bias

    # Select biases
    indx = metadata['frametype'] == 'bias'
    nbias = np.sum(indx)
    if nbias == 0:
        return None

    log.info(f'Found {nbias} bias frames to process and combine.')

    frames = dir_nav.frame_path(metadata, select=indx, path_key='rawdir')
    super_bias = stack_frames(frames, flag_cosmics=False)

    if rdxdir is None:
        return super_bias
    
    super_bias.header["OBJECT"] = "Bias"
    super_bias.write(ofile, overwrite=overwrite)
    log.info(f"Saved super bias to {ofile}")
    return super_bias


def get_super_flats(metadata, rdxdir=None, flattype='skyflat', bias=None, filter=None,
                    overwrite=False):
    """
    Create combined (super) flat frames (one per filter) from individual flat
    frames.

    Parameters
    ----------
    metadata : :class:`astropy.table.Table`
        Astropy Table containing file information.
    rdxdir : str, :class:`pathlib.Path`, optional
        Directory to save the super flat frames.
    flattype : str, optional
        Type of flats to use.  Must be 'domeflat' or 'skyflat'.
    bias : :class:`astropy.nddata.CCDData`, optional
        Super bias
    filter : list, optional
        List of strings with the filters to process.  If None, process all
        available filters.
    overwrite : bool, optional
        Overwrite any existing files.

    Returns
    -------
    dict
        Dictionary of super flat CCDData objects keyed by filter.
    """
    indx = metadata['frametype'] == flattype
    if np.sum(indx) == 0:
        log.warning(f'No {flattype} frames found!')
        return None
    
    _rdxdir = None if rdxdir is None else Path(rdxdir).absolute()

    filters = np.unique(metadata['filter'][indx]) if filter is None else filter

    title = 'SkyFlat' if flattype == 'skyflat' else 'DomeFlat'

    super_flat = {}

    for f in filters:
        _indx = indx & (metadata['filter'] == f)
        if np.sum(_indx) == 0:
            continue

        if _rdxdir is not None:
            ofile = _rdxdir / f'{title}_{f}.fits'
            if ofile.is_file() and not overwrite:
                log.info(f'{ofile} exists.  Loading existing data and continuing.')
                super_flat[f] = CCDData.read(ofile)
                continue
    
        log.info(f'Processing and combining {np.sum(_indx)} {f} flats')
        frames = dir_nav.frame_path(metadata, select=_indx, path_key='rawdir')
        super_flat[f] = stack_frames(frames, scale=True, bias=bias, flag_cosmics=False)
        if _rdxdir is None:
            continue

        super_flat[f].header["OBJECT"] = title
        super_flat[f].write(ofile, overwrite=True)
        log.info(f'Saving {f} super {title} to {ofile}')

    return super_flat


