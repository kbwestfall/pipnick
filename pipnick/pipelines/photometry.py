from pathlib import Path

from IPython import embed

import numpy as np

from astropy.table import Table

from pipnick.photometry.psf_photometry import psf_analysis, consolidate_groups
from pipnick.photometry.aperture_photometry import aperture_analysis
from pipnick.utils import dir_nav
from pipnick import log


def photometry_all(rdxdir=None, raw_table=None, phot_table=None, thresh=8.0, group=False,
                   mode='all', fittype='circ', overwrite=False):
    r"""
    Perform photometric analysis on reduced files in a specified directory.

    Procedure:
    
        - Extract image files from the provided directory.
        - Perform PSF photometry on each image using 
          :func:`~pipnick.photometry.psf_photometry.psf_analysis`.
        - Optionally consolidate detected sources into groups
          using :func:`~pipnick.photometry.psf_photometry.consolidate_groups`.
        - Conduct aperture photometry on the results with
          :func:`~pipnick.photometry.aperture_photometry.aperture_analysis`.
        - Save the source catalogs to CSV files in the output directory.

    Parameters
    ----------
    rdxdir : str, :class:`pathlib.Path`, optional
        Path to directory containing the reduced images for photometry.  Must
        provide this or ``raw_table``.
    raw_table : str, :class:`pathlib.Path`, optional
        File with the metadata table of the raw data.  This is the same as the
        ``table`` argument in :func:`~pipnick.pipelines.reduction.reduce_all`.
        Must provide this or ``rdxdir``.
    phot_table : str, :class:`pathlib.Path`, optional
        Name for the output photometry data table
    thresh : float, optional
        Threshold for source detection, defined as a multiple of the background 
        standard deviation.
    group : bool, optional
        If True, consolidate detected sources that belong to the same group 
        into a single source.
    mode : str, optional
        Mode to run PSF photometry. Use 'all' for a complete analysis or 'new'
        for only analyzing newly detected sources.
    fittype : str, optional
        Type of Moffat fit to use for PSF photometry. Use 'circ' for a circular
        function or 'ellip' for an elliptical function.
    overwrite : bool, optional
        Overwrite any existing files.

    Returns
    -------
    source_catalog_paths : list
        A list of file paths to the generated source catalogs (in CSV format).
        Each entry corresponds to an analyzed image file.
    """
    if rdxdir is None and raw_table is None:
        raise ValueError('Must provide the reduction directory or the raw data table.')

    if raw_table is None:
        # Find all the reduced files in the directory
        framemeta = dir_nav.build_camera_datatable(rdxdir, ext='_rdx.fits', path_key='rdxdir')
        # Change the name of the column
        framemeta['file'].name = 'rdxfile'
        # Set the reduction directory path
        _rdxdir = Path(framemeta.meta['rdxdir']).absolute()
    else:
        # Read the table with the raw metadata
        framemeta = Table.read(raw_table, format='ascii.ecsv')
        # Check that the reduction directory is defined
        if 'rdxdir' not in framemeta.meta:
            raise KeyError('Provided raw data table does not include rdxdir metadata keyword.')
        # Set the reduction directory path
        _rdxdir = Path(framemeta.meta['rdxdir']).absolute()
        # Check that the directory exists
        if not _rdxdir.is_dir():
            raise NotADirectoryError(f'Reduction directory in {raw_table} does not exist.  '
                                     f'Expected to find reduced data in: {_rdxdir}')
        # Define the paths to the raw files
        raw_frames = dir_nav.frame_path(framemeta)
        # Use those to build the names of the reduced files
        rdx_frames = [_rdxdir / dir_nav.default_reduced_file_name(f) for f in raw_frames]
        # Find which of the reduced files exist
        found_frames = [rf.is_file() for rf in rdx_frames]
        # Check if any of the files exist
        if not any(found_frames):
            return ValueError(f'None of the frames in {raw_table} have reduced counterparts in '
                              f'{_rdxdir}.')
        # Add the column with the expected name of the reduced files.
        framemeta['rdxfile'] = [r.name for r in rdx_frames]
        # Downselect to only the rows with frames that exist
        framemeta = framemeta[found_frames]

    # Define the file for the output photometry metadata
    if phot_table is None:
        _phot_table = Path(dir_nav.default_table_name(framemeta.meta['camera'], 'phot')).absolute()
    else:
        _phot_table = Path(phot_table).absolute()
    if _phot_table.is_file() and not overwrite:
        raise ValueError(f'{_phot_table} exists and overwrite=False!')

    log.info(f'Performing photometry on data:')
    log.info(f'    Reduction directory: {_rdxdir}')

    source_catalog_paths = []  # List to hold paths to source catalogs

    rdx_frames = dir_nav.frame_path(framemeta, col_key='rdxfile', path_key='rdxdir')
    
    for frame in rdx_frames:
        # Perform PSF photometry on the file
        log.info(f'Performing PSF analysis for {frame}.')
        psf_par, psf_data = psf_analysis(frame, thresh=thresh, mode=mode, fittype=fittype)
        
        embed()
        exit()
        
        # Consolidate groups of sources if enabled
        if group:
            psf_data = consolidate_groups(psf_data)
        
        # Perform aperture photometry
        all_data = aperture_analysis(frame, psf_data)
        
        # Save result
        all_data.write(output_file, format='csv', overwrite=True)
        source_catalog_paths.append(output_file)
    
    log.info('---- photometry_all() call ended')
    # Return the list of paths to the source catalogs
    return source_catalog_paths