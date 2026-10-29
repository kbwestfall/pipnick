"""
Locate well-isolated, well-detected point sources and extract image "stamps"
around them, used to model the point-spread function.
"""

import numpy as np
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure

from astropy import table
from astropy.stats import sigma_clipped_stats
from astropy.nddata import CCDData
from skimage import transform

from pipnick import log
from pipnick.image.detect import isolated
from pipnick.image.io import read_reduced


def as_masked_array(image):
    """
    Return the image data as a masked array.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`, :class:`astropy.nddata.CCDData`, :class:`numpy.ndarray`
        The image: a file, a :class:`~astropy.nddata.CCDData` object (whose
        mask is applied), or an array (masked or not).

    Returns
    -------
    :class:`numpy.ma.MaskedArray`
        The image data and mask.
    """
    if isinstance(image, np.ndarray):
        return np.ma.MaskedArray(image, mask=np.ma.getmaskarray(image))
    if isinstance(image, CCDData):
        mask = np.zeros(image.data.shape, dtype=bool) if image.mask is None else image.mask
        return np.ma.MaskedArray(np.asarray(image.data, dtype=float), mask=mask)
    data, mask, _, _ = read_reduced(image)
    return np.ma.MaskedArray(data, mask=mask)


def plot_stamps(stamps, pdf_file, vmin=None, vmax=None):
    """
    Plot images of the source stamps in a PDF file.

    The figures are constructed without :mod:`matplotlib.pyplot`, so that they
    can be written in any environment.

    Parameters
    ----------
    stamps : :class:`numpy.ndarray`
        Array of stamp images.
    pdf_file : str
        Path to the output PDF file.
    vmin : float, optional
        Minimum data value to map to colormap.
    vmax : float, optional
        Maximum data value to map to colormap.
    """
    dim = 8  # Dimension of the grid for displaying stamps
    start = 0.01  # Starting position for the grid
    buff = 0.01  # Buffer between stamps
    dw = (1-2*start-(dim-1)*buff)/dim  # Width of each stamp
    nstamps = stamps.shape[0]  # Total number of stamps

    _vmin = np.ma.amin(stamps) if vmin is None else vmin
    _vmax = np.ma.amax(stamps) if vmax is None else vmax

    # Create a PDF file to save the stamps
    with PdfPages(pdf_file) as pdf:
        i = 0
        while i < nstamps:
            fig = Figure(figsize=(7.2, 7.2))
            for j in range(dim):
                for k in range(dim):
                    if i < nstamps:
                        # Add a subplot for each stamp
                        ax = fig.add_axes(
                            [start+k*(dw+buff), start+(dim-j-1)*dw+(dim-j-1)*buff, dw, dw]
                        )
                        ax.xaxis.set_major_formatter(ticker.NullFormatter())
                        ax.xaxis.set_visible(False)
                        ax.yaxis.set_major_formatter(ticker.NullFormatter())
                        ax.yaxis.set_visible(False)
                        ax.imshow(
                            stamps[i], origin='lower', interpolation='nearest', vmin=_vmin,
                            vmax=_vmax
                        )
                        ax.text(
                            0.1, 0.9, str(i), ha='center', va='center', transform=ax.transAxes,
                            color='C1'
                        )
                        ax.scatter(
                            [stamps.shape[1]//2], [stamps.shape[2]//2], marker='+', s=50,
                            color='C3', zorder=10, lw=0.5
                        )
                        i += 1
            pdf.savefig(fig)  # Save the figure to the PDF file


def select_psf_stars(
    image, sources, fwhm=7., inner_bkg=3., outer_bkg=5., min_dist_fwhm=5., min_flux=1e3
):
    """
    Measure the Kron radii and fluxes of detected sources, and select the
    well-separated sources with good fluxes, used to model the PSF.

    Parameters
    ----------
    image : str, :class:`~pathlib.Path`, :class:`~astropy.nddata.CCDData`, array-like
        The image; see :func:`as_masked_array`.
    sources : :class:`astropy.table.Table`
        The detected sources; see :func:`~pipnick.image.detect.find_sources`.
    fwhm : float, optional
        Rough FWHM of the point-spread function in pixels
    inner_bkg : float, optional
        Inner radius used for determination of the local background near
        sources, in units of the provided FWHM.
    outer_bkg : float, optional
        Outer radius used for determination of the local background near
        sources, in units of the provided FWHM.
    min_dist_fwhm : float, optional
        Minimum distance, in units of the FWHM, for "isolated" sources; see
        :func:`~pipnick.image.detect.isolated`.
    min_flux : float, optional
        Minimum Kron flux for a well-detected source.

    Returns
    -------
    :class:`astropy.table.Table`

        Table with the sources and their Kron measurements.  The table
        metadata provides the detection threshold (``threshold``, in the
        units of the image), the other input parameters, and
        'mean_kron_radius', the sigma-clipped mean of the Kron radius measured
        for each source.  The table columns are:

            - 'id': running ID number for the sources
            - 'xcen': source X centroid
            - 'ycen': source Y centroid
            - 'kron_bkg_raw': masked array with the background flux for the
              "raw" Kron measurements
            - 'kron_radius': masked array with the Kron radius for source
            - 'kron_flux_raw': masked array with the Kron flux of each source
              measured using "kron_radius"
            - 'kron_bkg': masked array with the background flux for the Kron
              measurements based on the "mean_kron_radius"
            - 'kron_flux': masked array with the Kron flux of each source
              measured using "mean_kron_radius"
            - 'good': a boolean flag selecting well-separated sources with good
              fluxes

    Raises
    ------
    ValueError
        Raised if there are no sources.
    """
    _image = as_masked_array(image)
    # Create images for coordinates
    nx, ny = _image.shape
    x, y = np.meshgrid(np.arange(ny, dtype=float), np.arange(nx, dtype=float))

    nsources = len(sources)
    log.debug(f'Detected {nsources} sources.')
    if nsources == 0:
        raise ValueError('No sources detected.')

    # Get the size and flux of each source via Kron measurements
    ap_centers = np.column_stack((sources['x'], sources['y']))
    kron_bkg_raw, kron_radius, kron_flux_raw = kron_stats(
        x, y, _image, ap_centers, fwhm, inner_bkg=inner_bkg, outer_bkg=outer_bkg
    )
    mean_kron_radius = sigma_clipped_stats(kron_radius, sigma=3., maxiters=None)[0]
    log.debug(f'Mean Kron radius of sources: {mean_kron_radius:.1f} pixels')
    kron_bkg, _, kron_flux = kron_stats(
        x, y, _image, ap_centers, fwhm, inner_bkg=inner_bkg, outer_bkg=outer_bkg,
        source_radius=mean_kron_radius
    )

    # Keep the brightest of sources that are too close to each other, and
    # remove sources with low flux
    good = (
        isolated(
            np.asarray(sources['x']), np.asarray(sources['y']), kron_flux, min_dist_fwhm * fwhm
        )
        & (kron_flux.data > min_flux)
        & np.logical_not(np.ma.getmaskarray(kron_flux))
    )
    log.debug(f'Selected {np.sum(good)} well-separated sources with good flux')

    meta = {
        'threshold': sources.meta.get('threshold'), 'fwhm': fwhm, 'inner_bkg': inner_bkg,
        'outer_bkg': outer_bkg, 'min_dist_fwhm': min_dist_fwhm, 'min_flux': min_flux,
        'mean_kron_radius': mean_kron_radius
    }

    # Construct and return the table
    return table.Table(
        [
            table.Column(name='id', data=np.arange(nsources, dtype=int)),
            table.Column(name='xcen', data=np.asarray(sources['x'])),
            table.Column(name='ycen', data=np.asarray(sources['y'])),
            table.MaskedColumn(
                name='kron_bkg_raw', data=kron_bkg_raw.data, mask=np.ma.getmaskarray(kron_bkg_raw)
            ),
            table.MaskedColumn(
                name='kron_radius', data=kron_radius.data, mask=np.ma.getmaskarray(kron_radius)
            ),
            table.MaskedColumn(
                name='kron_flux_raw', data=kron_flux_raw.data,
                mask=np.ma.getmaskarray(kron_flux_raw)
            ),
            table.MaskedColumn(
                name='kron_bkg', data=kron_bkg.data, mask=np.ma.getmaskarray(kron_bkg)
            ),
            table.MaskedColumn(
                name='kron_flux', data=kron_flux.data, mask=np.ma.getmaskarray(kron_flux)
            ),
            table.Column(name='good', data=good)
        ],
        meta=meta
    )


def kron_stats(x, y, img, cen, fwhm, inner_bkg=3., outer_bkg=5., source_radius=None):
    """
    Measure the Kron radius and flux for each source.

    Parameters
    ----------
    x, y : :class:`numpy.ndarray`
        X and Y coordinates in pixels.
    img : :class:`numpy.ndarray`, :class:`numpy.ma.MaskedArray`
        Image data
    cen : :class:`numpy.ndarray`
        Array with the center coordinates of the detected sources.  Shape must
        be (N,2) for N sources.
    inner_bkg : float, optional
        Inner radius used for determination of the local background near
        sources, in units of the provided FWHM.
    outer_bkg : float, optional
        Outer radius used for determination of the local background near
        sources, in units of the provided FWHM.
    source_radius : bool, optional
        If provided, use this radius as the Kron radius when measuring the flux,
        not the value that is directly measured.

    Returns
    -------
    kron_bkg : :class:`numpy.ndarray`
        The local background near each source.
    kron_rad : :class:`numpy.ndarray`
        The measured Kron radius in pixels for each source.
    kron_flx : :class:`numpy.ndarray`
        The measured Kron flux for each source.
    """

    stamp_width = int(np.ceil(outer_bkg * fwhm))
    nx, ny = img.shape

    nsrc = cen.shape[0]
    kron_bkg = np.ma.zeros(nsrc, dtype=float)
    kron_flx = np.ma.zeros(nsrc, dtype=float)
    kron_rad = np.ma.zeros(nsrc, dtype=float)

    for i in range(nsrc):
        sx = np.clip(int(cen[i,0] - stamp_width // 2), 0, nx)
        sy = np.clip(int(cen[i,1] - stamp_width // 2), 0, ny)
        slc = np.s_[sy:sy+stamp_width, sx:sx+stamp_width]
        if np.all(np.ma.getmaskarray(img)[slc]):
            kron_bkg[i] = np.ma.masked
            kron_flx[i] = np.ma.masked
            kron_rad[i] = np.ma.masked
            continue
        r = np.sqrt((x[slc] - cen[i,0])**2 + (y[slc] - cen[i,1])**2)
        if np.all(np.ma.getmaskarray(img)[slc][r > inner_bkg * fwhm]):
            kron_bkg[i] = np.ma.masked
            kron_flx[i] = np.ma.masked
            kron_rad[i] = np.ma.masked
            continue
        kron_bkg[i] = np.ma.median(img[slc][r > inner_bkg * fwhm])
        kron_rad[i] = np.ma.sum(r * (img[slc]-kron_bkg[i])) / np.ma.sum(img[slc]-kron_bkg[i])
        indx = r < kron_rad[i] if source_radius is None else r < source_radius
        if not np.any(indx):
            kron_rad[i] = np.ma.masked
            kron_flx[i] = np.ma.masked
            continue
        kron_flx[i] = np.ma.sum(img[slc][indx]-kron_bkg[i])
    return kron_bkg, kron_rad, kron_flx


def extract_stamps(image, source_data, inner_bkg=5., window=10., radius_unit='kron'):
    """
    Extract image stamps around the detected sources.

    Parameters
    ----------
    image : str, :class:`~pathlib.Path`, :class:`~astropy.nddata.CCDData`, array-like
        Image with sources; see :func:`as_masked_array`.
    source_data : :class:`astropy.table.Table`
        Table with the source data (see :func:`select_psf_stars`)
    inner_bkg : float, optional
        Inner radius for background determination
    window : float, optional
        Size of the stamp window.
    radius_unit : str, optional
        The unit for ``inner_bkg`` and ``window``.  Must be either ``pixel`` or
        ``kron``.

    Returns
    -------
    bkg : :class:`numpy.ndarray`
        Array with measured background for each stamp
    stamps : :class:`numpy.ndarray`
        A 3D array with the stamps organized along the first axis.

    """
    if radius_unit == 'kron':
        if 'mean_kron_radius' not in source_data.meta:
            raise ValueError(
                'Provided source data table does not provide mean_kron_radius.  '
                'Unit for stamp width cannot be "kron".'
            )
        stamp_width = int(np.ceil(window * source_data.meta['mean_kron_radius']))
        bkg_radius = inner_bkg * source_data.meta['mean_kron_radius']
    elif radius_unit == 'pixel':
        stamp_width = int(window)
        bkg_radius = inner_bkg
    else:
        raise ValueError('window_unit must be "kron" or "pixel"')

    # Force an odd shape for the stamps such that the source center lands in the
    # center of the center pixel.
    if stamp_width % 2 == 0:
        stamp_width += 1
    stamp_shape = (stamp_width, stamp_width)

    # Radius image
    x, y = np.meshgrid(
        np.arange(stamp_width, dtype=float) - stamp_width // 2,
        np.arange(stamp_width, dtype=float) - stamp_width // 2
    )
    is_bkg = np.sqrt(x**2 + y**2) > bkg_radius

    # Read the image
    _image = as_masked_array(image)
    img_median = sigma_clipped_stats(_image, sigma=3.)[1]
    _img_sub = (_image - img_median).filled(0.0)

    nsources = len(source_data)
    ap_centers = np.column_stack((source_data['xcen'], source_data['ycen']))

    # Stamp setup
    # NOTE: Scikit-Image uses the bottom left corner to define the coordinate
    # system.  Photutils uses the center of the pixel, leading to a 0.5 pixel
    # offset, which is why the transformation below has that additional offset.

    center = np.atleast_1d(stamp_shape) / 2  # Center of the stamp
    stamps = np.ma.zeros((nsources,) + stamp_shape, dtype=float)  # Array to hold the stamps
    bkg = np.zeros(nsources, dtype=float)  # Array to hold background values

    # Process each detected source and add to stamp
    for i in range(nsources):
        t = transform.EuclideanTransform(translation=-ap_centers[i][::-1] - 0.5 + center)
        stamps[i] = transform.warp(
            _img_sub.T, t.inverse, output_shape=stamp_shape, cval=0, order=0
        ).T
        stamps[i, stamps[i] == 0.] = np.ma.masked  # Mask the background
        if np.all(np.ma.getmaskarray(stamps[i, is_bkg])):
            stamps[i] = np.ma.masked
            continue
        bkg[i] = np.ma.median(stamps[i, is_bkg])
        stamps[i] -= bkg[i]  # Subtract the background
        if np.all(np.ma.getmaskarray(stamps[i])):
            stamps[i] = np.ma.masked

    return bkg+img_median, stamps
