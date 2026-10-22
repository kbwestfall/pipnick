
#################################################################
########  Locate sources and create "stamps" to analyze  ########
#################################################################

from IPython import embed

import numpy as np
from matplotlib import pyplot, ticker
from matplotlib.backends.backend_pdf import PdfPages

from astropy import table
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.nddata import CCDData
from photutils.detection import DAOStarFinder
from skimage import transform

from pipnick.utils.fits_class import Fits_Simple
from pipnick import log


def plot_stamps(stamps, pdf_file, vmin=None, vmax=None):
    """
    Plot images of the source stamps in a PDF file.

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
            w, h = pyplot.figaspect(1)
            fig = pyplot.figure(figsize=(1.5*w, 1.5*h))
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
            pdf.savefig()  # Save the figure to the PDF file
            fig.clear()
            pyplot.close()


#def generate_stamps(images, output_base, thresh=15.):
#    """
#    Generate and analyze stamps from the given images.
#
#    Args:
#        images (list): List of image file paths.
#        output_base (Path): Base of path to files for storing stamp data
#        thresh (float): Star detection threshold (threshold = thresh*std)
#    """
#    num_images = len(images)
#    
#    # Load images and create masked arrays
#    images = [Fits_Simple(image) for image in images]
#    masked_images = np.ma.array([image.masked_array for image in images])
#    
#    ofits = output_base.with_suffix('.rdx.fits')
#    src_ofile = output_base.with_suffix('.src.db')
#
#    # Setting up stamp & source detection parameters
#    default_fwhm = 7.0
#    stamp_width = int(np.ceil(5 * default_fwhm))
#    if stamp_width % 2 == 0:
#        stamp_width += 1
#    stamp_shape = (stamp_width, stamp_width)
#    x, y = np.meshgrid(np.arange(stamp_width, dtype=float) - stamp_width // 2, 
#                       np.arange(stamp_width, dtype=float) - stamp_width // 2)
#    r = np.sqrt(x**2 + y**2)
#
#    mean_rkron = np.zeros(num_images, dtype=float)  # Mean Kron radius for each image
#    source_data = np.zeros((0, 9), dtype=float)  # Array to hold source data
#    all_stamps = np.empty(num_images, dtype=object)  # Array to hold all stamps
#
#    # Detect sources and generate stamps
#    for index in range(num_images):
#        # Source detection
#        mean, median, std = sigma_clipped_stats(masked_images[index], sigma=3.)
#        starfind = DAOStarFinder(fwhm=default_fwhm, threshold=thresh*std)
#        sources = starfind(masked_images[index].filled(0.0) - median,
#                           mask=np.ma.getmaskarray(masked_images[index]))
#        if sources is None:
#            continue
#        nsources = len(sources)
#        # Source positions
#        ap_centers = np.column_stack((sources['xcentroid'], sources['ycentroid']))
#
#        # Stamp setup
#            # Scikit-Image uses the bottom left corner to define the coordinate system.  
#            # Photutils uses the center of the pixel. ...  I think.  In any case, there 
#            # seems to be an 0.5 pixel offset, which is why the transformation below 
#            # has that additional offset
#        center = np.atleast_1d(stamp_shape) / 2  # Center of the stamp
#        stamps = np.ma.zeros((nsources,) + stamp_shape, dtype=float)  # Array to hold the stamps
#        rkron = np.zeros(nsources, dtype=float)  # Array to hold Kron radii
#        kflux_raw = np.zeros(nsources, dtype=float)  # Array to hold raw flux
#        bkg = np.zeros(nsources, dtype=float)  # Array to hold background values
#
#        # Process each detected source and add to stamp
#        for i in range(nsources):
#            t = transform.EuclideanTransform(translation=-ap_centers[i][::-1] - 0.5 + center)
#            stamps[i] = transform.warp(masked_images[index].filled(0.0).T, t.inverse,
#                                       output_shape=stamp_shape, cval=0, order=0).T
#            stamps[i, stamps[i] == 0.] = np.ma.masked  # Mask the background
#            if np.all(np.ma.getmaskarray(stamps[i, r > 3 * default_fwhm])):
#                continue
#            bkg[i] = np.ma.median(stamps[i, r > 3 * default_fwhm])
#            stamps[i] -= bkg[i]  # Subtract the background
#            if np.all(np.ma.getmaskarray(stamps[i])):
#                continue
#            _rkron = np.ma.sum(r * stamps[i]) / np.ma.sum(stamps[i])
#            if _rkron == np.ma.masked or _rkron < 0:
#                continue
#            rkron[i] = _rkron  # Calculate the Kron radius
#            indx = r < rkron[i]
#            if not np.any(indx):
#                continue
#            kflux_raw[i] = np.ma.sum(stamps[i, indx])  # Calculate the raw flux
#
#        mean_rkron[index] = sigma_clipped_stats(rkron, sigma=3., maxiters=None)[0]
#        kflux = np.array([np.sum(s[r < mean_rkron[index]]) for s in stamps])
#
#        # Remove sources that are too close to each other
#        i, j = np.triu_indices(nsources, k=1)
#        dist = np.sqrt((sources['xcentroid'][j] - sources['xcentroid'][i])**2 + 
#                       (sources['ycentroid'][j] - sources['ycentroid'][i])**2)
#        too_close = np.identity(nsources).astype(bool)
#        min_dist_fwhm = 5
#        too_close[i, j] = dist < min_dist_fwhm * default_fwhm
#        too_close[j, i] = dist < min_dist_fwhm * default_fwhm
#
#        # Remove sources with low flux
#        min_flux = 1e3
#        _kflux = np.ma.MaskedArray(np.tile(kflux, (nsources, 1)), mask=np.logical_not(too_close))
#        keep = ((np.ma.argmax(_kflux, axis=1) == np.arange(nsources)) & 
#                (kflux > min_flux) & 
#                (np.sum(stamps.mask, axis=(1, 2)) == 0))
#        nkeep = np.sum(keep)
#        if nkeep == 0:
#            all_stamps[index] = None
#            continue
#
#        srt = np.argsort(kflux[keep])[::-1]
#        all_stamps[index] = stamps[keep][srt]
#        source_data = np.append(source_data,
#                                np.column_stack((np.full(nkeep, index+1, dtype=float),
#                                                    np.arange(nkeep, dtype=float),
#                                                    sources['xcentroid'][keep][srt],
#                                                    sources['ycentroid'][keep][srt],
#                                                    bkg[keep][srt],
#                                                    rkron[keep][srt],
#                                                    kflux_raw[keep][srt],
#                                                    kflux[keep][srt],
#                                                    np.ones(nkeep, dtype=float))),
#                                axis=0)
#        
#        log.info(f'Working on image {images[index]} ({len(stamps[keep])} stamps)')
#
#        # Plot and save the stamps to a PDF file
#        pdf_file = output_base.with_suffix(f".{images[index].path.stem.split('_')[0]}.stamps.pdf")
#        show_stamps(all_stamps[index], pdf_file)
#
#    # Save the source data to a text file
#    header_base = 'Mean per-index Kron Radii:\n'
#    for index in range(num_images):
#        header_base += f'Image {images[index]}: {mean_rkron[index]:.2f}\n'
#    np.savetxt(str(src_ofile), source_data,
#                  fmt='%6.0f %4.0f %7.2f %7.2f %8.2f %7.2f %11.4e %11.4e %2.0f',
#                  header=header_base + f"{'CHIP':>4} {'ID':>4} {'X':>7} {'Y':>7} {'BKG':>8} "
#                                       f"{'R':>7} {'C_RAW':>11} {'C':>11} {'F':>2}")
#
#    # Save the image data and stamps to a FITS file
#    hdr = fits.Header()
#    stamp_hdus = []
#    for index in range(num_images):
#        hdr[f'MRK_{index + 1:>02}'] = mean_rkron[index]
#        if all_stamps[index] is None:
#            continue
#        stamp_hdus += [fits.ImageHDU(data=all_stamps[index].data, name=f'STAMPS_{index + 1:>02}')]
#    fits.HDUList([fits.PrimaryHDU(header=hdr),
#                  fits.ImageHDU(data=masked_images.data, name='IMAGES'),
#                  fits.ImageHDU(data=np.ma.getmaskarray(masked_images).astype(np.uint8),
#                                name='MASKS')]
#                  + stamp_hdus).writeto(str(ofits), overwrite=True)
#    
#    return source_data


# TODO: We might want to instead use photutils.segmentation:
# https://photutils.readthedocs.io/en/latest/reference/segmentation_api.html
def detect_sources(
    image, thresh=15., fwhm=7., inner_bkg=3., outer_bkg=5., min_dist_fwhm=5., min_flux=1e3
):
    """
    Perform point-source detection on the provided image.

    Parameters
    ----------
    image : str, :class:`pathlib.Path`, :class:`astropy.nddata.CCDData`
        Image to process
    thresh : float, optional
        Detection threshold in units of the standard deviation in the
        background.
    fwhm : float, optional
        Rough FWHM of the point-spread function in pixels
    inner_bkg : float, optional
        Inner radius used for determination of the local background near
        sources, in units of the provided FWHM.
    outer_bkg : float, optional
        Outer radius used for determination of the local background near
        sources, in units of the provided FWHM.
    min_dist_fwhm : float, optional
        Minimum distance in pixels for "isolated" sources.
    min_flux : float, optional
        Minimum flux for a well-detected source (takes precedence over ``thresh``)

    Returns
    -------
    :class:`astropy.table.Table`

        Table with detected sources and initial metrics. The table metadata
        provides all of the input parameters (except for the image), as well as
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
    """
    # Read the image, and convert CCDData to a masked array
#    _image = np.asanyarray(image if isinstance(image, CCDData) else CCDData.read(image))
    _image = image if isinstance(image, CCDData) else CCDData.read(image)
    # Create images for coordinates
    nx, ny = _image.shape
    x, y = np.meshgrid(np.arange(ny, dtype=float), np.arange(nx, dtype=float))

    # Determine the clipped statistics.
    # TODO:
    #   - Add clipping parameters to function call?
    #   - Preferrentially clip high values?
    # NOTE: sigma_clipped_stats can work directly on _image (at least for
    # astropy 7.0.1)
    mean, median, std = sigma_clipped_stats(_image, sigma=3.)
    log.info(f'Median and standard deviation: {median:.1f} {std:.2f}')

    # Source detection
    starfind = DAOStarFinder(fwhm=fwhm, threshold=thresh*std)
    sources = starfind((_image - median).filled(0.0), mask=np.ma.getmaskarray(_image))
    nsources = len(sources)
    log.info(f'Detected {nsources} sources.')

    # Get the size and flux of each source via Kron measurements
    ap_centers = np.column_stack((sources['xcentroid'], sources['ycentroid']))  # Source positions
    kron_bkg_raw, kron_radius, kron_flux_raw = kron_stats(
        x, y, np.asanyarray(_image), ap_centers, fwhm, inner_bkg=inner_bkg, outer_bkg=outer_bkg
    )
    mean_kron_radius = sigma_clipped_stats(kron_radius, sigma=3., maxiters=None)[0]
    log.info(f'Mean Kron radius of sources: {mean_kron_radius:.1f}')
    kron_bkg, _, kron_flux = kron_stats(
        x, y, np.asanyarray(_image), ap_centers, fwhm, inner_bkg=inner_bkg, outer_bkg=outer_bkg,
        source_radius=mean_kron_radius
    )

    # Remove sources that are too close to each other
    i, j = np.triu_indices(nsources, k=1)
    dist = np.sqrt(
        (sources['xcentroid'][j] - sources['xcentroid'][i])**2
        + (sources['ycentroid'][j] - sources['ycentroid'][i])**2
    )
    too_close = np.identity(nsources).astype(bool)
    too_close[i,j] = dist < min_dist_fwhm * fwhm
#    too_close[j, i] = dist < min_dist_fwhm * fwhm
    too_close[j,i] = too_close[i,j]

    # Remove sources with low flux
    _kflux = np.tile(kron_flux, (nsources, 1))
    _kflux[np.logical_not(too_close)] = np.ma.masked
    good = (
        (np.ma.argmax(_kflux, axis=1) == np.arange(nsources))
        & (kron_flux.data > min_flux)
        & np.logical_not(np.ma.getmaskarray(kron_flux))
    )
    log.info(f'Selected {np.sum(good)} well-separated sources with good flux')

    meta = {
        'thresh': thresh, 'fwhm': fwhm, 'inner_bkg': inner_bkg, 'outer_bkg': outer_bkg,
        'min_dist_fwhm': min_dist_fwhm, 'min_flux': min_flux, 'mean_kron_radius': mean_kron_radius
    }

    # Construct and return the table
    return table.Table(
        [
            table.Column(name='id', data=np.arange(nsources, dtype=int)),
            table.Column(name='xcen', data=sources['xcentroid']),
            table.Column(name='ycen', data=sources['ycentroid']),
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
    image : :class:`astropy.nddata.CCDData`
        Image with sources
    source_data : :class:`astropy.table.Table`
        Table with the source data (see :func:`detect_sources`)
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
    _image = np.asanyarray(image if isinstance(image, CCDData) else CCDData.read(image))
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

