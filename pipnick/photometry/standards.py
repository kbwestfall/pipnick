"""
The standard stars in the images of a night.

The sources of each image are matched to the Landolt standards (see
:func:`~pipnick.photometry.landolt.match_landolt`), and the measurements of
the standards are collected in a table used to fit the calibration; see
:mod:`pipnick.photometry.nightfit`.
"""

import numpy as np

from astropy.table import MaskedColumn, Table
from astropy.time import Time

from pipnick.photometry.landolt import catalog_color, landolt_band, match_landolt


def match_standards(images, catalog, settings):
    """
    Find the Landolt standards in the source tables of a night.

    Parameters
    ----------
    images : list
        The images, as returned by :class:`~pipnick.datamodel.image.ReducedImage`.
    catalog : :class:`astropy.table.Table`
        The Landolt catalog; see :func:`~pipnick.photometry.landolt.load_landolt`.
    settings : dict
        The calibration settings, including the color used to calibrate each
        band (``colors``); see :func:`~pipnick.pipelines.calibrate.calibrate_all`.

    Returns
    -------
    :class:`astropy.table.Table`
        One row per matched standard and image, with the image, band, star,
        its number of Landolt observations, the match separation (arcsec), the
        airmass and MJD, the catalog magnitude and color and their errors, the
        instrumental magnitude and its error, the number of residual
        detections near the standard (``n_residual``; see
        :meth:`~pipnick.datamodel.catalog.SourceCatalog.neighbors`),
        the measurement error used to weight the fit (``err``), and the reason
        the measurement is not used (``'-'`` if it is usable).
    """
    rows = []
    for image in images:
        phot = image.photometry()
        ra = np.ma.asarray(phot['ra'], dtype=float)
        dec = np.ma.asarray(phot['dec'], dtype=float)
        good = np.logical_not(np.ma.getmaskarray(ra) | np.ma.getmaskarray(dec))
        epoch = Time(image.mjd, format='mjd', scale='utc').jyear
        isrc, icat, sep = match_landolt(
            catalog, ra[good].filled(0.), dec[good].filled(0.), epoch,
            radius=settings['match_radius']
        )
        isrc = np.where(good)[0][isrc]
        band = landolt_band(image.filter)
        mag_i, err_i = phot.instrumental_mag(settings['mag'])
        cat = catalog[icat]
        mag_c = np.ma.asarray(cat[band], dtype=float)
        err_c = np.ma.asarray(cat[f'e_{band}'], dtype=float)
        color, color_err = catalog_color(cat, settings['colors'][band])
        flags = np.asarray(phot['flags'])
        aper_masked = np.ma.getmaskarray(np.ma.asarray(phot['flux_aper']))
        for j, (s, c) in enumerate(zip(isrc, icat)):
            n_residual, n_contaminating = phot.neighbors(
                s, residual_radius=settings['residual_radius'],
                max_neighbor_flux=settings['max_neighbor_flux']
            )
            if np.ma.getmaskarray(mag_c)[j]:
                reason = 'no catalog magnitude'
            elif np.ma.getmaskarray(color)[j]:
                reason = 'no catalog color'
            elif np.ma.getmaskarray(mag_i)[s]:
                reason = 'masked magnitude'
            elif flags[s] != 0:
                reason = 'flagged'
            elif aper_masked[s]:
                reason = 'saturated or otherwise masked'
            elif n_contaminating > 0:
                reason = 'contaminated'
            elif (
                not np.ma.getmaskarray(err_c)[j] and err_c[j] > settings['max_cat_err']
            ):
                reason = 'catalog error'
            elif cat['nobs'][j] < settings['min_cat_obs']:
                reason = 'few catalog observations'
            else:
                reason = '-'
            rows.append(
                {
                    'image': image.path.stem, 'band': band, 'star': str(cat['name'][j]),
                    'nobs': int(cat['nobs'][j]), 'source_id': int(phot['id'][s]),
                    'sep': float(sep[j]), 'airmass': float(image.airmass),
                    'mjd': float(image.mjd), 'mag_cat': mag_c[j], 'mag_cat_err': err_c[j],
                    'color': color[j], 'color_err': color_err[j], 'mag_inst': mag_i[s],
                    'mag_inst_err': err_i[s], 'n_residual': n_residual, 'reason': reason
                }
            )
    return _measurement_table(rows, settings['sys_err'])



def _measurement_table(rows, sys_err):
    """
    Construct the table of standard measurements; see
    :func:`~pipnick.photometry.standards.match_standards`.
    """
    names = {
        'image': str, 'band': str, 'star': str, 'nobs': int, 'source_id': int, 'sep': float,
        'airmass': float, 'mjd': float, 'mag_cat': float, 'mag_cat_err': float,
        'color': float, 'color_err': float, 'mag_inst': float, 'mag_inst_err': float,
        'n_residual': int, 'reason': str
    }
    tbl = Table()
    for name, dtype in names.items():
        values = [r[name] for r in rows]
        if dtype is float:
            mask = [v is np.ma.masked for v in values]
            data = np.array([0. if m else float(v) for v, m in zip(values, mask)], dtype=float)
            tbl[name] = MaskedColumn(data, mask=np.array(mask, dtype=bool))
        else:
            tbl[name] = (
                np.array(values, dtype=dtype) if len(values) > 0
                else np.array([], dtype=dtype if dtype is int else 'U1')
            )
    # The error of each measurement of m - M: the instrumental and catalog
    # errors (0 if unknown) and the systematic error, in quadrature
    tbl['err'] = np.sqrt(
        np.ma.filled(tbl['mag_inst_err'], 0.)**2 + np.ma.filled(tbl['mag_cat_err'], 0.)**2
        + sys_err**2
    )
    formatted = [
        'sep', 'airmass', 'mag_cat', 'mag_cat_err', 'color', 'color_err', 'mag_inst',
        'mag_inst_err', 'err'
    ]
    for name in formatted:
        tbl[name].info.format = '.4f'
    tbl['sep'].unit = 'arcsec'
    return tbl

