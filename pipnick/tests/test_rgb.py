"""
Tests for the three-color images.
"""

import matplotlib.image
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError
from pipnick.scripts.rgb import RGB
from pipnick.utils import rgb

SHAPE = (60, 80)


def image_wcs(dx=0., cd=True):
    """
    Return a WCS, offset by ``dx`` pixels, written with a CD matrix or with
    CDELT and PC.
    """
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [282.77, -6.27]
    wcs.wcs.crpix = [40. - dx, 30.]
    if cd:
        wcs.wcs.cd = 0.333 / 3600 * np.array([[-1., 0.], [0., 1.]])
    else:
        wcs.wcs.cdelt = [-0.333 / 3600, 0.333 / 3600]
        wcs.wcs.pc = np.eye(2)
    return wcs


def write_image(path, level, wcs, shape=SHAPE, masked=None):
    """
    Write an image with a star at (20, 10) on a uniform background, with a
    mask.
    """
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    data = 10. + level * np.exp(-0.5 * ((xx - 20)**2 + (yy - 10)**2) / 2.**2)
    mask = np.zeros(shape, dtype=np.uint8)
    if masked is not None:
        mask[masked] = 1
    hdr = fits.Header() if wcs is None else wcs.to_header()
    fits.HDUList(
        [fits.PrimaryHDU(data=data, header=hdr), fits.ImageHDU(data=mask, name='MASK')]
    ).writeto(path)
    return path


def test_color_image(tmp_path):
    paths = [
        write_image(tmp_path / 'r.fits', 100., image_wcs(), masked=(slice(50, 55), slice(60, 70))),
        write_image(tmp_path / 'g.fits', 200., image_wcs(cd=False)),
        write_image(tmp_path / 'b.fits', 50., image_wcs())
    ]
    images = [rgb.read_image(p) for p in paths]
    rgb.check_same_wcs(paths, [i[2] for i in images], [i[0].shape for i in images])
    for method in ['lupton', 'linear', 'log']:
        with np.errstate(all='raise'):
            color = rgb.color_image(*[i[:2] for i in images], method=method)
        assert color.shape == SHAPE + (3,) and color.dtype == np.uint8, (
            'The color image should be an 8-bit RGB array with the shape of the images'
        )
        assert np.all(color[10,20] > 0), 'The star should be visible in all channels'
        assert color[52,65,0] == 0, 'Masked pixels should be black in their channel'
    with pytest.raises(PipnickParameterError):
        rgb.color_image(*[i[:2] for i in images], method='bad')

    ofile = tmp_path / 'rgb.png'
    RGB.main(RGB.parse_args([str(p) for p in paths] + ['-o', str(ofile)]))
    written = matplotlib.image.imread(ofile)
    assert written.shape[:2] == SHAPE, 'The written image should have the shape of the images'
    star = np.unravel_index(np.argmax(written[...,1]), SHAPE)
    assert star == (SHAPE[0] - 1 - 10, 20), (
        'The first row of the images should be written at the bottom (north up)'
    )


def test_wcs_mismatch(tmp_path):
    good = [write_image(tmp_path / f'{c}.fits', 100., image_wcs()) for c in 'rg']
    shifted = write_image(tmp_path / 'b.fits', 100., image_wcs(dx=1.))
    paths = good + [shifted]
    images = [rgb.read_image(p) for p in paths]
    with pytest.raises(PipnickDataModelError, match='pipnick_mosaic'):
        rgb.check_same_wcs(paths, [i[2] for i in images], [i[0].shape for i in images])
    with pytest.raises(PipnickDataModelError, match='pipnick_mosaic'):
        RGB.main(RGB.parse_args([str(p) for p in paths] + ['-o', str(tmp_path / 'x.png')]))
    assert not (tmp_path / 'x.png').exists(), 'No image should be written'

    other = write_image(tmp_path / 'small.fits', 100., image_wcs(), shape=(50, 80))
    images = [rgb.read_image(p) for p in good + [other]]
    with pytest.raises(PipnickDataModelError, match='different shapes'):
        rgb.check_same_wcs(good + [other], [i[2] for i in images], [i[0].shape for i in images])

    nowcs = write_image(tmp_path / 'nowcs.fits', 100., None)
    images = [rgb.read_image(p) for p in good + [nowcs]]
    with pytest.raises(PipnickDataModelError, match='have no WCS'):
        rgb.check_same_wcs(good + [nowcs], [i[2] for i in images], [i[0].shape for i in images])
