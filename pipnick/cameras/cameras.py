
from pathlib import Path
from IPython import embed

import numpy as np
from astropy.io import fits

from pipnick.utils import masking

__all__ = ['NickelC2', 'NickelAndor', 'identify_camera']

class NickelCamera:
    """
    Class with parameters generic to both Nickel cameras
    """

    object_keys = {
        'bias': 'Bias',
        'domeflat': 'Dome flat',
        'skyflat': ('Sky flat', 'Flat'),
        'dark': 'dark',
        'focus': 'focus',
    }
    """
    Header keyword parameters for 'OBJTYPE' used to identify different frames.
    """

    meta_keys = {
        'object': 'OBJECT',
        'filter': 'FILTNAM',
        'exptime': 'EXPTIME',
        'airmass': 'AIRMASS',
        'ra': 'RA',
        'dec': 'DEC',
        'xbin': 'CDELT1U',
        'ybin': 'CDELT2U',
        'mode': 'READ-SPD'
    }
    """
    Link between generalized metadata terms and header keywords
    """

    ra_dec_units = ('hour', 'degree')
    """
    Format of RA / Dec values in header, for astroquery
    """

    def __init__(self):
        pass

    @property
    def name(self):
        return self.__class__.__name__
    
    @classmethod
    def metadata_cols(cls):
        """
        Return the metadata columns for this class.  See :func:`parse_metadata`.
        """
        return ['file', 'frametype'] + list(cls.meta_keys.keys())
    
    @classmethod
    def parse_metadata(cls, rawfile, keys=None):
        """
        Parse the relevant metadata needed for reduction from a raw file header.

        Parameters
        ----------
        rawfile : str, Path
            Raw file to reduce
        keys : list of str
            A list of metadata keys to (try to) extract from the header.  The
            valid keys are listed by ``cls.meta_keys.keys()``.

        Returns
        -------
        list
            List of the metadata extracted from the file.  See
            :func:`metadata_cols`.  The first two elements are always the file
            name (str) and the automatically determined frame type (str).  The
            remaining metadata can be specific to each camera.
        """
        _file = Path(rawfile).absolute()
        hdr = get_header(_file)

        frametype = 'None'
        for k,v in cls.object_keys.items():
            if (isinstance(v, tuple) and hdr[cls.meta_keys['object']] in v) \
                    or hdr[cls.meta_keys['object']] == v:
                frametype = k
                break

        metadata = [hdr[cls.meta_keys[key]] if cls.meta_keys[key] in hdr else None
                    for key in cls.meta_keys.keys()]
        return [_file.name, frametype] + metadata
    
    @classmethod
    def bpm(cls):
        """
        Use the class internals to construct the bad-pixel mask.
        """
        return masking.create_mask(cls.raw_shape, bad_cols=cls.mask_columns,
                                   bad_slices=cls.mask_slices, bad_regions=cls.mask_regions)



# TODO: Convert these to unbinned values and then account for the binning.
class NickelC2(NickelCamera):
    """
    Class providing metadata for the Nickel UCAM camera, C2
    """

    pixelsize = 15
    """
    Size of the pixel in microns.
    """

    plate_scale_approx = 0.37
    """
    Approximate plate scale in arcsec per pixel
    """

    # For slow readout speed and 2x2 binning
    gain = 1.8
    readnoise = 10.7
    saturation = 65535      # Saturation in ADU
    # nonlinear

    raw_shape = (1024,1056)
    mask_columns = [255, 256, 783, 784, 1002]
    mask_slices = None
    mask_regions = [np.array([[-0.5, 33.5], [34.5, -0.5], [-0.5, -0.5]]),
                    np.array([[-0.5, 960.5], [64.5, 1024.5], [-0.5, 1024.5]])]
    trim_shape = (1024,1024)

    @staticmethod
    def valid_file(hdr):
        """
        Check that the provided file is from this camera.
        """
        return hdr['INSTRUME'] == 'Nickel Direct Camera'
    
    @staticmethod
    def oscansec(hdr):
        nc = hdr['NAXIS1']
        no = hdr['COVER']
        nr = hdr['NAXIS2']
        return f'[{nc-no+1}:{nc},1:{nr}]'
    

class NickelAndor(NickelCamera):
    """
    Class providing metadata for the Nickel Andor Camera
    """

    pixelsize = 13.5
    """
    Size of the pixel in microns.
    """

    plate_scale_approx = 0.33
    """
    Approximate plate scale in arcsec per pixel
    """

    # For slow readout speed and 2x2 binning
    gain = 1.8
    readnoise = 10.7
    saturation = 65535      # Saturation in ADU
    # nonlinear

    raw_shape = (1024,1056)
    mask_columns = [255, 256, 783, 784, 1002]
    mask_slices = None
    mask_regions = [np.array([[-0.5, 33.5], [34.5, -0.5], [-0.5, -0.5]]),
                    np.array([[-0.5, 960.5], [64.5, 1024.5], [-0.5, 1024.5]])]
    trim_shape = (1024,1024)

    @staticmethod
    def valid_file(hdr):
        """
        Check that the provided file is from this camera.
        """
        return 'INSTRUME' not in hdr

    @staticmethod
    def oscansec(hdr):
        nc = hdr['NAXIS1']
        no = hdr['COVER']
        nr = hdr['NAXIS2']
        return f'[{nc-no+1}:{nc},1:{nr}]'

def identify_camera(inp):
    """
    Identify the camera based on the header information.
    """
    hdr = get_header(inp)
    if NickelC2.valid_file(hdr):
        return 'NickelC2'
    if NickelAndor.valid_file(hdr):
        return 'NickelAndor'
    raise ValueError('Unable to determine camera used to obtain data.')


def get_header(inp):
    """
    Convenience function used to get a header from the input.

    Parameters
    ----------
    inp : str, Path, fits.Header
        File or header to read.  If a Header object, the object is simply
        returned.

    Returns
    -------
    hdr : fits.Header
        Fits header
    """
    if isinstance(inp, fits.Header):
        return inp
    return fits.getheader(Path(inp).absolute())

