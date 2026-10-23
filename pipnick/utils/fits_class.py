import numpy.ma as ma
import re
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from pathlib import Path
from typing import Union

from pipnick import log
from pipnick import cameras
from pipnick.pkg.exceptions import PipnickDataModelError

class Fits_Simple:
    """
    A simple class to handle FITS files.

    Attributes
    ----------
    path : :class:`pathlib.Path`
        The file path to the FITS image.
    filename : str
        The name of the FITS file.
    filenamebase : str
        The base name of the FITS file (before any underscores).
    _data : :class:`numpy.ndarray`
        Cached image data from the FITS file.
    _mask : :class:`numpy.ndarray`
        Cached mask data from the FITS file.
    """

    def __new__(cls, *args, **kwargs):
        """
        Returns an existing instance if the first argument is a Fits_Simple instance,
        otherwise creates a new instance.

        Parameters
        ----------
        *args : tuple
            Variable length argument list.
        **kwargs : dict
            Arbitrary keyword arguments.

        Returns
        -------
        :class:`~pipnick.utils.fits_class.Fits_Simple`
            An instance of the Fits_Simple class.
        """
        if len(args) > 0 and isinstance(args[0], cls):
            return args[0]
        return super(Fits_Simple, cls).__new__(cls)

    def __init__(self, image_path: Union[str, Path]):
        """
        Initializes the Fits_Simple object with the given image path.

        Parameters
        ----------
        image_path : str, :class:`pathlib.Path`
            The file path to the FITS image. Can be a string or a Path object.
        """
        if isinstance(image_path, Fits_Simple):
            return
        else:
            image_path = Path(image_path)
            self.path = image_path
            self.filename = image_path.name
            self.filenamebase = image_path.name.split('_')[0]
            self._data = None
            self._mask = None

    @property
    def header(self):
        """
        The header information of the FITS file.

        Returns
        -------
        header : :class:`astropy.io.fits.Header`
            The header information from the FITS file.
        """
        with fits.open(self.path) as hdul:
            return hdul[0].header

    @property
    def data(self):
        """
        The data contained in the FITS file.

        Returns
        -------
        data : :class:`numpy.ndarray`
            The image data from the FITS file.
        """
        if self._data is None:
            with fits.open(self.path) as hdul:
                self._data = hdul[0].data
        return self._data

    @data.setter
    def data(self, new_data):
        """
        Setter for the data attribute.

        Parameters
        ----------
        new_data : :class:`numpy.ndarray`
            The new data to replace the existing data in the FITS file.
        """
        self._data = new_data

    @property
    def mask(self):
        """
        Mask for the FITS data. If no mask is present in the FITS file,
        returns a default mask based on the image shape.

        Returns
        -------
        mask : :class:`numpy.ndarray`
            The mask data from the FITS file or a default mask.
        """
        if self._mask is None:
            with fits.open(self.path) as hdul:
                try:
                    self._mask = hdul['MASK'].data
                except KeyError:
                    log.debug(
                        f'No mask in FITS file {self.path.name}; returning the camera '
                        'bad-pixel mask'
                    )
                    self._mask = self.camera_bpm()
        return self._mask

    @property
    def camera(self):
        """
        The camera used to obtain the image, identified from its header.

        Returns
        -------
        camera : :class:`~pipnick.cameras.cameras.NickelCamera` or None
            The camera, or None if it cannot be identified.
        """
        try:
            return cameras.identify_camera(self.header)
        except PipnickDataModelError:
            return None

    def camera_bpm(self):
        """
        The archival bad-pixel mask of the camera used to obtain the image.

        Returns
        -------
        bpm : :class:`numpy.ndarray` or None
            The bad-pixel mask, or None if the camera cannot be identified or
            the image is not a raw frame (e.g., it has been trimmed).
        """
        camera = self.camera
        if camera is None:
            return None
        try:
            return camera.bpm(self.header)
        except PipnickDataModelError:
            return None

    def _header_value(self, quantity, keyword):
        """
        Get a header quantity, using the camera keywords if the camera can be
        identified, or the provided keyword otherwise.
        """
        camera = self.camera
        if camera is not None:
            return camera.header_value(self.header, quantity)
        return self.header.get(keyword)

    @mask.setter
    def mask(self, new_mask):
        """
        Setter for the mask attribute.

        Parameters
        ----------
        new_mask : :class:`numpy.ndarray`
            The new mask to replace the existing mask in the FITS file.

        Raises
        ------
        :exc:`ValueError`
            If the shape of the new mask does not match the shape of the data.
        """
        if new_mask.shape != self.shape:
            raise ValueError(
                f"new_mask must have same shape as data {self.data.shape}. "
                f"new_mask.shape = {new_mask.shape}"
            )
        self._mask = new_mask

    @property
    def masked_array(self):
        """
        Returns a masked array combining data and mask.

        Returns
        -------
        masked_array : :class:`numpy.ma.MaskedArray`
            The masked array combining the FITS data and mask.
        """
        return ma.masked_array(self.data, self.mask)

    @property
    def shape(self):
        """
        The shape (dimensions) of the FITS image.

        Returns
        -------
        shape : tuple
            The shape of the FITS image as a tuple (NAXIS1, NAXIS2).
        """
        header = self.header
        shape = (int(header['NAXIS1']), int(header['NAXIS2']))
        return shape

    @property
    def image_num(self):
        """
        Extracts and returns the image number from the filename.

        Returns
        -------
        image_num : int
            The image number parsed from the filename.
        """
        numbers = re.findall(r'\d+', self.filename)
        return int(''.join(numbers))

    @property
    def object(self):
        """
        The object name from the FITS header.

        Returns
        -------
        object : str or None
            The object name from the FITS header, or None if not found.
        """
        return self._header_value('object', 'OBJECT')

    @property
    def filtnam(self):
        """
        The filter name from the FITS header.

        Returns
        -------
        filtnam : str or None
            The filter name from the FITS header, or None if not found.
        """
        return self._header_value('filter', 'FILTNAM')

    @property
    def exptime(self):
        """
        The exposure time from the FITS header.

        Returns
        -------
        exptime : float or None
            The exposure time from the FITS header, or None if not found.
        """
        return self._header_value('exptime', 'EXPTIME')

    @property
    def airmass(self):
        """
        The airmass value from the FITS header.

        Returns
        -------
        airmass : float or None
            The airmass value from the FITS header, or None if not found.
        """
        return self._header_value('airmass', 'AIRMASS')

    def __str__(self) -> str:
        """
        Returns a string representation of the Fits_Simple object.

        Returns
        -------
        str
            A string containing the filename, object name, and filter.
        """
        return f"{self.filename} ({self.object} - {self.filtnam})"

    def display(self, header: bool = False, data: bool = False):
        """
        Displays the FITS image and optionally prints the header and data.

        Parameters
        ----------
        header : bool, optional
            Whether to print the header information (default is False).
        data : bool, optional
            Whether to print the data (default is False).
        """
        with fits.open(self.path) as hdul:
            if header:
                print("HDU List Info:")
                print(hdul.info())
                print("\nHDU Header")
                print(repr(hdul[0].header))
            if data:
                print(hdul[0].data)

            # Display the image with ZScale interval normalization
            plt.figure(figsize=(8, 6))
            interval = ZScaleInterval()
            vmin, vmax = interval.get_limits(hdul[0].data)
            plt.imshow(hdul[0].data, origin='lower', vmin=vmin, vmax=vmax)
            plt.gcf().set_dpi(300)
            plt.colorbar()
            plt.show()


def unzip_directories(
    path_list: list, output_format: str = 'Fits_Simple', allow_exceptions: bool = False
) -> list:
    """
    Extract all files in a list of directories or files (or their subdirectories) and
    return a list of image objects.

    Parameters
    ----------
    path_list : list
        A list of paths to directories or files that need to be unzipped.
    output_format : str, optional
        The format of the output objects. Can be 'Fits_Simple' (default) or 'Path'.
    allow_exceptions : bool, optional
        Whether to allow and handle exceptions for file-related errors such
        as `KeyError` or `OSError`. Defaults to False.

    Returns
    -------
    list
        A list of image objects in the specified output format. The objects can be
        instances of `Fits_Simple` or `Path`, depending on input `output_format`.
    """
    if output_format == 'Path':
        output = Path
    elif output_format == 'Fits_Simple':
        output = Fits_Simple
    
    images = []
    for elem_path in path_list:
        elem_path = Path(elem_path)
        if elem_path.is_dir():
            for sub_elem_path in sorted(elem_path.iterdir()):
                images += unzip_directories([sub_elem_path], output_format, allow_exceptions)
        elif elem_path.is_file():
            if allow_exceptions:
                try:
                    images.append(output(elem_path))
                except (KeyError, OSError):
                    log.debug(f"unzip_directories() threw KeyError or OSError on {elem_path}")
            else:
                images.append(output(elem_path))
    return images
