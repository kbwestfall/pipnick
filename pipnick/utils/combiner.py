"""
Combine images while propagating their uncertainties.

:class:`ccdproc.Combiner` sets the uncertainty of a combined image from the
scatter of the combined values at each pixel (e.g., their standard deviation
divided by the square root of their number), and ignores the uncertainties of
the input images.  For small numbers of images, this estimate is very noisy;
and if the images are not identical (e.g., unregistered images of a field, or
images with different exposure times), it does not represent the noise.

:class:`Combiner` extends :class:`ccdproc.Combiner` to keep the uncertainties
of the input images, and adds a ``propagate_uncertainty`` option to
:meth:`Combiner.average_combine` that propagates them instead, using the
uncertainty propagation of :mod:`astropy.nddata` where possible.

.. note::

    This module depends only on ``numpy``, ``astropy``, and ``ccdproc``, and is
    written to be contributed to ``ccdproc``, where the uncertainties would be
    kept by :class:`ccdproc.Combiner` itself and ``propagate_uncertainty``
    would become an option of its ``average_combine`` method.
"""

import numpy as np
from astropy.nddata import CCDData, StdDevUncertainty
import ccdproc


class Combiner(ccdproc.Combiner):
    """
    A :class:`ccdproc.Combiner` that keeps the uncertainties of the images,
    so that they can be propagated to the combined image.

    Parameters
    ----------
    ccd_iter : list or generator
        The :class:`astropy.nddata.CCDData` objects to combine.  If any of them
        has no uncertainty, the uncertainties cannot be propagated.
    dtype : str or `numpy.dtype` or None, optional
        See :class:`ccdproc.Combiner`.

    Attributes
    ----------
    uncertainty_arr : :class:`numpy.ndarray`, None
        The standard deviation of each pixel of each image, with the same shape
        as ``data_arr``, or None if any image has no uncertainty.
    """
    def __init__(self, ccd_iter, dtype=None):
        ccd_list = list(ccd_iter)
        super().__init__(ccd_list, dtype=dtype)
        if any(ccd.uncertainty is None for ccd in ccd_list):
            self.uncertainty_arr = None
        else:
            self.uncertainty_arr = np.array(
                [
                    np.asarray(ccd.uncertainty.represent_as(StdDevUncertainty).array)
                    for ccd in ccd_list
                ],
                dtype=float
            )

    def _propagated_average_uncertainty(self, scale_to=None):
        r"""
        Propagate the uncertainties of the images to their (weighted) average.

        At each pixel, the uncertainty is

        .. math::

            \sigma = \frac{\sqrt{\sum_i (w_i s_i \sigma_i)^2}}{\sum_i w_i},

        where the sums are over the images that are not masked (or rejected)
        at that pixel, and :math:`w_i`, :math:`s_i`, and :math:`\sigma_i` are
        the weight, scaling, and uncertainty of image :math:`i`.  The
        uncertainty is 0 where all images are masked.  As for all
        uncertainty propagation in :mod:`astropy.nddata`, the images are
        assumed to be uncorrelated.

        Without weights, the uncertainty is propagated by ``astropy``: the
        scaled images and their uncertainties are stacked in a single
        :class:`~astropy.nddata.CCDData` object, masked where the images are
        masked or rejected, and averaged along the stacking axis with
        :meth:`~astropy.nddata.NDArithmeticMixin.mean`.  ``astropy`` does not
        provide a weighted mean, so the formula above is computed directly
        when the images are weighted.

        Parameters
        ----------
        scale_to : float, :class:`numpy.ndarray`, optional
            The scaling used for the combination; see
            :meth:`average_combine`.

        Returns
        -------
        :class:`numpy.ndarray`
            The uncertainty of the average.
        """
        nimg = len(self.data_arr)
        shape = (nimg,) + (1,) * (self.data_arr.ndim - 1)
        if scale_to is not None:
            scale = np.broadcast_to(np.asarray(scale_to, dtype=float), self.data_arr.shape)
        elif self.scaling is not None:
            scale = np.asarray(self.scaling, dtype=float).reshape(shape)
        else:
            scale = np.ones(shape, dtype=float)
        mask = np.ma.getmaskarray(self.data_arr)

        if self.weights is None:
            stack = CCDData(
                np.asarray(self.data_arr.data, dtype=float) * scale, unit=self.unit,
                uncertainty=StdDevUncertainty(np.absolute(scale) * self.uncertainty_arr),
                mask=mask
            )
            uncertainty = stack.mean(axis=0).uncertainty.array
            # Pixels where all images are masked have no uncertainty
            return np.ma.filled(np.ma.asarray(uncertainty, dtype=float), 0.)

        weights = np.asarray(self.weights, dtype=float)
        if weights.ndim == 1:
            weights = weights.reshape(shape)
        weights = np.broadcast_to(weights, self.data_arr.shape) * np.logical_not(mask)
        variance = np.sum((weights * scale * self.uncertainty_arr)**2, axis=0)
        norm = np.sum(weights, axis=0)
        return np.sqrt(variance) / np.where(norm > 0, norm, 1.)

    def average_combine(
        self, scale_func=None, scale_to=None, uncertainty_func=None, sum_func=None,
        propagate_uncertainty=False
    ):
        """
        Average combine the images; see :meth:`ccdproc.Combiner.average_combine`.

        Parameters
        ----------
        scale_func, scale_to, uncertainty_func, sum_func : optional
            See :meth:`ccdproc.Combiner.average_combine`.  If None, the
            ``ccdproc`` defaults are used.
        propagate_uncertainty : bool, optional
            Set the uncertainty of the combined image by propagating the
            uncertainties of the images (see
            :meth:`_propagated_average_uncertainty`), instead of from the
            scatter of the combined values (``uncertainty_func``).

        Returns
        -------
        :class:`astropy.nddata.CCDData`
            The combined image.

        Raises
        ------
        ValueError
            Raised if ``propagate_uncertainty`` is True and the images do not
            all have uncertainties.
        """
        if propagate_uncertainty and self.uncertainty_arr is None:
            raise ValueError('Cannot propagate uncertainties; not all images have one.')
        kwargs = {'scale_func': scale_func, 'scale_to': scale_to}
        if uncertainty_func is not None:
            kwargs['uncertainty_func'] = uncertainty_func
        if sum_func is not None:
            kwargs['sum_func'] = sum_func
        combined = super().average_combine(**kwargs)
        if propagate_uncertainty:
            combined.uncertainty = StdDevUncertainty(
                self._propagated_average_uncertainty(scale_to=scale_to)
            )
        return combined
