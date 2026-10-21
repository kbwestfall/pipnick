
#############################################
########  Fit objects to model PSFs  ########
#############################################

from IPython import embed

import numpy as np
from scipy import optimize

from astropy.modeling.functional_models import Moffat2D, Fittable2DModel
from astropy.modeling.parameters import Parameter
from astropy.units import UnitsError
from astropy.stats import sigma_clipped_stats

class FitMoffat2D:
    """
    Fit a 2D Moffat model to the given data.

    There are six model parameters:

        0. ``x0``: x-coordinate of the center
        1. ``y0``: y-coordinate of the center
        2. ``amplitude``: peak value
        3. ``gamma``: width
        4. ``alpha``: fall-off slope
        5. ``background``: background

    Parameters
    ----------
    c : :class:`numpy.ndarray`
        2D array of data to fit.
    """
    def __init__(self, c):
        self.c = c
        self.shape = self.c.shape
        self.x, self.y = np.meshgrid(np.arange(self.shape[0]).astype(float),
                                        np.arange(self.shape[1]).astype(float))
        self.moff = Moffat2D()
        self.par = self.guess_par()
        self.npar = self.par.size
        self.nfree = self.npar
        self.free = np.ones(self.npar, dtype=bool)

    def guess_par(self, image_moments=True, alpha=3.5):
        """
        Guess initial parameters for the Moffat fit based on data.

        Parameters
        ----------
        image_moments : bool, optional
            Use the image moments to set the center and amplitude.  If False,
            these are set using the center of the image.
        alpha : float, optional
            Moffat fall-off slope.

        Returns
        -------
        :class:`numpy.ndarray`
            Initial parameter guesses: (x0, y0, amplitude, gamma, alpha,
            background).
        """
        par = np.zeros(6, dtype=float)
        # Clip high values to guess the background
        par[5] = sigma_clipped_stats(self.c, sigma_lower=3., sigma_upper=1.)[1]
        c_bkg = self.c - par[5]
        if image_moments:
            # Use moments to calculate the center, peak value, and width
            par[0] = np.ma.sum(self.x * c_bkg)/np.ma.sum(c_bkg)
            par[1] = np.ma.sum(self.y * c_bkg)/np.ma.sum(c_bkg)
            # Set the amplitude to be the value at the image center
            par[2] = c_bkg[int(np.round(par[0])), int(np.round(par[1]))]
        else:
            # Set the parameters based on the image center
            par[0] = self.x[self.shape[0]//2, self.shape[1]//2]
            par[1] = self.y[self.shape[0]//2, self.shape[1]//2]
            par[2] = self.c[self.shape[0]//2, self.shape[1]//2]
        # Always set alpha explicitly
        par[4] = alpha
        # Set the width based on the first moment of the radius
        r = np.sqrt((self.x-par[0])**2 + (self.y-par[1])**2)
        fwhm = 2 * np.ma.sum(r * c_bkg)/np.ma.sum(c_bkg)
        par[3] = self.to_gamma(fwhm, par[4])
        return par

    def default_bounds(self):
        """
        Define default bounds for the parameters.

        Returns
        -------
        lb : :class:`numpy.ndarray`
            Lower bounds, ordered as the model parameters.
        ub : :class:`numpy.ndarray`
            Upper bounds, ordered as the model parameters.
        """
        lb = np.zeros(6, dtype=float)
        ub = np.zeros(6, dtype=float)
        lb[0], ub[0] = 0., float(self.shape[0])
        lb[1], ub[1] = 0., float(self.shape[1])
        lb[2], ub[2] = 0., 1.2*np.amax(self.c)
        lb[3], ub[3] = 1., max(self.shape)/2.
        lb[4], ub[4] = 1., 6.
        lb[5] = np.amin(self.c)-0.1*np.absolute(np.amin(self.c))
        ub[5] = 0.1*np.amax(self.c)
        return lb, ub

    def _set_par(self, par):
        """
        Set the parameter values for the model.

        Parameters
        ----------
        par : :class:`numpy.ndarray`
            Array of parameter values: (x0, y0, amplitude, gamma, alpha,
            background).  Can also be only the free parameters.
        """
        if par.ndim != 1:
            raise ValueError('Parameter array must be a 1D vector.')
        if par.size == self.npar:
            self.par = par.copy()
            return
        if par.size != self.nfree:
            raise ValueError(f'Must provide {self.npar} or {self.nfree} parameters.')
        self.par[self.free] = par.copy()

    def model(self, par=None):
        """
        Evaluate the Moffat model with the object's parameters or
        the given parameters.

        Parameters
        ----------
        par : :class:`numpy.ndarray`, optional
            Parameters for the model: (x0, y0, amplitude, gamma, alpha,
            background).  If None, use the current parameters.

        Returns
        -------
        :class:`numpy.ndarray`
            Evaluated model.
        """
        if par is not None:
            self._set_par(par)
        return self.moff.evaluate(self.x, self.y, self.par[2], self.par[0], self.par[1],
                                  self.par[3], self.par[4]) + self.par[5]

    def resid(self, par):
        """
        Calculate the residuals between the data and the model.
        """
        return self.c.ravel() - self.model(par).ravel()

    def deriv_resid(self, par):
        """
        Calculate the derivative of the residuals with respect to the parameters.
        """
        self._set_par(par)
        dmoff = self.moff.fit_deriv(self.x.ravel(), self.y.ravel(), self.par[2], self.par[0],
                                    self.par[1], self.par[3], self.par[4])
        return np.array([-d for d in dmoff]
                            + [np.full(np.prod(self.shape), -1., dtype=float)]).T

    def fit(self, p0=None, fix=None, lb=None, ub=None):
        """
        Fit the Moffat model to the data and modify parameters accordingly.

        Parameters
        ----------
        p0 : :class:`numpy.ndarray`, optional
            Initial parameter guesses.  If None, use :meth:`guess_par`.
        fix : :class:`numpy.ndarray`, optional
            Boolean array indicating which parameters to fix.
        lb : :class:`numpy.ndarray`, optional
            Lower bounds for the parameters.  If None, use
            :meth:`default_bounds`.
        ub : :class:`numpy.ndarray`, optional
            Upper bounds for the parameters.  If None, use
            :meth:`default_bounds`.
        """
        if p0 is None:
            p0 = self.guess_par()
        _p0 = np.atleast_1d(p0)
        if _p0.size != self.npar:
            raise ValueError('Incorrect number of model parameters.')
        self.par = _p0.copy()
        _free = np.ones(self.npar, dtype=bool) if fix is None else np.logical_not(fix)
        if _free.size != self.npar:
            raise ValueError('Incorrect number of model parameter fitting flags.')
        self.free = _free.copy()
        self.nfree = np.sum(self.free)
        
        _lb, _ub = self.default_bounds()
        if lb is None:
            lb = _lb
        if ub is None:
            ub = _ub
        if len(lb) != self.npar or len(ub) != self.npar:
            raise ValueError('Length of one or both of the bounds vectors is incorrect.')

        p = self.par[self.free]
        result = optimize.least_squares(self.resid, p, method='trf', xtol=1e-12,
                                        bounds=(lb[self.free], ub[self.free]),
                                        verbose=0) #, jac=self.deriv_resid)
        self._set_par(result.x)

    @staticmethod
    def to_gamma(fwhm, alpha):
        """
        Convert full-width half-maximum (FWHM) to gamma.
        """
        return fwhm / 2 / np.sqrt(2**(1/alpha)-1)
    
    @staticmethod
    def to_fwhm(gamma, alpha):
        """
        Convert gamma to full-width half-maximum (FWHM).
        """
        return 2 * gamma * np.sqrt(2**(1/alpha)-1)



class FitEllipticalMoffat2D:
    """
    Fit an elliptical 2D Moffat model to the given data.

    There are eight model parameters:

        0. ``x0``: x-coordinate of the center
        1. ``y0``: y-coordinate of the center
        2. ``amplitude``: peak value
        3. ``gamma1``: width in the 1st arbitrary(?) direction
        4. ``gamma2``: width in the 2nd arbitrary(?) direction
        5. ``phi``: rotation angle in radians
        6. ``alpha``: fall-off slope
        7. ``background``: background

    Parameters
    ----------
    c : :class:`numpy.ndarray`
        2D array of data to fit.
    """
    def __init__(self, c):
        self.c = c
        self.shape = self.c.shape
        self.x, self.y = np.meshgrid(np.arange(self.shape[0]).astype(float),
                                        np.arange(self.shape[1]).astype(float))
        self.par = self.guess_par()
        self.npar = self.par.size
        self.nfree = self.npar
        self.free = np.ones(self.npar, dtype=bool)

#    def guess_par(self):
#        """
#        Guess initial parameters for the elliptical Moffat fit based on data
#
#        Returns:
#            ndarray: Initial parameter guesses.
#                     (x0, y0, amplitude, gamma1, gamma2, phi, alpha, background)
#        """
#        par = np.zeros(8, dtype=float)
#        par[0] = self.x[self.shape[0]//2, self.shape[1]//2]
#        par[1] = self.y[self.shape[0]//2, self.shape[1]//2]
#        par[2] = self.c[self.shape[0]//2, self.shape[1]//2]
#        r = np.sqrt((self.x-par[0])**2 + (self.y-par[1])**2)
#        par[3] = np.sum(r*self.c)/np.sum(self.c)/2
#        par[4] = np.sum(r*self.c)/np.sum(self.c)/2
#        par[5] = 0.
#        par[6] = 3.5
#        par[7] = 0.
#        return par

    def guess_par(self, image_moments=True, alpha=3.5):
        """
        Guess initial parameters for the Moffat fit based on data.

        Parameters
        ----------
        image_moments : bool, optional
            Use the image moments to set the center and amplitude.  If False,
            these are set using the center of the image.
        alpha : float, optional
            Moffat fall-off slope.

        Returns
        -------
        :class:`numpy.ndarray`
            Initial parameter guesses: (x0, y0, amplitude, gamma1, gamma2, phi,
            alpha, background).
        """
        par = np.zeros(6, dtype=float)
        # Clip high values to guess the background
        par[7] = sigma_clipped_stats(self.c, sigma_lower=3., sigma_upper=1.)[1]
        c_bkg = self.c - par[7]
        if image_moments:
            # Use moments to calculate the center, peak value, and width
            par[0] = np.ma.sum(self.x * c_bkg)/np.ma.sum(c_bkg)
            par[1] = np.ma.sum(self.y * c_bkg)/np.ma.sum(c_bkg)
            # Set the amplitude to be the value at the image center
            par[2] = c_bkg[int(np.round(par[0])), int(np.round(par[1]))]
        else:
            # Set the parameters based on the image center
            par[0] = self.x[self.shape[0]//2, self.shape[1]//2]
            par[1] = self.y[self.shape[0]//2, self.shape[1]//2]
            par[2] = self.c[self.shape[0]//2, self.shape[1]//2]
        # Always set alpha explicitly
        par[6] = alpha
        # Set the width based on the first moment of the radius
        r = np.sqrt((self.x-par[0])**2 + (self.y-par[1])**2)
        fwhm = 2 * np.ma.sum(r * c_bkg)/np.ma.sum(c_bkg)
        par[3] = self.to_gamma(fwhm, par[4])
        # Guess it's circular
        par[4] = par[3]
        # And that the rotation angle is 0.
        par[5] = 0.
        return par

    def default_bounds(self):
        """
        Define default bounds for the parameters.

        Returns
        -------
        lb : :class:`numpy.ndarray`
            Lower bounds, ordered as the model parameters.
        ub : :class:`numpy.ndarray`
            Upper bounds, ordered as the model parameters.
        """
        lb = np.zeros(8, dtype=float)
        ub = np.zeros(8, dtype=float)
        lb[0], ub[0] = 0., float(self.shape[0])
        lb[1], ub[1] = 0., float(self.shape[1])
        lb[2], ub[2] = 0., 1.2*np.amax(self.c)
        lb[3], ub[3] = 1., max(self.shape)/2.
        lb[4], ub[4] = 1., max(self.shape)/2.
        lb[5], ub[5] = -1*np.pi/2, np.pi/2
        lb[6], ub[6] = 1., 10.
        lb[7] = np.amin(self.c)-0.1*np.absolute(np.amin(self.c))
        ub[7] = 0.1*np.amax(self.c)
        return lb, ub

    def _set_par(self, par):
        """
        Set the parameter values for the model.

        Parameters
        ----------
        par : :class:`numpy.ndarray`
            Array of parameter values: (x0, y0, amplitude, gamma1, gamma2, phi,
            alpha, background).  Can also be only the free parameters.
        """
        if par.ndim != 1:
            raise ValueError('Parameter array must be a 1D vector.')
        if par.size == self.npar:
            self.par = par.copy()
            return
        if par.size != self.nfree:
            raise ValueError(f'Must provide {self.npar} or {self.nfree} parameters.')
        self.par[self.free] = par.copy()

    @staticmethod
    def _eval_moffat(par, x, y):
        """
        Evaluate the Moffat function given a set of Moffat fit parameters at a
        set of x, y coordinates.

        Parameters
        ----------
        par : :class:`numpy.ndarray`
            Model parameters: (x0, y0, amplitude, gamma1, gamma2, phi, alpha,
            background).
        x : array-like
            x-coordinate values.
        y : array-like
            y-coordinate values.

        Returns
        -------
        :class:`numpy.ndarray`
            Evaluated Moffat function values at the given x, y coordinates.
        """
        cosr = np.cos(par[5])
        sinr = np.sin(par[5])
        A = (cosr/par[3])**2 + (sinr/par[4])**2
        B = (sinr/par[3])**2 + (cosr/par[4])**2
        C = 2 * sinr * cosr * (par[3]**-2 - par[4]**-2)
        dx = (x - par[0])
        dy = (y - par[1])
        D = 1 + A*dx**2 + B*dy**2 + C*dx*dy
        return par[2] / D**par[6]

    @staticmethod
    def _eval_moffat_deriv(par, x, y):
        """
        Calculate the derivatives of the Moffat function with respect to a
        set of given parameters at a set of x, y coordinates.

        Parameters
        ----------
        par : :class:`numpy.ndarray`
            Model parameters: (x0, y0, amplitude, gamma1, gamma2, phi, alpha,
            background).
        x : array-like
            x-coordinate values.
        y : array-like
            y-coordinate values.

        Returns
        -------
        list
            Derivatives of the Moffat function with respect to each parameter.
        """
        cosr = np.cos(par[5])
        sinr = np.sin(par[5])

        a1 = (cosr/par[3])**2
        a2 = (sinr/par[4])**2
        A = a1 + a2
        dA_dg1 = -2*a1/par[3]
        dA_dg2 = -2*a2/par[4]

        b1 = (sinr/par[3])**2
        b2 = (cosr/par[4])**2
        B = b1 + b2
        dB_dg1 = -2*b1/par[3]
        dB_dg2 = -2*b2/par[4]

        C = 2 * sinr * cosr * (par[3]**-2 - par[4]**-2)
        dC_dg1 = -4 * sinr * cosr / par[3]**3
        dC_dg2 = 4 * sinr * cosr / par[4]**3

        dA_dphi = -C
        dB_dphi = C
        dC_dphi = C * (cosr / sinr - sinr / cosr)  # = 2 * C / tan(2 * phi)

        dx = (x.ravel() - par[0])
        dy = (y.ravel() - par[1])

        D = 1 + A*dx**2 + B*dy**2 + C*dx*dy
        dD_dx0 = -2*A*dx - C
        dD_dy0 = -2*B*dy - C
        dD_dg1 = dA_dg1*dx**2 + dB_dg1*dy**2 + dC_dg1*dx*dy
        dD_dg2 = dA_dg2*dx**2 + dB_dg2*dy**2 + dC_dg2*dx*dy
        dD_dphi = dA_dphi*dx**2 + dB_dphi*dy**2 + dC_dphi*dx*dy

        I = par[2] / D**par[6]
        f = -par[6] / D

        return [f * I * dD_dx0,
                f * I * dD_dy0,
                I/par[2],
                f * I * dD_dg1,
                f * I * dD_dg2,
                f * I * dD_dphi,
                -I * np.log(D)]

    def model(self, par=None, x=None, y=None):
        """
        Evaluate the Moffat model with the given or stored parameters.

        Parameters
        ----------
        par : :class:`numpy.ndarray`, optional
            Moffat fit parameters: (x0, y0, amplitude, gamma1, gamma2, phi,
            alpha, background).  If None, use the current parameters.
        x : array-like, optional
            x-coordinate values.  If None, use the data coordinates.
        y : array-like, optional
            y-coordinate values.  If None, use the data coordinates.

        Returns
        -------
        :class:`numpy.ndarray`
            Evaluated Moffat model values at the given or stored x, y
            coordinates.
        """
        if par is not None:
            self._set_par(par)
        return self._eval_moffat(self.par, self.x if x is None else x, self.y if y is None else y) \
                    + self.par[7]

    def resid(self, par):
        """
        Calculate the residuals between the data and the Moffat model.
        """
        return self.c.ravel() - self.model(par).ravel()

    def deriv_resid(self, par):
        """
        Calculate the derivatives of the residuals with respect to the Moffat fit parameters.
        """
        self._set_par(par)
        dmoff = self._eval_moffat_deriv(self.par, self.x, self.y)
        return np.array([-d for d in dmoff]
                            + [np.full(np.prod(self.shape), -1., dtype=float)]).T

    def fit(self, p0=None, fix=None, lb=None, ub=None):
        """
        Fit the Moffat model to the data and modify parameters accordingly.

        Parameters
        ----------
        p0 : :class:`numpy.ndarray`, optional
            Initial parameter guesses: (x0, y0, amplitude, gamma1, gamma2, phi,
            alpha, background).  If None, use :meth:`guess_par`.
        fix : :class:`numpy.ndarray`, optional
            Boolean array indicating which parameters to fix.
        lb : :class:`numpy.ndarray`, optional
            Lower bounds for the parameters.  If None, use
            :meth:`default_bounds`.
        ub : :class:`numpy.ndarray`, optional
            Upper bounds for the parameters.  If None, use
            :meth:`default_bounds`.
        """
        # Guess parameters if no initial guess given
        if p0 is None:
            p0 = self.guess_par()
        _p0 = np.atleast_1d(p0)
        if _p0.size != self.npar:
            raise ValueError('Incorrect number of model parameters.')
        self.par = _p0.copy()
        # Fix designated parameters
        _free = np.ones(self.npar, dtype=bool) if fix is None else np.logical_not(fix)
        if _free.size != self.npar:
            raise ValueError('Incorrect number of model parameter fitting flags.')
        self.free = _free.copy()
        self.nfree = np.sum(self.free)
        
        # Set bounds for fit
        _lb, _ub = self.default_bounds()
        if lb is None:
            lb = _lb
        if ub is None:
            ub = _ub
        if len(lb) != self.npar or len(ub) != self.npar:
            raise ValueError('Length of one or both of the bounds vectors is incorrect.')

        # Perform fit
        p = self.par[self.free]
        result = optimize.least_squares(self.resid, p, method='trf', xtol=1e-12,
                                        bounds=(lb[self.free], ub[self.free]))
        # Sets phi to 0 if fit is too circular for phi to be accurate
        new_par = result.x
        if 0.93 < new_par[3]/new_par[4] < 1.07:
            new_par[5] = 0
        self._set_par(new_par)

    @staticmethod
    def to_gamma(fwhm, alpha):
        """
        Convert full-width half-maximum (FWHM) to gamma.
        """
        return fwhm / 2 / np.sqrt(2**(1/alpha)-1)

    @staticmethod
    def to_fwhm(gamma, alpha):
        """
        Convert gamma to full-width half-maximum (FWHM).
        """
        return 2 * gamma * np.sqrt(2**(1/alpha)-1)

    @staticmethod
    def get_nice_phi(par):
        """
        Convert phi to angle between semi-major axis and +x-axis in degrees.

        Parameters
        ----------
        par : array-like
            Fit parameters: (x0, y0, amplitude, gamma1, gamma2, phi, alpha,
            background).

        Returns
        -------
        float
            Nice phi in degrees.
        """
        nice_phi = np.rad2deg(par[5]) + 0.000000000001
        if par[3] < par[4]:
            if nice_phi > 0.001:
                return nice_phi - 90
            else:
                return nice_phi + 90
        else:
            return nice_phi

    @staticmethod
    def get_orig_phi(gamma1, gamma2, nice_phi):
        """
        Convert nice phi back to original phi (rad).

        Parameters
        ----------
        gamma1 : float
            Gamma value in x-direction.
        gamma2 : float
            Gamma value in y-direction.
        nice_phi : float
            Nice phi in degrees.

        Returns
        -------
        float
            Original phi in radians.
        """
        if gamma1 < gamma2:
            if nice_phi >= -89.999:  # This accounts for the previous condition nice_phi - 90 > 0.001
                original_phi = nice_phi + 90
            else:
                original_phi = nice_phi - 90
        else:
            original_phi = nice_phi
        
        original_rad = np.deg2rad(original_phi - 0.000000000001)
        return original_rad


def make_ellipse(a, b, phi, n=500):
    """
    Generate coordinates for an ellipse.

    Parameters
    ----------
    a : float
        Semi-major axis length.
    b : float
        Semi-minor axis length.
    phi : float
        Rotation angle in radians.
    n : int, optional
        Number of points to generate.

    Returns
    -------
    x : :class:`numpy.ndarray`
        x coordinates of the ellipse.
    y : :class:`numpy.ndarray`
        y coordinates of the ellipse.
    """
    cosr = np.cos(phi)
    sinr = np.sin(phi)
    theta = np.linspace(0., 2*np.pi, num=n)
    x = a * np.cos(theta)
    y = b * np.sin(theta)
    return x*cosr - y*sinr, y*cosr + x*sinr


class EllipticalMoffat2DModel(Fittable2DModel):
    """
    Two dimensional Moffat model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the model.
    x_0 : float
        x position of the maximum of the Moffat model.
    y_0 : float
        y position of the maximum of the Moffat model.
    gamma1 : float
        Core width in 1st arbitrary(?) direction of the Moffat model
    gamma2 : float
        Core width in 2nd arbitrary(?) direction of the Moffat model
    phi : float
        Rotation angle of the Moffat model (rad)
    alpha : float
        Power index of the Moffat model.

    See Also
    --------
    Gaussian2D, Box2D

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2} +
        \\left(y - y_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter(default=1, description="Amplitude (peak value) of the model")
    x_0 = Parameter(
        default=0, description="X position of the maximum of the Moffat model"
    )
    y_0 = Parameter(
        default=0, description="Y position of the maximum of the Moffat model"
    )
    gamma1 = Parameter(default=1, description="Core width in 1st arbitrary(?) direction of the Moffat model")
    gamma2 = Parameter(default=1, description="Core width in 1st arbitrary(?) direction of the Moffat model")
    phi = Parameter(default=0, description="Rotation angle of the Moffat model (rad)")
    alpha = Parameter(default=1, description="Power index of the Moffat model")

    @property
    def fwhm(self):
        """
        Moffat full width at half maximum.
        Derivation of the formula is available in
        `this notebook by Yoonsoo Bach
        <https://nbviewer.jupyter.org/github/ysbach/AO_2017/blob/master/04_Ground_Based_Concept.ipynb#1.2.-Moffat>`_.
        """
        return 2.0 * np.abs(self.gamma) * np.sqrt(2.0 ** (1.0 / self.alpha) - 1.0)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, gamma1, gamma2, phi, alpha):
        """Two dimensional Moffat model function."""
        cosr = np.cos(phi)
        sinr = np.sin(phi)
        A = (cosr/gamma1)**2 + (sinr/gamma2)**2
        B = (sinr/gamma1)**2 + (cosr/gamma2)**2
        C = 2 * sinr * cosr * (gamma1**-2 - gamma2**-2)
        dx = (x - x_0)
        dy = (y - y_0)
        D = 1 + A*dx**2 + B*dy**2 + C*dx*dy
        return amplitude / D**alpha

    @staticmethod
    def fit_deriv(x, y, amplitude, x_0, y_0, gamma1, gamma2, phi, alpha):
        """Two dimensional Moffat model derivative with respect to parameters."""
        cosr = np.cos(phi)
        sinr = np.sin(phi)

        a1 = (cosr/gamma1)**2
        a2 = (sinr/gamma2)**2
        A = a1 + a2
        dA_dg1 = -2*a1/gamma1
        dA_dg2 = -2*a2/gamma2

        b1 = (sinr/gamma1)**2
        b2 = (cosr/gamma2)**2
        B = b1 + b2
        dB_dg1 = -2*b1/gamma1
        dB_dg2 = -2*b2/gamma2

        C = 2 * sinr * cosr * (gamma1**-2 - gamma2**-2)
        dC_dg1 = -4 * sinr * cosr / gamma1**3
        dC_dg2 = 4 * sinr * cosr / gamma2**3

        dA_dphi = -C
        dB_dphi = C
        dC_dphi = C * (cosr / sinr - sinr / cosr)  # = 2 * C / tan(2 * phi)

        dx = (x.ravel() - x_0)
        dy = (y.ravel() - y_0)

        D = 1 + A*dx**2 + B*dy**2 + C*dx*dy
        dD_dx0 = -2*A*dx - C
        dD_dy0 = -2*B*dy - C
        dD_dg1 = dA_dg1*dx**2 + dB_dg1*dy**2 + dC_dg1*dx*dy
        dD_dg2 = dA_dg2*dx**2 + dB_dg2*dy**2 + dC_dg2*dx*dy
        dD_dphi = dA_dphi*dx**2 + dB_dphi*dy**2 + dC_dphi*dx*dy

        I = amplitude / D**alpha
        f = -alpha / D

        return [f * I * dD_dx0,
                f * I * dD_dy0,
                I/amplitude,
                f * I * dD_dg1,
                f * I * dD_dg2,
                f * I * dD_dphi,
                -I * np.log(D)]

    @property
    def input_units(self):
        if self.x_0.input_unit is None:
            return None
        else:
            return {
                self.inputs[0]: self.x_0.input_unit,
                self.inputs[1]: self.y_0.input_unit,
            }

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit[self.inputs[0]] != inputs_unit[self.inputs[1]]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {
            "x_0": inputs_unit[self.inputs[0]],
            "y_0": inputs_unit[self.inputs[0]],
            "gamma1": inputs_unit[self.inputs[0]],
            "gamma2": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }