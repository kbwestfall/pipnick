import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from loess.loess_2d import loess_2d

from photutils.aperture import CircularAperture
from astropy.visualization import ZScaleInterval

from pipnick import log
from pipnick.utils.nickel_data import ccd_shape
from pipnick.utils.fits_class import Fits_Simple

def smooth_contour(data_x, data_y, data_vals, color_range, backgrd_ax=None, 
                   frac=0.3, title=None, category_str=None):
    """
    Smooths and plots a contour map of the input data using LOESS smoothing.

    Parameters
    ----------
    data_x : :class:`numpy.ndarray`
        Array of x-coordinates of the data points.
    data_y : :class:`numpy.ndarray`
        Array of y-coordinates of the data points.
    data_vals : :class:`numpy.ndarray`
        Array of data values at the corresponding (x, y) points.
    color_range : tuple
        Tuple specifying the range of colors for the contour levels (min, max).
    backgrd_ax : :class:`matplotlib.axes.Axes`, optional
        Background axes to plot on. If None, a new figure and axes are created.
    frac : float, optional
        Fraction of data used to compute locally weighted regression (default is 0.3).
    title : str, optional
        Title of the plot.
    category_str : str, optional
        Category string to include in the plot title.

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
        Axes object with the contour plot.
    cp : :class:`matplotlib.contour.QuadContourSet`, None
        Contour set from the plot. None if the plot could not be generated.
    """
    
    if backgrd_ax is None:
        # Create new figure and axes if not provided
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.set_xlabel('X (pixels)')
        ax.set_ylabel('Y (pixels)')
    else:
        ax = backgrd_ax
    
    # Create a grid on which to sample the smoothed parameters
    subplot_size = 15
    border = int(subplot_size / 2)
    grid_x, grid_y = np.mgrid[border:ccd_shape[0]-border:subplot_size, 
                              border:ccd_shape[1]-border:subplot_size]
    
    try:
        # Perform LOESS smoothing on the data
        param_list, _ = loess_2d(data_x, data_y, data_vals, xnew=grid_x.flatten(),
                                 ynew=grid_y.flatten(), frac=frac)
    except np.linalg.LinAlgError:
        log.warning("LinAlgError: SVD did not converge in Linear Least Squares")
        log.warning("Skipping this contour plot")
        return ax, None
    
    # Reshape the smoothed parameters to match the grid shape
    param_list = param_list.reshape(grid_x.shape)
    
    # Define colors for the contour plot
    colors = ["#cd0000", "#cb4000", "#c97f00", "#c7bc00", "#91c400", "#52c200", 
              "#00bc62", "#00ba9c", "#009cb8", "#0061b6", "#1100b1", "#4800af", "#7e00ad"]
    levels = np.linspace(color_range[0], color_range[1], len(colors))

    # Create the contour plot
    cp = ax.contourf(grid_x, grid_y, param_list, levels=levels, colors=colors)
    ax.set_title(f'{title} Graph - {category_str}')
    
    return ax, cp


def scatter_sources(data_x, data_y, data_vals, color_range, backgrd_ax=None, 
                    title=None, category_str=None):
    """
    Creates a scatter plot of the input data, with points color-coded by value.

    Parameters
    ----------
    data_x : :class:`numpy.ndarray`
        Array of x-coordinates of the data points.
    data_y : :class:`numpy.ndarray`
        Array of y-coordinates of the data points.
    data_vals : :class:`numpy.ndarray`
        Array of data values at the corresponding (x, y) points.
    color_range : tuple
        Tuple specifying the range of colors for the scatter points (min, max).
    backgrd_ax : :class:`matplotlib.axes.Axes`, optional
        Background axes to plot on. If None, a new figure and axes are created.
    title : str, optional
        Title of the plot.
    category_str : str, optional
        Category string to include in the plot title.

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
        Axes object with the scatter plot.
    cmap_custom : :class:`matplotlib.colors.ListedColormap`
        Colormap used for the scatter plot.
    """
    
    if backgrd_ax is None:
        # Create new figure and axes if not provided
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.set_xlabel('X (pixels)')
        ax.set_ylabel('Y (pixels)')
    else:
        ax = backgrd_ax
    
    # Define the colors and range of contour levels
    colors = ["#cd0000", "#cb4000", "#c97f00", "#c7bc00", "#91c400", "#52c200", 
              "#00bc62", "#00ba9c", "#009cb8", "#0061b6", "#1100b1", "#4800af"]
    levels = np.linspace(color_range[0], color_range[1], len(colors))
    cmap_custom = ListedColormap(colors)
    
    # Add a random coordinate offset to show overlapping points
    jitter_x = np.random.normal(scale=7, size=len(data_x))
    jitter_y = np.random.normal(scale=7, size=len(data_y))
    
    # Create a scatter plot, color-coded by data values
    ax.set_title(f'{title} Graph - {category_str}')
    ax.scatter(data_x + jitter_x, data_y + jitter_y, s=15, c=data_vals, cmap=cmap_custom, 
               vmin=levels[0], vmax=levels[-1], alpha=1.0,
               linewidths=0.7, edgecolors='k')
    
    return ax, cmap_custom


def plot_sources(phot_data, given_fwhm, image=None, flux_name='flux_fit',
                 x_name='x_fit', y_name='y_fit', label_name='group_id',
                 scale=1):
    """
    Plots sources from a photometric data table on a corresponding image, highlighting 
    the grouped/ungrouped sources.

    Parameters
    ----------
    phot_data : :class:`astropy.table.Table`
        Photometric data containing positions and fluxes of sources.
    given_fwhm : float
        Full-width half-maximum (FWHM) of the sources to set aperture sizes.
    image : :class:`~pipnick.utils.fits_class.Fits_Simple`, optional
        Image to plot the sources on. If None, it is loaded from the metadata.
    flux_name : str, optional
        Name of the flux column in phot_data (default is 'flux_fit').
    x_name : str, optional
        Name of the x-coordinate column in phot_data (default is 'x_fit').
    y_name : str, optional
        Name of the y-coordinate column in phot_data (default is 'y_fit').
    label_name : str, optional
        Name of the label column in phot_data to use for annotating sources (default is 'group_id').
    scale : float, optional
        Scaling factor for the aperture sizes and annotation text (default is 1).

    Returns
    -------
    None
    """
    
    if image is None:
        # Load the image from the metadata if not provided
        image = Fits_Simple(phot_data.meta['image_path'])
    log.info(f'Plotting image {image}')
    
    if flux_name == 'flux_fit' and 'flux_fit' not in phot_data.colnames:
        flux_name = 'flux_psf'
    
    def get_apertures(phot_data):
        """Create circular apertures for the sources based on their positions."""
        x = phot_data[x_name]
        y = phot_data[y_name]
        positions = np.transpose((x, y))
        return CircularAperture(positions, r=2 * given_fwhm * scale)
    
    # Separate good and bad photometric data based on group size
    good_phot_data = phot_data[phot_data['group_size'] <= 1]
    bad_phot_data = phot_data[phot_data['group_size'] > 1]
    bad_apertures = get_apertures(bad_phot_data)
    good_apertures = get_apertures(good_phot_data)
    
    # Determine image display limits using ZScale
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(image.masked_array)
    
    # Set colormap and mask bad pixels with red
    cmap = plt.get_cmap()
    cmap.set_bad('r', alpha=0.5)
    
    # Plot the image and the good/bad sources
    plt.figure(figsize=(12,10))
    plt.title(image)
    plt.imshow(image.masked_array, origin='lower', vmin=vmin, vmax=vmax,
               cmap=cmap, interpolation='nearest')
    plt.colorbar()
    good_apertures.plot(color='purple', lw=1.5*scale, alpha=1)
    bad_apertures.plot(color='r', lw=1.5*scale, alpha=1)
    
    # Annotate singular sources with label_name and flux_name values
    y_offset = 3.5*given_fwhm*scale
    for i in range(len(good_phot_data)):
        plt.text(good_phot_data[x_name][i], good_phot_data[y_name][i]+y_offset, 
                 f'{good_phot_data[label_name][i]}: {good_phot_data[flux_name][i]:.0f}',
                 color='white', fontsize=8*scale, ha='center', va='center')
    
    # Annotate grouped sources with label_name and flux_name values in one large stack
    group_ids = set(bad_phot_data[label_name])
    for id in group_ids:
        group = bad_phot_data[bad_phot_data[label_name] == id]
        group_x = np.median(group[x_name])
        group_y = np.median(group[y_name]) + y_offset
        for i in range(len(group)):
            plt.text(group_x, group_y+i*20*scale, 
                     f'{id}: {group[flux_name][i]:.0f}',
                     color='red', fontsize=8*scale, ha='center', va='center')
    
    # Show plot
    plt.gcf().set_dpi(300)
    plt.show()
