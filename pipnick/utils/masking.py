from pathlib import Path

import numpy as np

def polygon_winding_number(polygon, point):
    """
    Determine the winding number of a 2D polygon about a point.
    
    The code does **not** check if the polygon is simple (no interesecting line
    segments).  Algorithm taken from Numerical Recipes Section 21.4.

    Parameters
    ----------
    polygon : np.ndarray
        An Nx2 array containing the x,y coordinates of a polygon.  The points
        should be ordered either counter-clockwise or clockwise.
    point : np.ndarray
        One or more points for the winding number calculation.  Must be either a
        2-element array for a single (x,y) pair, or an Nx2 array with N (x,y)
        points.

    Returns
    -------

    :obj:`int`, np.ndarray
        The winding number of each point with respect to the provided polygon.
        Points inside the polygon have winding numbers of 1 or -1; see
        :func:`point_inside_polygon`.

    Raises
    ------
    ValueError
        Raised if ``polygon`` is not 2D, if ``polygon`` does not have two
        columns, or if the last axis of ``point`` does not have 2 and only 2
        elements.
    """
    # Check input shape is for 2D only
    if len(polygon.shape) != 2:
        raise ValueError('Polygon must be an Nx2 array.')
    if polygon.shape[1] != 2:
        raise ValueError('Polygon must be in two dimensions.')
    _point = np.atleast_2d(point)
    if _point.shape[1] != 2:
        raise ValueError('Point must contain two elements.')

    # Get the winding number
    nvert = polygon.shape[0]
    npnt = _point.shape[0]

    dl = np.roll(polygon, 1, axis=0)[None,:,:] - _point[:,None,:]
    dr = polygon[None,:,:] - point[:,None,:]
    dx = dl[...,0]*dr[...,1] - dl[...,1]*dr[...,0]

    indx_l = dl[...,1] > 0
    indx_r = dr[...,1] > 0

    wind = np.zeros((npnt, nvert), dtype=int)
    wind[indx_l & np.logical_not(indx_r) & (dx < 0)] = -1
    wind[np.logical_not(indx_l) & indx_r & (dx > 0)] = 1

    return np.sum(wind, axis=1)[0] if point.ndim == 1 else np.sum(wind, axis=1)


def point_inside_polygon(polygon, point):
    """
    Determine if one or more points is inside the provided polygon.

    Primarily a wrapper for :func:`polygon_winding_number`, that
    returns True for each point that is inside the polygon.

    Parameters
    ----------
    polygon : np.ndarray
        An Nx2 array containing the x,y coordinates of a polygon.  The points
        should be ordered either counter-clockwise or clockwise.
    point : np.ndarray
        One or more points for the winding number calculation.  Must be either a
        2-element array for a single (x,y) pair, or an Nx2 array with N (x,y)
        points.

    Returns
    -------
    :obj:`bool`, `np.ndarray`
        Boolean indicating whether or not each point is within the polygon.
    """
    return np.absolute(polygon_winding_number(polygon, point)) == 1


def mask_polygons(shape, polygons):
    """
    Mask a set of polygons in an 2D image.

    Coordinates of the polygon vertices must be in x,y coordinates (in pixel
    indices) following the numpy.  I.e., the vertex coordinates should follow
    the x,y convention when displaying the array using matplotlib (x is the
    second axis, y is the first axis).  Note the pixel coordinates are at the
    center of the pixel, and the vertices of the first pixel are at :math:`x,y =
    (-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)`.

    Parameters
    ----------
    shape : tuple
        Two-tuple with the image shape.
    polygons : list
        A list of Nx2 arrays providing the vertices of a set of polygons.  The
        number of vertices in each polygon can be different.

    Returns
    -------
    np.ndarray
        A boolean numpy array selecting pixels that should be masked
        (masked=True; unmasked=False).
    """
    # Check input
    if len(shape) != 2:
        raise ValueError('Shape must be 2 integers.')
    if not isinstance(polygons, list):
        raise TypeError('Entered polygons object must be a list.  It can be a list of one array.')
    # Build the array coordinates
    m, n = shape
    coo = np.array(list(map(lambda x : x.ravel(), np.meshgrid(np.arange(n), np.arange(m))))).T
    # Iteratively mask each polygon
    mask = np.zeros(shape, dtype=bool)
    for p in polygons:
        mask |= point_inside_polygon(p, coo).reshape(m,n)
    return mask


def create_mask(shape, bad_cols=None, bad_slices=None, bad_regions=None):
    """
    Apply masks to an image by masking specified columns, slices, or regions.

    Parameters
    ----------
    shape : tuple
        Two-tuple with the shape of the 2D data array.
    bad_cols : list, optional
        A list of column indices to mask.
    bad_slices : list, optional
        A list of array slices to mask.
    bad_regions : list, optional
        A list of Nx2 np.ndarray objects defining polygons that identify regions
        to mask.

    Returns
    -------
    mask : np.ndarray
        Boolean mask (masked=True, unmasked=False)
    """
    mask = np.zeros(shape, dtype=bool) if bad_regions is None \
            else mask_polygons(shape, bad_regions)
    if bad_slices is not None:
        for s in bad_slices:
            mask[s] = True
    if bad_cols is not None:
        mask[:,bad_cols] = True
    return mask


#def get_masks_from_file(mode: str) -> np.ndarray:
#    """
#    Load and retrieve specific masks from a .npz file based on the given mode.
#
#    Parameters
#    ----------
#    mode : str
#        Specifies which mask to retrieve. Possible values are:
#        - 'mask': Returns the standard nickel mask.
#        - 'fov_mask': Returns the standard nickel mask with overscan columns.
#        - 'mask_cols_only': Returns the nickel mask with only bad columns masked.
#        - 'fov_mask_cols_only': Returns the nickel mask with overscan columns and only bad columns masked.
#
#    Returns
#    -------
#    np.ndarray
#        The requested mask array corresponding to the specified mode.
#
#    Raises
#    ------
#    ValueError
#        If the mode is not one of the specified options.
#    """
#    loaded_masks = np.load(mask_file)
#    
#    if mode == 'mask':
#        return loaded_masks['nickel_mask']
#    elif mode == 'fov_mask':
#        return loaded_masks['nickel_fov_mask']
#    elif mode == 'mask_cols_only':
#        return loaded_masks['nickel_mask_cols_only']
#    elif mode == 'fov_mask_cols_only':
#        return loaded_masks['nickel_fov_mask_cols_only']
#    else:
#        raise ValueError("Invalid mode. Expected one of: 'mask', 'fov_mask', 'mask_cols_only', 'fov_mask_cols_only'.")
 
 
#def generate_masks() -> tuple:
#    """
#    Generate and save the masks for Nickel images, including masks for bad columns, 
#    blind corners, and field of view (FOV) padding.
#
#    Returns
#    -------
#    tuple
#        A tuple containing four masks:
#        - nickel_mask_cols_only: Mask for Nickel images with only bad columns masked.
#        - nickel_mask: Mask for Nickel images with bad columns and blind corners masked.
#        - nickel_fov_mask_cols_only: Mask for Nickel images with FOV padding and only bad columns masked.
#        - nickel_fov_mask: Mask for Nickel images with FOV padding, bad columns, and blind corners masked.
#    """
#    # Mask for Nickel images (masking bad columns and blind corners)
#    nickel_mask_cols_only = add_mask(np.zeros(ccd_shape), bad_columns, [], []).mask
#    nickel_mask = add_mask(np.zeros(ccd_shape), bad_columns, bad_triangles, bad_rectangles).mask
#
#    # Calculate the padding needed
#    pad_height = fov_shape[0] - ccd_shape[0]
#    pad_width = fov_shape[1] - ccd_shape[1]
#    # Apply padding
#    nickel_fov_mask_cols_only = np.pad(nickel_mask_cols_only, ((0, pad_height), (0, pad_width)), mode='constant', constant_values=0)
#    nickel_fov_mask = np.pad(nickel_mask, ((0, pad_height), (0, pad_width)), mode='constant', constant_values=0)
#    
#    np.savez(mask_file, nickel_mask_cols_only=nickel_mask_cols_only,
#                        nickel_mask=nickel_mask,
#                        nickel_fov_mask_cols_only=nickel_fov_mask_cols_only,
#                        nickel_fov_mask=nickel_fov_mask)
#    
#    return nickel_mask_cols_only, nickel_mask, nickel_fov_mask_cols_only, nickel_fov_mask


#def add_mask(data: np.ndarray, cols_to_mask: list, tris_to_mask: list, rects_to_mask: list) -> np.ma.MaskedArray:
#    """
#    Apply masks to an image by masking specified columns, triangles, and rectangles.
#
#    Parameters
#    ----------
#    data : np.ndarray
#        A 2D numpy array representing the image to mask.
#    cols_to_mask : list
#        A list of column indices to mask.
#    tris_to_mask : list
#        A list of tuples, where each tuple contains 3 coordinates representing a triangle to mask.
#    rects_to_mask : list
#        A list of tuples, where each tuple contains 4 coordinates representing a rectangle to mask.
#
#    Returns
#    -------
#    np.ma.MaskedArray
#        A masked array with the specified regions masked.
#    """
#    rows, cols = data.shape
#    mask = np.zeros((rows, cols), dtype=bool)
#    
#    # Mask the rectangles
#    for rectangle in rects_to_mask:
#        mask[rectangle[0][0]:rectangle[1][0], rectangle[0][1]:rectangle[1][1]] = True
#    # Transpose mask so that correct areas are masked (FITS indexing is odd)
#    mask = mask.T
#
#    for triangle in tris_to_mask:
#        # Create a path object for the triangle
#        path = matPath(triangle)
#        
#        # Determine which points are inside the triangle
#        y, x = np.mgrid[:rows, :cols]
#        points = np.vstack((x.flatten(), y.flatten())).T
#        mask = np.logical_or(mask, path.contains_points(points).reshape(rows, cols))
#
#    # Mask the specified columns
#    for col in cols_to_mask:
#        mask[:, col] = True
#    
#    # Create the masked array
#    masked_data = np.ma.masked_array(data, mask)
#    return masked_data

