import numpy as np

# For the Nickel Telescope original science camera
plate_scale_approx = 0.37
ccd_shape = np.array([1024, 1024])
fov_shape = np.array([1024, 1056])
fov_size_approx = 6.3   # in arcmin

sat_columns = [255]

bad_columns = [255, 256, 783, 784, 1002]
bad_photometry_columns = [252, 253, 254, 255, 256, 257, 258, 259, 260, 783, 784, 1002]
bad_triangles = [((0, 960), (64, 1024), (0, 1024)), ((0, 33), (34, 0), (0, 0))]
bad_rectangles = []

# Header keys
ra_key = 'RA'
dec_key = 'DEC'

# Format of RA / Dec values in header, for astroquery
ra_dec_units = ('hour', 'degree')

# Expected values for the Nickel Camera
fwhm_approx = 5.0   # pixels

# Labels for different types of astronomical frames
bias_label = 'Bias'
dome_flat_label = 'Dome flat'
sky_flat_label = 'Flat'
sky_flat_label_alt = 'Sky flat'
dark_label = 'dark'
focus_label = 'focus'

# For slow readout speed and 2x2 binning
gain = 1.8
read_noise = 10.7
