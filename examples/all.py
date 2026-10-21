from pipnick.pipelines.reduction import reduce_all
from pipnick.pipelines.astrometry import astrometry_all
from pipnick.pipelines.photometry import photometry_all
from pipnick.pipelines.final_calib import final_calib_all
from pipnick.photometry.fit import plot_sources

import logging
from pipnick.utils.log import adjust_global_logger
adjust_global_logger('INFO', __name__)
logger = logging.getLogger(__name__)

# Paste your own API Key
api_key = 'exampleapikey'

# Define directory containing raw images
maindir = 'data_example'

# Basic reduction
red_files = reduce_all(maindir)

# Astrometric calibration
astro_calib_files = astrometry_all(maindir, api_key=api_key)

# Photometric calibration
src_catalog_paths = photometry_all(maindir, group=True, plot_final=True,
                                   plot_inters=False)

# Final calibration (convert pixel coordinates -> RA/Dec)
astrophot_data_tables = final_calib_all(maindir)

# Display images & annotate sources
for object, src_table_dict in astrophot_data_tables.items():
    for file_key, src_table in src_table_dict.items():
        plot_sources(src_table, given_fwhm=8.0, flux_name='flux_psf', scale=1.5)
