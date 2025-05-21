"""
Perform photometric calibration on reduced images
"""

from pipnick.scripts import scriptbase
from pipnick.utils import log

class PhotometryPipeline(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Performs photometry on reduced images',
                                    width=width)

        # Indicate where to find the reduced files to be processed.  Must either
        # provide the directory with the reduced images or the table used,
        # and/or created, by the data-reduction script (or an edited version of
        # it).  Cannot provide both.
        inp_group = parser.add_mutually_exclusive_group(required=True)
        inp_group.add_argument('-r', '--rdxdir', type=str, default=None,
                              help='Path to the directory with the reduced image data.')
        inp_group.add_argument('-t', '--raw_table', default=None,
                               help='File with metadata for the raw files processed by the data '
                                    'reduction script.')

#        parser.add_argument('-o', '--photdir', type=str,
#                            help='Path for the photometry data products.  Default is the '
#                                 'directory with the reductions.')
        parser.add_argument('-p', '--phot_table', default=None,
                            help='Summary file with the metadata for the photometry results.')

#        parser.add_argument('maindir', type=str,
#                            help='Path to main directory containing reduced directory with the files to be photometrically analyzed.')
#        parser.add_argument('-t', '--use_table', action='store_true',
#                            help='Whether to use the table file to automatically exclude files that have been commented-out')
#        parser.add_argument('--excl_files', default=[], type=list,
#                            help='List of file stems substrings to exclude (exact match not necessary).')
#        parser.add_argument('--excl_objs', default=[], type=list,
#                            help='List of object substrings to exclude (exact match not necessary).')
#        parser.add_argument('--excl_filts', default=[], type=list,
#                            help='List of filter substrings to exclude (exact match not necessary).')

        parser.add_argument('-s', '--sigma_thresh', default=8.0, type=float,
                            help='Threshold for source detection = background std * thresh.')
        parser.add_argument('-g', '--group', action='store_true', 
                            help='Consolidates groups of sources detected together into one source')
        parser.add_argument('-m', '--mode', default='all', type=str,
                            help='Mode to run photutils PSFPhotometry. `all` recommended.',
                            choices=['all', 'new'])
        parser.add_argument('-f', '--fittype', default='circ', type=str,
                            help='Which type of Moffat fit to use for PSF photometry',
                            choices=['circ', 'ellip'])

#        parser.add_argument('-pf', '--plot_final', action='store_true', 
#                            help="Displays images with sources & flux labelled")
#        parser.add_argument('-pi', '--plot_inters', action='store_true', 
#                            help="Displays images with initial sources & source groups for inspection")
        parser.add_argument('-l', '--log', type=str, default=log.default_log_name(__name__),
                            help='Log file')
        parser.add_argument('--verbose', '-v', action='count', default=0,
                            help='Level of verbosity to display (5=highest; "-vv" sets level=2)')
        parser.add_argument('--overwrite', '-o', default=False, action='store_true',
                            help='Overwrite any existing files')
        return parser

    @staticmethod
    def main(args):
        
        from pipnick.pipelines.photometry import photometry_all

        # Initialize logger
        log.adjust_global_logger(log_level='INFO' if args.verbose < 1 else args.verbose,
                                 log_file=args.log)
        
        src_catalogs = photometry_all(rdxdir=args.rdxdir, raw_table=args.raw_table,
                                      phot_table=args.phot_table,
                                      thresh=args.sigma_thresh, group=args.group, 
                                      mode=args.mode, fittype=args.fittype,
                                      overwrite=args.overwrite)
        
#        return src_catalogs
 