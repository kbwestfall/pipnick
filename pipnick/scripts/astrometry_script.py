"""
Perform astrometric calibration on reduced images
""" 

from pipnick.scripts import scriptbase

class AstrometryPipeline(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Performs astrometric calibration on reduced images', width=width
        )
        parser.add_argument('maindir', type=str,
                            help='Path to the main directory containing the reduced directory '
                                 'with the FITS files to be astrometrically calibrated.')
        parser.add_argument('apikey', type=str,
                            help='API key from https://nova.astrometry.net account')
        parser.add_argument('-t', '--use_table', action='store_true',
                            help='Whether to use the table file to automatically exclude '
                                 'files that have been commented-out')
        parser.add_argument('-r', '--resolve', action='store_true', 
                            help="re-solves images with previously generated local solves")
        parser.add_argument('--excl_files', default=[], type=list,
                            help='List of file stems substrings to exclude (exact match not '
                                 'necessary).')
        parser.add_argument('--excl_objs', default=[], type=list,
                            help='List of object substrings to exclude (exact match not '
                                 'necessary).')
        parser.add_argument('--excl_filts', default=[], type=list,
                            help='List of filter substrings to exclude (exact match not '
                                 'necessary).')
        return parser

    @classmethod
    def main(cls, args):

        from pipnick.pipelines.astrometry import astrometry_all

        # Initialize logger
        cls.init_log(args)

        calib_files = astrometry_all(args.maindir, args.apikey, args.use_table, 
                                     args.resolve, args.excl_files, args.excl_objs, 
                                     args.excl_filts)
        
        return calib_files
        
        