"""
Perform reduction of raw astronomical data frames (overscan subtraction,
bias subtraction, flat division, cosmic ray masking)
""" 

from pipnick.scripts import scriptbase

class MetaData(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Construct a table with metadata pulled from data file headers.',
            width=width
        )
        parser.add_argument('-d', '--rawdir', type=str, default=None,
                            help='Path to parent directory with the raw data to be reduced.')
        parser.add_argument('-t', '--raw_table', default=None,
                            help='File for the metadata')
        parser.add_argument('--verbose', '-v', action='count', default=0,
                            help='Level of verbosity to display (5=highest; "-vv" sets level=2)')
        parser.add_argument('--overwrite', '-o', default=False, action='store_true',
                            help='Overwrite any existing files')
        parser.add_argument('--append', '-a', default=False, action='store_true',
                            help='Append new files to existing reduction.  This only appends new '
                                 'science observations!  New calibration frames are ignored, '
                                 'unless overwrite is set, which is identical to re-running the '
                                 'entire script again from scratch.')
        return parser

    @staticmethod
    def main(args):
        from pipnick.utils.log import adjust_global_logger
        from pipnick.utils.dir_nav import build_metadata

        # Initialize logger
        adjust_global_logger(log_level='INFO' if args.verbose < 1 else args.verbose, name=__name__)

        # Organize raw files based on input directory or table
        build_metadata(
            rawdir=args.rawdir,
            filename=args.raw_table,
            overwrite=args.overwrite or args.append
        )

