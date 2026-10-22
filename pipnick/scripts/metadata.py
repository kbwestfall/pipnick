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
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files')
        parser.add_argument('-a', '--append', default=False, action='store_true',
                            help='Append new files to existing reduction.  This only appends new '
                                 'science observations!  New calibration frames are ignored, '
                                 'unless overwrite is set, which is identical to re-running the '
                                 'entire script again from scratch.')
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.metadata import build_metadata

        # Initialize logger
        cls.init_log(args)

        # Organize raw files based on input directory or table
        build_metadata(
            rawdir=args.rawdir,
            filename=args.raw_table,
            overwrite=args.overwrite or args.append
        )

