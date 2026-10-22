"""
Perform reduction of raw astronomical data frames (overscan subtraction,
bias subtraction, flat division, cosmic ray masking)
""" 

from pipnick.scripts import scriptbase

class Reduce(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Perform basic image processing: (1) subtract & '
                        'trim overscan, (2) subtract bias, and (3) '
                        'divide flat',
            width=width
        )
        parser.add_argument(
            '-d', '--rawdir', type=str, default=None,
            help='Path to parent directory with the raw data to be reduced.'
        )
        parser.add_argument(
            '-r', '--rdxdir', type=str,
            help='Path for the reduced data.  If not provided, set to the parent '
                 'directory of rawdir.'
        )
        parser.add_argument(
            '-t', '--raw_table', default=None,
            help='File with or for relevant metadata for the raw files to be '
                 'processed by the data reduction script.'
        )
#        parser.add_argument('-d', '--display', action='store_true', help="Display reduced images")
        parser.add_argument(
            '-o', '--overwrite', default=False, action='store_true',
            help='Overwrite any existing files'
        )
        parser.add_argument(
            '-a', '--append', default=False, action='store_true',
            help='Append new files to existing reduction.  This only appends new '
                 'science observations!  New calibration frames are ignored, '
                 'unless overwrite is set, which is identical to re-running the '
                 'entire reduction again from scratch.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.pipelines.reduction import reduce_all
        from pipnick.utils.display_fits import display_many_nickel

        # Initialize logger
        cls.init_log(args)

        # Reduce the images
        metadata = reduce_all(
            rawdir=args.rawdir, table=args.raw_table, rdxdir=args.rdxdir, overwrite=args.overwrite,
            append=args.append
        )

#        # Display reduced images
#        if args.display:
#            display_many_nickel(red_files)
