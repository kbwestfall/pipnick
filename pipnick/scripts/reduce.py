"""
Perform reduction of raw astronomical data frames (overscan subtraction,
bias subtraction, flat division, cosmic ray masking)
"""

from pipnick.scripts import scriptbase

class Reduce(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Perform basic image processing: (1) subtract & trim overscan, (2) '
                        'subtract bias, and (3) divide flat.  Frames are combined in the groups '
                        'set by the comb_id column of a per-setup metadata file.  Rerunning the '
                        'reduction only reprocesses the outputs whose raw frames or '
                        'calibrations have changed.',
            width=width
        )
        parser.add_argument(
            '-t', '--raw_table', required=True,
            help='Per-setup metadata file for the raw frames to process, written by '
                 'pipnick_metadata --split.  To add new frames, use pipnick_metadata --split '
                 '--append and rerun the reduction.'
        )
        parser.add_argument(
            '-r', '--rdxdir', type=str,
            help='Path for the reduced data, which can only hold the products of one setup.  '
                 'If not provided, the path saved in the metadata file is used or, if there is '
                 'none, the parent directory of the raw data directory.'
        )
        parser.add_argument(
            '-o', '--overwrite', default=False, action='store_true',
            help='Reprocess all outputs, even if they are up to date.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.datamodel.dataset import reduction_dir
        from pipnick.pipelines.reduction import reduce_all
        from pipnick.pkg.report import StepReport

        # Initialize logger
        cls.init_log(args, log_dir=cls.log_dir(reduction_dir(args.raw_table, args.rdxdir)))

        # Reduce the images
        report = StepReport('REDUCTION')
        reduce_all(args.raw_table, rdxdir=args.rdxdir, overwrite=args.overwrite, report=report)
        return report.exit_status
