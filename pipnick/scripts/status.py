"""
Report the status of the processing of the science images of a data set.
"""

from pipnick.scripts import scriptbase


class Status(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Report the status of each processing step for each science image of '
                        'a per-setup metadata file, and the next command to run.  Nothing is '
                        'changed.  Each step is ok, stale (it must be run again), missing, '
                        'failed, or n/a.',
            width=width
        )
        parser.add_argument(
            '-t', '--raw_table', required=True, help='Per-setup metadata file.'
        )
        parser.add_argument(
            '-r', '--rdxdir', type=str,
            help='Path for the reduced data.  If not provided, the path saved in the metadata '
                 'file is used.'
        )
        parser.add_argument(
            '--mosaic_dir', type=str, default=None,
            help='Directory with the mosaics.  If not provided, the mosaic subdirectory of the '
                 'reduction directory is used.'
        )
        parser.add_argument(
            '--ofile', type=str, default=None,
            help='Also write the status table to this ECSV file.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.datamodel.dataset import reduction_dir
        from pipnick.pipelines.status import status_all

        # Initialize logger
        cls.init_log(args, log_dir=cls.log_dir(reduction_dir(args.raw_table, args.rdxdir)))

        status_all(
            args.raw_table, rdxdir=args.rdxdir, mosaic_dir=args.mosaic_dir, ofile=args.ofile
        )
