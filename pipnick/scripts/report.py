"""
Write a static HTML report of the processing of a data set.
"""

from pipnick.scripts import scriptbase


class Report(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Write a static HTML report of the processing of the science images of '
                        'a per-setup metadata file to the report directory of the reduction: '
                        'the status of every step for every image, the provenance and QA '
                        'figures of each image, the summary tables and QA figures of every '
                        'step, and the logs.  Nothing is changed or rebuilt.  The links are '
                        'relative, so the report works from a copy of the whole reduction '
                        'directory.',
            width=width
        )
        parser.add_argument(
            'raw_table', help='Per-setup metadata file.'
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
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.datamodel.dataset import reduction_dir
        from pipnick.pipelines.report import write_report
        from pipnick.pkg.report import StepReport

        # Initialize logger
        cls.init_log(args, log_dir=cls.log_dir(reduction_dir(args.raw_table, args.rdxdir)))

        report = StepReport('REPORT')
        write_report(
            args.raw_table, rdxdir=args.rdxdir, mosaic_dir=args.mosaic_dir, report=report
        )
        return report.exit_status
