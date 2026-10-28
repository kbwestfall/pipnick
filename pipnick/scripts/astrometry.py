"""
Perform astrometric calibration on reduced images
"""

from pipnick.scripts import scriptbase


class Astrometry(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Find the astrometric solutions of the reduced science images, using as '
                        'few astrometry.net submissions as possible: images are grouped by '
                        'pointing, one image per group is solved, and the others are solved by '
                        'matching their sources to it or to the Gaia catalog.  The solutions are '
                        'written to the image headers.  The astrometry.net API key is read from '
                        'the ASTROMETRY_NET_API_KEY environment variable or the astroquery '
                        'configuration file.',
            width=width
        )
        parser.add_argument(
            '-t', '--raw_table', required=True,
            help='Per-setup metadata file used to reduce the images (see pipnick_reduce).'
        )
        parser.add_argument(
            '-r', '--rdxdir', type=str,
            help='Path to the reduced data.  If not provided, the path saved in the metadata '
                 'file is used.'
        )
        parser.add_argument(
            '--images', nargs='+', default=None,
            help='Only solve these images, given by file name (with or without the extension) '
                 'or combination group (comb_id).  The other images are not changed, but their '
                 'solutions are used as references.'
        )
        parser.add_argument(
            '--resolve', default=False, action='store_true',
            help='Ignore existing and cached solutions and solve the images (only those given '
                 'by --images, if provided) again.  This can submit jobs to astrometry.net for '
                 'images that were solved before.'
        )
        parser.add_argument(
            '--upload_image', default=False, action='store_true',
            help='Upload the images to astrometry.net (with masked pixels set to 0), instead '
                 'of the lists of their detected sources.'
        )
        parser.add_argument(
            '--dry_run', default=False, action='store_true',
            help='Detect the sources, group the images, and report what would be submitted to '
                 'astrometry.net, without solving or submitting anything.'
        )
        parser.add_argument(
            '--qa', default=False, action='store_true',
            help='Write a QA figure for each image solved in this run (for every image, with '
                 '--dry_run) to the astrometry/qa directory, showing the detected sources and '
                 'the matched sources with their residuals.'
        )
        parser.add_argument(
            '--max_sources', type=int, default=100,
            help='Maximum number of (brightest) sources detected in each image, and submitted '
                 'to astrometry.net.'
        )
        parser.add_argument(
            '--detect_thresh', type=float, default=5.,
            help='Source detection threshold, as a multiple of the standard deviation of the '
                 'background.'
        )
        parser.add_argument(
            '--fwhm', type=float, default=None,
            help='FWHM of the sources in arcsec, for their detection.  If not provided, the '
                 'typical value for the camera is used.'
        )
        parser.add_argument(
            '--max_submissions', type=int, default=None,
            help='Maximum number of submissions to astrometry.net.  If not provided, there is '
                 'no limit.'
        )
        parser.add_argument(
            '--group_radius', type=float, default=1.,
            help='Maximum distance in arcmin between the pointings of images in a group.'
        )
        parser.add_argument(
            '--min_matches', type=int, default=4,
            help='Minimum number of matched sources for a solution found from another image or '
                 'the Gaia catalog.'
        )
        parser.add_argument(
            '--max_rms', type=float, default=0.3,
            help='Maximum RMS residual in arcsec for a solution found from another image or the '
                 'Gaia catalog.'
        )
        parser.add_argument(
            '--search', type=float, default=60.,
            help='Maximum error in arcsec in the expected offset between matched images, or '
                 'between an image and the Gaia catalog, e.g., from the error in the telescope '
                 'pointing.'
        )
        parser.add_argument(
            '--wide_search', type=float, default=10.,
            help='For images that cannot be solved otherwise, the maximum error in arcmin in '
                 'the pointing in the metadata table allowed when matching to the Gaia '
                 'catalog.  Set to 0 to skip this search.'
        )
        parser.add_argument(
            '--wide_min_matches', type=int, default=10,
            help='Minimum number of matched sources for a solution found by the wide search.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.datamodel.dataset import reduction_dir
        from pipnick.pipelines.astrometry import astrometry_all
        from pipnick.pkg.report import StepReport

        # Initialize logger
        cls.init_log(args, log_dir=cls.log_dir(reduction_dir(args.raw_table, args.rdxdir)))

        report = StepReport('ASTROMETRY')
        astrometry_all(
            args.raw_table, rdxdir=args.rdxdir, resolve=args.resolve,
            upload_image=args.upload_image, dry_run=args.dry_run,
            max_submissions=args.max_submissions, group_radius=args.group_radius,
            min_matches=args.min_matches, max_rms=args.max_rms, search=args.search,
            wide_search=args.wide_search,
            wide_min_matches=args.wide_min_matches, images=args.images, qa=args.qa,
            max_sources=args.max_sources, detect_thresh=args.detect_thresh, fwhm=args.fwhm,
            report=report
        )
        return report.exit_status
