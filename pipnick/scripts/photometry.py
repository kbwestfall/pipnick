"""
Perform photometry on reduced images
"""

from pipnick.scripts import scriptbase


class Photometry(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Measure the sources in the reduced science images with PSF and aperture '
                        'photometry, writing a source table for each image to the photometry '
                        'directory.  Pixel coordinates are always provided; celestial '
                        'coordinates are provided for images with an astrometric solution (see '
                        'pipnick_astrometry), and masked otherwise.  Rerunning only measures '
                        'images that, or whose settings, have changed, and updates the celestial '
                        'coordinates of images whose astrometric solution has changed.',
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
            help='Only measure these images, given by file name (with or without the extension) '
                 'or combination group (comb_id).'
        )
        parser.add_argument(
            '-o', '--overwrite', default=False, action='store_true',
            help='Measure the images again, even if they are up to date.'
        )
        parser.add_argument(
            '--no_qa', dest='qa', default=True, action='store_false',
            help='Do not write the QA figures.  By default, the figures of the PSF model and '
                 'the measured sources are written to the photometry/qa directory for each '
                 'measured image, and for any image that lacks them.'
        )
        parser.add_argument(
            '--thresh', default=8.0, type=float,
            help='Source detection threshold, as a multiple of the background noise.'
        )
        parser.add_argument(
            '--psf_thresh', default=15.0, type=float,
            help='Detection threshold, as a multiple of the background noise, for the stars used '
                 'to model the PSF.'
        )
        parser.add_argument(
            '--fwhm', default=None, type=float,
            help='Initial guess for the FWHM of the PSF in arcsec.  If not provided, the typical '
                 'value for the camera is used.'
        )
        parser.add_argument(
            '--fittype', default='circ', type=str, choices=['circ', 'ellip'],
            help='The PSF model: a circular or elliptical Moffat function.'
        )
        parser.add_argument(
            '--stack', default='mean', type=str, choices=['mean', 'median'],
            help='How the stars used to model the PSF are stacked.'
        )
        parser.add_argument(
            '--mode', default='all', type=str, choices=['all', 'new'],
            help='Mode of the iterative PSF photometry (see photutils IterativePSFPhotometry).'
        )
        parser.add_argument(
            '--min_sep', default=1.0, type=float,
            help='Minimum separation of detected sources, in units of the FWHM.'
        )
        parser.add_argument(
            '--group_sep', default=2.0, type=float,
            help='Sources closer than this, in units of the FWHM, are fit together.'
        )
        parser.add_argument(
            '--group', default=False, action='store_true',
            help='Consolidate the sources fit together into single sources.'
        )
        parser.add_argument(
            '--aper_radius', default=2.0, type=float,
            help='Radius of the photometric apertures, in units of the FWHM.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.datamodel.dataset import reduction_dir
        from pipnick.pipelines.photometry import photometry_all
        from pipnick.pkg.report import StepReport

        # Initialize logger
        cls.init_log(args, log_dir=cls.log_dir(reduction_dir(args.raw_table, args.rdxdir)))

        report = StepReport('PHOTOMETRY')
        photometry_all(
            args.raw_table, rdxdir=args.rdxdir, images=args.images, overwrite=args.overwrite,
            thresh=args.thresh, fwhm=args.fwhm, fittype=args.fittype, mode=args.mode,
            stack=args.stack, group=args.group, min_sep=args.min_sep, group_sep=args.group_sep,
            aper_radius=args.aper_radius, psf_thresh=args.psf_thresh, qa=args.qa, report=report
        )
        return report.exit_status
