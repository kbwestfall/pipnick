"""
Construct mosaics of the reduced science images
"""

from pipnick.scripts import scriptbase


class Mosaic(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Construct the mosaics defined by the mosaic and wcs_group columns of '
                        'the per-setup metadata tables (see pipnick_metadata --assign_mosaics).  '
                        'Images with the same mosaic label are combined, and the mosaics with '
                        'the same wcs_group label share an output WCS.  The images are matched '
                        'to a common PSF and scaled to a common flux scale before they are '
                        'combined.  Run pipnick_astrometry and pipnick_photometry first.',
            width=width
        )
        parser.add_argument(
            'raw_table', nargs='+',
            help='Per-setup metadata files used to reduce the images (see pipnick_reduce).  '
                 'The images in all of the tables with the same labels are combined.'
        )
        parser.add_argument(
            '--outdir', type=str, default=None,
            help='Output directory.  If not provided, the mosaic subdirectory of the reduction '
                 'directory is used, which requires a single table.'
        )
        parser.add_argument(
            '--mosaics', nargs='+', default=None, help='Only construct these mosaics (labels).'
        )
        parser.add_argument(
            '--psf_match', default='mosaic', choices=['mosaic', 'group', 'none'],
            help='Convolve the images to the worst seeing of their mosaic, of their WCS group '
                 '(e.g., for color images), or not at all.'
        )
        parser.add_argument(
            '--kernel', default='moffat', choices=['moffat', 'gaussian'],
            help='The PSF matching kernel: the kernel matching the Moffat PSF models (with a '
                 'Gaussian kernel as the fallback), or a Gaussian kernel.'
        )
        parser.add_argument(
            '--psf_margin', default=1.05, type=float,
            help='The target FWHM of the PSF matching, in units of the worst FWHM.'
        )
        parser.add_argument(
            '--scale', default='sources', choices=['sources', 'calib', 'none'],
            help='Scale the images to the flux scale of the reference image using the sources '
                 'common to overlapping images, using their photometric calibrations '
                 '(pipnick_calibrate --apply image), or not at all.'
        )
        parser.add_argument(
            '--reference', default=None,
            help='The reference image of the flux scale (file name) of the mosaic that '
                 'includes it.  The other mosaics use the image with the lowest airmass.'
        )
        parser.add_argument(
            '--match_radius', default=1., type=float,
            help='Radius (arcsec) for matching sources between images.'
        )
        parser.add_argument(
            '--min_common', default=3, type=int,
            help='Minimum number of sources common to a pair of images used for the scaling.'
        )
        parser.add_argument(
            '--min_snr', default=20., type=float,
            help='Minimum signal-to-noise ratio of the sources used for the scaling.'
        )
        parser.add_argument(
            '--background', default='constant', choices=['constant', '2d', 'none'],
            help='Background subtraction: a constant (the median of the pixels without '
                 'sources), a smooth 2D background (256-pixel boxes), or none.  Use constant '
                 'for objects that fill much of the field.'
        )
        parser.add_argument(
            '--match_background', default=False, action='store_true',
            help='Match the background levels of the images, using the median differences in '
                 'their overlaps.'
        )
        parser.add_argument(
            '--reproject', default='adaptive', choices=['adaptive', 'interp', 'exact'],
            help='The reprojection method of the reproject package.'
        )
        parser.add_argument(
            '--pixel_scale', default=None, type=float,
            help='Output pixel scale (arcsec).  If not provided, the finest of the images.'
        )
        parser.add_argument(
            '--rej', default=3., type=float, help='Outlier rejection threshold (sigma).'
        )
        parser.add_argument(
            '--rej_frac', default=0.1, type=float,
            help='Fractional error added in quadrature to the rejection threshold, so that small '
                 'differences between the images (e.g., from the PSF matching and scaling) do '
                 'not reject the cores of bright sources.'
        )
        parser.add_argument(
            '--no_reject', default=False, action='store_true', help='Do not reject outliers.'
        )
        parser.add_argument(
            '-o', '--overwrite', default=False, action='store_true',
            help='Rebuild the mosaics, even if they are up to date.'
        )
        parser.add_argument(
            '--no_qa', dest='qa', default=True, action='store_false',
            help='Do not write the QA figures.  By default, a figure of each mosaic is written '
                 'to the qa subdirectory of the output directory when the mosaic is '
                 'constructed, or if it lacks one.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.datamodel.dataset import reduction_dir
        from pipnick.pipelines.mosaic import mosaic_all
        from pipnick.pkg.report import StepReport

        # Initialize logger: in the output directory, or the reduction
        # directory of the (first) table
        root = reduction_dir(args.raw_table[0]) if args.outdir is None else args.outdir
        cls.init_log(args, log_dir=cls.log_dir(root))

        report = StepReport('MOSAICS')
        mosaic_all(
            args.raw_table, outdir=args.outdir, mosaics=args.mosaics, psf_match=args.psf_match,
            kernel=args.kernel, psf_margin=args.psf_margin, scale=args.scale,
            reference=args.reference, match_radius=args.match_radius,
            min_common=args.min_common, min_snr=args.min_snr, background=args.background,
            match_background=args.match_background,
            reproject=args.reproject, pixel_scale=args.pixel_scale,
            rej=None if args.no_reject else args.rej, rej_frac=args.rej_frac,
            overwrite=args.overwrite, qa=args.qa, report=report
        )
        return report.exit_status
