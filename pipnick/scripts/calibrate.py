"""
Photometric calibration with the Landolt standard stars
"""

from pipnick.scripts import scriptbase


class Calibrate(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Find the Landolt standard stars in the photometry source tables and '
                        'fit the photometric calibration of each band (U, B, V, R, I) of each '
                        'night.  Run pipnick_astrometry and pipnick_photometry first: the '
                        'standards are identified by their celestial coordinates.  The '
                        'calibration is written to an editable TOML file in the calibrate '
                        'directory, with a table of the residuals of the standards.  If the '
                        'standards of a band do not cover a sufficient range of airmass, the '
                        'nominal extinction is used instead of a measurement, and this is '
                        'reported in all of the products (see --no_fallback).',
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
            help='Only search these images for standards, given by file name (with or without '
                 'the extension) or combination group (comb_id).  This does not limit the '
                 'images the calibration is applied to.'
        )
        parser.add_argument(
            '--terms', nargs='+', default=['zp', 'airmass', 'color'],
            choices=['zp', 'airmass', 'color', 'color2', 'airmass_color'],
            help='Terms of the calibration model: the zero point (required), the extinction '
                 '(airmass), the linear and quadratic color terms, and the color dependence '
                 'of the extinction (airmass_color).'
        )
        parser.add_argument(
            '--fix', nargs='+', default=None,
            help='Fix terms at given values, as <term>=<value>.  The value is a number used '
                 'for all bands, a comma-separated list of <band>:<number> pairs (e.g., '
                 'airmass=B:0.28,V:0.17), or, for the airmass term only, "nominal" for the '
                 'nominal extinction of the observatory.'
        )
        parser.add_argument(
            '--color', nargs='+', default=None,
            help='The color used to calibrate a band, as <band>=<color>, e.g. R=R-I.  The '
                 'defaults are U-B for U, B-V for B and V, V-R for R, and V-I for I.'
        )
        parser.add_argument(
            '--mag', default='psf', type=str, choices=['psf', 'aper'],
            help='The instrumental magnitudes to calibrate: PSF or aperture.  No aperture '
                 'correction is applied, so a calibration of aperture magnitudes is valid '
                 'only for magnitudes measured in the same aperture.'
        )
        parser.add_argument(
            '--match_radius', default=2., type=float,
            help='Radius (arcsec) for matching sources to the standards.'
        )
        parser.add_argument(
            '--max_neighbor_flux', default=0.01, type=float,
            help='Standards with a neighbor within their aperture brighter than this fraction '
                 'of their flux are not used.'
        )
        parser.add_argument(
            '--residual_radius', default=1.5, type=float,
            help='Sources detected in the residuals of a standard (in a later iteration of the '
                 'PSF photometry) within this radius, in FWHM, do not exclude it.'
        )
        parser.add_argument(
            '--max_cat_err', default=0.05, type=float,
            help='Standards with a larger catalog error (mag) in the band are not used.'
        )
        parser.add_argument(
            '--min_cat_obs', default=1, type=int,
            help='Standards observed fewer times by Landolt are not used.'
        )
        parser.add_argument(
            '--sys_err', default=0.01, type=float,
            help='Systematic error (mag) added in quadrature to each measurement.'
        )
        parser.add_argument(
            '--rej', default=3., type=float,
            help='Outlier rejection threshold, in standard deviations.'
        )
        parser.add_argument(
            '--maxiter', default=5, type=int, help='Maximum number of rejection iterations.'
        )
        parser.add_argument(
            '--min_stars', default=3, type=int,
            help='Minimum number of standard stars in each band.'
        )
        parser.add_argument(
            '--min_airmass_span', default=0.2, type=float,
            help='Minimum airmass span required to fit the extinction.'
        )
        parser.add_argument(
            '--min_color_span', default=0.3, type=float,
            help='Minimum color span (mag) required to fit the color terms.'
        )
        parser.add_argument(
            '--force', default=False, action='store_true',
            help='Downgrade the minimum-stars and color-span errors to warnings.'
        )
        parser.add_argument(
            '--no_fallback', default=False, action='store_true',
            help='If the airmass span is too small to fit the extinction, fail instead of '
                 'using the nominal extinction.'
        )
        parser.add_argument(
            '--drift_warn', default=0.02, type=float,
            help='Warn if the residuals drift by more than this (mag/hour), and the drift is '
                 'significant.  The drift is not fit.'
        )
        parser.add_argument(
            '--no_qa', dest='qa', default=True, action='store_false',
            help='Do not write the QA figures.  By default, a figure of the residuals of each '
                 'night is written to the calibrate/qa directory when the night is fit, or if '
                 'it lacks one.'
        )
        parser.add_argument(
            '-o', '--overwrite', default=False, action='store_true',
            help='Refit, even if the inputs have not changed or the calibration file has been '
                 'edited.'
        )
        parser.add_argument(
            '--apply', default=None, type=str, choices=['catalog', 'image', 'both'],
            help='Apply the calibration to all science images of each night whose band is '
                 'calibrated: write calibrated source tables to the calibrate directory '
                 '(catalog), record the calibration in the image headers (image), or both.'
        )
        parser.add_argument(
            '--default_color', nargs='+', default=None,
            help='Color assumed for sources without a measured color, as <color>=<value>, '
                 'e.g. B-V=0.6.  Colors not given default to 0.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.datamodel.dataset import reduction_dir
        from pipnick.pipelines.calibrate import calibrate_all, parse_assignments, parse_fix
        from pipnick.pkg.report import StepReport

        # Initialize logger
        cls.init_log(args, log_dir=cls.log_dir(reduction_dir(args.raw_table, args.rdxdir)))
        report = StepReport('CALIBRATION')

        default_color = {
            k: float(v) for k, v in parse_assignments(args.default_color, 'color').items()
        }
        calibrate_all(
            args.raw_table, rdxdir=args.rdxdir, images=args.images, terms=args.terms,
            fix=parse_fix(args.fix), colors=parse_assignments(args.color, 'color'),
            mag=args.mag, match_radius=args.match_radius,
            max_neighbor_flux=args.max_neighbor_flux, residual_radius=args.residual_radius,
            max_cat_err=args.max_cat_err,
            min_cat_obs=args.min_cat_obs, sys_err=args.sys_err, rej=args.rej,
            maxiter=args.maxiter, min_stars=args.min_stars,
            min_airmass_span=args.min_airmass_span, min_color_span=args.min_color_span,
            force=args.force, fallback=not args.no_fallback, drift_warn=args.drift_warn,
            qa=args.qa, overwrite=args.overwrite, apply=args.apply,
            default_color=default_color, report=report
        )
        return report.exit_status
