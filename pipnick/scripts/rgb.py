"""
Construct a three-colour image from three images with the same WCS
"""

from pipnick.scripts import scriptbase


class RGB(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Construct a three-colour (RGB) image from three images, e.g., the '
                        'mosaics of one WCS group made by pipnick_mosaic.  The images must have '
                        'the same WCS and shape; if they do not, run pipnick_mosaic with the '
                        'images in the same WCS group (a mosaic can have a single image).',
            width=width
        )
        parser.add_argument('red', type=str, help='The image for the red channel.')
        parser.add_argument('green', type=str, help='The image for the green channel.')
        parser.add_argument('blue', type=str, help='The image for the blue channel.')
        parser.add_argument(
            '-o', '--ofile', type=str, default='rgb.png',
            help='The output image; its format (e.g., png or jpg) is set by its extension.'
        )
        parser.add_argument(
            '--method', default='lupton', choices=['lupton', 'linear', 'log'],
            help='The stretch: the arcsinh scheme of Lupton et al. (2004), linear, or '
                 'logarithmic.'
        )
        parser.add_argument(
            '--stretch', default=0.5, type=float,
            help='The linear stretch of the Lupton scheme, in units of the normalized images.'
        )
        parser.add_argument(
            '--Q', default=8., type=float, help='The asinh softening of the Lupton scheme.'
        )
        parser.add_argument(
            '--log_a', default=1000., type=float, help='The parameter of the log stretch.'
        )
        parser.add_argument(
            '--percentile', default=99.5, type=float,
            help='Each image is divided by this percentile of its pixels, after subtracting '
                 'its background.'
        )
        parser.add_argument(
            '--common', default=False, action='store_true',
            help='Divide all of the images by the largest of their percentiles, to keep their '
                 'relative fluxes, instead of normalizing each separately.'
        )
        parser.add_argument(
            '--no_subtract', default=False, action='store_true',
            help='Do not subtract the background (median) of each image.'
        )
        parser.add_argument(
            '--scales', nargs=3, default=[1., 1., 1.], type=float,
            help='Factors applied to the normalized red, green, and blue images, to adjust the '
                 'colour balance.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick import log
        from pipnick.utils.rgb import (
            check_same_wcs, color_image, read_image, write_color_image
        )

        # Initialize logger
        cls.init_log(args)

        paths = [args.red, args.green, args.blue]
        images = [read_image(p) for p in paths]
        check_same_wcs(paths, [i[2] for i in images], [i[0].shape for i in images])
        rgb = color_image(
            images[0][:2], images[1][:2], images[2][:2], method=args.method,
            stretch=args.stretch, Q=args.Q, log_a=args.log_a, scales=args.scales,
            percentile=args.percentile, common=args.common, subtract=not args.no_subtract
        )
        write_color_image(args.ofile, rgb)
        log.info(f'Wrote {args.ofile}')
