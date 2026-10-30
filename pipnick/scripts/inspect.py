"""
Inspect a reduced science image and its sources.
"""

from pipnick.scripts import scriptbase


class Inspect(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Print the provenance of a reduced science image: its reduction, '
                        'astrometric solution, photometry, photometric calibration, and mosaic, '
                        'with its QA figures.  Optionally, trace one of its sources through the '
                        'source table, the removed sources, the calibrated source tables, and '
                        'the standard-star measurements, and show the image with its sources.  '
                        'Nothing is changed.',
            width=width
        )
        parser.add_argument(
            'raw_table', help='Per-setup metadata file.'
        )
        parser.add_argument(
            'image',
            help='The image: its file name (with or without the extension) or its combination '
                 'group (comb_id).'
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
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            '--source', type=str, default=None,
            help='Trace the source at these coordinates, given as RA,Dec in degrees, or '
                 'sexagesimal with the RA in hours (e.g., 18:51:05.0,-06:16:12).  The image '
                 'must have an astrometric solution.'
        )
        group.add_argument(
            '--id', type=int, default=None,
            help='Trace the source with this id in the source table of the image.'
        )
        parser.add_argument(
            '--radius', type=float, default=1.,
            help='Match radius in arcsec for --source.'
        )
        parser.add_argument(
            '--display', default=False, action='store_true',
            help='Show the image with its sources (labeled by id), the removed sources, and '
                 'the traced source in an interactive window.'
        )
        parser.add_argument(
            '--ofile', type=str, default=None,
            help='Write the image with its sources to this file instead of showing it.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pipnick.datamodel.dataset import reduction_dir
        from pipnick.pipelines.inspect import inspect_image, parse_coordinates

        # Initialize logger
        cls.init_log(args, log_dir=cls.log_dir(reduction_dir(args.raw_table, args.rdxdir)))

        ra, dec = (None, None) if args.source is None else parse_coordinates(args.source)
        inspect_image(
            args.raw_table, args.image, rdxdir=args.rdxdir, mosaic_dir=args.mosaic_dir, ra=ra,
            dec=dec, source_id=args.id, radius=args.radius, display=args.display,
            ofile=args.ofile
        )
