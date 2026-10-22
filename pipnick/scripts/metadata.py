"""
Construct the table of metadata for a set of raw observations, which defines
how the data are reduced.
"""

from pipnick.scripts import scriptbase


class MetaData(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Construct a table with metadata pulled from data file headers.  Frames '
                        'that fail the quality checks are commented out of the table; comment '
                        'out (or restore) any other rows, or edit the table, before reducing '
                        'the data.',
            width=width
        )
        parser.add_argument(
            '-d', '--rawdir', type=str, default=None,
            help='Path to the directory with the raw data from a single night.'
        )
        parser.add_argument(
            '-t', '--raw_table', default=None,
            help='File for the metadata, or the base name for the per-setup files if --split '
                 'is set.  If not provided, a default name is used.'
        )
        parser.add_argument(
            '--split', default=False, action='store_true',
            help='Write one metadata file per instrument setup (read speed, binning, and '
                 'readout window), named <base>_<xbin>x<ybin>_<readmode>_<window>, where the '
                 'window is "full" for the full readout or win00, win01, etc. for custom '
                 'windows.  Per-setup files also have a comb_id column that groups the frames '
                 'to combine, and pointing and focus frames are commented out by default.'
        )
        parser.add_argument(
            '--camera', default=None,
            help='The camera used to obtain the data (e.g., NickelC2 or '
                 'NickelAndor).  If not provided, the camera is identified from '
                 'the file headers, and all files must be from the same camera.  '
                 'If provided, files from other cameras are skipped.'
        )
        parser.add_argument(
            '--ext', default='.fits', help='Extension of the raw files.'
        )
        parser.add_argument(
            '--classify', default='header', choices=['header', 'stats'],
            help='How to classify and check the frames.  "header" uses only the '
                 'file headers.  "stats" also uses statistics of the image data, '
                 'which is slower, but catches problems like empty images and '
                 'flats with too little light.'
        )
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            '-o', '--overwrite', default=False, action='store_true',
            help='Rebuild the table from scratch, overwriting an existing file '
                 '(and discarding any edits to it).'
        )
        group.add_argument(
            '-a', '--append', default=False, action='store_true',
            help='Add new files in the raw directory to an existing table, '
                 'preserving the existing rows, including any edits and '
                 'commented-out rows.'
        )
        return parser

    @classmethod
    def main(cls, args):
        from pathlib import Path
        from pipnick import log
        from pipnick.cameras import available_cameras
        from pipnick.metadata import build_metadata, find_setup_files

        # Initialize logger
        cls.init_log(args)

        if args.rawdir is None:
            log.error('Must provide the directory with the raw data (--rawdir).')
            return
        if args.camera is not None and args.camera not in available_cameras():
            log.error(
                f'Unknown camera {args.camera}.  Available cameras are: '
                f'{", ".join(available_cameras())}.'
            )
            return
        if args.raw_table is not None and not args.overwrite and not args.append:
            if args.split and len(find_setup_files(args.raw_table)) > 0:
                log.error(
                    f'Per-setup files for {args.raw_table} already exist.  Use --append to add '
                    'new files to them, or --overwrite to rebuild them.'
                )
                return
            if not args.split and Path(args.raw_table).is_file():
                log.error(
                    f'{args.raw_table} already exists.  Use --append to add new files to it, '
                    'or --overwrite to rebuild it.'
                )
                return

        build_metadata(
            rawdir=args.rawdir, filename=args.raw_table, ext=args.ext, overwrite=args.overwrite,
            append=args.append, camera=args.camera, classify=args.classify, split=args.split
        )
