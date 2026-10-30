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
            'rawdir', type=str,
            help='Path to the directory with the raw data from a single night.'
        )
        parser.add_argument(
            'raw_table', type=str,
            help='Name of the metadata file, with or without the .ecsv extension (e.g., n3 or '
                 'n3.ecsv, both written to n3.ecsv), or the base name for the per-setup files '
                 'if --split is set.'
        )
        parser.add_argument(
            '--split', default=False, action='store_true',
            help='Write one metadata file per instrument setup (read speed, binning, and '
                 'readout window), named <base>_<xbin>x<ybin>_<readmode>_<window>.ecsv, where '
                 '<base> is raw_table without the .ecsv extension (e.g., '
                 'n3_2x2_fast_full.ecsv for n3), and the window is "full" for the full readout '
                 'or win00, win01, etc. for custom windows.  Per-setup files also have a '
                 'comb_id column that groups the frames to combine, and pointing and focus '
                 'frames are commented out by default.'
        )
        parser.add_argument(
            '--assign_mosaics', default=False, action='store_true',
            help='Assign default mosaics (requires --split): the science frames are grouped by '
                 'pointing and, within each pointing, by filter, and labeled in the mosaic and '
                 'wcs_group columns.  Only frames without labels are assigned.  Without this '
                 'option, the columns are empty, i.e., no mosaics are defined.  With --append, '
                 'the labels are assigned even if there are no new frames.'
        )
        parser.add_argument(
            '--mosaic_radius', default=None, type=float,
            help='Linking length (arcmin) of the pointing groups used by --assign_mosaics.  '
                 'If not provided, half the field of view is used.'
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
        from pipnick.metadata.build import build_metadata
        from pipnick.metadata.names import table_file_name
        from pipnick.metadata.setups import find_setup_files
        from pipnick.pkg.report import StepReport

        # Initialize logger, in the default reduction directory (the parent of
        # the raw directory)
        cls.init_log(args, log_dir=cls.log_dir(Path(args.rawdir).absolute().parent))

        if args.assign_mosaics and not args.split:
            log.error('--assign_mosaics requires --split.')
            return 1
        if args.camera is not None and args.camera not in available_cameras():
            log.error(
                f'Unknown camera {args.camera}.  Available cameras are: '
                f'{", ".join(available_cameras())}.'
            )
            return 1
        table = table_file_name(args.raw_table)
        if not args.overwrite and not args.append:
            if args.split and len(find_setup_files(table)) > 0:
                log.error(
                    f'Per-setup files for {table} already exist.  Use --append to add new files '
                    'to them, or --overwrite to rebuild them.'
                )
                return 1
            if not args.split and table.is_file():
                log.error(
                    f'{table} already exists.  Use --append to add new files to it, or '
                    '--overwrite to rebuild it.'
                )
                return 1

        report = StepReport('METADATA')
        files, tbl = build_metadata(
            rawdir=args.rawdir, filename=table, ext=args.ext, overwrite=args.overwrite,
            append=args.append, camera=args.camera, classify=args.classify, split=args.split,
            assign=args.assign_mosaics, mosaic_radius=args.mosaic_radius, report=report
        )
        if args.classify == 'header':
            report.note(
                'The checks that need the pixel data (bias level, flat illumination, '
                'saturation) were skipped; use --classify stats to apply them.'
            )
        files = files if isinstance(files, list) else [files]
        names = ', '.join(Path(f).name for f in files)
        following = (
            f'review {names} (comment out bad frames with "#", and check the object names and '
            'comb_id groups), then run: '
        )
        if args.split:
            following += ' or '.join(f'pipnick_reduce {f}' for f in files)
        else:
            following += (
                f'pipnick_metadata {args.rawdir} {args.raw_table} --split to write a table for '
                'each setup'
            )
        report.finish(f'{len(tbl)} frames in {len(files)} table(s).', next_step=following)
        return report.exit_status
