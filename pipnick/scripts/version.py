"""
Report the pipnick version and the versions of its key dependencies.
"""

from pipnick.scripts import scriptbase


class Version(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Report the pipnick version and the versions of its key dependencies.',
            width=width, include_log_options=False
        )
        parser.add_argument('--dependencies', default=False, action='store_true',
                            help='Also report the versions of the packages pipnick depends on.')
        return parser

    @staticmethod
    def main(args):
        import importlib.metadata

        import pipnick

        print(f'pipnick {pipnick.__version__}')
        if not args.dependencies:
            return

        packages = ['numpy', 'scipy', 'astropy', 'photutils', 'ccdproc', 'scikit-image',
                    'matplotlib']
        for package in packages:
            try:
                print(f'  {package:<14} {importlib.metadata.version(package)}')
            except importlib.metadata.PackageNotFoundError:
                print(f'  {package:<14} not installed')
