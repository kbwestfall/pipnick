"""
Smoke tests of the command-line scripts: each script's ``main`` runs on
synthetic data, with options that must reach its pipeline.
"""

import matplotlib.image
from astropy.table import Table

from pipnick import log
from pipnick.datamodel.dataset import DataSet
from pipnick.scripts.astrometry import Astrometry
from pipnick.scripts.calibrate import Calibrate
from pipnick.scripts.inspect import Inspect
from pipnick.scripts.metadata import MetaData
from pipnick.scripts.mosaic import Mosaic
from pipnick.scripts.photometry import Photometry
from pipnick.scripts.reduce import Reduce
from pipnick.scripts.report import Report
from pipnick.scripts.status import Status
from pipnick.scripts.version import Version
from pipnick.tests import test_astrometry, test_calibrate, test_photometry, test_reduction
from pipnick.tests.test_mosaic import three_dithers, write_dithers


def run(script, options):
    """
    Run a script as from the command line, and close its log file.
    """
    status = script.main(script.parse_args([str(o) for o in options]))
    log.close_file()
    return status


def test_metadata_reduce_status_report(tmp_path):
    rawdir = tmp_path / 'raw'
    test_reduction.write_night(rawdir)
    assert run(MetaData, [rawdir, tmp_path / 'n5', '--split']) == 0, (
        'pipnick_metadata should succeed'
    )
    tables = sorted(tmp_path.glob('n5_*.ecsv'))
    assert len(tables) == 1, 'pipnick_metadata --split should write one table per setup'
    assert run(MetaData, [rawdir, tmp_path / 'n5.ecsv', '--split']) == 1, (
        'pipnick_metadata should refuse to overwrite existing tables without --overwrite, '
        'whether or not the base name has the .ecsv extension'
    )
    assert run(MetaData, [rawdir, tmp_path / 'n6', '--assign_mosaics']) == 1, (
        'pipnick_metadata should refuse --assign_mosaics without --split'
    )

    rdxdir = tmp_path / 'rdx'
    assert run(Reduce, [tables[0], '-r', rdxdir, '--no_qa']) == 0, (
        'pipnick_reduce should succeed'
    )
    assert len(sorted(rdxdir.glob('*.fits'))) > 0, 'pipnick_reduce should write the products'
    assert not (rdxdir / 'reduce' / 'qa').exists(), '--no_qa should write no QA figures'
    assert len(sorted((rdxdir / 'log').glob('pipnick_reduce_*.log'))) == 1, (
        'pipnick_reduce should write a log to the reduction directory'
    )

    ofile = tmp_path / 'status.ecsv'
    run(Status, [tables[0], '-r', rdxdir, '--ofile', ofile])
    status = Table.read(ofile)
    assert list(status['reduced']) == ['ok'] * len(status), (
        'pipnick_status should report the reduced images'
    )
    assert run(Report, [tables[0], '-r', rdxdir]) == 0, 'pipnick_report should succeed'
    assert (rdxdir / 'report' / 'index.html').is_file(), 'pipnick_report should write the report'


def test_astrometry_dry_run(tmp_path):
    table, rdxdir, _, _, _ = test_astrometry.write_field(tmp_path)
    # A dry run never submits anything to astrometry.net
    assert run(
        Astrometry, [table, '-r', rdxdir, '--dry_run', '--no_qa', '--max_sources', 20]
    ) == 0, 'pipnick_astrometry --dry_run should succeed'
    lists = sorted((rdxdir / 'astrometry' / 'dry_run').glob('*_sources.ecsv'))
    assert len(lists) == len(test_astrometry.IMAGES), (
        'A dry run should write the source list of each image'
    )
    assert all(len(Table.read(f)) == 20 for f in lists), '--max_sources should be applied'
    assert not (rdxdir / 'astrometry' / 'qa').exists(), '--no_qa should write no QA figures'


def test_photometry_script(tmp_path):
    table, rdxdir, _ = test_photometry.write_field(tmp_path)
    assert run(Photometry, [table, '-r', rdxdir, '--no_qa', '--thresh', 10]) == 0, (
        'pipnick_photometry should succeed'
    )
    tables = sorted((rdxdir / 'photometry').glob('*_phot.ecsv'))
    assert len(tables) == 1, 'pipnick_photometry should write the source table'
    assert Table.read(tables[0]).meta['settings']['thresh'] == 10., (
        'The options should reach the photometry'
    )
    assert not (rdxdir / 'photometry' / 'qa').exists(), '--no_qa should write no QA figures'


def test_calibrate_script(tmp_path):
    table, rdxdir, _ = test_calibrate.write_night(tmp_path, test_calibrate.wide_design())
    assert run(
        Calibrate, [table, '-r', rdxdir, '--apply', 'both', '--no_qa', '--fix', 'color=0.']
    ) == 0, 'pipnick_calibrate should succeed'
    content, _ = test_calibrate.read_calibration(test_calibrate.calibration_file(rdxdir))
    assert content['band']['V']['source'][2] == 'fixed', 'The options should reach the fit'
    assert len(sorted((rdxdir / 'calibrate').glob('*_cal.ecsv'))) > 0, (
        '--apply both should write the calibrated source tables'
    )
    assert not (rdxdir / 'calibrate' / 'qa').exists(), '--no_qa should write no QA figures'


def test_mosaic_and_inspect(tmp_path):
    table, rdxdir, _ = write_dithers(tmp_path, three_dithers())
    assert run(Mosaic, [table, '--no_qa', '--scale', 'none']) == 0, (
        'pipnick_mosaic should succeed'
    )
    mosaics = sorted((rdxdir / 'mosaic').glob('*.fits'))
    assert [f.name for f in mosaics] == ['M11_V.fits', 'M11_wcs.fits'], (
        'pipnick_mosaic should write the mosaic and the WCS of its group'
    )
    assert Table.read(mosaics[0], hdu='INPUTS')['scaled'].sum() == 3, (
        'The options should reach the mosaics'
    )
    assert not (rdxdir / 'mosaic' / 'qa').exists(), '--no_qa should write no QA figures'

    image = DataSet.from_file(table).science_images()[0]
    ofile = tmp_path / 'inspect.png'
    run(Inspect, [table, image.path.stem, '--id', 1, '--ofile', ofile])
    assert matplotlib.image.imread(ofile).ndim == 3, (
        'pipnick_inspect --ofile should write the image of the sources'
    )


def test_version(capsys):
    Version.main(Version.parse_args(['--dependencies']))
    out = capsys.readouterr().out
    assert out.startswith('pipnick ') and 'photutils' in out, (
        'pipnick_version should report the versions of pipnick and its dependencies'
    )


def test_help_defaults():
    text = Reduce.get_parser().format_help()
    assert '(default: True)' not in text, (
        'The help of --no_qa should not report the default of the setting it turns off'
    )
    assert '(default: False)' in text, (
        'The help of other options should report their defaults'
    )
