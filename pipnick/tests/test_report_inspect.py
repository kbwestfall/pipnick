"""
Tests for the report and the inspection of the products of a data set.
"""

import re

import matplotlib.image
import numpy as np
import pytest
from astropy.table import MaskedColumn, Table

from pipnick.datamodel.dataset import DataSet
from pipnick.pipelines import mosaic
from pipnick.pipelines.calibrate import calibrate_all
from pipnick.pipelines.inspect import inspect_image, parse_coordinates
from pipnick.pipelines.report import write_report
from pipnick.pipelines.status import status_all
from pipnick.pkg.exceptions import PipnickDataModelError, PipnickParameterError
from pipnick.pkg.paths import removed_table
from pipnick.report.html import relative, status_class
from pipnick.report.records import format_record, image_record
from pipnick.tests.test_calibrate import narrow_design, wide_design, write_night
from pipnick.tests.test_mosaic import three_dithers, write_dithers


def broken_links(root):
    """
    Return the links in the pages of a report that do not resolve to a file.
    """
    broken = []
    for page in sorted(root.rglob('*.html')):
        for link in re.findall(r'(?:href|src)="([^"]+)"', page.read_text()):
            if not (page.parent / link).exists():
                broken.append(f'{page.name}: {link}')
    return broken


def test_status_class():
    assert status_class('ok') == 'good' and status_class('up to date') == 'good', (
        'Successful states should be good'
    )
    assert status_class('stale') == 'warn' and status_class('failed') == 'bad', (
        'Stale and failed states should be flagged'
    )
    assert status_class('M11') == '', 'Other values should not be colored'


def test_relative(tmp_path):
    assert relative(tmp_path / 'a' / 'b.png', tmp_path / 'report' / 'images') == '../../a/b.png', (
        'Links should be relative to the page'
    )


def test_parse_coordinates():
    assert parse_coordinates('282.77,-6.27') == pytest.approx((282.77, -6.27)), (
        'Coordinates in degrees should be parsed'
    )
    ra, dec = parse_coordinates('18:51:04.8, -06:16:12')
    assert ra == pytest.approx(282.77) and dec == pytest.approx(-6.27), (
        'Sexagesimal coordinates, with the RA in hours, should be parsed'
    )
    with pytest.raises(PipnickParameterError):
        parse_coordinates('282.77')


def test_report_and_inspect_mosaics(tmp_path):
    table, rdxdir, _ = write_dithers(tmp_path, three_dithers())
    mosaic.mosaic_all(table)
    # A source removed by the photometry, away from the measured sources
    removed = Table()
    removed['x_init'] = [5.]
    removed['y_init'] = [5.]
    removed['x_fit'] = MaskedColumn([0.], mask=[True])
    removed['y_fit'] = MaskedColumn([0.], mask=[True])
    removed['flux_fit'] = MaskedColumn([0.], mask=[True])
    removed['reason'] = ['failed']
    images = DataSet.from_file(table).science_images()
    removed.write(removed_table(images[0].path), format='ascii.ecsv')

    index = write_report(table)
    root = index.parent
    assert root == rdxdir.absolute() / 'report', 'The report should be in the reduction directory'
    pages = sorted(f.name for f in root.rglob('*.html'))
    assert len(pages) == 6 + len(images), 'There should be a page per step and per image'
    assert broken_links(root) == [], 'Every link of the report should resolve'
    text = index.read_text()
    assert all(f'images/{i.path.stem}.html' in text for i in images), (
        'The overview should link to the page of each image'
    )
    assert 'class="good"' in text, 'The states should be colored'
    assert 'M11_V.png' in (root / 'mosaic.html').read_text(), (
        'The mosaic page should show the QA figure of the mosaic'
    )

    # The record of an image collects its provenance
    record = image_record(images[0])
    assert record.sections['Astrometry'][1] == ('method', 'gaia'), (
        'The record should give the astrometric method'
    )
    mosaic_entries = dict(record.sections['Mosaic'])
    assert mosaic_entries['status'] == 'ok' and mosaic_entries['mosaic'] == 'M11_V', (
        'The record should give the mosaic of the image'
    )
    assert 'flux scale' in mosaic_entries, 'The record should give the scale in the mosaic'
    assert 'QA figures' in format_record(record), 'The text record should list the figures'

    # Trace a source by id, by coordinates, and a removed source
    _, trace = inspect_image(table, images[0].name, source_id=3)
    phot = images[0].photometry()
    assert dict(trace.sections['Source table'])['id'] == '3', 'The source should be found'
    ra, dec = float(phot['ra'][2]), float(phot['dec'][2])
    _, trace = inspect_image(table, images[0].path.stem, ra=ra, dec=dec)
    assert dict(trace.sections['Source table'])['id'] == '3', (
        'The source should be found by its coordinates'
    )
    coo = images[0].wcs.pixel_to_world(5., 5.)
    _, trace = inspect_image(
        table, images[0].name, ra=float(coo.ra.deg), dec=float(coo.dec.deg)
    )
    assert dict(trace.sections['Removed sources'])['reason'] == 'failed', (
        'A removed source should be reported with its reason'
    )
    ofile = tmp_path / 'inspect.png'
    inspect_image(table, images[0].name, source_id=3, ofile=ofile)
    assert matplotlib.image.imread(ofile).ndim == 3, 'The display should be written to a file'
    with pytest.raises(PipnickParameterError):
        inspect_image(table, 'nothing')

    # The status, with the mosaics in another directory
    tbl, _ = status_all(table, mosaic_dir=tmp_path / 'elsewhere')
    assert list(tbl['mosaic']) == ['missing'] * 3, (
        'Mosaics should be looked for in the given directory'
    )


def test_report_calibration(tmp_path):
    table, rdxdir, _ = write_night(tmp_path, wide_design(), science=True)
    calibrate_all(table, rdxdir=rdxdir, apply='both')
    index = write_report(table, rdxdir=rdxdir)
    root = index.parent
    assert broken_links(root) == [], 'Every link of the report should resolve'
    assert '_residuals.png' in (root / 'calibrate.html').read_text(), (
        'The calibration page should show the residuals'
    )

    # A standard used by the calibration
    meas = Table.read(sorted((rdxdir / 'calibrate').glob('*_residuals.ecsv'))[0])
    row = meas[np.asarray(meas['used'])][0]
    _, trace = inspect_image(table, str(row['image']), rdxdir=rdxdir, source_id=row['source_id'])
    standard = dict(trace.sections['Standard star'])
    assert standard['star'] == row['star'] and standard['measurement'] == 'used', (
        'A standard should be traced to its measurement'
    )
    assert 'Calibrated source table' in trace.sections, (
        'The calibrated magnitude of a source should be traced'
    )

    # A science source, linked to the other images of its object
    science = [
        i for i in DataSet.from_file(table, rdxdir=rdxdir).science_images()
        if str(i.values['object']).startswith('M11')
    ]
    _, trace = inspect_image(table, science[0].name, rdxdir=rdxdir, source_id=1)
    linked = [k for k, _ in trace.sections['Calibrated source table'] if k.startswith('in ')]
    assert len(linked) == len(science) - 1, (
        'A science source should be found in the other images of its object'
    )
    with pytest.raises(PipnickDataModelError, match='coordinates'):
        # The synthetic images have no WCS
        inspect_image(table, science[0].name, rdxdir=rdxdir, ra=282.77, dec=-6.27)


def test_report_fallback_banner(tmp_path):
    table, rdxdir, _ = write_night(tmp_path, narrow_design())
    calibrate_all(table, rdxdir=rdxdir, qa=False)
    root = write_report(table, rdxdir=rdxdir).parent
    assert 'Extinction not measured' in (root / 'calibrate.html').read_text(), (
        'The calibration page should list the bands whose extinction was not measured'
    )
