"""
Tests for the closing verdicts of the steps.
"""

from pipnick.pkg.report import StepReport
from pipnick.pkg.version import version


def test_status():
    report = StepReport('TEST')
    assert report.status == 'OK' and report.exit_status == 0, 'A new report should be OK'
    report.note('something to know')
    assert report.status == 'OK', 'A note should not need attention'
    report.attention('something to check')
    assert report.status == 'NEEDS ATTENTION' and report.exit_status == 0, (
        'A problem should need attention, but not fail the step'
    )
    report.fail('a product failed')
    assert report.status == 'FAILED' and report.exit_status == 1, (
        'A failure should fail the step, with a non-zero exit status'
    )


def test_check_range():
    report = StepReport('TEST')
    assert report.check_range('FWHM', 1.5, (0.8, 3.5), ' arcsec'), (
        'A value inside the range should pass'
    )
    assert report.check_range('RMS', 0.2, None), 'A value without a range should pass'
    assert report.check_range('RMS', 0.4, (None, 0.5)), 'A value below a maximum should pass'
    assert not report.check_range('RMS', 0.7, (None, 0.5), ' arcsec', advice='Resolve it.'), (
        'A value above a maximum should fail'
    )
    assert report.problems == [
        'RMS is 0.70 arcsec, outside the expected range (at most 0.50 arcsec).  Resolve it.'
    ], 'The problem should give the value, the expected range, and the advice'


def test_versions(caplog):
    report = StepReport('TEST')
    report.check_version('current.fits', version)
    report.check_version('old.fits', '0.0.1')
    report.check_version('none.fits', None)
    assert report.mismatched == {'old.fits': '0.0.1', 'none.fits': None}, (
        'Products made by another version, or with none recorded, should be mismatched'
    )
    status = report.finish('3 products.', next_step='pipnick_next')
    assert status == 'OK', 'A version mismatch should not change the status'
    assert caplog.text.count('made with a different version of pipnick') == 2, (
        'The mismatch should be warned once, and repeated in the verdict'
    )
    assert 'an earlier version (none recorded)' in caplog.text, (
        'Products without a version should be described as made by an earlier version'
    )
    assert 'Next: pipnick_next' in caplog.text, 'The next step should be given'


def test_finish_failed(caplog):
    report = StepReport('TEST')
    report.fail('X failed')
    report.finish('1 product.', next_step='pipnick_next')
    assert 'TEST FAILED: 1 product.' in caplog.text, 'The verdict should give the status'
    assert 'FAILED: X failed' in caplog.text, 'The verdict should repeat the failures'
    assert 'Next:' not in caplog.text, 'No next step should be given for a failed step'
