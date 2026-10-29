"""
Tests for the bookkeeping shared by the loops of the processing steps.
"""

import logging

import pytest

from pipnick.pkg.report import StepReport
from pipnick.pkg.step import Attempt, count_statuses, needs_qa, progress, reuse
from pipnick.pkg.version import version


def test_progress_and_counts():
    assert progress(2, 17) == '[3/17]', 'The progress label should be 1-based'
    rows = [{'status': 'measured'}, {'status': 'up to date'}, {'status': 'measured'}]
    assert count_statuses(rows, ['measured', 'up to date', 'failed']) == (
        '2 measured, 1 up to date'
    ), 'Statuses without items should be omitted by default'
    assert count_statuses(rows, ['measured', 'failed'], keep_zero=True) == (
        '2 measured, 0 failed'
    ), 'Statuses without items should be listed if requested'


def test_reuse(caplog):
    report = StepReport('TEST')
    with caplog.at_level(logging.INFO):
        reuse(report, 'a.fits', version, label='[1/2]')
    assert '[1/2] a.fits is up to date.' in caplog.text, 'The reuse should be logged'
    reuse(report, 'b.fits', '0.0.1')
    assert len(report.mismatched) == 1, 'A product made by another version should be recorded'


def test_attempt(caplog):
    report = StepReport('TEST')
    row = {'reason': '-'}
    attempt = Attempt(report, 'Photometry of a.fits', row=row, then='the other images')
    with attempt:
        raise ValueError('no stars')
    assert attempt.failed and row['reason'] == 'ValueError: no stars', (
        'A failure should be recorded in the summary row'
    )
    assert report.failures == [
        'Photometry of a.fits failed: ValueError: no stars.  Continuing with the other images.'
    ], 'A failure should be recorded in the report'

    attempt = Attempt(report, 'The QA figure of a.fits', problem=True)
    with attempt:
        raise RuntimeError('no display')
    assert attempt.failed and len(report.problems) == 1 and len(report.failures) == 1, (
        'A failed QA figure should be a problem, not a failure'
    )

    attempt = Attempt(report, 'Nothing')
    with attempt:
        pass
    assert not attempt.failed, 'An attempt without an exception should not fail'
    with pytest.raises(KeyboardInterrupt):
        with Attempt(report, 'Interrupted'):
            raise KeyboardInterrupt


def test_needs_qa(tmp_path):
    figure = tmp_path / 'a.png'
    assert not needs_qa(False, True, [figure]), 'No figures should be drawn without QA'
    assert needs_qa(True, False, [figure]), 'A missing figure should be drawn'
    figure.write_text('')
    assert not needs_qa(True, False, [figure]), 'An existing figure should not be redrawn'
    assert needs_qa(True, True, [figure]), 'The figure of a rebuilt product should be redrawn'
