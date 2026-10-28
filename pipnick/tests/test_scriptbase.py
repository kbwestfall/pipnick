"""
Tests for the logs and exit status of the scripts.
"""

import datetime
from argparse import Namespace

import pytest

from pipnick import log
from pipnick.pkg.exceptions import PipnickPathError
from pipnick.pkg.version import version
from pipnick.scripts import scriptbase
from pipnick.scripts.scriptbase import ScriptBase


class Failing(ScriptBase):
    """A script whose step fails."""

    @staticmethod
    def main(args):
        return 1


class Passing(ScriptBase):
    """A script whose step succeeds."""

    @staticmethod
    def main(args):
        return 0


def args(log_file=None):
    """The parsed logging arguments of a script."""
    return Namespace(verbosity=1, log_file=log_file, log_level=None)


class Fixed(datetime.datetime):
    """A clock stopped at 2026-10-28 15:36:14 UT."""

    @classmethod
    def now(cls, tz=None):
        return cls(2026, 10, 28, 15, 36, 14, 123456, tzinfo=tz)


def test_default_log_file(tmp_path, monkeypatch):
    monkeypatch.setattr(scriptbase.datetime, 'datetime', Fixed)
    path = ScriptBase.default_log_file(tmp_path / 'log')
    assert path.name == 'pipnick_scriptbase_20261028T153614.log', (
        'The default log should be named for the script and the time, to the second'
    )
    path.parent.mkdir()
    path.write_text('')
    assert ScriptBase.default_log_file(tmp_path / 'log').name == (
        'pipnick_scriptbase_20261028T153614_2.log'
    ), 'An existing log should never be reused; a counter should be appended'
    assert ScriptBase.log_dir(tmp_path) == tmp_path / 'log', (
        'The default log directory should be the log subdirectory'
    )


def test_init_log(tmp_path):
    ScriptBase.init_log(args(), log_dir=tmp_path / 'log')
    files = sorted((tmp_path / 'log').glob('*.log'))
    log.init()
    assert len(files) == 1, 'A log should be written by default'
    text = files[0].read_text()
    assert text.splitlines()[0].endswith(f'pipnick {version}'), (
        'The first line of the log should be the full pipnick version'
    )

    bespoke = tmp_path / 'mine.log'
    ScriptBase.init_log(args(str(bespoke)), log_dir=tmp_path / 'log')
    log.init()
    assert bespoke.is_file() and len(sorted((tmp_path / 'log').glob('*.log'))) == 1, (
        'A bespoke log file should be used instead of the default'
    )
    with pytest.raises(PipnickPathError, match='never overwritten'):
        ScriptBase.init_log(args(str(bespoke)))
    with pytest.raises(FileExistsError):
        log.init(log_file=bespoke)
    log.init()


def test_exit_status(monkeypatch):
    monkeypatch.setattr(Failing, 'parse_args', classmethod(lambda cls: None))
    monkeypatch.setattr(Passing, 'parse_args', classmethod(lambda cls: None))
    with pytest.raises(SystemExit) as e:
        Failing.entry_point()
    assert e.value.code == 1, 'A failed step should give a non-zero exit status'
    Passing.entry_point()
