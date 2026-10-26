"""
Test configuration.

No test may submit a job to astrometry.net or query the Gaia archive.  In
addition to ``pytest-remotedata``, which blocks internet access for the test
session (``remote_data_strict = true`` in ``pyproject.toml``) unless ``pytest``
is run with ``--remote-data``, every method that would contact either service
is replaced by one that raises an error.
"""

import pytest
from astroquery.astrometry_net.core import AstrometryNetClass
from astroquery.utils.tap.core import TapPlus

from pipnick.astrometry.services import GaiaCatalog


def _blocked(*args, **kwargs):
    """
    Replacement for the methods that contact a remote service.
    """
    raise RuntimeError('Tests must not contact astrometry.net or the Gaia archive.')


@pytest.fixture(autouse=True)
def block_remote_services(monkeypatch):
    for name in ['solve_from_source_list', 'solve_from_image', 'monitor_submission', '_login']:
        monkeypatch.setattr(AstrometryNetClass, name, _blocked)
    monkeypatch.setattr(GaiaCatalog, '_query', _blocked)
    monkeypatch.setattr(TapPlus, 'launch_job', _blocked)
