# pipnick

<p align="center">
  <img src="docs/figures/pipnick-logo-with-text.png" alt="The pipnick logo" width="240">
</p>

[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue)](LICENSE.rst)

**Photometric pipeline for the 1-meter Nickel Telescope at Lick Observatory.**

`pipnick` reduces and calibrates data taken with the science cameras on the
Nickel Telescope at Lick Observatory on Mt. Hamilton, CA (the previous C2 camera
and the new Andor camera), from the raw images to calibrated source catalogs and
image mosaics:

- **Reduction:** bad-pixel and saturation masking, overscan subtraction and
  trimming, bias subtraction, flat-field division, cosmic-ray masking, the
  combination of repeated frames, and the propagation of the uncertainty of
  each pixel.
- **Astrometry:** astrometric solutions, found by matching the sources to the
  other images of a night or to the Gaia catalog, and only when needed by
  submitting an image to astrometry.net
  ([nova.astrometry.net](https://nova.astrometry.net/)).
- **Photometry:** PSF modeling with a circular or elliptical 2D Moffat
  function, and PSF and aperture photometry.
- **Photometric calibration:** zero points, extinction, and color terms from
  the Landolt standards.
- **Mosaics:** PSF-matched, flux-scaled mosaics of dithered images, sharing one
  WCS across filters, and three-color images.

`pipnick` was originally developed for an Akamai 2024 internship project.

## Status

`pipnick` is under active development and has not yet been released.
Interfaces may change without notice.

## Installation

`pipnick` requires Python 3.12 or later.  We strongly recommend installing it
in its own virtual environment:

```console
python3.13 -m venv ~/.virtualenvs/pipnick
source ~/.virtualenvs/pipnick/bin/activate
```

`pipnick` is not yet released on PyPI, so install the latest version directly
from GitHub:

```console
pip install "git+https://github.com/UCObservatories/pipnick.git"
```

To check the installation, run `pipnick_version --dependencies`, which prints
the installed versions of `pipnick` and its dependencies.

`pipnick_astrometry` submits some images to astrometry.net, which requires a
free account at [nova.astrometry.net](https://nova.astrometry.net/).  Keep your
API key secret: set it in the `ASTROMETRY_NET_API_KEY` environment variable
(or the `api_key` entry of the `[astrometry_net]` section of
`~/.astropy/config/astroquery.cfg`), never on the command line.

### Developer installation

To contribute, fork the repository, then clone and install your fork with the
testing and documentation dependencies:

```console
git clone https://github.com/<your_github_handle>/pipnick.git
cd pipnick
pip install -e ".[dev]"
```

New contributions must pass all the tests, which use synthetic data and never
contact astrometry.net or the Gaia archive:

```console
pytest --pyargs pipnick
```

and the documentation must build without warnings (see below).

## Usage

Processing is run with a sequence of command-line scripts, each a wrapper of
a Python function that can also be called directly (e.g., from a Jupyter
notebook).  Most scripts take the table of raw frames written by the first step
and edited by you.  A typical night, with one instrument setup, is processed
as follows:

```console
pipnick_metadata raw n3 --split --assign_mosaics
# Edit the n3_2x2_fast_full.ecsv file!
pipnick_reduce n3_2x2_fast_full.ecsv
pipnick_astrometry n3_2x2_fast_full.ecsv
pipnick_photometry n3_2x2_fast_full.ecsv
pipnick_calibrate n3_2x2_fast_full.ecsv  # If Landolt standards were observed
pipnick_mosaic n3_2x2_fast_full.ecsv
```

Each step reports whether it succeeded and suggests the next command, and a
rerun redoes only what has changed.  Use `-h` with any script for the full list
of options, and see the quick-start tutorial in the documentation.

| Script | Purpose |
| --- | --- |
| `pipnick_metadata` | Build the tables of the raw frames of a night, one per instrument setup |
| `pipnick_reduce` | Reduce the raw frames |
| `pipnick_astrometry` | Find the astrometric solutions of the reduced images |
| `pipnick_photometry` | Detect and measure the sources in the reduced images |
| `pipnick_calibrate` | Calibrate the photometry with the Landolt standards |
| `pipnick_mosaic` | Combine the reduced images into mosaics |
| `pipnick_rgb` | Combine three images with the same WCS into a three-color image |
| `pipnick_status` | Report the state of every step for every image, and the next command |
| `pipnick_report` | Write an HTML report of the processing, with every summary and QA figure |
| `pipnick_inspect` | Report the provenance of one image, and trace one of its sources |
| `pipnick_version` | Report the installed version of `pipnick` and, optionally, its dependencies |

## Documentation

The documentation is built with Sphinx:

```console
cd docs
make clean ; make html
```

The HTML pages are written to `docs/_build/html`.  Start with the installation
instructions and the quick-start tutorial under "Getting started".

Please [submit an issue](https://github.com/UCObservatories/pipnick/issues) if
you run into bugs or have feature requests.

## Citation

If you use `pipnick` in your research, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff).

## License

`pipnick` is distributed under the BSD 3-Clause License; see
[`LICENSE.rst`](LICENSE.rst).  See [`licenses/README.rst`](licenses/README.rst)
for the terms governing third-party code that `pipnick` adapts or depends on.
