# pipnick

[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue)](LICENSE.rst)

**Photometric pipeline for the 1-meter Nickel Telescope at Lick Observatory.**

`pipnick` performs basic image processing and photometric analysis of data
taken with the science imaging camera on the Nickel Telescope at Lick
Observatory on Mt. Hamilton, CA:

- **Reduction:** overscan subtraction and trimming, bias subtraction,
  flat-field division, and cosmic-ray masking.
- **Photometry:** source detection, PSF modeling with a 2D Moffat function,
  and PSF and aperture photometry.
- **Astrometry:** astrometric calibration using
  [astrometry.net](https://nova.astrometry.net/).

`pipnick` was originally developed for an Akamai 2024 internship project.

## Status

`pipnick` is under active development and has not yet been released.
Interfaces may change without notice, and the astrometry pipeline and the
usage documentation are currently being updated.

## Installation

`pipnick` requires Python 3.12 or later.  Install it from a clone of the
repository:

```console
git clone https://github.com/allisoneto/pipnick.git
cd pipnick
pip install -e .
```

To also install the testing and documentation dependencies, use
`pip install -e ".[dev]"`.

## Usage

Processing is run with a set of command-line scripts.  Use `-h` with any of
them for the full list of options.

| Script | Purpose |
| --- | --- |
| `pipnick_metadata` | Build a metadata table from the headers of the raw data files |
| `pipnick_reduce` | Reduce the raw images |
| `pipnick_photometry` | Perform PSF and aperture photometry on the reduced images |
| `pipnick_astrometry` | Astrometrically calibrate the reduced images |
| `pipnick_version` | Report the installed version of `pipnick` and its key dependencies |

## Documentation

The documentation is built with Sphinx:

```console
cd docs
make html
```

The HTML pages are written to `docs/_build/html`.

## Citation

If you use `pipnick` in your research, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff).

## License

`pipnick` is distributed under the BSD 3-Clause License; see
[`LICENSE.rst`](LICENSE.rst).  See [`licenses/README.rst`](licenses/README.rst)
for the terms governing third-party code that `pipnick` adapts or depends on.
