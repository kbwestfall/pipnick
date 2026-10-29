# pipnick

<p align="center">
  <img src="docs/figures/pipnick-logo-with-text.png" alt="The pipnick logo" width="240">
</p>

[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue)](LICENSE.rst)

**Photometric pipeline for the 1-meter Nickel Telescope at Lick Observatory.**

`pipnick` reduces and calibrates data taken with the science imaging camera on
the Nickel Telescope at Lick Observatory on Mt. Hamilton, CA, from the raw
images to calibrated source catalogs and image mosaics:

- **Reduction:** overscan subtraction and trimming, bias subtraction,
  flat-field division, cosmic-ray masking, and the propagation of the
  uncertainty of each pixel.
- **Astrometry:** astrometric solutions, from astrometry.net
  ([nova.astrometry.net](https://nova.astrometry.net/)), the Gaia catalog,
  and the other images of a night.
- **Photometry:** PSF modeling with a 2D Moffat function, and PSF and aperture
  photometry.
- **Photometric calibration:** zero points, extinction, and color terms from
  the Landolt standards.
- **Mosaics:** PSF-matched, flux-scaled mosaics of dithered images, and color
  images.

`pipnick` was originally developed for an Akamai 2024 internship project.

## Status

`pipnick` is under active development and has not yet been released.
Interfaces may change without notice.

## Installation

`pipnick` requires Python 3.12 or later.  Install it from a clone of the
repository:

```console
git clone https://github.com/UCObservatories/pipnick.git
cd pipnick
pip install -e .
```

To also install the testing and documentation dependencies, use
`pip install -e ".[dev]"`.

## Usage

Processing is run with a set of command-line scripts.  Use `-h` with any of
them for the full list of options, and see the documentation for a quick-start
tutorial.

| Script | Purpose |
| --- | --- |
| `pipnick_metadata` | Build the table of metadata of the raw frames of a night |
| `pipnick_reduce` | Reduce the raw frames |
| `pipnick_astrometry` | Find the astrometric solutions of the reduced images |
| `pipnick_photometry` | Measure the sources in the reduced images |
| `pipnick_calibrate` | Calibrate the photometry with the Landolt standards |
| `pipnick_mosaic` | Combine the reduced images into mosaics |
| `pipnick_rgb` | Combine three mosaics into a color image |
| `pipnick_status` | Report the state of every step for every image |
| `pipnick_report` | Write an HTML report of the processing |
| `pipnick_inspect` | Report the provenance of one image, and trace its sources |
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
