"""
The objects shared by the pipnick steps: the data set of a setup and its
products (:mod:`~pipnick.datamodel.dataset`), the reduced images
(:mod:`~pipnick.datamodel.image`), and the PSF models
(:mod:`~pipnick.datamodel.psf`).

The objects wrap the files that pipnick writes (header keywords, ECSV tables
and their metadata), whose formats they do not change.  This package never
imports :mod:`pipnick.pipelines`.
"""
