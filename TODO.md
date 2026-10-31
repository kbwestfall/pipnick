
# TODO items

## Checks

- Performance benchmarks, particularly for the photometry.  I expect the
  slow-down is driven by photutils, but maybe there are some Moffat
  function evaluations we can speed up...

## Docs

- Review all docstrings and particularly any method/function
  cross-references.

- Provide pseudo-examples that essentially strip most of the
  infrastructure of the pipeline and boil the processing down to the
  critical functionality used.

- Construct a notebook that effectively performs the basic reductions
  and astrometric solution.


## Refactor

- Are there any objects that have significant overlap with astropy
  objects?  Can these be reconciled?


## Enhancements

- Investigate the details regarding how source groups are treated in the PSF
  photometry.  This may need refinement.  Based on the QA figures for M11, it
  seems like the code is finding too many clustered sources and missing obvious
  faint sources.

- Apart from the toml files in the data/cameras directory, is there
  anything in the code that is specific to the Nickel telescope?  I.e.,
  can we provide an interface for people that want to use pipnick with
  data from different telescopes?

- Restructure the scripts so that they run via a super script; e.g.,
  `pipnick reduce`, etc?

- Create a single toml configuration file that is used throughout the
  reductions; i.e., all settings are defined in one place.

- The photometry is currently not deterministic.

- Enable cross-matching between different sets of detected source.
  Basically, I want to be able to isolate and compare the photometry for
  the stars used to construct the PSF from the rest.

- Replace the settings classes with dc3-style ParSet subclasses.

- suppress the FITSFixedWarning warning

- Add a scatter plot of residuals to the astrometry QA plots

