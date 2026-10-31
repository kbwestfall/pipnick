
# TODO items

## Checks

- How are pre/over-scans treated for the new Andor camera?

- Make sure that every step can force a re-reduction; if not, docs
  should indicate that previous reductions need to be deleted.  Maybe
  provide a simple `pipnick_clean` convenience script?

- How often and where are globals used?  Minimize this.

- Performance benchmarks, particularly for the photometry.  I expect the
  slow-down is driven by photutils, but maybe there are some Moffat
  function evaluations we can speed up...

- Status of a failed band: when a band's calibration fails (here, I),
  pipnick_status reports its images as calibrated: missing and keeps
  suggesting pipnick_calibrate, though rerunning it can't help.

- Photometry: some sources removed as failed were detected far outside
  the image by the photometry's second iteration (4 of 39 in SA 111 1965
  V C021).

- Does the photometric calibration accept data from different nights?
  The documentation seems to indicate it does, but it only takes a
  single metadata table for a fixed setup.


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

- How are source groups treated in the PSF photometry?  (When are)/Can
  source groups combined into a single source for photometry?

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

