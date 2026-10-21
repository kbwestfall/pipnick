Known Issues
============

The code is still in active development.  Here are a few known issues with the
code.  If you find more, please submit an issue on our GitHub page.

- Basic reduction does not always fully eliminate inconsistencies, likely
  due to flat-fielding issues. Occasionally, a dust particle is preserved
  or sensitivity still varies slightly across the pixel field.
- Photometry code currently detects multiple sources (a 'group' of sources)
  where there likely is only one source. Sometimes sources in a group have
  negative flux.
- Astrometry code is sometimes unable to solve images with very few stars.
