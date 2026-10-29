"""
The status, provenance, and report of the processing of a data set.

The package collects what every step did, from the products and their
provenance, without changing them: the state of each step for each science
image (:mod:`~pipnick.report.status`), the provenance of an image and of its
sources (:mod:`~pipnick.report.records`), and the rendering of the HTML
report (:mod:`~pipnick.report.html`).  It is used by ``pipnick_status``,
``pipnick_report``, and ``pipnick_inspect``.
"""
