# Makefile for the pipnick documentation, structured after PypeIt's.
#
#   make html          regenerate the API pages, then build the HTML
#   make htmlonly      build the HTML from the files as they stand
#   make picky         as html, but treat every warning as an error and report
#                      every broken cross-reference
#   make apirst        regenerate the API pages only
#   make clean         remove the build and the generated API pages

BUILDDIR      = _build
SPHINXBUILD   = sphinx-build
SPHINXAPI     = sphinx-apidoc
SPHINXAPIOPT  = members,private-members,undoc-members,show-inheritance
SPHINXOPTS    = -b html -d $(BUILDDIR)/doctrees -aE -w ./sphinx_warnings.out
LOCALFILES    = $(BUILDDIR)/* api/*.rst sphinx_warnings.out

.PHONY: help clean apirst html htmlonly picky

help:
	@echo "Options are 'make <target>' where <target> is one of"
	@echo "  html       regenerate the API pages, then build the HTML"
	@echo "  htmlonly   build the HTML from the files as they stand"
	@echo "  picky      as html, but fail on any warning and report broken cross-references"
	@echo "  apirst     regenerate the API pages only"
	@echo "  clean      remove the build and the generated API pages"

clean:
	rm -rf $(LOCALFILES)

apirst:
	SPHINX_APIDOC_OPTIONS=$(SPHINXAPIOPT) $(SPHINXAPI) --separate --force --no-toc -o ./api \
		../pipnick ../pipnick/tests ../pipnick/pkg/version.py

html:
	$(MAKE) apirst
	$(SPHINXBUILD) $(SPHINXOPTS) . $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

htmlonly:
	$(SPHINXBUILD) $(SPHINXOPTS) . $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

picky:
	$(MAKE) apirst
	$(SPHINXBUILD) $(SPHINXOPTS) -n -W --keep-going . $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."
