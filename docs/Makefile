# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
# SPHINXOPTS    = "-W"  # This flag turns warnings into errors.
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = edrixs
SOURCEDIR     = source
BUILDDIR      = build
DOXYFILE      = ../Doxyfile

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile html doxygen

# Run Doxygen before building HTML so Breathe has fresh XML input.
html: doxygen
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

doxygen:
	cd .. && doxygen Doxyfile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
