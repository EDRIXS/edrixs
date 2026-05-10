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

.PHONY: help Makefile html doxygen clean

# Run Doxygen before building HTML so Breathe has fresh XML input.
# After the build, sphinx-gallery's intermediate files (copied .py/.ipynb/.rst
# under source/auto_examples and source/sg_execution_times.rst) are removed:
# the rendered HTML/zip artefacts already live in build/html/ and the source-
# tree copies just clutter the working directory.  Cleanup runs whether the
# build succeeded or failed, and the original sphinx-build exit code is
# preserved.
html: doxygen
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O); ec=$$?; \
	rm -rf "$(SOURCEDIR)/auto_examples"; \
	rm -f  "$(SOURCEDIR)/sg_execution_times.rst"; \
	exit $$ec

doxygen:
	cd .. && doxygen Doxyfile

# Clean removes the Sphinx build output, the sphinx-gallery thumbnail cache
# under source/auto_examples/images/thumb, and the runtime artefacts that
# sphinx-gallery leaves behind in examples/sphinx/ when it executes the
# example scripts (Fock bases, Hamiltonian inputs, eigenvectors, Lanczos
# pole files, density matrices, spectra).  The .py example scripts and
# README.rst are preserved.
clean:
	@$(SPHINXBUILD) -M clean "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@rm -rf "$(SOURCEDIR)/auto_examples/images/thumb"
	@rm -f ../examples/sphinx/*.in
	@rm -f ../examples/sphinx/*.dat
	@rm -f ../examples/sphinx/eigvec.*
	@rm -f ../examples/sphinx/xas_poles.*

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
