# @HEADER

# @HEADER

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = Doxyfile.core Doxyfile.utils Doxyfile.trilinos \
	$(srcdir)/includes/footer.html \
	$(srcdir)/includes/header.html \
	$(srcdir)/includes/stylesheet.css

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

dox: 
	doxygen Doxyfile.utils
	doxygen Doxyfile.core

trilinos-dox:
	doxygen Doxyfile.trilinos



