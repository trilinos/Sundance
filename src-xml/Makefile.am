# @HEADER
# @HEADER

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
#np# There is, at present, no distinction between most of the categories, 
#np# other than those which start with EXTRA_.  All categories which do not
#np# start with EXTRA_ are included in all builds.  Categories which start
#np# with EXTRA_ are included in some builds (typically controlled by 
#np# --with-xyz or --without-xyz)
#
#np# Make sure to list all source files in one of the following categories.

noinst_PROGRAMS = Sundance

XML_SOURCES = \
$(top_srcdir)/src-xml/SundanceExprParser.cpp \
$(top_srcdir)/src-xml/SundanceExprScanner.cpp \
$(top_srcdir)/src-xml/SundanceToken.cpp 

EXEEXT = .exe


Sundance_SOURCES = $(srcdir)/SundanceInterpreter.cpp $(XML_SOURCES)

Sundance_DEPENDENCIES = \
$(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a

Sundance_CXXFLAGS = $(STD_FWK_INCLUDES)

Sundance_LDADD = $(STD_FWK_LIBS) 

include $(top_srcdir)/Makefile.common







