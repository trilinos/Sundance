#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = FunctorTest

LIBDEPS = $(top_builddir)/src-utils/libsundanceUtils.a 

FunctorTest_SOURCES = $(srcdir)/FunctorTest.cpp 

FunctorTest_DEPENDENCIES = $(top_builddir)/src-utils/libsundanceUtils.a 

FunctorTest_CXXFLAGS = $(UTIL_INCLUDES)

FunctorTest_LDADD = $(UTIL_LIBS) 





include $(top_srcdir)/Makefile.common
