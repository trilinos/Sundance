#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = SimpleEval StdTransformation ListEval Integrals \
EvalCache

LIBDEPS = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

SimpleEval_SOURCES = $(srcdir)/SimpleEval.cpp 

SimpleEval_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

SimpleEval_CXXFLAGS = $(CORE_INCLUDES)

SimpleEval_LDADD = $(CORE_LIBS) 



StdTransformation_SOURCES = $(srcdir)/StdTransformation.cpp 

StdTransformation_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

StdTransformation_CXXFLAGS = $(CORE_INCLUDES)

StdTransformation_LDADD = $(CORE_LIBS) 


ListEval_SOURCES = $(srcdir)/ListEval.cpp 

ListEval_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

ListEval_CXXFLAGS = $(CORE_INCLUDES)

ListEval_LDADD = $(CORE_LIBS) 


Integrals_SOURCES = $(srcdir)/Integrals.cpp 

Integrals_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

Integrals_CXXFLAGS = $(CORE_INCLUDES)

Integrals_LDADD = $(CORE_LIBS)


EvalCache_SOURCES = $(srcdir)/EvalCache.cpp 

EvalCache_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

EvalCache_CXXFLAGS = $(CORE_INCLUDES)

EvalCache_LDADD = $(CORE_LIBS) 



include $(top_srcdir)/Makefile.common
