#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = SimpleEval TestEval StdTransformation VecOpTest \
ScalabilityTest DepthScalabilityTest VecTimings


LIBDEPS = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

SimpleEval_SOURCES = $(srcdir)/SimpleEval.cpp 

SimpleEval_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

SimpleEval_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

SimpleEval_LDADD = $(CORE_LIBS) 


VecOpTest_SOURCES = $(srcdir)/VecOpTest.cpp 

VecOpTest_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

VecOpTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

VecOpTest_LDADD = $(CORE_LIBS) 




VecTimings_SOURCES = $(srcdir)/VecTimings.cpp 

VecTimings_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

VecTimings_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

VecTimings_LDADD = $(CORE_LIBS) 






ScalabilityTest_SOURCES = $(srcdir)/ScalabilityTest.cpp 

ScalabilityTest_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

ScalabilityTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

ScalabilityTest_LDADD = $(CORE_LIBS) 



DepthScalabilityTest_SOURCES = $(srcdir)/DepthScalabilityTest.cpp 

DepthScalabilityTest_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

DepthScalabilityTest_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

DepthScalabilityTest_LDADD = $(CORE_LIBS) 




TestEval_SOURCES = $(srcdir)/TestEval.cpp 

TestEval_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

TestEval_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

TestEval_LDADD = $(CORE_LIBS) 



StdTransformation_SOURCES = $(srcdir)/StdTransformation.cpp 

StdTransformation_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-utils/libsundanceUtils.a 

StdTransformation_CXXFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

StdTransformation_LDADD = $(CORE_LIBS) 


#ListEval_SOURCES = $(srcdir)/ListEval.cpp 
#
#ListEval_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
#$(top_builddir)/src-utils/libsundanceUtils.a 
#
#ListEval_CXXFLAGS = $(CORE_INCLUDES)
#
#ListEval_LDADD = $(CORE_LIBS) 
#
#
#Integrals_SOURCES = $(srcdir)/Integrals.cpp 
#
#Integrals_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
#$(top_builddir)/src-utils/libsundanceUtils.a 
#
#Integrals_CXXFLAGS = $(CORE_INCLUDES)
#
#Integrals_LDADD = $(CORE_LIBS)
#
#
#EvalCache_SOURCES = $(srcdir)/EvalCache.cpp 
#
#EvalCache_DEPENDENCIES = $(top_builddir)/src-core/libsundanceCore.a \
#$(top_builddir)/src-utils/libsundanceUtils.a 
#
#EvalCache_CXXFLAGS = $(CORE_INCLUDES)
#
#EvalCache_LDADD = $(CORE_LIBS) 



include $(top_srcdir)/Makefile.common
