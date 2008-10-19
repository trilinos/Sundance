# @HEADER
# ************************************************************************
# 
#                              Sundance
#                 Copyright (2005) Sandia Corporation
# 
# Copyright (year first published) Sandia Corporation.  Under the terms 
# of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
# retains certain rights in this software.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#                                                                                 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA                                                                                
# Questions? Contact Kevin Long (krlong@sandia.gov), 
# Sandia National Laboratories, Livermore, California, USA
# 
# ************************************************************************
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




SUNDANCE_ELEMENTS_H = \
$(top_srcdir)/src-std-framework/Elements/SundanceBasisDOFTopologyBase.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceBasisFamilyBase.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceBasisFamily.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceBasisReferenceEvaluationBase.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceClosedNewtonCotes.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceClosedNewtonCotesType.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceDiscreteFunction.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceDiscreteFunctionData.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceDiscreteSpace.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceFIATExpansion.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceFIATLagrange.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceFIATP1NC.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceFuncWithBasis.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceGauss1D.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceGaussianQuadrature.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceGaussianQuadratureType.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceLagrange.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceNedelec.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceP1NC.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceQuadratureFamily.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceQuadratureFamilyBase.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceQuadratureType.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceQuadratureTypeBase.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceTestFunction.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceTestFunctionData.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceTetQuadrature.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceTriangleQuadrature.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceUnknownFunction.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceUnknownFunctionData.hpp \
$(top_srcdir)/src-std-framework/Elements/SundanceVectorBasisComponent.hpp

SUNDANCE_ELEMENTS = \
$(top_srcdir)/src-std-framework/Elements/SundanceBasisFamily.cpp\
$(top_srcdir)/src-std-framework/Elements/SundanceBasisFamilyBase.cpp\
$(top_srcdir)/src-std-framework/Elements/SundanceClosedNewtonCotes.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceClosedNewtonCotesType.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceDiscreteFunction.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceDiscreteFunctionData.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceDiscreteSpace.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceFIATExpansion.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceFIATLagrange.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceFIATP1NC.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceGauss1D.cpp\
$(top_srcdir)/src-std-framework/Elements/SundanceGaussianQuadrature.cpp\
$(top_srcdir)/src-std-framework/Elements/SundanceGaussianQuadratureType.cpp\
$(top_srcdir)/src-std-framework/Elements/SundanceLagrange.cpp\
$(top_srcdir)/src-std-framework/Elements/SundanceNedelec.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceP1NC.cpp\
$(top_srcdir)/src-std-framework/Elements/SundanceQuadratureFamily.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceQuadratureFamilyBase.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceQuadratureType.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceTestFunction.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceTestFunctionData.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceTetQuadrature.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceTriangleQuadrature.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceUnknownFunction.cpp \
$(top_srcdir)/src-std-framework/Elements/SundanceUnknownFunctionData.cpp  \
$(top_srcdir)/src-std-framework/Elements/SundanceVectorBasisComponent.cpp

SUNDANCE_DOFMAPPING_H = \
$(top_srcdir)/src-std-framework/DOFMapping/SundanceBlock.hpp\
$(top_srcdir)/src-std-framework/DOFMapping/SundanceDOFMapBase.hpp\
$(top_srcdir)/src-std-framework/DOFMapping/SundanceDOFMapBuilder.hpp\
$(top_srcdir)/src-std-framework/DOFMapping/SundanceInhomogeneousNodalDOFMap.hpp \
$(top_srcdir)/src-std-framework/DOFMapping/SundanceMapStructure.hpp\
$(top_srcdir)/src-std-framework/DOFMapping/SundanceMixedDOFMap.hpp \
$(top_srcdir)/src-std-framework/DOFMapping/SundanceNodalDOFMap.hpp \
$(top_srcdir)/src-std-framework/DOFMapping/SundancePartialElementDOFMap.hpp \
$(top_srcdir)/src-std-framework/DOFMapping/SundanceSpatiallyHomogeneousDOFMapBase.hpp 

SUNDANCE_DOFMAPPING = \
$(top_srcdir)/src-std-framework/DOFMapping/SundanceDOFMapBase.cpp\
$(top_srcdir)/src-std-framework/DOFMapping/SundanceDOFMapBuilder.cpp\
$(top_srcdir)/src-std-framework/DOFMapping/SundanceInhomogeneousNodalDOFMap.cpp \
$(top_srcdir)/src-std-framework/DOFMapping/SundanceMapStructure.cpp\
$(top_srcdir)/src-std-framework/DOFMapping/SundanceMixedDOFMap.cpp \
$(top_srcdir)/src-std-framework/DOFMapping/SundanceNodalDOFMap.cpp \
$(top_srcdir)/src-std-framework/DOFMapping/SundancePartialElementDOFMap.cpp \
$(top_srcdir)/src-std-framework/DOFMapping/SundanceSpatiallyHomogeneousDOFMapBase.cpp 

SUNDANCE_ASSEMBLY_H = \
$(top_srcdir)/src-std-framework/Assembly/SundanceAssembler.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceAssemblyKernelBase.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceElementIntegral.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceFunctionalAssemblyKernel.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceFunctionalGradientAssemblyKernel.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceFunctionalEvaluator.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceGrouperBase.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceIntegralGroup.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceLocalDOFMap.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceLocalMatrixContainer.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceMapBundle.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceMatrixVectorAssemblyKernel.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceQuadratureEvalMediator.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceQuadratureIntegral.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceRefIntegral.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceStdFwkEvalMediator.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceTrivialGrouper.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceVectorAssemblyKernel.hpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceVectorFillingAssemblyKernel.hpp 

SUNDANCE_ASSEMBLY = \
$(top_srcdir)/src-std-framework/Assembly/SundanceAssembler.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceElementIntegral.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceFunctionalAssemblyKernel.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceFunctionalEvaluator.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceGrouperBase.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceIntegralGroup.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceLocalDOFMap.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceLocalMatrixContainer.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceMapBundle.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceMatrixVectorAssemblyKernel.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceQuadratureEvalMediator.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceQuadratureIntegral.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceRefIntegral.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceStdFwkEvalMediator.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceTrivialGrouper.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceVectorAssemblyKernel.cpp \
$(top_srcdir)/src-std-framework/Assembly/SundanceVectorFillingAssemblyKernel.cpp 

SUNDANCE_PROBLEM_H = \
$(top_srcdir)/src-std-framework/Problem/Sundance.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceAToCDensitySampler.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceAToCPointLocator.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceCToAInterpolator.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceEigensolution.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceExprFieldWrapper.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceFunctional.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceL2Projector.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceLinearEigenproblem.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceLinearProblem.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceModelEvaluatorBase.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceNLPModelEvaluator.hpp \
$(top_srcdir)/src-std-framework/Problem/SundanceNonlinearProblem.hpp 


SUNDANCE_PROBLEM = \
$(top_srcdir)/src-std-framework/Problem/Sundance.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceAToCDensitySampler.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceAToCPointLocator.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceCToAInterpolator.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceExprFieldWrapper.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceFunctional.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceL2Projector.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceLinearEigenproblem.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceLinearProblem.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceModelEvaluatorBase.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceNLPModelEvaluator.cpp \
$(top_srcdir)/src-std-framework/Problem/SundanceNonlinearProblem.cpp 



STD_FWK_SOURCES = \
$(SUNDANCE_ELEMENTS) \
$(SUNDANCE_DOFMAPPING) \
$(SUNDANCE_ASSEMBLY)  \
$(SUNDANCE_PROBLEM)  

STD_FWK_SOURCES_H = \
$(SUNDANCE_ELEMENTS_H) \
$(SUNDANCE_DOFMAPPING_H)  \
$(SUNDANCE_ASSEMBLY_H) \
$(SUNDANCE_PROBLEM_H)  


include $(top_srcdir)/Makefile.common

#
#np#  Files which are included in select builds are typically have
#np#  names that start with EXTRA
#
#np#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#np#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif

AM_CPPFLAGS = $(STD_FWK_INCLUDES) $(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Sundance library specifications
# ------------------------------------------------------------------------

STD_FWK_LIB = libsundanceStdFwk.a


#np# replace new_package with the name of the package being autotool'ed here
libsundanceStdFwk_a_SOURCES = $(STD_FWK_SOURCES)

#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libsundanceStdFwk_a_SOURCES = \
$(EXTRA_MPI) \
$(EXTRA_MPI_H) 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libsundanceStdFwk_a_AR = $(ALTERNATE_AR)
else

libsundanceStdFwk_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libsundanceStdFwk_a_LIBADD = $(SUNDANCE_TEMPLATE_OBJS) $(XTRALDADD) 

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(STD_FWK_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(STD_FWK_SOURCES_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SundanceConfig.h was installed when the utiities were built, so skip it.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 








