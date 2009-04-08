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


SUNDANCE_USEREXPRS_H = \
$(top_srcdir)/src-core/UserExprs/SundanceCellDiameterExpr.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceCellVectorExpr.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceComplexExpr.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceCoordExpr.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceDerivative.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceEssentialBC.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceExpr.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceFunctionalDerivative.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceIntegral.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceParameter.hpp \
$(top_srcdir)/src-core/UserExprs/SundancePointwiseUserDefFunctor.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceStdMathOps.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceUnknownParameter.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceUserDefFunctor.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceUserDefOp.hpp 


SUNDANCE_USEREXPRS = \
$(top_srcdir)/src-core/UserExprs/SundanceCellDiameterExpr.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceCellVectorExpr.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceComplexExpr.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceCoordExpr.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceDerivative.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceEssentialBC.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceExpr.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceFunctionalDerivative.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceIntegral.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceParameter.cpp \
$(top_srcdir)/src-core/UserExprs/SundancePointwiseUserDefFunctor.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceStdMathOps.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceUnknownParameter.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceUserDefFunctor.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceUserDefOp.cpp 


SUNDANCE_INTERNALEXPRS_H = \
$(top_srcdir)/src-core/InternalExprs/SundanceBinaryExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceConstantExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceCoordDeriv.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDerivOfSymbFunc.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDiffOp.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDiscreteFuncElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceEvaluatableExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExplicitFunctionalDerivativeElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExprBase.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExprWithChildren.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceFuncElementBase.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceFunctionalPolynomial.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceListExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceNonlinearExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceNonlinearUnaryOp.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceParameterData.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceProductExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceScalarExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSparsitySubset.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSparsitySuperset.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSpatiallyConstantExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumOfBCs.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumOfIntegrals.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbPreprocessor.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbolicFuncElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbolicFunc.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceTestFuncElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryMinus.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnknownFuncElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnknownParameterElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUserDefFunctorElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUserDefOpElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceZeroExpr.hpp 

SUNDANCE_INTERNALEXPRS = \
$(top_srcdir)/src-core/InternalExprs/SundanceBinaryExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceConstantExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceCoordDeriv.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDerivOfSymbFunc.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDiffOp.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDiscreteFuncElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceEvaluatableExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExplicitFunctionalDerivativeElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExprBase.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExprWithChildren.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceFuncElementBase.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceFunctionalPolynomial.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceListExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceNonlinearExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceNonlinearUnaryOp.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceProductExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceScalarExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSparsitySubset.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSparsitySuperset.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSpatiallyConstantExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumOfBCs.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumOfIntegrals.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbPreprocessor.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbolicFuncElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbolicFunc.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceTestFuncElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryMinus.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnknownFuncElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnknownParameterElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUserDefFunctorElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUserDefOpElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceZeroExpr.cpp 

SUNDANCE_EVALUATION_H = \
$(top_srcdir)/src-core/Evaluation/SundanceBinaryEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceCellDiameterEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceCellVectorEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceChainRuleEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceChainRuleSum.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceConstantEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceCoordExprEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceDerivOfSymbFuncEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceDiscreteFuncEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceDiffOpEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceEFDEEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvalManager.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvaluatorFactory.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceNonlinearUnaryOpEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceNullEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceProductEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceStringEvalMediator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceSubtypeEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceSumEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceSymbolicFuncEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceTempStack.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceUnaryEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceUnaryMinusEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceUserDefOpCommonEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceUserDefOpEvaluator.hpp 


SUNDANCE_EVALUATION = \
$(top_srcdir)/src-core/Evaluation/SundanceCellDiameterEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceCellVectorEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceChainRuleEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceChainRuleSum.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceConstantEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceCoordExprEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceDerivOfSymbFuncEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceDiscreteFuncEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceDiffOpEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceEFDEEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvalManager.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceNonlinearUnaryOpEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceNullEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceProductEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceStringEvalMediator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceSumEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceSymbolicFuncEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceTempStack.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceUnaryMinusEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceUserDefOpCommonEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceUserDefOpEvaluator.cpp 


SUNDANCE_FRAMEWORKINTERFACE_H = \
$(top_srcdir)/src-core/FrameworkInterface/SundanceAbstractEvalMediator.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceCommonFuncDataStub.hpp\
$(top_srcdir)/src-core/FrameworkInterface/SundanceDerivBase.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDeriv.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDerivSet.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDiscreteFuncDataStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDiscreteFunctionStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEquationSet.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalContext.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalVector.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceFunctionalDeriv.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceFunctionSupportResolver.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceLoadableVector.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceMultiIndex.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceMultipleDeriv.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceQuadratureFamilyStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceRegionQuadCombo.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceTestFuncDataStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceTestFunctionStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceUnknownFuncDataStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceUnknownFunctionStub.hpp 

SUNDANCE_FRAMEWORKINTERFACE = \
$(top_srcdir)/src-core/FrameworkInterface/SundanceAbstractEvalMediator.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDerivBase.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDeriv.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDerivSet.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDiscreteFunctionStub.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEquationSet.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalVector.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceFunctionalDeriv.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceFunctionSupportResolver.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceLoadableVector.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceMultiIndex.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceMultipleDeriv.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceQuadratureFamilyStub.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceRegionQuadCombo.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceTestFunctionStub.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceUnknownFunctionStub.cpp 


SUNDANCE_TRANSFORMATIONS_H = \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceSymbolicTransformation.hpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceSumTransformation.hpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceSumTransformationSequence.hpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceStdSumTransformations.hpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceProductTransformation.hpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceProductTransformationSequence.hpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceStdProductTransformations.hpp

SUNDANCE_TRANSFORMATIONS = \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceSymbolicTransformation.cpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceSumTransformation.cpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceSumTransformationSequence.cpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceStdSumTransformations.cpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceProductTransformation.cpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceProductTransformationSequence.cpp \
$(top_srcdir)/src-core/SymbolicTransformations/SundanceStdProductTransformations.cpp

SUNDANCE_SPECTRAL_H = \
$(top_srcdir)/src-core/Spectral/chaos.h \
$(top_srcdir)/src-core/Spectral/VECMAT.h \
$(top_srcdir)/src-core/Spectral/cijk.h \
$(top_srcdir)/src-core/Spectral/SundanceHermiteSpectralBasis.hpp \
$(top_srcdir)/src-core/Spectral/SundanceSpectralBasis.hpp \
$(top_srcdir)/src-core/Spectral/SundanceSpectralBasisBase.hpp \
$(top_srcdir)/src-core/Spectral/SundanceSpectralExpr.hpp 

SUNDANCE_SPECTRAL = \
$(top_srcdir)/src-core/Spectral/SundanceHermiteSpectralBasis.cpp \
$(top_srcdir)/src-core/Spectral/SundanceSpectralExpr.cpp 

SUNDANCE_TESTUTILITIES_H = \
$(top_srcdir)/src-core/TestingUtilities/SundanceADBasis.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceADCoord.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceADDerivative.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceADField.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceEvaluationTester.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestDiscreteFuncData.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestDiscreteFunction.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestEvalMediator.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestUnknownFuncData.hpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestUnknownFunction.hpp 

SUNDANCE_TESTUTILITIES = \
$(top_srcdir)/src-core/TestingUtilities/SundanceADBasis.cpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceADCoord.cpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceADDerivative.cpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceADField.cpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceEvaluationTester.cpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestDiscreteFunction.cpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestEvalMediator.cpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestUnknownFuncData.cpp \
$(top_srcdir)/src-core/TestingUtilities/SundanceTestUnknownFunction.cpp 




CORE_SOURCES = \
$(SUNDANCE_USEREXPRS) \
$(SUNDANCE_INTERNALEXPRS) \
$(SUNDANCE_EVALUATION) \
$(SUNDANCE_FRAMEWORKINTERFACE) \
$(SUNDANCE_TRANSFORMATIONS) \
$(SUNDANCE_SPECTRAL) \
$(SUNDANCE_TESTUTILITIES)

CORE_SOURCES_H = \
$(SUNDANCE_USEREXPRS_H) \
$(SUNDANCE_INTERNALEXPRS_H) \
$(SUNDANCE_EVALUATION_H) \
$(SUNDANCE_FRAMEWORKINTERFACE_H) \
$(SUNDANCE_TRANSFORMATIONS_H) \
$(SUNDANCE_SPECTRAL_H) \
$(SUNDANCE_TESTUTILITIES_H) 

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

AM_CPPFLAGS = $(CORE_INCLUDES) $(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Sundance library specifications
# ------------------------------------------------------------------------

CORE_LIB = libsundanceCore.a


#np# replace new_package with the name of the package being autotool'ed here
libsundanceCore_a_SOURCES = $(CORE_SOURCES)

#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libsundanceCore_a_SOURCES = \
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

libsundanceCore_a_AR = $(ALTERNATE_AR)
else

libsundanceCore_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libsundanceCore_a_LIBADD = $(SUNDANCE_TEMPLATE_OBJS) $(XTRALDADD) 

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(CORE_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(CORE_SOURCES_H) 

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






