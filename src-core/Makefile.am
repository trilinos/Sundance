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


SUNDANCE_USEREXPRS_H = \
$(top_srcdir)/src-core/UserExprs/SundanceCoordExpr.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceDerivative.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceEssentialBC.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceExpr.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceIntegral.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceParameter.hpp \
$(top_srcdir)/src-core/UserExprs/SundanceStdMathOps.hpp 


SUNDANCE_USEREXPRS = \
$(top_srcdir)/src-core/UserExprs/SundanceCoordExpr.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceDerivative.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceEssentialBC.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceExpr.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceIntegral.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceParameter.cpp \
$(top_srcdir)/src-core/UserExprs/SundanceStdMathOps.cpp 


SUNDANCE_INTERNALEXPRS_H = \
$(top_srcdir)/src-core/InternalExprs/SundanceBinaryExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceConstantExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceCoordDeriv.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDiffOp.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDiscreteFuncElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceEvaluatableExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExprBase.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExprWithChildren.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceFuncElementBase.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceLeafExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceListExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceNonlinearUnaryOp.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceProductExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceScalarExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSparsityPattern.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSpatiallyConstantExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumOfBCs.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumOfIntegrals.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbPreprocessor.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbolicFuncElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbolicFunc.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceTestFuncElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryExpr.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryFunctor.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryMinus.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnknownFuncElement.hpp \
$(top_srcdir)/src-core/InternalExprs/SundanceZeroExpr.hpp 

SUNDANCE_INTERNALEXPRS = \
$(top_srcdir)/src-core/InternalExprs/SundanceBinaryExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceConstantExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceCoordDeriv.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDiffOp.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceDiscreteFuncElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceEvaluatableExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExprBase.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceExprWithChildren.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceFuncElementBase.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceLeafExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceListExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceNonlinearUnaryOp.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceProductExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceScalarExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSparsityPattern.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSpatiallyConstantExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumOfBCs.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSumOfIntegrals.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbPreprocessor.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbolicFuncElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceSymbolicFunc.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceTestFuncElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryExpr.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryFunctor.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnaryMinus.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceUnknownFuncElement.cpp \
$(top_srcdir)/src-core/InternalExprs/SundanceZeroExpr.cpp 

SUNDANCE_EVALUATION_H = \
$(top_srcdir)/src-core/Evaluation/SundanceBruteForceEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvalManager.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvaluatorFactory.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceInstructionCachingEvaluator.hpp \
$(top_srcdir)/src-core/Evaluation/SundanceTempStack.hpp 

SUNDANCE_EVALUATION = \
$(top_srcdir)/src-core/Evaluation/SundanceBruteForceEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvalManager.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceEvaluatorFactory.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceInstructionCachingEvaluator.cpp \
$(top_srcdir)/src-core/Evaluation/SundanceTempStack.cpp 


SUNDANCE_FRAMEWORKINTERFACE_H = \
$(top_srcdir)/src-core/FrameworkInterface/SundanceAbstractEvalMediator.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceCellFilterStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDerivBase.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDeriv.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDerivSet.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDiscreteFunctionStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEquationSet.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalContextStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalVector.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalVectorArray.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceFunctionalDeriv.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceLoadableVector.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceMultiIndex.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceMultipleDeriv.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceNullCellFilterStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceQuadratureFamilyStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceRegionQuadCombo.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceTestFunctionStub.hpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceUnknownFunctionStub.hpp 

SUNDANCE_FRAMEWORKINTERFACE = \
$(top_srcdir)/src-core/FrameworkInterface/SundanceAbstractEvalMediator.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceCellFilterStub.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDerivBase.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDeriv.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDerivSet.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceDiscreteFunctionStub.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEquationSet.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalContextStub.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalVector.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceEvalVectorArray.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceFunctionalDeriv.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceLoadableVector.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceMultiIndex.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceMultipleDeriv.cpp \
$(top_srcdir)/src-core/FrameworkInterface/SundanceNullCellFilterStub.cpp \
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





CORE_SOURCES = \
$(SUNDANCE_USEREXPRS) \
$(SUNDANCE_INTERNALEXPRS) \
$(SUNDANCE_EVALUATION) \
$(SUNDANCE_FRAMEWORKINTERFACE) \
$(SUNDANCE_TRANSFORMATIONS) 

CORE_SOURCES_H = \
$(SUNDANCE_USEREXPRS_H) \
$(SUNDANCE_INTERNALEXPRS_H) \
$(SUNDANCE_EVALUATION_H) \
$(SUNDANCE_FRAMEWORKINTERFACE_H) \
$(SUNDANCE_TRANSFORMATIONS_H)

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

AM_CPPFLAGS = $(CORE_INCLUDES)

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






