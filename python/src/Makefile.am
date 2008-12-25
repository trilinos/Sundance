# @HEADER
# ************************************************************************
# 
#              PyTrilinos.Epetra: Python Interface to Epetra
#                   Copyright (2005) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER


# Include swig rules
include $(top_srcdir)/config/swig.am
@HAVE_MPI_TRUE@MPI_FLAG = -DHAVE_MPI

SWIG_CXXFLAGS = -I$(top_srcdir)/src-std-mesh \
                 -I$(srcdir)                       \
                -I$(top_builddir)/src-std-mesh  \
                 $(MPI_FLAG)

EXTRA_DIST = PySundance.i setup.py \
PyTeuchos_Utils.cpp PyTeuchos_Utils.hpp \
PySundanceNOXSolverHandle.cpp PySundanceNOXSolverHandle.hpp \
PySundanceCallback.cpp PySundanceCallback.hpp \
PySundanceCellPredicate.cpp PySundanceCellPredicate.hpp \
PySundanceLinearSolver.cpp PySundanceLinearSolver.hpp \
PySundanceFIATScalarAdapter.cpp PySundanceFIATScalarAdapter.hpp \
PySundanceFIATQuadratureAdapter.cpp PySundanceFIATQuadratureAdapter.hpp

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','temp.%s-%s'%(get_platform(),sys.version[0:3]))")

# path to shared library
solibdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','lib.%s-%s'%(get_platform(),sys.version[0:3]),'PySundance')")

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Define the SWIG interface files and the resulting wrapper files
INTERFACES = PySundance.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cpp)
PROXIES    = $(INTERFACES:.i=.py      )

LIBDEPENDENCIES = $(top_builddir)/src-std-framework/libsundanceStdFwk.a \
$(top_builddir)/src-std-mesh/libsundanceStdMesh.a \
$(top_builddir)/src-core/libsundanceCore.a \
$(top_builddir)/src-interop/libsundanceInterop.a \
$(top_builddir)/src-utils/libsundanceUtils.a $(top_builddir)/src-solvers/libTSFExtended.a

include_HEADERS = AToC.i \
Basis.i \
CellFilter.i\
CoordinateSystem.i\
Discrete.i\
Integral.i\
LinearProblem.i\
Mesh.i\
NonlinearProblem.i\
ParameterList.i\
PySundance.i\
Quadrature.i\
Spectral.i\
Symbolics.i\
TSF.i\
Utils.i\
Viz.i \
PyTeuchos_Utils.hpp \
PySundanceNOXSolverHandle.hpp \
PySundanceCallback.hpp \
PySundanceCellPredicate.hpp \
PySundanceLinearSolver.hpp \
PySundanceFIATScalarAdapter.hpp \
PySundanceFIATQuadratureAdapter.hpp



build: $(WRAPPERS) $(tempdir) $(LIBDEPENDENCIES)
	rm -f $(solibdir)/_PySundance.so ../../python/python/src/PySundance_wrap.o
	$(PYTHON) $(srcdir)/setup.py build

install: build install-includeHEADERS
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(PYTHON_PREFIX)

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean
	$(RM) -r ../../packages

CLEANFILES += $(WRAPPERS) $(PROXIES)
