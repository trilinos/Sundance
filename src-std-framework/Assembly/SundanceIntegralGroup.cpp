/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceIntegralGroup.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

static Time& integrationTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("integration"); 
  return *rtn;
}


IntegralGroup
::IntegralGroup(const Array<RefCountPtr<ElementIntegral> >& integrals,
                const Array<int>& resultIndices)
  : order_(0),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(),
    unkID_(),
    integrals_(integrals),
    resultIndices_(resultIndices)
{
  verbosity() = classVerbosity();
}

IntegralGroup
::IntegralGroup(const Array<int>& testID,
                const Array<RefCountPtr<ElementIntegral> >& integrals,
                const Array<int>& resultIndices)
  : order_(1),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(testID),
    unkID_(),
    integrals_(integrals),
    resultIndices_(resultIndices)
{
  verbosity() = classVerbosity();

  for (int i=0; i<integrals.size(); i++)
    {
      int nt = integrals[i]->nNodesTest();
      if (i > 0) 
        {
          TEST_FOR_EXCEPTION(nt != nTestNodes_, InternalError,
                             "IntegralGroup ctor detected integrals with inconsistent numbers of test nodes");
        }
      nTestNodes_ = nt;
    }
}

IntegralGroup
::IntegralGroup(const Array<int>& testID,
                const Array<int>& unkID,
                const Array<RefCountPtr<ElementIntegral> >& integrals,
                const Array<int>& resultIndices)
  : order_(2),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(testID),
    unkID_(unkID),
    integrals_(integrals),
    resultIndices_(resultIndices)
{
  verbosity() = classVerbosity();

  for (int i=0; i<integrals.size(); i++)
    {
      int nt = integrals[i]->nNodesTest();
      int nu = integrals[i]->nNodesUnk();
      if (i > 0) 
        {
          TEST_FOR_EXCEPTION(nt != nTestNodes_, InternalError,
                             "IntegralGroup ctor detected integrals with inconsistent numbers of test nodes");
          TEST_FOR_EXCEPTION(nu != nUnkNodes_, InternalError,
                             "IntegralGroup ctor detected integrals with inconsistent numbers of unk nodes");
        }
      nTestNodes_ = nt;
      nUnkNodes_ = nu;
    }
}

bool IntegralGroup
::evaluate(const CellJacobianBatch& J,
           const Array<RefCountPtr<EvalVector> >& vectorCoeffs,
           const Array<double>& constantCoeffs,
           RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(integrationTimer());
  Tabs tab0;


  SUNDANCE_OUT(verbosity() > VerbSilent,
               tab0 << "evaluating an integral group of size "
               << integrals_.size());

  /* initialize the return vector */
  if (integrals_[0]->nNodes() == -1) A->resize(1);
  else A->resize(J.numCells() * integrals_[0]->nNodes());
  double* aPtr = &((*A)[0]);
  int n = A->size();
  for (int i=0; i<n; i++) aPtr[i] = 0.0;

  SUNDANCE_OUT(verbosity() > VerbHigh, tab0 << "begin A=" << *A);

  /* do the integrals */
  for (int i=0; i<integrals_.size(); i++)
    {
      Tabs tab;
      const RefIntegral* ref 
        = dynamic_cast<const RefIntegral*>(integrals_[i].get());
      const QuadratureIntegral* quad 
        = dynamic_cast<const QuadratureIntegral*>(integrals_[i].get());

      if (ref!=0)
        {
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       tab << "Integrating term group " << i 
                       << " by transformed reference integral");
          double f = constantCoeffs[resultIndices_[i]];
          SUNDANCE_OUT(verbosity() > VerbSilent,
                       tab << "Coefficient is " << f);
          ref->transform(J, f, A);
        }
      else 
        {
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       tab << "Integrating term group " << i 
                       << " by quadrature");
          TEST_FOR_EXCEPTION(vectorCoeffs[resultIndices_[i]]->length()==0,
                             InternalError,
                             "zero-length coeff vector detected in "
                             "quadrature integration branch of "
                             "IntegralGroup::evaluate(). String value is ["
                             << vectorCoeffs[resultIndices_[i]]->str()
                             << "]");
          const double* const f = vectorCoeffs[resultIndices_[i]]->start();
          quad->transform(J, f, A);
        }
      SUNDANCE_OUT(verbosity() > VerbHigh, tab << "i=" << i << " A=" << *A);
    }
  SUNDANCE_OUT(verbosity() > VerbSilent, tab0 << "done");

  return true;
}


