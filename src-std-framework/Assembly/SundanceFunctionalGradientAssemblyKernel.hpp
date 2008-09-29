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

#ifndef SUNDANCE_FUNCTIONALGRADIENTASSEMBLYKERNEL_H
#define SUNDANCE_FUNCTIONALGRADIENTASSEMBLYKERNEL_H

#include "SundanceDefs.hpp"
#include "SundanceVectorFillingAssemblyKernel.hpp"
#include "SundanceFunctionalAssemblyKernel.hpp"

namespace SundanceStdFwk
{
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;

namespace Internal
{
using namespace Teuchos;

/** 
 * FunctionalGradientAssemblyKernel does assembly of a functional and
 * its gradient. 
 */
class FunctionalGradientAssemblyKernel : public AssemblyKernelBase
{
public:
  /** */
  FunctionalGradientAssemblyKernel(const MPIComm& comm,
    const Array<RefCountPtr<DOFMapBase> >& dofMap,
    const Array<RefCountPtr<Array<int> > >& isBCIndex,
    const Array<int>& lowestLocalIndex,
    Vector<double>& grad,
    bool partitionBCs,
    double* value, 
    int verb)
    : AssemblyKernelBase(verb),
      funcKernel_(rcp(new FunctionalAssemblyKernel(comm, value, verb))),
      vecKernel_(rcp(new VectorAssemblyKernel(dofMap, isBCIndex,
            lowestLocalIndex, grad, partitionBCs, verb)))
    {}

  /** */
  void prepareForWorkSet(
    const Array<Set<int> >& requiredTests,
    const Array<Set<int> >& requiredUnks,
    RefCountPtr<StdFwkEvalMediator> mediator) 
    {
      funcKernel_->prepareForWorkSet(requiredTests, requiredUnks, mediator);
      vecKernel_->prepareForWorkSet(requiredTests, requiredUnks, mediator);
    }

  /** */
  void fill(bool isBC,
    const IntegralGroup& group,
    const RefCountPtr<Array<double> >& localValues) 
    {
      if (group.isOneForm())
      {
        vecKernel_->fill(isBC, group, localValues);
      }
      else if (group.isZeroForm())
      {
        funcKernel_->fill(isBC, group, localValues);
      }
      else
      {
        TEST_FOR_EXCEPT(group.isTwoForm());
      }
    }

  /** */
  void postLoopFinalization()
    {
      funcKernel_->postLoopFinalization();
      vecKernel_->postLoopFinalization();
    }

private:
  RefCountPtr<FunctionalAssemblyKernel> funcKernel_;
  RefCountPtr<VectorAssemblyKernel> vecKernel_;
};

}
}



#endif
