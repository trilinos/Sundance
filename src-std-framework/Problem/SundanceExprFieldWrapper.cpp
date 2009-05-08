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

#include "SundanceExprFieldWrapper.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


ExprFieldWrapper::ExprFieldWrapper(const Expr& expr)
  : expr_(expr),
    df_(),
    discreteSpace_(),
    map_(),
    indices_(),
    isPointData_(true)
{
  if (expr.size()==1)
    {
      const DiscreteFunction* df 
        = dynamic_cast<const DiscreteFunction*>(expr[0].ptr().get());
      const DiscreteFuncElement* dfe 
        = dynamic_cast<const DiscreteFuncElement*>(expr[0].ptr().get());
      if (df != 0)
        {
          discreteSpace_ = df->discreteSpace();
          map_ = df->map();
          indices_ = tuple(0);
          BasisFamily basis = discreteSpace_.basis()[0];
          const Lagrange* lagr = dynamic_cast<const Lagrange*>(basis.ptr().get());
          if (lagr != 0 && lagr->order()==0) isPointData_ = false;
          df_ = df->data();
        }
      else if (dfe != 0)
        {
          const DiscreteFunctionData* f = DiscreteFunctionData::getData(dfe);

          TEST_FOR_EXCEPTION(f == 0, RuntimeError,
                             "ExprFieldWrapper ctor argument "
                             << expr << " is not a discrete function");
          discreteSpace_ = f->discreteSpace();
          map_ = f->map();
          indices_ = tuple(dfe->myIndex());
          BasisFamily basis = discreteSpace_.basis()[indices_[0]];
          const Lagrange* lagr = dynamic_cast<const Lagrange*>(basis.ptr().get());
          if (lagr != 0 && lagr->order()==0) isPointData_ = false;
          df_ = f;
          
        }
      else
        {
          TEST_FOR_EXCEPTION(df == 0 && dfe == 0, RuntimeError,
                             "ExprFieldWrapper ctor argument is not a discrete "
                             "function");
        }
    }
  else
    {
      TEST_FOR_EXCEPTION(expr.size() != 1, RuntimeError,
                         "non-scalar expr given to ExprFieldWrapper ctor");
    }
}


double ExprFieldWrapper::getData(int cellDim, int cellID, int elem) const
{
  Array<int> dofs;

  map_->getDOFsForCell(cellDim, cellID, indices_[elem], dofs);

  TEST_FOR_EXCEPTION(dofs.size() > 1, RuntimeError,
                     "too many DOFs found in ExprFieldWrapper::getData()");

  TEST_FOR_EXCEPTION(dofs.size() == 0, RuntimeError,
                     "no DOFs found in ExprFieldWrapper::getData()");

  return df_->ghostView()->getElement(dofs[0]);
}

/*
void ExprFieldWrapper::getData(int cellDim, const Array<int>& cellLID, 
                               const double& undefinedVal,
                               std::vector<double>& vals) const 
{
  Array<Array<int> > dofs;
  Array<int> nNodes;
  Set<int> reqFuncs;
  for (int i=0; i<indices_.size(); i++) funcs.put(indices_[i]);

  RefCountPtr<const Set<int> > allowedFuncs 
    = map_->allowedFuncsOnCellBatch(cellDim, cellLID);

  if (reqFuncs.set

  

  vals.resize(cellLID.size() * indices_.size());

  RefCountPtr<const MapStructure> ms 
    = map_->getDOFsForCellBatch(cellDim, cellLID, funcs, dofs, nNodes);

  const int* globalIndices = &(dofs[0][0]);
  df_->ghostView()->getElements(globalIndices, 

  return df_->ghostView()->getElement(dofs[0]);
}
*/
    
bool ExprFieldWrapper::isDefined(int cellDim, int cellID, int elem) const
{
  RefCountPtr<const Set<int> > allowedFuncs 
    = map_->allowedFuncsOnCellBatch(cellDim, tuple(cellID));
  return allowedFuncs->contains(indices_[elem]);
}
