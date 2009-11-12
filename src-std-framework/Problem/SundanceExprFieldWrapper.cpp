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
    //map_(),
    indices_(),
    Expr_size_(1),
    isPointData_(true)
{
	int index = 0;
	Expr_size_ = expr.size();
	// Now it is independent of the size of the size of the Expression
	for(index = 0 ; index < Expr_size_ ; index++)
	{
	  const DiscreteFunction* df
	        = dynamic_cast<const DiscreteFunction*>(expr[index].ptr().get());
	  const DiscreteFuncElement* dfe
	        = dynamic_cast<const DiscreteFuncElement*>(expr[index].ptr().get());
      if (df != 0)
        {
          discreteSpace_ = df->discreteSpace();
          //map_ = df->map();
          indices_.append(tuple(0));
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
          //map_ = f->map();
          indices_.append(tuple(dfe->myIndex()));
          BasisFamily basis = discreteSpace_.basis()[indices_[index][0]];
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
}


double ExprFieldWrapper::getData(int cellDim, int cellID, int elem) const
{
  Array<int> dofs;

  discreteSpace_.map()->getDOFsForCell(cellDim, cellID, indices_[elem][0] , dofs); //indecies[elem][0] should be OK!

  //cout << "Arguments ExprFieldWrapper::getData " << cellDim << "," << cellID << "," << elem << std::endl;

  TEST_FOR_EXCEPTION(dofs.size() > 1, RuntimeError,
                     "too many DOFs found in ExprFieldWrapper::getData()");

  return df_->ghostView()->getElement(dofs[0]);
}


bool ExprFieldWrapper::isDefined(int cellDim, int cellID, int elem) const
{
  // this works only for the first
  RefCountPtr<const Set<int> > allowedFuncs 
    = discreteSpace_.map()->allowedFuncsOnCellBatch(cellDim, tuple(cellID));

  //cout << "Arguments ExprFieldWrapper::isDefined" << cellDim << "," << cellID << "," << elem << std::endl;

  return allowedFuncs->contains(indices_[elem][0]);
}
