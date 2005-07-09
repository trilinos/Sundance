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

#include "SundanceDiscreteFunction.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceSubtypeEvaluator.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;


static Time& getLocalValsTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("DF getLocalValues()"); 
  return *rtn;
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(space_.createVector()),
    ghostView_(),
    ghostsAreValid_(false)
{}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const double& constantValue,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(space_.createVector()),
    ghostView_(),
    ghostsAreValid_(false)
{
  vector_.setToConstant(constantValue);
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const Vector<double>& vector,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc()), 
    FuncWithBasis(space.basis()),
    space_(space),
    vector_(vector),
    ghostView_(),
    ghostsAreValid_(false)
{}

void DiscreteFunction::setVector(const Vector<double>& vec) 
{
  ghostsAreValid_ = false;
  vector_ = vec;
}

void DiscreteFunction::updateGhosts() const
{
  if (!ghostsAreValid_)
    {
      space_.importGhosts(vector_, ghostView_);
      ghostsAreValid_ = true;
    }
}


void DiscreteFunction::getLocalValues(int cellDim, 
                        const Array<int>& cellLID,
                        Array<double>& localValues) const 
{
  TimeMonitor timer(getLocalValsTimer());
  Tabs tab;

  if (Evaluator::classVerbosity() > VerbHigh)
    {
      cerr << tab << "getting DF local values" << endl;
    }
  updateGhosts();

  const RefCountPtr<DOFMapBase>& map = space_.map();
  static Array<int> dofs;
  static Array<int> indices;
  int nNodes;
  map->getDOFsForCellBatch(cellDim, cellLID, dofs, nNodes);
  int nFunc = space_.nFunc();
  int nCells = cellLID.size();
  localValues.resize(nFunc*cellLID.size()*nNodes);
  indices.resize(dofs.size());
  
  /* store the values in fortran order with node number running fastest,
   * then function ID, then cell ID */
  for (int c=0; c<cellLID.size(); c++)
    {
      for (int f=0; f<nFunc; f++)
        {
          for (int n=0; n<nNodes; n++)
            {
              indices[(c*nFunc + f)*nNodes + n]
                = dofs[(f*nCells + c)*nNodes + n];
            }
        }
    }
  ghostView_->getElements(&(indices[0]), indices.size(), localValues);

  if (Evaluator::classVerbosity() > VerbHigh)
    {
      cerr << tab << "local values are " << localValues << endl;
    }
}


const DiscreteFunction* DiscreteFunction::discFunc(const Expr& expr)
{
  const DiscreteFunction* df 
    = dynamic_cast<const DiscreteFunction*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(df==0, RuntimeError,
                     "failed to cast " << expr << " to a discrete function");

  return df;
}



DiscreteFunction* DiscreteFunction::discFunc(Expr& expr)
{
  DiscreteFunction* df 
    = dynamic_cast<DiscreteFunction*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(df==0, RuntimeError,
                     "failed to cast " << expr << " to a discrete function");

  return df;
}



