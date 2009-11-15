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

#include "SundanceDiscreteFunctionData.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceSubtypeEvaluator.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#endif


using namespace SundanceStdMesh;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace Teuchos;




DiscreteFunctionData::DiscreteFunctionData(const DiscreteSpace& space)
  : DiscreteFuncDataStub(), 
    space_(space),
    vector_(space_.createVector()),
    ghostView_(),
    ghostsAreValid_(false)
{}

DiscreteFunctionData::DiscreteFunctionData(const DiscreteSpace& space, 
  const double& constantValue)
  : DiscreteFuncDataStub(), 
    space_(space),
    vector_(space_.createVector()),
    ghostView_(),
    ghostsAreValid_(false)
{
  vector_.setToConstant(constantValue);
}

DiscreteFunctionData::DiscreteFunctionData(const DiscreteSpace& space, 
  const Vector<double>& vector)
  : DiscreteFuncDataStub(), 
    space_(space),
    vector_(vector),
    ghostView_(),
    ghostsAreValid_(false)
{}

const DiscreteFunctionData* DiscreteFunctionData::getData(const DiscreteFuncElement* dfe)
{
  TEST_FOR_EXCEPTION(dfe==0, RuntimeError, "null argument to DiscreteFunctionData::getData()");
  RCP<const DiscreteFunctionData> rtn 
    = rcp_dynamic_cast<const DiscreteFunctionData>(dfe->commonData());
  TEST_FOR_EXCEPTION(rtn.get()==0, RuntimeError, 
    "cast to DiscreteFunctionData* failed for "
    "discrete function element " << dfe->toXML());
  return rtn.get();
}

void DiscreteFunctionData::setVector(const Vector<double>& vec) 
{
  ghostsAreValid_ = false;
  vector_ = vec;
}

void DiscreteFunctionData::updateGhosts() const
{
  if (!ghostsAreValid_)
  {
    space_.importGhosts(vector_, ghostView_);
    ghostsAreValid_ = true;
  }
}


RefCountPtr<const MapStructure> DiscreteFunctionData
::getLocalValues(int cellDim, 
  const Array<int>& cellLID,
  Array<Array<double> >& localValues) const 
{
  Tabs tab;

  if (Evaluator::classVerbosity() > 3)
  {
    cerr << tab << "getting DF local values" << endl;
  }
  updateGhosts();

  const RefCountPtr<DOFMapBase>& map = space_.map();
  static Array<Array<int> > dofs;
  Array<int> nNodes;

  RefCountPtr<const Set<int> > requestedFuncs = map->allowedFuncsOnCellBatch(cellDim,
    cellLID);

  RefCountPtr<const MapStructure> s = map->getDOFsForCellBatch(cellDim, cellLID,
    *requestedFuncs,
    dofs, nNodes, Evaluator::classVerbosity());
  localValues.resize(s->numBasisChunks());
  for (unsigned int b=0; b<nNodes.size(); b++)
  {
    int nFuncs = s->numFuncs(b);
    localValues[b].resize(nFuncs*nNodes[b]*cellLID.size());
    ghostView_->getElements(&(dofs[b][0]), dofs[b].size(), localValues[b]);
  }

  return s;
}



