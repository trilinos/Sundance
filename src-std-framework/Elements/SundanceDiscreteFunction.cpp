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
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "PlayaDefaultBlockVectorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaDefaultBlockVectorImpl.hpp"
#endif

namespace Sundance
{
using namespace Teuchos;
using std::runtime_error;

static Time& getLocalValsTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("DF getLocalValues"); 
  return *rtn;
}
static Time& dfCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("DF ctor"); 
  return *rtn;
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const std::string& name)
  : DiscreteFunctionStub(tuple(name), space.dimStructure(),
    getRCP(new DiscreteFunctionData(space))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const Array<string>& name)
  : DiscreteFunctionStub(name, space.dimStructure(),
    getRCP(new DiscreteFunctionData(space))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const double& constantValue,
  const std::string& name)
  : DiscreteFunctionStub(tuple(name), space.dimStructure(),
    getRCP(new DiscreteFunctionData(space, constantValue))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
  Vector<double> vec = data_->getVector();
  vec.setToConstant(constantValue);
  data_->setVector(vec);
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const double& constantValue,
  const Array<string>& name)
  : DiscreteFunctionStub(name, space.dimStructure(),
    getRCP(new DiscreteFunctionData(space, constantValue))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
  Vector<double> vec = data_->getVector();
  vec.setToConstant(constantValue);
  data_->setVector(vec);
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const Vector<double>& vec,
  const std::string& name)
  : DiscreteFunctionStub(tuple(name), space.dimStructure(),
    getRCP(new DiscreteFunctionData(space, vec))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const Vector<double>& vec,
  const Array<string>& name)
  : DiscreteFunctionStub(name, space.dimStructure(),
    getRCP(new DiscreteFunctionData(space, vec))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

void DiscreteFunction::setVector(const Vector<double>& vec) 
{
  data_->setVector(vec);
}

void DiscreteFunction::updateGhosts() const
{
  data_->updateGhosts();
}


RCP<const MapStructure> DiscreteFunction::getLocalValues(int cellDim, 
  const Array<int>& cellLID,
  Array<Array<double> >& localValues) const 
{
  TimeMonitor timer(getLocalValsTimer());
  return data_->getLocalValues(cellDim, cellLID, localValues);
}


const DiscreteFunction* DiscreteFunction::discFunc(const Expr& expr)
{
  const ExprBase* e = expr.ptr().get();
  const DiscreteFunction* df 
    = dynamic_cast<const DiscreteFunction*>(e);

  TEST_FOR_EXCEPTION(df==0, std::runtime_error,
    "failed to cast " << expr << " to a discrete function. "
    "It appears to be of type " << e->typeName());

  return df;
}



DiscreteFunction* DiscreteFunction::discFunc(Expr& expr)
{
  DiscreteFunction* df 
    = dynamic_cast<DiscreteFunction*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(df==0, std::runtime_error,
    "failed to cast " << expr << " to a discrete function. "
    "It appears to be of type " << expr.ptr()->typeName());

  return df;
}


RCP<DiscreteFuncDataStub> DiscreteFunction::getRCP(DiscreteFunctionData* ptr)
{
  return rcp_dynamic_cast<DiscreteFuncDataStub>(rcp(ptr));
}



void updateDiscreteFunction(const Expr& newVals, Expr old)
{
  const DiscreteFunction* in = DiscreteFunction::discFunc(newVals);
  TEST_FOR_EXCEPTION(in==0, std::runtime_error,
    "input argument " << newVals 
    << " is not a discrete function in updateDiscreteFunction()");

  DiscreteFunction* out = DiscreteFunction::discFunc(old);
  TEST_FOR_EXCEPTION(out==0, std::runtime_error,
    "output argument " << old 
    << " is not a discrete function in updateDiscreteFunction()");

  TEST_FOR_EXCEPTION(
    in->getVector().space() != out->getVector().space(),
    std::runtime_error,
    "incompatible spaces " << in->getVector().space()
    << " and " << out->getVector().space()
    << " in updateDiscreteFunction()");

  Vector<double> vec = in->getVector();
  out->setVector(vec);
}

Expr copyDiscreteFunction(const Expr& u0, const string& name)
{
  const DiscreteFunction* in = DiscreteFunction::discFunc(u0);  

  TEST_FOR_EXCEPTION(in==0, std::runtime_error,
    "input argument " << u0 
    << " is not a discrete function in copyDiscreteFunction()");

  const DiscreteSpace& space = in->discreteSpace();
  const Vector<double>& vec = in->getVector();
  Vector<double> vecCopy = vec.copy();
  return new DiscreteFunction(space, vec, name);
}

void addVecToDiscreteFunction(Expr u, const Vector<double>& v)
{
  DiscreteFunction* in = DiscreteFunction::discFunc(u);

  TEST_FOR_EXCEPTION(in==0, std::runtime_error,
    "input argument " << u
    << " is not a discrete function in addVecToDiscreteFunction()");
  
  Vector<double> vec = in->getVector();

  TEST_FOR_EXCEPTION(
    in->getVector().space() != v.space(),
    std::runtime_error,
    "incompatible spaces " << in->getVector().space()
    << " and " << v.space()
    << " in addVecToDiscreteFunction()");
  
  vec.update(1.0, v);
}

Vector<double> getDiscreteFunctionVector(const Expr& u)
{
  const DiscreteFunction* df = DiscreteFunction::discFunc(u);  
  if (df != 0)
  {
    return df->getVector();
  }
  else
  {
    TEST_FOR_EXCEPTION(df==0 && u.size()==1, runtime_error,
      "non-block vector should be a discrete function in getDiscreteFunctionVector()");
    Array<Vector<double> > vec(u.size());
    for (int b=0; b<u.size(); b++)
    {
      vec[b] = getDiscreteFunctionVector(u[b]);
    }
    return blockVector(vec);
  }
}


}

