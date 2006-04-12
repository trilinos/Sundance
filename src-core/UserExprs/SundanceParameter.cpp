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

#include "SundanceParameter.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

Parameter::Parameter(const double& value, const string& name)
	: FuncElementBase(name, ""),
    DiscreteFuncElement(rcp(new ParameterData(value)), name, "", 0), 
    SpatiallyConstantExpr()
{;}

void Parameter::setValue(const double& value)
{
  
  data()->setValue(value);
}

const double& Parameter::value() const {return data()->value();}


void Parameter::findNonzeros(const EvalContext& context,
                             const Set<MultiIndex>& multiIndices,
                             const Set<MultiSet<int> >& activeFuncIDs,
                             bool regardFuncsAsConstant) const 
{
  SpatiallyConstantExpr::findNonzeros(context, multiIndices, activeFuncIDs,
                                      regardFuncsAsConstant);
}


Set<MultipleDeriv> Parameter::internalFindW(int order, const EvalContext& context) const 
{
  return SpatiallyConstantExpr::internalFindW(order, context);
}

Evaluator* Parameter::createEvaluator(const EvaluatableExpr* expr,
                                      const EvalContext& context) const 
{
  return SpatiallyConstantExpr::createEvaluator(expr, context);
}


XMLObject Parameter::toXML() const 
{
	XMLObject rtn("Parameter");
	rtn.addAttribute("name", name());
	rtn.addAttribute("value", Teuchos::toString(value()));
	return rtn;
}


const ParameterData* Parameter::data() const
{
  const ParameterData* pd = dynamic_cast<const ParameterData*>(commonData());

  TEST_FOR_EXCEPTION(pd==0, InternalError, "cast failed");

  return pd;
}



ParameterData* Parameter::data()
{
  ParameterData* pd = dynamic_cast<ParameterData*>(commonData());

  TEST_FOR_EXCEPTION(pd==0, InternalError, "cast failed");

  return pd;
}
