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

#include "SundanceXMLFunctionBuilder.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceXML;


XMLFunctionBuilder
::XMLFunctionBuilder(const RefCountPtr<XMLBasisBuilder>& basis,
                     const RefCountPtr<XMLDiscreteSpaceBuilder>& space)
  : XMLObjectBuilder<Expr>(),
    basis_(basis),
    space_(space)
{;}





Expr XMLFunctionBuilder::create(const XMLObject& xml) const 
{
  static Set<string> valid = makeSet("UnknownFunction",
                                     "TestFunction",
                                     "DiscreteFunction");

  checkOptionValidity(xml.getTag(), valid);

  string name = xml.getRequired("name");

  Expr rtn;

  if (xml.getTag()=="DiscreteFunction")
    {
      string spaceName = xml.getRequired("space");
      double value = xml.getRequiredDouble("value");
      DiscreteSpace space = space_->get(spaceName);
      rtn = new DiscreteFunction(space, value, name);
    }
  else
    {
      string basisName = xml.getRequired("basis");
      BasisFamily basis = basis_->get(basisName);
      if (xml.getTag()=="UnknownFunction")
        {
          rtn = new UnknownFunction(basis, name);
        }
      else
        {
          rtn = new TestFunction(basis, name);
        }
    }
  

  addToMap(name, rtn);

  return rtn;
}


