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

#include "SundanceXMLExprBuilder.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;


XMLExprBuilder::XMLExprBuilder()
  : Repository<ExprFamily>()
{;}



Expr XMLExprBuilder::assembleFromXML(const XMLObject& xml)
{
  const string& tag = xml.getTag();

  if (binaryOps.contains(tag))
    {
      checkNumChildren(xml, 2);
      const XMLObject& left = xml.getChild(0);
      const XMLObject& right = xml.getChild(1);
      Expr leftExpr = assembleFromXML(left);
      Expr rightExpr = assembleFromXML(right);
      if (tag=="Sum")
        {
          return leftExpr + rightExpr;
        }
      else if (tag=="Product")
        {
          return leftExpr * rightExpr; 
        }
    }
  else if (unaryOps.contains(tag))
    {
      checkNumChildren(xml, 1);
      const XMLObject& arg = xml.getChild(0);
      Expr argExpr = assembleFromXML(arg);
      if (tag=="UnaryMinus")
        {
          return -argExpr;
        }
    }
  else if (tag=="Var")
    {
      const string& name = xml.getRequired("name");
      if (contains(name))
        {
          return get(name);
        }
      if (functions_->contains(name))
        {
          return functions_->get(name);
        }
    }

}

