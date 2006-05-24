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

#include "SundanceSumExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"



using namespace SundanceCore::Internal;
using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;


SumExpr::SumExpr(const RefCountPtr<ScalarExpr>& left,
                 const RefCountPtr<ScalarExpr>& right, int sign)
	: BinaryExpr(left, right, sign)
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "forming SumExpr " << toString());
  typedef Set<int>::const_iterator setIter;

  if (isEvaluatable(left.get()) && isEvaluatable(right.get()))
    {

      for (int d=0; d<MultiIndex::maxDim(); d++) 
        {
          int lod = leftEvaluatable()->orderOfSpatialDependency(d);
          int rod = rightEvaluatable()->orderOfSpatialDependency(d);
          if (lod < 0 || rod < 0) setOrderOfDependency(d, -1);
          else setOrderOfDependency(d, max(lod, rod));
        }

      Set<MultiSet<int> > tmp = leftEvaluatable()->funcIDSet();
      tmp.merge(rightEvaluatable()->funcIDSet());
      setFuncIDSet(tmp);
      Tabs tab1;
      SUNDANCE_VERB_HIGH(tab1 << "dependencies are " << tmp);
    }
}

bool SumExpr::isHungryDiffOp() const
{
  return leftScalar()->isHungryDiffOp() || rightScalar()->isHungryDiffOp();
}


const string& SumExpr::xmlTag() const 
{
	static string plusStr = "Plus";
	static string minusStr = "Minus";
	if (sign() < 0) return minusStr;
	return plusStr;
}

const string& SumExpr::opChar() const 
{
	static string plusStr = "+";
	static string minusStr = "-";
	if (sign() < 0) return minusStr;
	return plusStr;
}


bool SumExpr::allTermsHaveTestFunctions() const
{
  return leftEvaluatable()->allTermsHaveTestFunctions()
    && rightEvaluatable()->allTermsHaveTestFunctions();
}

