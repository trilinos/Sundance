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

#include "SundanceTrivialGrouper.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceMap.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

void TrivialGrouper::findGroups(const EquationSet& eqn,
  const CellType& maxCellType,
  int spatialDim,
  const CellType& cellType,
  int cellDim,
  const QuadratureFamily& quad,
  const RefCountPtr<SparsitySuperset>& sparsity,
  Array<IntegralGroup>& groups) const
{
  Tabs tab;
  const ParameterList& verbParams = params();

  SUNDANCE_LEVEL1("find groups",
    tab << "in TrivialGrouper::findGroups(), num derivs = " 
    << sparsity->numDerivs());
  SUNDANCE_LEVEL2("find groups",
    tab << "cell type = " << cellType);
  SUNDANCE_LEVEL2("find groups",
    tab << "sparsity = " << std::endl << *sparsity << std::endl);

  int vecCount=0;
  int constCount=0;

  /* turn off grouping for submaximal cells. This works around a bug detected by Rob Kirby that
   * shows up with Nitsche BCs in mixed-element discretizations */
  bool doGroups = true;
  if (cellType != maxCellType) doGroups = false;

  typedef SundanceUtils::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<RefCountPtr<ElementIntegral> > > twoFormMap;
  typedef SundanceUtils::Map<OrderedPair<int,BasisFamily>, Array<RefCountPtr<ElementIntegral> > > oneFormMap;
  SundanceUtils::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<RefCountPtr<ElementIntegral> > > twoForms;
  SundanceUtils::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<int> > twoFormResultIndices;
  SundanceUtils::Map<OrderedPair<int,BasisFamily>, Array<RefCountPtr<ElementIntegral> > > oneForms;
  SundanceUtils::Map<OrderedPair<int,BasisFamily>, Array<int> > oneFormResultIndices;

  for (int i=0; i<sparsity->numDerivs(); i++)
  {
    Tabs tab1;
    const MultipleDeriv& d = sparsity->deriv(i);
    SUNDANCE_LEVEL3("find groups",
      tab1 << "defining integration policy for " << d);
      
    if (d.order()==0) 
    {
      Tabs tab2;
      RefCountPtr<ElementIntegral> integral ;
      int resultIndex;
      if (sparsity->isConstant(i))
      {
        integral = rcp(new RefIntegral(spatialDim, maxCellType, 
            cellDim, cellType, verbParams));
        resultIndex = constCount++;
      }
      else
      {
        integral = rcp(new QuadratureIntegral(spatialDim, maxCellType, 
            cellDim, cellType, quad,
            verbParams));
        resultIndex = vecCount++;
      }
      SUNDANCE_LEVEL3("find groups", tab2 << "is zero-form");
      groups.append(IntegralGroup(tuple(integral),
          tuple(resultIndex), verbParams));
    }
    else
    {
      Tabs tab2;
      BasisFamily testBasis;
      BasisFamily unkBasis;
      MultiIndex miTest;
      MultiIndex miUnk;
      int rawTestID;
      int rawUnkID;
      int testID;
      int unkID;
      int testBlock;
      int unkBlock;
      bool isOneForm;
      extractWeakForm(eqn, d, testBasis, unkBasis, miTest, miUnk, 
        rawTestID, rawUnkID,
        testID, unkID,
        testBlock, unkBlock,
        isOneForm);

      /* In variational problems we might have (u,v) and (v,u). Because the derivative
       * is stored as an unordered multiset it can't distinguish between the two cases.
       * We need to check the equation set to see if the two functions show up as
       * variations and unknowns. If so, then we need to produce the transposed integral.
       */
      bool transposeNeeded = false;
      if (!isOneForm && rawTestID!=rawUnkID 
        && eqn.hasVarID(rawUnkID) && eqn.hasUnkID(rawTestID))
      {
        transposeNeeded = true;
      }


      if (isOneForm)
      {
        SUNDANCE_LEVEL3("find groups", tab2 << "is one-form");
      }
      else
      {
        SUNDANCE_LEVEL3("find groups", tab2 << "is two-form");
      }

      SUNDANCE_LEVEL3("find groups", 
        tab2 << "test ID: " << testID << " block=" << testBlock);

      if (!isOneForm)
      {
        SUNDANCE_LEVEL3("find groups", tab2 << "unk funcID: " << unkID << " block=" << unkBlock);
      }
                   
      SUNDANCE_LEVEL3("find groups", tab2 << "deriv = " << d);
      if (sparsity->isConstant(i))
      {
        SUNDANCE_LEVEL3("find groups", tab2 << "coeff is constant");
      }
      else
      {
        SUNDANCE_LEVEL3("find groups", tab2 << "coeff is non-constant");
      }

      RefCountPtr<ElementIntegral> integral;
      RefCountPtr<ElementIntegral> transposedIntegral;
      int resultIndex;
      if (sparsity->isConstant(i))
      {
        Tabs tab3;
        if (isOneForm)
        {
          int alpha=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          SUNDANCE_LEVEL3("find groups",
            tab3 << "creating reference integral for one-form");
          integral = rcp(new RefIntegral(spatialDim, maxCellType, 
              cellDim, cellType,
              testBasis, alpha, 
              miTest.order(), verbParams));
        }
        else
        {
          int alpha=0;
          int beta=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          if (miUnk.order()==1)
          {
            beta = miUnk.firstOrderDirection();
          }
          SUNDANCE_LEVEL3("find groups",
            tab3 << "creating reference integral for two-form");
          integral = rcp(new RefIntegral(spatialDim, maxCellType,
              cellDim, cellType,
              testBasis, alpha, miTest.order(),
              unkBasis, beta, miUnk.order(), verbParams));
          if (transposeNeeded)
          {
            transposedIntegral = rcp(new RefIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                unkBasis, beta, miUnk.order(),
                testBasis, alpha, miTest.order(), verbParams));
          }
        }
        resultIndex = constCount++;
      }
      else
      {
        Tabs tab3;
        if (isOneForm)
        {
          int alpha=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          SUNDANCE_LEVEL3("find groups",
            tab3 << "creating quadrature integral for one-form");
          integral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
              cellDim, cellType,
              testBasis, alpha, 
              miTest.order(), quad, verbParams));
        }
        else
        {
          int alpha=0;
          int beta=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          if (miUnk.order()==1)
          {
            beta = miUnk.firstOrderDirection();
          }
          SUNDANCE_LEVEL3("find groups",
            tab3 << "creating quadrature integral for two-form");
          integral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
              cellDim, cellType,
              testBasis, alpha, 
              miTest.order(),
              unkBasis, beta, 
              miUnk.order(), quad, verbParams));
          if (transposeNeeded)
          {
            transposedIntegral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                unkBasis, beta, miUnk.order(),
                testBasis, alpha, miTest.order(), quad, verbParams));
          }
        }
        resultIndex = vecCount++;
      }
          
      if (isOneForm)
      {
        if (doGroups)
        {
          OrderedPair<int,BasisFamily> testKey(rawTestID, testBasis);
          if (!oneForms.containsKey(testKey))
          {
            oneForms.put(testKey, tuple(integral));
            oneFormResultIndices.put(testKey, tuple(resultIndex));
          }
          else
          {
            oneForms[testKey].append(integral);
            oneFormResultIndices[testKey].append(resultIndex);
          }
        }
        else
        {
          groups.append(IntegralGroup(tuple(testID), tuple(testBlock),
              tuple(integral),
              tuple(resultIndex), tuple(d), verbParams));
        }
      }
      else
      {
        if (!doGroups)
        {
          groups.append(IntegralGroup(tuple(testID), tuple(testBlock),
              tuple(unkID), tuple(unkBlock),
              tuple(integral),
              tuple(resultIndex), tuple(d), verbParams));
          if (transposeNeeded)
          {
            groups.append(IntegralGroup(tuple(unkID), tuple(unkBlock),
                tuple(testID), tuple(testBlock),
                tuple(transposedIntegral),
                tuple(resultIndex), tuple(d), verbParams));
          }
        }
        else
        {
          Tabs tab3;
          OrderedQuartet<int, BasisFamily, int, BasisFamily> testUnkKey(rawTestID, testBasis, rawUnkID, unkBasis);


          SUNDANCE_LEVEL1("find groups", tab3 << "key=" << testUnkKey);
          if (!twoForms.containsKey(testUnkKey))
          {
            Tabs tab4;
            SUNDANCE_LEVEL1("find groups", tab4 << "key not found");
            twoForms.put(testUnkKey, tuple(integral));
            twoFormResultIndices.put(testUnkKey, tuple(resultIndex));
          }
          else
          {
            Tabs tab4;
            SUNDANCE_LEVEL1("find groups", tab4 << "key found");
            twoForms[testUnkKey].append(integral);
            twoFormResultIndices[testUnkKey].append(resultIndex);
          }
          if (transposeNeeded)
          {
            OrderedQuartet<int, BasisFamily, int, BasisFamily> unkTestKey(rawUnkID, unkBasis, rawTestID, testBasis);
            
            if (!twoForms.containsKey(unkTestKey))
            {
              Tabs tab4;
              SUNDANCE_LEVEL1("find groups", tab4 << "key not found");
              twoForms.put(unkTestKey, tuple(transposedIntegral));
              twoFormResultIndices.put(unkTestKey, tuple(resultIndex));
            }
            else
            {
              Tabs tab4;
              SUNDANCE_LEVEL1("find groups", tab4 << "key found");
              twoForms[unkTestKey].append(transposedIntegral);
              twoFormResultIndices[unkTestKey].append(resultIndex);
            }
          }
        }
      }
    }
  }

  if (doGroups)
  {
    Tabs tab2;
    SUNDANCE_LEVEL1("find groups", tab2 << "creating integral groups");
    for (twoFormMap::const_iterator i=twoForms.begin(); i!=twoForms.end(); i++)
    {
      Tabs tab3;
      SUNDANCE_LEVEL1("find groups", tab3 << "integral group number="
        << groups.size());
      int rawTestID = i->first.a();
      BasisFamily testBasis = i->first.b();
      int rawUnkID = i->first.c();
      BasisFamily unkBasis = i->first.d();
      int testID = eqn.reducedVarID(rawTestID);
      int unkID = eqn.reducedUnkID(rawUnkID);
      int testBlock = eqn.blockForVarID(rawTestID);
      int unkBlock = eqn.blockForUnkID(rawUnkID);
      const Array<RefCountPtr<ElementIntegral> >& integrals = i->second;
      const Array<int>& resultIndices 
        = twoFormResultIndices.get(i->first);
      SUNDANCE_LEVEL1("find groups", tab3 << "creating two-form integral group" << std::endl
        << tab3 << "testID=" << rawTestID << std::endl
        << tab3 << "unkID=" << rawUnkID << std::endl
        << tab3 << "testBlock=" << testBlock << std::endl
        << tab3 << "unkBlock=" << unkBlock << std::endl
        << tab3 << "testBasis=" << testBasis << std::endl
        << tab3 << "unkBasis=" << unkBasis << std::endl
        << tab3 << "resultIndices=" << resultIndices);
      Array<MultipleDeriv> grpDerivs;
      for (unsigned int j=0; j<resultIndices.size(); j++)
      {
        MultipleDeriv d = sparsity->deriv(resultIndices[j]);
        SUNDANCE_LEVEL1("find groups", tab3 << "deriv " << j << " " 
          << d);
        grpDerivs.append(d);
      }
      groups.append(IntegralGroup(tuple(testID), tuple(testBlock), tuple(unkID), 
          tuple(unkBlock),
          integrals, resultIndices, grpDerivs, verbParams));
    }

    for (oneFormMap::const_iterator i=oneForms.begin(); i!=oneForms.end(); i++)
    {
      Tabs tab3;
      SUNDANCE_LEVEL1("find groups", tab3 << "integral group number="
        << groups.size());
      int rawTestID = i->first.first();
      int testID = eqn.reducedVarID(rawTestID);
      int testBlock = eqn.blockForVarID(rawTestID);
      const Array<RefCountPtr<ElementIntegral> >& integrals = i->second;
      const Array<int>& resultIndices 
        = oneFormResultIndices.get(i->first);
      SUNDANCE_LEVEL1("find groups", tab3 << "creating one-form integral group" << std::endl
        << tab3 << "testID=" << testID << std::endl
        << tab3 << "resultIndices=" << resultIndices);
      Array<MultipleDeriv> grpDerivs;
      for (unsigned int j=0; j<resultIndices.size(); j++)
      {
        MultipleDeriv d = sparsity->deriv(resultIndices[j]);
        SUNDANCE_LEVEL1("find groups", tab3 << "deriv " << j << " " 
          << d);
        grpDerivs.append(d);
      }
      groups.append(IntegralGroup(tuple(testID), tuple(testBlock),
          integrals, resultIndices, grpDerivs, verbParams));
    }
  }
  
  
}

