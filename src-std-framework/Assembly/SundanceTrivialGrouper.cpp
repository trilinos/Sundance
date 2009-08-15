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
#include "SundanceEquationSet.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceQuadratureFamily.hpp"
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


void TrivialGrouper::findGroups(const EquationSet& eqn,
  const CellType& maxCellType,
  int spatialDim,
  const CellType& cellType,
  int cellDim,
  const QuadratureFamily& quad,
  const RefCountPtr<SparsitySuperset>& sparsity,
  bool isInternalBdry,
  Array<RCP<IntegralGroup> >& groups) const
{
  Tabs tab(0);

  SUNDANCE_MSG1(setupVerb(),
    tab << "in TrivialGrouper::findGroups(), num derivs = " 
    << sparsity->numDerivs());
  SUNDANCE_MSG1(setupVerb(), 
    tab << "cell type = " << cellType);
  SUNDANCE_MSG1(setupVerb(), 
    tab << "sparsity = " << std::endl << *sparsity << std::endl);

  int vecCount=0;
  int constCount=0;

  /* turn off grouping for submaximal cells. This works around 
   * a bug detected by Rob Kirby that
   * shows up with Nitsche BCs in mixed-element discretizations */
  bool doGroups = true;
  if (cellType != maxCellType) doGroups = false;

  typedef SundanceUtils::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<RefCountPtr<ElementIntegral> > > twoFormMap;
  typedef SundanceUtils::Map<OrderedTriple<int,int,BasisFamily>, Array<RefCountPtr<ElementIntegral> > > oneFormMap;
  SundanceUtils::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<RefCountPtr<ElementIntegral> > > twoForms;
  SundanceUtils::Map<OrderedQuartet<int, BasisFamily, int, BasisFamily>, Array<int> > twoFormResultIndices;
  SundanceUtils::Map<OrderedTriple<int,int,BasisFamily>, Array<RefCountPtr<ElementIntegral> > > oneForms;
  SundanceUtils::Map<OrderedTriple<int,int,BasisFamily>, Array<int> > oneFormResultIndices;

  for (int i=0; i<sparsity->numDerivs(); i++)
  {
    const MultipleDeriv& d = sparsity->deriv(i);
    SUNDANCE_MSG3(setupVerb(),
      tab << "--------------------------------------------------");
    SUNDANCE_MSG3(setupVerb(),
      tab << "defining integration policy for " << d);
    SUNDANCE_MSG3(setupVerb(),
      tab << "--------------------------------------------------");
      
    if (d.order()==0) 
    {
      RefCountPtr<ElementIntegral> integral ;
      int resultIndex;
      if (sparsity->isConstant(i))
      {
        integral = rcp(new RefIntegral(spatialDim, maxCellType, 
            cellDim, cellType, isInternalBdry, setupVerb()));
        resultIndex = constCount++;
      }
      else
      {
        integral = rcp(new QuadratureIntegral(spatialDim, maxCellType, 
            cellDim, cellType, quad, isInternalBdry,
            setupVerb()));
        resultIndex = vecCount++;
      }
      integral->setVerbosity(integrationVerb(), transformVerb());
      SUNDANCE_MSG3(setupVerb(), tab << "is zero-form");
      groups.append(rcp(new IntegralGroup(tuple(integral),
            tuple(resultIndex), setupVerb())));
    }
    else
    {
      BasisFamily testBasis;
      BasisFamily unkBasis;
      MultiIndex miTest;
      MultiIndex miUnk;
      int rawTestID = -1;
      int rawUnkID = -1;
      int rawParamID = -1;
      int testID = -1;
      int unkID = -1;
      int paramID = -1;
      int testBlock = -1;
      int unkBlock = -1;
      bool isOneForm;
      bool hasParam;
      extractWeakForm(eqn, d, testBasis, unkBasis, miTest, miUnk, 
        rawTestID, rawUnkID,
        testID, unkID,
        testBlock, unkBlock,
        rawParamID, paramID,
        isOneForm, hasParam);
      
      TEST_FOR_EXCEPT(hasParam && !isOneForm);

      /* The parameter index acts as an index into a multivector. If
       * this one-form is not a parametric derivative, we use zero as
       * the multivector index */
      int mvIndex = 0;
      if (hasParam) mvIndex = paramID; 
      
      /* In variational problems we might have (u,v) and (v,u). Because 
       * the derivative is stored as an unordered multiset it can't 
       * distinguish between the two cases. We need to check the equation 
       * set to see if the two functions show up as variations and 
       * unknowns. If so, then we need to produce the transposed integral.
       */
      bool transposeNeeded = false;
      if (!isOneForm && rawTestID!=rawUnkID 
        && eqn.hasVarID(rawUnkID) && eqn.hasUnkID(rawTestID))
      {
        transposeNeeded = true;
      }


      if (isOneForm)
      {
        SUNDANCE_MSG3(setupVerb(), tab << "is one-form");
      }
      else
      {
        SUNDANCE_MSG3(setupVerb(), tab << "is two-form");
      }


      if (hasParam)
      {
        SUNDANCE_MSG3(setupVerb(), tab << "is a parametric derivative");
      }
      else
      {
        SUNDANCE_MSG3(setupVerb(), tab << "is not a parametric derivative");
      }

      SUNDANCE_MSG3(setupVerb(), 
        tab << "test ID: " << testID << " block=" << testBlock);

      if (!isOneForm)
      {
        SUNDANCE_MSG3(setupVerb(), tab << "unk funcID: " << unkID << " block=" << unkBlock);
      }

      if (hasParam)
      {
        SUNDANCE_MSG3(setupVerb(), tab << "paramID=" << paramID);
      }
                   
      SUNDANCE_MSG3(setupVerb(), tab << "deriv = " << d);
      if (sparsity->isConstant(i))
      {
        SUNDANCE_MSG3(setupVerb(), tab << "coeff is constant");
      }
      else
      {
        SUNDANCE_MSG3(setupVerb(), tab << "coeff is non-constant");
      }

      RefCountPtr<ElementIntegral> integral;
      RefCountPtr<ElementIntegral> transposedIntegral;
      int resultIndex;
      if (sparsity->isConstant(i))
      {
        if (isOneForm)
        {
          int alpha=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          integral = rcp(new RefIntegral(spatialDim, maxCellType, 
              cellDim, cellType,
              testBasis, alpha, 
              miTest.order(), isInternalBdry, setupVerb()));
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
          integral = rcp(new RefIntegral(spatialDim, maxCellType,
              cellDim, cellType,
              testBasis, alpha, miTest.order(),
              unkBasis, beta, miUnk.order(), isInternalBdry, setupVerb()));
          if (transposeNeeded)
          {
            transposedIntegral = rcp(new RefIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                unkBasis, beta, miUnk.order(),
                testBasis, alpha, miTest.order(), isInternalBdry, setupVerb()));
          }
        }
        resultIndex = constCount++;
      }
      else /* sparsity->isVector(i) */
      {
        if (isOneForm)
        {
          int alpha=0;
          if (miTest.order()==1)
          {
            alpha = miTest.firstOrderDirection();
          }
          integral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
              cellDim, cellType,
              testBasis, alpha, 
              miTest.order(), quad, isInternalBdry, setupVerb()));
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
          integral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
              cellDim, cellType,
              testBasis, alpha, 
              miTest.order(),
              unkBasis, beta, 
              miUnk.order(), quad, isInternalBdry, setupVerb()));
          if (transposeNeeded)
          {
            transposedIntegral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
                cellDim, cellType,
                unkBasis, beta, miUnk.order(),
                testBasis, alpha, miTest.order(), quad, isInternalBdry, setupVerb()));
          }
        }
        resultIndex = vecCount++;
      }

      /* Set the verbosity for the integrals */
      integral->setVerbosity(integrationVerb(), transformVerb());
      if (transposeNeeded)
      {
        transposedIntegral->setVerbosity(integrationVerb(), transformVerb());
      }
          
      
      if (isOneForm)
      {
        if (doGroups)
        {
          OrderedTriple<int,int,BasisFamily> testKey(rawTestID, mvIndex, testBasis);
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
          groups.append(rcp(new IntegralGroup(tuple(testID), tuple(testBlock),
                tuple(mvIndex),
                tuple(integral),
                tuple(resultIndex), tuple(d), setupVerb())));
        }
      }
      else
      {
        if (!doGroups)
        {
          groups.append(rcp(new IntegralGroup(tuple(testID), tuple(testBlock),
                tuple(unkID), tuple(unkBlock),
                tuple(integral),
                tuple(resultIndex), tuple(d), setupVerb())));
          if (transposeNeeded)
          {
            groups.append(rcp(new IntegralGroup(tuple(unkID), tuple(unkBlock),
                  tuple(testID), tuple(testBlock),
                  tuple(transposedIntegral),
                  tuple(resultIndex), tuple(d), setupVerb())));
          }
        }
        else
        {
          Tabs tab3;
          OrderedQuartet<int, BasisFamily, int, BasisFamily> testUnkKey(rawTestID, testBasis, rawUnkID, unkBasis);


          SUNDANCE_MSG2(setupVerb(), tab3 << "key=" << testUnkKey);
          if (!twoForms.containsKey(testUnkKey))
          {
            Tabs tab4;
            SUNDANCE_MSG2(setupVerb(), tab4 << "key not found");
            twoForms.put(testUnkKey, tuple(integral));
            twoFormResultIndices.put(testUnkKey, tuple(resultIndex));
          }
          else
          {
            Tabs tab4;
            SUNDANCE_MSG2(setupVerb(), tab4 << "key found");
            twoForms[testUnkKey].append(integral);
            twoFormResultIndices[testUnkKey].append(resultIndex);
          }
          if (transposeNeeded)
          {
            OrderedQuartet<int, BasisFamily, int, BasisFamily> unkTestKey(rawUnkID, unkBasis, rawTestID, testBasis);
            
            if (!twoForms.containsKey(unkTestKey))
            {
              Tabs tab4;
              SUNDANCE_MSG2(setupVerb(), tab4 << "key not found");
              twoForms.put(unkTestKey, tuple(transposedIntegral));
              twoFormResultIndices.put(unkTestKey, tuple(resultIndex));
            }
            else
            {
              Tabs tab4;
              SUNDANCE_MSG2(setupVerb(), tab4 << "key found");
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
    Tabs tab;
    SUNDANCE_MSG2(setupVerb(), tab << "creating integral groups");
    for (twoFormMap::const_iterator i=twoForms.begin(); i!=twoForms.end(); i++)
    {
      Tabs tab3;
      SUNDANCE_MSG2(setupVerb(), tab3 << "integral group number="
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
      SUNDANCE_MSG2(setupVerb(), tab3 << "creating two-form integral group" << std::endl
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
        SUNDANCE_MSG2(setupVerb(), tab3 << "deriv " << j << " " 
          << d);
        grpDerivs.append(d);
      }
      groups.append(rcp(new IntegralGroup(tuple(testID), tuple(testBlock), 
            tuple(unkID), tuple(unkBlock),
            integrals, resultIndices, grpDerivs, setupVerb())));
    }

    for (oneFormMap::const_iterator i=oneForms.begin(); i!=oneForms.end(); i++)
    {
      Tabs tab3;
      SUNDANCE_MSG2(setupVerb(), tab3 << "integral group number="
        << groups.size());
      int rawTestID = i->first.a();
      int mvIndex = i->first.b();
      int testID = eqn.reducedVarID(rawTestID);
      int testBlock = eqn.blockForVarID(rawTestID);
      const Array<RefCountPtr<ElementIntegral> >& integrals = i->second;
      const Array<int>& resultIndices 
        = oneFormResultIndices.get(i->first);
      SUNDANCE_MSG2(setupVerb(), tab3 << "creating one-form integral group" << std::endl
        << tab3 << "testID=" << testID << std::endl
        << tab3 << "resultIndices=" << resultIndices);
      Array<MultipleDeriv> grpDerivs;
      for (unsigned int j=0; j<resultIndices.size(); j++)
      {
        MultipleDeriv d = sparsity->deriv(resultIndices[j]);
        SUNDANCE_MSG2(setupVerb(), tab3 << "deriv " << j << " " 
          << d);
        grpDerivs.append(d);
      }
      groups.append(rcp(new IntegralGroup(tuple(testID), tuple(testBlock),
            tuple(mvIndex),
            integrals, resultIndices, grpDerivs, setupVerb())));
    }
  }
  
  
}

