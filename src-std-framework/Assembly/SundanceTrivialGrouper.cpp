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
using namespace SundanceCore::Internal;
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
  VerbositySetting verb = GrouperBase::classVerbosity();
  SUNDANCE_OUT(verb > VerbLow, 
               tab << "trivial grouper num derivs = " << sparsity->numDerivs() << endl);
  SUNDANCE_OUT(verb > VerbLow, 
               tab << "cell type = " << cellType);

  SUNDANCE_OUT(verb > VerbMedium,  
               tab << "sparsity = " << endl << *sparsity << endl);

  int vecCount=0;
  int constCount=0;

  bool doGroups = true;

  typedef SundanceUtils::Map<OrderedPair<int, int>, Array<RefCountPtr<ElementIntegral> > > twoFormMap;
  typedef SundanceUtils::Map<int, Array<RefCountPtr<ElementIntegral> > > oneFormMap;
  SundanceUtils::Map<OrderedPair<int, int>, Array<RefCountPtr<ElementIntegral> > > twoForms;
  SundanceUtils::Map<OrderedPair<int, int>, Array<int> > twoFormResultIndices;
  SundanceUtils::Map<int, Array<RefCountPtr<ElementIntegral> > > oneForms;
  SundanceUtils::Map<int, Array<int> > oneFormResultIndices;

  for (int i=0; i<sparsity->numDerivs(); i++)
    {
      Tabs tab1;

      const MultipleDeriv& d = sparsity->deriv(i);
      if (d.order()==0) 
        {
          RefCountPtr<ElementIntegral> integral ;
          int resultIndex;
          if (sparsity->isConstant(i))
            {
              integral = rcp(new RefIntegral(spatialDim, maxCellType, 
                                             cellDim, cellType));
              resultIndex = constCount++;
            }
          else
            {
              integral = rcp(new QuadratureIntegral(spatialDim, maxCellType, 
                                                    cellDim, cellType, quad));
              resultIndex = vecCount++;
            }
          groups.append(IntegralGroup(tuple(integral),
                                      tuple(resultIndex)));
        }
      else
        {

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


          SUNDANCE_VERB_MEDIUM(tab1 << "test ID: " << testID << " block=" << testBlock);

          if (!isOneForm)
            {
              SUNDANCE_VERB_MEDIUM(tab1 << "unk funcID: " << unkID << " block=" << unkBlock);
            }
                   
          SUNDANCE_OUT(verb > VerbMedium, tab1 << "deriv = " << d);
          SUNDANCE_OUT(verb > VerbMedium && sparsity->isConstant(i), 
                       tab1 << "coeff is constant");
          SUNDANCE_OUT(verb > VerbMedium && !sparsity->isConstant(i), 
                       tab1 << "coeff is non-constant");

          RefCountPtr<ElementIntegral> integral;
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
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab1 << "creating reference integral for one-form");
                  integral = rcp(new RefIntegral(spatialDim, maxCellType, 
                                                 cellDim, cellType,
                                                 testBasis, alpha, 
                                                 miTest.order()));
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
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab1 << "creating reference integral for two-form");
                  integral = rcp(new RefIntegral(spatialDim, maxCellType,
                                                 cellDim, cellType,
                                                 testBasis, alpha, miTest.order(),
                                                 unkBasis, beta, miUnk.order()));
                }
              resultIndex = constCount++;
            }
          else
            {
              if (isOneForm)
                {
                  int alpha=0;
                  if (miTest.order()==1)
                    {
                      alpha = miTest.firstOrderDirection();
                    }
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab1 << "creating quadrature integral for two-form");
                  integral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
                                                        cellDim, cellType,
                                                        testBasis, alpha, 
                                                        miTest.order(), quad));
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
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab1 << "creating quadrature integral for two-form");
                  integral = rcp(new QuadratureIntegral(spatialDim, maxCellType,
                                                        cellDim, cellType,
                                                        testBasis, alpha, 
                                                        miTest.order(),
                                                        unkBasis, beta, 
                                                        miUnk.order(), quad));
                }
              resultIndex = vecCount++;
            }
          
          if (isOneForm)
            {
              if (doGroups)
                {
                  if (!oneForms.containsKey(rawTestID))
                    {
                      oneForms.put(rawTestID, tuple(integral));
                      oneFormResultIndices.put(rawTestID, tuple(resultIndex));
                    }
                  else
                    {
                      oneForms[rawTestID].append(integral);
                      oneFormResultIndices[rawTestID].append(resultIndex);
                    }
                }
              else
                {
                  groups.append(IntegralGroup(tuple(testID), tuple(testBlock),
                                              tuple(integral),
                                              tuple(resultIndex)));
                }
            }
          else
            {
              if (!doGroups)
                {
                  groups.append(IntegralGroup(tuple(testID), tuple(testBlock),
                                              tuple(unkID), tuple(unkBlock),
                                              tuple(integral),
                                              tuple(resultIndex)));
                }
              else
                {
                  OrderedPair<int, int> testUnk(rawTestID, rawUnkID);
                  if (!twoForms.containsKey(testUnk))
                    {
                      twoForms.put(testUnk, tuple(integral));
                      twoFormResultIndices.put(testUnk, tuple(resultIndex));
                    }
                  else
                    {
                      twoForms[testUnk].append(integral);
                      twoFormResultIndices[testUnk].append(resultIndex);
                    }
                }
            }
        }
    }

  if (doGroups)
    {
      for (twoFormMap::const_iterator i=twoForms.begin(); i!=twoForms.end(); i++)
        {
          int rawTestID = i->first.first();
          int rawUnkID = i->first.second();
          int testID = eqn.reducedVarID(rawTestID);
          int unkID = eqn.reducedUnkID(rawUnkID);
          int testBlock = eqn.blockForVarID(rawTestID);
          int unkBlock = eqn.blockForUnkID(rawUnkID);
          const Array<RefCountPtr<ElementIntegral> >& integrals = i->second;
          const Array<int>& resultIndices 
            = twoFormResultIndices.get(i->first);
          SUNDANCE_OUT(verb > VerbLow, tab << "creating integral group" << endl
                               << tab << "testID=" << testID << endl
                               << tab << "unkID=" << unkID << endl
                               << tab << "resultIndices=" << resultIndices);
          for (unsigned int j=0; j<resultIndices.size(); j++)
            {
              SUNDANCE_OUT(verb > VerbLow, tab << "deriv " << j << " " 
                           << sparsity->deriv(resultIndices[j]));
            }
          groups.append(IntegralGroup(tuple(testID), tuple(testBlock), tuple(unkID), 
                                      tuple(unkBlock),
                                      integrals, resultIndices));
        }

      for (oneFormMap::const_iterator i=oneForms.begin(); i!=oneForms.end(); i++)
        {
          int rawTestID = i->first;
          int testID = eqn.reducedVarID(rawTestID);
          int testBlock = eqn.blockForVarID(rawTestID);
          const Array<RefCountPtr<ElementIntegral> >& integrals = i->second;
          const Array<int>& resultIndices 
            = oneFormResultIndices.get(i->first);
          SUNDANCE_OUT(verb > VerbLow, tab << "creating integral group" << endl
                               << tab << "testID=" << testID << endl
                               << tab << "resultIndices=" << resultIndices);
          for (unsigned int j=0; j<resultIndices.size(); j++)
            {
              SUNDANCE_OUT(verb > VerbLow, tab << "deriv " << j << " " 
                           << sparsity->deriv(resultIndices[j]));
            }
          groups.append(IntegralGroup(tuple(testID), tuple(testBlock),
                                      integrals, resultIndices));
        }
    }
  
  
}

