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
               tab << "trivial grouper num derivs = " << sparsity->numDerivs() << std::endl);
  SUNDANCE_OUT(verb > VerbLow, 
               tab << "cell type = " << cellType);

  SUNDANCE_OUT(verb > VerbMedium,  
               tab << "sparsity = " << std::endl << *sparsity << std::endl);

  int vecCount=0;
  int constCount=0;

  /* turn off grouping for BCs. This works around a bug detected by Rob Kirby that
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
      SUNDANCE_VERB_MEDIUM(tab1 << "defining integration policy for "
                           << d);
      
      if (d.order()==0) 
        {
          Tabs tab2;
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
          SUNDANCE_VERB_MEDIUM(tab2 << "is zero-form");
          groups.append(IntegralGroup(tuple(integral),
                                      tuple(resultIndex)));
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


          if (isOneForm)
            {
              SUNDANCE_VERB_MEDIUM(tab2 << "is one-form");
            }
          else
            {
              SUNDANCE_VERB_MEDIUM(tab2 << "is two-form");
            }

          SUNDANCE_VERB_MEDIUM(tab2 << "test ID: " << testID << " block=" << testBlock);

          if (!isOneForm)
            {
              SUNDANCE_VERB_MEDIUM(tab2 << "unk funcID: " << unkID << " block=" << unkBlock);
            }
                   
          SUNDANCE_OUT(verb > VerbMedium, tab2 << "deriv = " << d);
          SUNDANCE_OUT(verb > VerbMedium && sparsity->isConstant(i), 
                       tab2 << "coeff is constant");
          SUNDANCE_OUT(verb > VerbMedium && !sparsity->isConstant(i), 
                       tab2 << "coeff is non-constant");

          RefCountPtr<ElementIntegral> integral;
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
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab3 << "creating reference integral for one-form");
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
                               tab3 << "creating reference integral for two-form");
                  integral = rcp(new RefIntegral(spatialDim, maxCellType,
                                                 cellDim, cellType,
                                                 testBasis, alpha, miTest.order(),
                                                 unkBasis, beta, miUnk.order()));
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
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab3 << "creating quadrature integral for one-form");
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
                               tab3 << "creating quadrature integral for two-form");
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
                  Tabs tab3;
                  OrderedQuartet<int, BasisFamily, int, BasisFamily> testUnkKey(rawTestID, testBasis, rawUnkID, unkBasis);

                  SUNDANCE_OUT(verb > VerbLow, tab3 << "key=" << testUnkKey);
                  if (!twoForms.containsKey(testUnkKey))
                    {
                      Tabs tab4;
                      SUNDANCE_OUT(verb > VerbLow, tab4 << "key not found");
                      twoForms.put(testUnkKey, tuple(integral));
                      twoFormResultIndices.put(testUnkKey, tuple(resultIndex));
                    }
                  else
                    {
                      Tabs tab4;
                      SUNDANCE_OUT(verb > VerbLow, tab4 << "key found");
                      twoForms[testUnkKey].append(integral);
                      twoFormResultIndices[testUnkKey].append(resultIndex);
                    }
                }
            }
        }
    }

  if (doGroups)
    {
      Tabs tab2;
      SUNDANCE_OUT(verb > VerbLow, tab2 << "creating integral groups");
      for (twoFormMap::const_iterator i=twoForms.begin(); i!=twoForms.end(); i++)
        {
          Tabs tab3;
          SUNDANCE_OUT(verb > VerbLow, tab3 << "integral group number="
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
          SUNDANCE_OUT(verb > VerbLow, tab3 << "creating two-form integral group" << std::endl
                               << tab3 << "testID=" << rawTestID << std::endl
                               << tab3 << "unkID=" << rawUnkID << std::endl
                               << tab3 << "testBlock=" << testBlock << std::endl
                               << tab3 << "unkBlock=" << unkBlock << std::endl
                               << tab3 << "testBasis=" << testBasis << std::endl
                               << tab3 << "unkBasis=" << unkBasis << std::endl
                               << tab3 << "resultIndices=" << resultIndices);
          for (unsigned int j=0; j<resultIndices.size(); j++)
            {
              SUNDANCE_OUT(verb > VerbLow, tab3 << "deriv " << j << " " 
                           << sparsity->deriv(resultIndices[j]));
            }
          groups.append(IntegralGroup(tuple(testID), tuple(testBlock), tuple(unkID), 
                                      tuple(unkBlock),
                                      integrals, resultIndices));
        }

      for (oneFormMap::const_iterator i=oneForms.begin(); i!=oneForms.end(); i++)
        {
          Tabs tab3;
          SUNDANCE_OUT(verb > VerbLow, tab3 << "integral group number="
                       << groups.size());
          int rawTestID = i->first.first();
          int testID = eqn.reducedVarID(rawTestID);
          int testBlock = eqn.blockForVarID(rawTestID);
          const Array<RefCountPtr<ElementIntegral> >& integrals = i->second;
          const Array<int>& resultIndices 
            = oneFormResultIndices.get(i->first);
          SUNDANCE_OUT(verb > VerbLow, tab3 << "creating one-form integral group" << std::endl
                               << tab3 << "testID=" << testID << std::endl
                               << tab3 << "resultIndices=" << resultIndices);
          for (unsigned int j=0; j<resultIndices.size(); j++)
            {
              SUNDANCE_OUT(verb > VerbLow, tab3 << "deriv " << j << " " 
                           << sparsity->deriv(resultIndices[j]));
            }
          groups.append(IntegralGroup(tuple(testID), tuple(testBlock),
                                      integrals, resultIndices));
        }
    }
  
  
}

