/* @HEADER@ */
/* @HEADER@ */

#include "SundanceTrivialGrouper.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceQuadratureIntegral.hpp"
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
                                const CellType& cellType,
                                int cellDim,
                                const QuadratureFamily& quad,
                                const RefCountPtr<SparsityPattern>& sparsity,
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

  for (int i=0; i<sparsity->numDerivs(); i++)
    {
      Tabs tab1;

      TEST_FOR_EXCEPTION(sparsity->isZero(i), InternalError,
                         "Structurally zero functional derivative detected "
                         "at the root level");

      const MultipleDeriv& d = sparsity->deriv(i);
      if (d.order()==0) continue;

      BasisFamily testBasis;
      BasisFamily unkBasis;
      MultiIndex miTest;
      MultiIndex miUnk;
      int testID;
      int unkID;
      bool isOneForm;
      Array<int> alpha;
      Array<int> beta;
      extractWeakForm(eqn, d, testBasis, unkBasis, miTest, miUnk, testID, unkID,
                      isOneForm);

      SUNDANCE_OUT(verb > VerbMedium, 
                   tab1 << "test ID: " << testID);

      SUNDANCE_OUT(!isOneForm && verb > VerbMedium,
                   tab1 << "unk funcID: " << unkID << endl);
                   
      SUNDANCE_OUT(verb > VerbMedium, tab1 << "deriv = " << d);
      SUNDANCE_OUT(verb > VerbMedium && sparsity->isConstant(i), 
                   tab1 << "coeff is constant");
      SUNDANCE_OUT(verb > VerbMedium && !sparsity->isConstant(i), 
                   tab1 << "coeff is non-constant");

      RefCountPtr<ElementIntegral> integral;
      if (sparsity->isConstant(i))
        {
          if (isOneForm)
            {
              if (miTest.order()==1)
                {
                  alpha = tuple(miTest.firstOrderDirection());
                }
              SUNDANCE_OUT(verb > VerbMedium,
                   tab1 << "creating reference integral for one-form");
              integral = rcp(new RefIntegral(cellDim, cellType,
                                             testBasis, alpha, 
                                             miTest.order()));
            }
          else
            {
              if (miTest.order()==1)
                {
                  alpha = tuple(miTest.firstOrderDirection());
                }
              if (miUnk.order()==1)
                {
                  beta = tuple(miUnk.firstOrderDirection());
                }
              SUNDANCE_OUT(verb > VerbMedium,
                   tab1 << "creating reference integral for two-form");
              integral = rcp(new RefIntegral(cellDim, cellType,
                                             testBasis, alpha, miTest.order(),
                                             unkBasis, beta, miUnk.order()));
            }
        }
      else
        {
          if (isOneForm)
            {
              if (miTest.order()==1)
                {
                  alpha = tuple(miTest.firstOrderDirection());
                }
              SUNDANCE_OUT(verb > VerbMedium,
                   tab1 << "creating quadrature integral for two-form");
              integral = rcp(new QuadratureIntegral(cellDim, cellType,
                                                    testBasis, alpha, 
                                                    miTest.order(), quad));
            }
          else
            {
              if (miTest.order()==1)
                {
                  alpha = tuple(miTest.firstOrderDirection());
                }
              if (miUnk.order()==1)
                {
                  beta = tuple(miUnk.firstOrderDirection());
                }
              SUNDANCE_OUT(verb > VerbMedium,
                   tab1 << "creating quadrature integral for two-form");
              integral = rcp(new QuadratureIntegral(cellDim, cellType,
                                                    testBasis, alpha, 
                                                    miTest.order(),
                                                    unkBasis, beta, 
                                                    miUnk.order(), quad));
            }
        }

      if (isOneForm)
        {
          groups.append(IntegralGroup(tuple(testID), tuple(integral),
                                      tuple(tuple(i))));
        }
      else
        {
          groups.append(IntegralGroup(tuple(testID), tuple(unkID),
                                      tuple(integral),
                                      tuple(tuple(i))));
        }
    }
}

