/* @HEADER@ */
/* @HEADER@ */

#include "SundanceIntegralGroup.hpp"
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

IntegralGroup
::IntegralGroup(const Array<int>& testID,
                const Array<RefCountPtr<ElementIntegral> >& integrals,
                const Array<Array<int> >& resultIndices)
  : isTwoForm_(false),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(testID),
    unkID_(),
    integrals_(integrals),
    resultIndices_(resultIndices)
{
  for (int i=0; i<integrals.size(); i++)
    {
      int nt = integrals[i]->nNodesTest();
      if (i > 0) 
        {
          TEST_FOR_EXCEPTION(nt != nTestNodes_, InternalError,
                             "IntegralGroup ctor detected integrals with inconsistent numbers of test nodes");
        }
      nTestNodes_ = nt;
    }
}

IntegralGroup
::IntegralGroup(const Array<int>& testID,
                const Array<int>& unkID,
                const Array<RefCountPtr<ElementIntegral> >& integrals,
                const Array<Array<int> >& resultIndices)
  : isTwoForm_(true),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(testID),
    unkID_(unkID),
    integrals_(integrals),
    resultIndices_(resultIndices)
{
  for (int i=0; i<integrals.size(); i++)
    {
      int nt = integrals[i]->nNodesTest();
      int nu = integrals[i]->nNodesUnk();
      if (i > 0) 
        {
          TEST_FOR_EXCEPTION(nt != nTestNodes_, InternalError,
                             "IntegralGroup ctor detected integrals with inconsistent numbers of test nodes");
          TEST_FOR_EXCEPTION(nu != nUnkNodes_, InternalError,
                             "IntegralGroup ctor detected integrals with inconsistent numbers of unk nodes");
        }
      nTestNodes_ = nt;
      nUnkNodes_ = nu;
    }
}

bool IntegralGroup::evaluate(const CellJacobianBatch& J,
                             const RefCountPtr<EvalVectorArray>& coeffs,
                             RefCountPtr<Array<double> >& A) const
{
  Tabs tab0;
  bool nonzero = false;

  SUNDANCE_OUT(verbosity() > VerbSilent,
               tab0 << "evaluating an integral group of size "
               << integrals_.size());

  for (int i=0; i<integrals_.size(); i++)
    {
      Tabs tab;
      const RefIntegral* ref 
        = dynamic_cast<const RefIntegral*>(integrals_[i].get());
      const QuadratureIntegral* quad 
        = dynamic_cast<const QuadratureIntegral*>(integrals_[i].get());

      if (ref!=0)
        {
          SUNDANCE_OUT(verbosity() > VerbSilent,
                       tab << "Integrating term group " << i 
                       << " by transformed reference integral");
          Array<double> f(resultIndices_[i].size());
          for (int j=0; j<f.size(); j++)
            {
              f[j] = (*coeffs)[resultIndices_[i][j]]->getConstantValue();
            }
          SUNDANCE_OUT(verbosity() > VerbSilent,
                       tab << "Coefficients are " << f);
          ref->transform(J, f, A);
        }
      else 
        {
          if ((*coeffs)[resultIndices_[i][0]]->isZero()) continue;
          SUNDANCE_OUT(verbosity() > VerbSilent,
                       tab << "Integrating term group " << i 
                       << " by quadrature");
          const double* const f = (*coeffs)[resultIndices_[i][0]]->start();
          quad->transform(J, f, A);
        }
      nonzero = true;
    }
  SUNDANCE_OUT(verbosity() > VerbSilent, tab0 << "done");

  return nonzero;
}


