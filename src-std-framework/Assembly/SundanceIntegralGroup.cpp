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

static Time& integrationTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("integration"); 
  return *rtn;
}


IntegralGroup
::IntegralGroup(const Array<RefCountPtr<ElementIntegral> >& integrals,
                const Array<Array<int> >& resultIndices)
  : order_(0),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(),
    unkID_(),
    integrals_(integrals),
    resultIndices_(resultIndices)
{
  verbosity() = classVerbosity();
}

IntegralGroup
::IntegralGroup(const Array<int>& testID,
                const Array<RefCountPtr<ElementIntegral> >& integrals,
                const Array<Array<int> >& resultIndices)
  : order_(1),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(testID),
    unkID_(),
    integrals_(integrals),
    resultIndices_(resultIndices)
{
  verbosity() = classVerbosity();

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
  : order_(2),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(testID),
    unkID_(unkID),
    integrals_(integrals),
    resultIndices_(resultIndices)
{
  verbosity() = classVerbosity();

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

bool IntegralGroup
::evaluate(const CellJacobianBatch& J,
           const Array<RefCountPtr<EvalVector> >& vectorCoeffs,
           const Array<double>& constantCoeffs,
           RefCountPtr<Array<double> >& A) const
{
  TimeMonitor timer(integrationTimer());
  Tabs tab0;


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
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       tab << "Integrating term group " << i 
                       << " by transformed reference integral");
          Array<double> f(resultIndices_[i].size());
          for (int j=0; j<f.size(); j++)
            {
              f[j] = constantCoeffs[resultIndices_[i][j]];
            }
          SUNDANCE_OUT(verbosity() > VerbSilent,
                       tab << "Coefficients are " << f);
          ref->transform(J, f, A);
        }
      else 
        {
          SUNDANCE_OUT(verbosity() > VerbMedium,
                       tab << "Integrating term group " << i 
                       << " by quadrature");
          TEST_FOR_EXCEPTION(vectorCoeffs[resultIndices_[i][0]]->length()==0,
                             InternalError,
                             "zero-length coeff vector detected in "
                             "quadrature integration branch of "
                             "IntegralGroup::evaluate(). String value is ["
                             << vectorCoeffs[resultIndices_[i][0]]->str()
                             << "]");
          const double* const f = vectorCoeffs[resultIndices_[i][0]]->start();
          quad->transform(J, f, A);
        }
    }
  SUNDANCE_OUT(verbosity() > VerbSilent, tab0 << "done");

  return true;
}


