/* @HEADER@ */
/* @HEADER@ */

#include "SundanceIntegral.hpp"
#include "SundanceSumOfIntegrals.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;
using namespace TSFExtended;
using namespace SundanceCore::Internal;

Expr SundanceCore::Integral(const Handle<CellFilterStub>& domain,
              const Expr& integrand)
{
  RefCountPtr<QuadratureFamilyStub> quad 
    = QuadratureFamilyStub::defaultQuadrature();
  return new SumOfIntegrals(domain.ptr(), integrand, quad);
}

Expr SundanceCore::Integral(const Handle<CellFilterStub>& domain,
                        const Expr& integrand,
                        const Handle<QuadratureFamilyStub>& quad)
{
  return new SumOfIntegrals(domain.ptr(), integrand, quad.ptr());
}
