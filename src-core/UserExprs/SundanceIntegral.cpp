/* @HEADER@ */
/* @HEADER@ */

#include "SundanceIntegral.hpp"
#include "SundanceSumOfIntegrals.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;
using namespace TSFExtended;
using namespace SundanceCore::FrameworkInterface;

Expr SundanceCore::Integral(const Handle<CellFilterBase>& domain,
              const Expr& integrand)
{
  RefCountPtr<QuadratureFamilyBase> quad = rcp(new QuadratureFamilyBase(0));
  return new SumOfIntegrals(domain.ptr(), integrand, quad);
}

Expr SundanceCore::Integral(const Handle<CellFilterBase>& domain,
                        const Expr& integrand,
                        const Handle<QuadratureFamilyBase>& quad)
{
  return new SumOfIntegrals(domain.ptr(), integrand, quad.ptr());
}
