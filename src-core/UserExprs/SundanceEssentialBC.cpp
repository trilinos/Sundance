/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEssentialBC.hpp"
#include "SundanceSumOfBCs.hpp"
#include "SundanceZeroExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;
using namespace TSFExtended;
using namespace SundanceCore::Internal;

Expr SundanceCore::EssentialBC(const Handle<CellFilterStub>& domain,
                               const Expr& integrand)
{
  const ZeroExpr* z = dynamic_cast<const ZeroExpr*>(integrand.ptr().get());
  const ConstantExpr* c = dynamic_cast<const ConstantExpr*>(integrand.ptr().get());
  if (z != 0 || (c != 0 && c->value()==0.0))
    {
      return integrand;
    }
  RefCountPtr<QuadratureFamilyStub> quad = rcp(new QuadratureFamilyStub(0));
  return new SumOfBCs(domain.ptr(), integrand, quad);
}

Expr SundanceCore::EssentialBC(const Handle<CellFilterStub>& domain,
                               const Expr& integrand,
                               const Handle<QuadratureFamilyStub>& quad)
{
  const ZeroExpr* z = dynamic_cast<const ZeroExpr*>(integrand.ptr().get());
  const ConstantExpr* c = dynamic_cast<const ConstantExpr*>(integrand.ptr().get());
  if (z != 0 || (c != 0 && c->value()==0.0))
    {
      return integrand;
    }
  return new SumOfBCs(domain.ptr(), integrand, quad.ptr());
}
