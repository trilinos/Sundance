/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEssentialBC.hpp"
#include "SundanceSumOfBCs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;
using namespace TSFExtended;
using namespace SundanceCore::FrameworkInterface;

Expr SundanceCore::EssentialBC(const Handle<CellFilterBase>& domain,
                               const Expr& integrand)
{
  RefCountPtr<QuadratureFamilyBase> quad = rcp(new QuadratureFamilyBase(0));
  return new SumOfBCs(domain.ptr(), integrand, quad);
}

Expr SundanceCore::EssentialBC(const Handle<CellFilterBase>& domain,
                               const Expr& integrand,
                               const Handle<QuadratureFamilyBase>& quad)
{
  return new SumOfBCs(domain.ptr(), integrand, quad.ptr());
}
