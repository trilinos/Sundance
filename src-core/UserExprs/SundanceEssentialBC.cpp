/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEssentialBC.hpp"
#include "SundanceSumOfBCs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace Teuchos;
using namespace SundanceCore::Internal;
using namespace TSFExtended;
using namespace SundanceCore::Internal;

Expr SundanceCore::EssentialBC(const Handle<CellFilterStub>& domain,
                               const Expr& integrand)
{
  RefCountPtr<QuadratureFamilyStub> quad = rcp(new QuadratureFamilyStub(0));
  return new SumOfBCs(domain.ptr(), integrand, quad);
}

Expr SundanceCore::EssentialBC(const Handle<CellFilterStub>& domain,
                               const Expr& integrand,
                               const Handle<QuadratureFamilyStub>& quad)
{
  return new SumOfBCs(domain.ptr(), integrand, quad.ptr());
}
