/* @HEADER@ */
/* @HEADER@ */

#include "SundanceEvalRegion.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;

EvalRegion::EvalRegion()
  : id_(-1)
{;}

EvalRegion::EvalRegion(const RefCountPtr<CellFilterBase>& domain,
                       const RefCountPtr<QuadratureFamilyBase>& quad)
  : id_(getID(domain, quad))
{;}

int EvalRegion::getID(const RefCountPtr<CellFilterBase>& domain,
                      const RefCountPtr<QuadratureFamilyBase>& quad)
{
  RegPair p(domain, quad);

  if (!domainAndQuadToIDMap().containsKey(p))
    {
      int id = topID();
      domainAndQuadToIDMap().put(p, id);
    }
  return domainAndQuadToIDMap().get(p);
}

string EvalRegion::toString() const
{
  return "EvalRegion[id="
    + Teuchos::toString(id_) + "]";
}

Map<RegPair, int>& EvalRegion::domainAndQuadToIDMap()
{
  static Map<RegPair, int> rtn = Map<RegPair, int>();
  return rtn;
}





