/* @HEADER@ */
/* @HEADER@ */

#include "SundanceRegionQuadCombo.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

RegionQuadCombo::RegionQuadCombo()
  : id_(-1)
{;}

RegionQuadCombo::RegionQuadCombo(const RefCountPtr<CellFilterStub>& domain,
                       const RefCountPtr<QuadratureFamilyStub>& quad)
  : id_(getID(domain, quad))
{;}

int RegionQuadCombo::getID(const RefCountPtr<CellFilterStub>& domain,
                      const RefCountPtr<QuadratureFamilyStub>& quad)
{
  RegPair p(domain, quad);

  if (!domainAndQuadToIDMap().containsKey(p))
    {
      int id = topID();
      domainAndQuadToIDMap().put(p, id);
    }
  return domainAndQuadToIDMap().get(p);
}

string RegionQuadCombo::toString() const
{
  return "RegionQuadCombo[id="
    + Teuchos::toString(id_) + "]";
}

Map<RegPair, int>& RegionQuadCombo::domainAndQuadToIDMap()
{
  static Map<RegPair, int> rtn = Map<RegPair, int>();
  return rtn;
}





