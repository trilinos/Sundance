/* @HEADER@ */
/* @HEADER@ */

#include "SundanceRegionQuadCombo.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

RegionQuadCombo::RegionQuadCombo()
  : id_(-1), domain_(), quad_()
{;}

RegionQuadCombo::RegionQuadCombo(const RefCountPtr<CellFilterStub>& domain,
                       const RefCountPtr<QuadratureFamilyStub>& quad)
  : id_(getID(domain, quad)), domain_(domain), quad_(quad)
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
    + Teuchos::toString(id_) + "cell filter=" + domain_->describe() 
    + " quad=" + quad_->describe() + "]";
}

Map<RegPair, int>& RegionQuadCombo::domainAndQuadToIDMap()
{
  static Map<RegPair, int> rtn = Map<RegPair, int>();
  return rtn;
}





