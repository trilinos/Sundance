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

EvalRegion::EvalRegion(const string& name)
  : id_(getID(name))
{;}

int EvalRegion::getID(const string& name)
{
  if (!nameToIDMap().containsKey(name))
    {
      int id = topID();
      nameToIDMap().put(name, id);
      idToNameMap().put(id, name);
    }
  return nameToIDMap().get(name);
}

string EvalRegion::toString() const
{
  return "EvalRegion[" + idToNameMap().get(id_) + ", id="
    + Teuchos::toString(id_) + "]";
}

Map<string, int>& EvalRegion::nameToIDMap()
{
  static Map<string, int> rtn = Map<string, int>();
  return rtn;
}

Map<int, string>& EvalRegion::idToNameMap()
{
  static Map<int, string> rtn = Map<int, string>();
  return rtn;
}




