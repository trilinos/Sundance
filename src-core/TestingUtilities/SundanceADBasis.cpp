/* @HEADER@ */
/* @HEADER@ */


#include "SundanceADBasis.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceTesting;
using namespace TSFExtended;
using namespace Teuchos;
using namespace std;

ADBasis::ADBasis(int order)
  : m_(order+1)
{;}

ADReal ADBasis::evaluate(const Point& pt) const
{
  ADReal x(pt[0], 0, 3);
  ADReal y(pt[1], 1, 3);
  ADReal z(pt[2], 2, 3);

  
  ADReal rtn = 1.0 + sin(m_ * x) * sin(m_ * y) * sin(m_ * z) + cos(m_ * x) * cos(m_ * y) * cos(m_ * z);

  return rtn;
}


