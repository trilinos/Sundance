/* @HEADER@ */
/* @HEADER@ */

#include "SundanceClosedNewtonCotes.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

ClosedNewtonCotes::ClosedNewtonCotes()
  : QuadratureFamilyBase(2)
{;}

XMLObject ClosedNewtonCotes::toXML() const 
{
  XMLObject rtn("ClosedNewtonCotes");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}

void ClosedNewtonCotes::getLineRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const 
{
  Array<double> x = tuple(0.0, 1.0/3.0, 2.0/3.0, 1.0);
	quadWeights = tuple(1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0);
  quadPoints.resize(x.size());

	for (int i=0; i<x.length(); i++)
		{
			quadPoints[i] = Point(x[i]);
		}
}

void ClosedNewtonCotes::getTriangleRule(Array<Point>& quadPoints,
                                          Array<double>& quadWeights) const 
{
  quadPoints = tuple(Point(0.0, 0.0), Point(0.5, 0.0), Point(1.0, 0.0),
                     Point(0.5, 0.5), Point(0.0, 1.0), Point(0.0, 0.5),
                     Point(1.0/3.0, 1.0/3.0));
	quadWeights = tuple(3.0/60.0, 8.0/60.0, 3.0/60.0,
                      8.0/60.0, 3.0/60.0, 8.0/60.0,
                      27.0/60.0);

	for (int i=0; i<quadWeights.size(); i++)
		{
			quadWeights[i] *=  0.5;
		}
}

