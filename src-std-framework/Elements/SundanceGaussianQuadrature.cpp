/* @HEADER@ */
/* @HEADER@ */

#include "SundanceGaussianQuadrature.hpp"
#include "SundanceGauss1D.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

GaussianQuadrature::GaussianQuadrature(int order)
  : QuadratureFamilyBase(order)
{;}

XMLObject GaussianQuadrature::toXML() const 
{
  XMLObject rtn("GaussianQuadrature");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}

void GaussianQuadrature::getLineRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const 
{
  int p = order() + 1;
  p = p + (p%2);
  int n = p/2;
			
  quadPoints.resize(n);
  quadWeights.resize(n);
			
  Gauss1D q1(n, 0.0, 1.0);
			
  for (int i=0; i<n; i++)
    {
      quadWeights[i] = q1.weights()[i];
      quadPoints[i] = Point(q1.nodes()[i]);
    }
}

void GaussianQuadrature::getTriangleRule(Array<Point>& quadPoints,
                                          Array<double>& quadWeights) const 
{
  Array<double> x;
  Array<double> y;
  Array<double> w;
			
  TriangleQuadrature::getPoints(order(), w, x, y);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = 0.5*w[i];
      quadPoints[i] = Point(x[i], y[i]);
    }  
}

void GaussianQuadrature::getTetRule(Array<Point>& quadPoints,
                                    Array<double>& quadWeights) const 
{
  Array<double> x;
  Array<double> y;
  Array<double> z;
  Array<double> w;
			
  TetQuadrature::getPoints(order(), w, x, y, z);
  quadPoints.resize(w.length());
  quadWeights.resize(w.length());
  for (int i=0; i<w.length(); i++)
    {
      quadWeights[i] = 0.5*w[i];
      quadPoints[i] = Point(x[i], y[i], z[i]);
    }  
}

