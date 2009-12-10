#include "SundanceBrickQuadrature.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceGauss1D.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace Teuchos;

void BrickQuadrature::getPoints(int order, Array<double>& wgt,
                                               Array<double>& x,
                                               Array<double>& y,
                                               Array<double>& z)
{
  int nNodes = (order+3)/2;
  Gauss1D rule(nNodes, 0.0, 1.0);
  Array<double> d1 = rule.nodes();
  Array<double> d2 = d1;
  Array<double> d3 = d1;
  Array<double> w = rule.weights();
  int n = rule.nPoints();

  wgt.resize(n*n*n);
  x.resize(n*n*n);
  y.resize(n*n*n);
  z.resize(n*n*n);

  int k=0;
  for (int i=0; i<n; i++)
    {
      double p = d1[i];
      for (int j=0; j<n; j++)
        {
          double q = d2[j]; //similar to the p value, caz we have quad
          for (int l=0; l<n; l++, k++){
        	 double r = d3[l];
             x[k] = p;
             y[k] = q;
             z[k] = r;
             wgt[k] = w[i]*w[j]*w[l];
          }
        }
    }
}

bool BrickQuadrature::test(int p)
{
	// TODO: implement this function (but has low priority)
	return true;
}

// TODO: this function is needed in the testing (but has low priority)
double BrickQuadrature::exact(int a, int b, int c)
{
	return fact(a)*fact(b)*fact(c)/fact(a+b+c+2);
}
double BrickQuadrature::fact(int x)
{
	if (x==0) return 1.0;
	return x*fact(x-1);
}


