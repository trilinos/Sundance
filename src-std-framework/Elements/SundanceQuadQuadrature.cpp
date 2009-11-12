#include "SundanceQuadQuadrature.hpp"
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

void QuadQuadrature::getPoints(int order, Array<double>& wgt,
                                   Array<double>& x,
                                   Array<double>& y)
{
    if (!getSymmetricPoints(order, wgt, x, y))
    {
      getNonsymmetricPoints(order, wgt, x, y);
    }
}


bool QuadQuadrature::getSymmetricPoints(int order, Array<double>& wgt,
                                            Array<double>& x,
                                            Array<double>& y)
{
	int np;
	Array<double> w;
	Array<int> multiplicity;
	Array<Array<double> > q;

	if (order==1)
		{
			multiplicity = tuple(1);
			np = 1;
			w = tuple(1.0);
			q.resize(1);
			q[0] = tuple(1.0/2.0, 1.0/2.0 );
		}
	else if (order==2)
		{
			multiplicity = tuple(4);
			np = 4;
			w = tuple(1.0/4.0);
			q.resize(1);
			q[0] = tuple( 1.0/3.0 , 2.0/3.0 );
		}
    // TODO: implement higher order quadrature for quad (low priority)
	// but then higher orders can be done then with getNonsymmetricPoints() functions
	// this should work for higher orders as well
	else
		{
      return false;
		}

	for (int i=0; i<q.length(); i++)
		{
			Array<Array<double> > qPerm;
			permute(multiplicity[i], q[i], qPerm);
			for (int j=0; j<multiplicity[i]; j++)
				{
					x.append(qPerm[j][0]);
					y.append(qPerm[j][1]);
					wgt.append(w[i]);
				}
		}
	return true;
}


void QuadQuadrature::getNonsymmetricPoints(int order, Array<double>& wgt,
                                               Array<double>& x,
                                               Array<double>& y)
{
  int nNodes = (order+3)/2;
  Gauss1D rule(nNodes, 0.0, 1.0);
  Array<double> s = rule.nodes();
  Array<double> t = s;
  Array<double> w = rule.weights();
  int n = rule.nPoints();

  wgt.resize(n*n);
  x.resize(n*n);
  y.resize(n*n);

  int k=0;
  for (int i=0; i<n; i++)
    {
      double p = s[i];
      for (int j=0; j<n; j++, k++)
        {
          double q = t[j]; //similar to the p value, caz we have quad
          x[k] = p;
          y[k] = q;
          wgt[k] = w[i]*w[j]; //we don't need factor 0.5 caz this is a quad (not triangle)
        }
    }
}


void QuadQuadrature::permute(int m, const Array<double>& q,
																 Array<Array<double> >& qPerm)
{
	qPerm.resize(m);
	if (m==1)
		{
			qPerm[0] = q;
		}
	else if (m==4)
		{
			qPerm[0] = tuple(q[0], q[0]);
			qPerm[1] = tuple(q[1], q[0]);
			qPerm[2] = tuple(q[0], q[1]);
			qPerm[3] = tuple(q[1], q[1]);
		}
	// higher orders we can do with getNonsymmetricPoints() function
	// so further developement here is not necesary
	else
		{
#ifndef TRILINOS_7
			SUNDANCE_ERROR("invalid multiplicity " 
                     << m <<
                     " in TriangleQuadrature::permute()");
#else
			SUNDANCE_ERROR7("invalid multiplicity " 
                     << m <<
                     " in TriangleQuadrature::permute()");
#endif
		}
}

bool QuadQuadrature::test(int p)
{
	// TODO: implement this function (but has low priority)
	return true;
}

// TODO: this function is needed in the testing (but has low priority)
double QuadQuadrature::exact(int a, int b, int c)
{
	return fact(a)*fact(b)*fact(c)/fact(a+b+c+2);
}
double QuadQuadrature::fact(int x)
{
	if (x==0) return 1.0;
	return x*fact(x-1);
}


