#include "SundanceFeketeTriangleQuadrature.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;
using namespace Teuchos;

/**
 *
 *  Reference:
 *  Mark Taylor, Beth Wingate, Rachel Vincent,
 *    An Algorithm for Computing Fekete Points in the Triangle,
 *    SIAM Journal on Numerical Analysis,
 *    Volume 38, Number 5, 2000, pages 1707-1720.
 *
 *	L. P. Bos,
 *	  Bounding the Lebesgue function for Lagrange interpolation in a simplex,
 *	  Journal of Approximation Theory,
 *	  Volume 38, Issue 1, 1983, pages 43-59.
 *
 *	Caveats:
 *	 - Weights sum up to 2 in the tables below (as given in Taylor et al.)
 *	   We take the halves at the end (and another half in calling routine)
 *	 - Some coordinates of 'orbit=3'-points had to be permuted compared to the
 *	   reference so that the correct three permutations will be done in the
 *	   permute() method
 *
 */
void FeketeTriangleQuadrature::getPoints(int order, Array<double>& wgt, Array<
		double>& x, Array<double>& y)
{
	Array<double> w;
	Array<int> multiplicity;
	Array<Array<double> > q;

	if (order == 1)
	{
		// Delete it? One gets 2nd order for the same price...
		multiplicity = tuple(3);
		q.resize(1);
		w = tuple(0.6666666666);
		q[0] = tuple(1.0000000000, 0.0000000000, 0.0000000000);
	}
	else if (order == 2)
	{
		// Corners have weight == 0 -> not listed
		multiplicity = tuple(3);
		q.resize(1);
		w = tuple(0.6666666666);
		q[0] = tuple(0.0000000000, 0.5000000000, 0.5000000000);
	}
	else if (order == 3)
	{
		multiplicity = tuple(1, 3, 6);
		q.resize(3);
		w = tuple(0.9000000000, 0.0333333333, 0.1666666667);
		q[0] = tuple(0.3333333333, 0.3333333333, 0.3333333334);
		q[1] = tuple(1.0000000000, 0.0000000000, 0.0000000000);
		q[2] = tuple(0.0000000000, 0.2763932023, 0.7236067977);
	}
	else if (order == 6)
	{
		multiplicity = tuple(1, 3, 3, 3, 6, 6, 6);
		q.resize(7);
		w = tuple(0.2178563571, 0.1104193374, 0.0358939762, 0.0004021278,
				0.1771348660, 0.0272344079, 0.0192969460);
		q[0] = tuple(0.3333333333, 0.3333333333, 0.3333333334);
		q[1] = tuple(0.7873290632, 0.1063354684, 0.1063354684);
		q[2] = tuple(0.0000000000, 0.5000000000, 0.5000000000);
		q[3] = tuple(1.0000000000, 0.0000000000, 0.0000000000);
		q[4] = tuple(0.1171809171, 0.3162697959, 0.5665492870);
		q[5] = tuple(0.0000000000, 0.2655651402, 0.7344348598);
		q[6] = tuple(0.0000000000, 0.0848854223, 0.9151145777);
	}
	else if (order == 9)
	{
		multiplicity = tuple(1, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6);
		q.resize(12);
		w = tuple(0.1096011288, 0.0767491008, 0.0646677819, 0.0276211659,
				0.0013925011, 0.0933486453, 0.0619010169, 0.0437466450,
				0.0114553907, 0.0093115568, 0.0078421987, 0.0022457501);
		q[0] = tuple(0.3333333333, 0.3333333333, 0.3333333334);
		q[1] = tuple(0.6591363598, 0.1704318201, 0.1704318201);
		q[2] = tuple(0.0600824712, 0.4699587644, 0.4699587644);
		q[3] = tuple(0.9021308608, 0.0489345696, 0.0489345696);
		q[4] = tuple(1.0000000000, 0.0000000000, 0.0000000000);
		q[5] = tuple(0.1784337588, 0.3252434900, 0.4963227512);
		q[6] = tuple(0.0588564879, 0.3010242110, 0.6401193011);
		q[7] = tuple(0.0551758079, 0.1543901944, 0.7904339977);
		q[8] = tuple(0.0000000000, 0.4173602935, 0.5826397065);
		q[9] = tuple(0.0000000000, 0.2610371960, 0.7389628040);
		q[10] = tuple(0.0000000000, 0.1306129092, 0.8693870908);
		q[11] = tuple(0.0000000000, 0.0402330070, 0.9597669930);
	}
	else
	{
#ifndef TRILINOS_7
		SUNDANCE_ERROR("symmetric Fekete quadrature rule order "
				<< order <<
				" for triangles not available");
#else
		SUNDANCE_ERROR7("symmetric Fekete quadrature rule order "
				<< order <<
				" for triangles not available");
#endif
	}

	for (int i = 0; i < q.length(); i++)
	{
		Array<Array<double> > qPerm;
		permute(multiplicity[i], q[i], qPerm);
		for (int j = 0; j < multiplicity[i]; j++)
		{
			x.append(qPerm[j][0]);
			y.append(qPerm[j][1]);
			wgt.append(0.5 * w[i]);
		}
	}

}

bool FeketeTriangleQuadrature::supportsOrder(int order)
{
	if (order == 1 || order == 2 || order == 3 || order == 6 || order == 9)
		return true;
	return false;
}

void FeketeTriangleQuadrature::permute(int m, const Array<double>& q, Array<
		Array<double> >& qPerm)
{
	qPerm.resize(m);
	if (m == 1)
	{
		qPerm[0] = q;
	}
	else if (m == 3)
	{
		qPerm[0] = tuple(q[0], q[1], q[2]);
		qPerm[1] = tuple(q[1], q[0], q[2]);
		qPerm[2] = tuple(q[2], q[1], q[0]);
	}
	else if (m == 6)
	{
		qPerm[0] = tuple(q[0], q[1], q[2]);
		qPerm[1] = tuple(q[0], q[2], q[1]);
		qPerm[2] = tuple(q[1], q[0], q[2]);
		qPerm[3] = tuple(q[1], q[2], q[0]);
		qPerm[4] = tuple(q[2], q[1], q[0]);
		qPerm[5] = tuple(q[2], q[0], q[1]);
	}
	else
	{
#ifndef TRILINOS_7
		SUNDANCE_ERROR("invalid multiplicity "
				<< m <<
				" in FeketeTriangleQuadrature::permute()");
#else
		SUNDANCE_ERROR7("invalid multiplicity "
				<< m <<
				" in FeketeTriangleQuadrature::permute()");
#endif
	}
}

bool FeketeTriangleQuadrature::test(int p)
{
	Array<double> w;
	Array<double> x;
	Array<double> y;

	getPoints(p, w, x, y);
	bool pass = true;

	for (int a = 0; a <= p; a++)
	{
		for (int b = 0; b < p - a; b++)
		{
			int cMax = p - a - b;
			for (int c = 0; c <= cMax; c++)
			{
				double sum = 0.0;
				for (int q = 0; q < w.length(); q++)
				{
					sum += 0.5 * w[q] * pow(x[q], (double) a) * pow(y[q],
							(double) b) * pow(1.0 - x[q] - y[q], (double) c);
				}
				double err = fabs(sum - exact(a, b, c));
				bool localPass = err < 1.0e-14;
				pass = pass && localPass;
				if (!localPass)
				{
					fprintf(
							stderr,
							"order=%d m (%d, %d, %d) q=%22.15g exact=%22.15g\n",
							p, a, b, c, sum, exact(a, b, c));
					std::cerr << "error = " << err << std::endl;
				}
			}
		}
	}
	return pass;
}

double FeketeTriangleQuadrature::exact(int a, int b, int c)
{
	return fact(a) * fact(b) * fact(c) / fact(a + b + c + 2);
}

double FeketeTriangleQuadrature::fact(int x)
{
	if (x == 0)
		return 1.0;
	return x * fact(x - 1);
}

