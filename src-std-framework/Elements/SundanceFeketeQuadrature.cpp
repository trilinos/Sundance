/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceFeketeBrickQuadrature.hpp"
#include "SundanceFeketeQuadQuadrature.hpp"
#include "SundanceFeketeQuadrature.hpp"
#include "SundanceFeketeTriangleQuadrature.hpp"
#include "SundanceGaussLobatto1D.hpp"

using namespace Sundance;
using namespace Teuchos;

FeketeQuadrature::FeketeQuadrature(int order) :
	QuadratureFamilyBase(order)
{

}

XMLObject FeketeQuadrature::toXML() const
{
	XMLObject rtn("FeketeQuadrature");
	rtn.addAttribute("order", Teuchos::toString(order()));
	return rtn;
}

void FeketeQuadrature::getLineRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	int p = order() + 3;
	p = p + (p % 2);
	int n = p / 2;

	quadPoints.resize(n);
	quadWeights.resize(n);

	GaussLobatto1D q1(n, 0.0, 1.0);

	for (int i = 0; i < n; i++)
	{
		quadWeights[i] = q1.weights()[i];
		quadPoints[i] = Point(q1.nodes()[i]);
	}
}

void FeketeQuadrature::getTriangleRule(Array<Point>& quadPoints, Array<
		double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> w;

	FeketeTriangleQuadrature::getPoints(order(), w, x, y);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = 0.5 * w[i];
		quadPoints[i] = Point(x[i], y[i]);
	}
}

void FeketeQuadrature::getQuadRule(Array<Point>& quadPoints,
		Array<double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> w;

	FeketeQuadQuadrature::getPoints(order(), w, x, y);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = w[i];
		quadPoints[i] = Point(x[i], y[i]);
	}
}

void FeketeQuadrature::getBrickRule(Array<Point>& quadPoints, Array<
		double>& quadWeights) const
{
	Array<double> x;
	Array<double> y;
	Array<double> z;
	Array<double> w;

	FeketeBrickQuadrature::getPoints(order(), w, x, y, z);
	quadPoints.resize(w.length());
	quadWeights.resize(w.length());
	for (int i = 0; i < w.length(); i++)
	{
		quadWeights[i] = w[i];
		quadPoints[i] = Point(x[i], y[i], z[i]);
	}
}
/*
void FeketeQuadrature::getPoints(const CellType& cellType, int cellDim,
		int celLID, int facetIndex, const Mesh& mesh,
		const ParametrizedCurve& globalCurve, Array<Point>& quadPoints, Array<
				double>& quadWeights, bool &isCut) const
{

	//	Array<int> celLIDs(1);
	//	celLIDs[0] = celLID;
	//Array<Point> tmpPoints = quadPoints;
	//mesh.pushForward( cellDim, celLIDs, quadPoints, tmpPoints );

	isCut = false;

	Array<Point> vertices;
	Array<int> vertexBool;

	Array<int> nodeLIDs;
	Array<int> orient;

	// Get coordinates
	mesh.getFacetArray(cellDim, celLID, 0, nodeLIDs, orient);

	// Make a list of inner and outer nodes
	vertices.resize(nodeLIDs.size());
	vertexBool.resize(nodeLIDs.size());
	for (int i; i < vertexBool.size(); i++)
		vertexBool[i] = 0;
	int nrInnerNodes = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		vertices[i] = mesh.nodePosition(nodeLIDs[i]);
		if (globalCurve.curveEquation(vertices[i]) < 0)
		{
			nrInnerNodes++;
			vertexBool[i] = 1;
		}
	}

	// If there are inner and outer nodes
	if (nrInnerNodes > 0 && nrInnerNodes < vertices.size())
	{
		// ToDo: Unterscheidung nach alpha1/alpha2 ob innen oder aussen (oder beides?) integrieren.

		isCut = true;

		// ToDo: Other cells!
		if (cellType == TriangleCell)
		{
			double s1 = vertices[0].distance(vertices[1]);
			double s2 = vertices[0].distance(vertices[2]);
			double s3 = vertices[1].distance(vertices[2]);

			// Make [ab] longest side
			Point& pa;
			Point& pb;
			Point& pc;
			if (s1 >= s2 && s1 >= s3)
			{
				pa = vertices[0];
				pb = vertices[1];
				pc = vertices[2];
			}
			else if (s2 >= s1 && s2 >= s3)
			{
				pa = vertices[2];
				pb = vertices[0];
				pc = vertices[1];
			}
			else
			{
				pa = vertices[1];
				pb = vertices[2];
				pc = vertices[0];
			}

			if (globalCurve.curveEquation(pa) < 0 && globalCurve.curveEquation(
					pb) < 0)
			{
				// Integration von a bis b
				Array<Point> tmp;
				globalCurve.returnIntersect(pa, pb, j, tmp);
				if (j != 1)
					SUNDANCE_MSG2(2, "Verfeinern!");
				Point intersect_ab = tmp[0];
				globalCurve.returnIntersect(pb, pc, j, tmp);
				if (j != 1)
					SUNDANCE_MSG2(2, "Verfeinern!");
				Point intersect_bc = tmp[0];

			}
			else if (globalCurve.curveEquation(pa) < 0
					&& globalCurve.curveEquation(pc) < 0)
			{
				// Integration von a bis proj_bc (event. minus eck)
			}
			else if (globalCurve.curveEquation(pa) < 0)
			{
				// Integration von a bis proj_ac (event. plus eck)
			}
			else if (globalCurve.curveEquation(pc) < 0
					&& globalCurve.curveEquation(pb) < 0)
			{
				// Integration von proj_ac (+- eck) bis b
			}
			else if (globalCurve.curveEquation(pb) < 0)
			{
				// Integration von proj_bc (+-eck) bis b
			}
			else if (globalCurve.curveEquation(pc) < 0)
			{
				// Gegenteil vom ersten
			}
			else
			{
				// GIBTS NICHT! TEST!
			}

			//Array<double> invJ;
			 //invJ.resize(4);
			 //invJ[0] = pc[1] - pa[1];
			 //invJ[1] = pa[0] - pc[0];
			 //invJ[2] = pa[1] - pb[1];
			 //invJ[3] = pb[0] - pa[0];
			 //double detJ = invJ[3] * invJ[0] - invJ[1] * invJ[2];
			 //for (int i = 0; i < 4; i++)
			 //invJ[i] /= detJ;

			 // hole integrationspunkte ueber dreieck
			 //Array<Point> localquadPts;
			 //Array<double> localquadWeights;
			 //QuadratureFamily quad = new GaussianQuadrature(order);
			 //quad.getPoints(TriangleCell, localquadPts, localquadWeights);


			// teile zelle in 2 dreiecke und pseudo-trapez (mit einer kurvigen seite)
			// Point proj_ac = ((pc - pa) * (pb - pa)) * (pb - pa) / ab;
			// Point proj_bc = ((pc - pb) * (pa - pb)) * (pa - pb) / ab;

		}

	}

	// SUNDANCE_MSG2(2, "Treffer!");


	// Either complete inner or outer cell
	else
	{
		double alpha = globalCurve.integrationParameter(vertices[0]);
		if (alpha != 1.0)
		{
			isCut = true;
			for (int i = 0; i < quadWeights.size(); i++)
			{
				quadWeights[i] *= alpha;
			}
		}
	}
}


double FeketeQuadrature::quadBasisFunc2D(const Point& xstart, const Point& xstop, const ParametrizedCurve& Curve) {

}


double FeketeQuadrature::evaluateBasisAtRefPoint(int basis_nr, double const x,
		double const y)
{
	if (!hasInvV_)
		computeInvV();

	Array<double> polvals;
	polvals.resize(numQuad_);
	evaluatePKDpolynoms(order_, &polvals[0], x, y);

	double result = 0.0;
	for (int n = 0; n < numQuad_; n++)
	{
		result += invV_[basis_nr + n * numQuad_] * polvals[n];
	}

	return result;
}


void FeketeQuadrature::computeInvV()
{
	TimeMonitor timer(VInversionTimer());
	if (hasInvV_)
		return;

	Tabs tabs;
	SUNDANCE_OUT(this->verb() > 2, tabs << "calculating inverted generalized Vandermonde");

	// Get Fekete points
	Array<Point> quadPts;
	Array<double> quadWeights;
	QuadratureFamily quad = new FeketeQuadrature(order_);
	quad.getPoints(cellType_, quadPts, quadWeights);

	invV_.resize(numQuad_ * numQuad_);

	// Build Vandermonde matrix of PKD basis at Fekete points
	for (int q = 0; q < numQuad_; q++)
	{
		double* start = &(invV_[q * numQuad_]);
		evaluatePKDpolynoms(order_, start, quadPts[q][0], quadPts[q][1]);
	}

	double* invVPtr = &(invV_[0]);

	Array<int> iPiv_;
	iPiv_.resize(numQuad_);
	int* iPivPtr = &(iPiv_[0]);

	// LAPACK error flag
	int info = 0;

	// LU factorization
	::dgetrf_(&numQuad_, &numQuad_, invVPtr, &numQuad_, iPivPtr, &info);

	TEST_FOR_EXCEPTION(info != 0, RuntimeError,
			"CellWeightsBatch::computeInvV(): factorization of V failed");

	// Inversion of the matrix
	Array<double> work;
	work.resize(1);

	int lwork = -1;
	::dgetri_(&numQuad_, invVPtr, &numQuad_, iPivPtr, &(work[0]), &lwork, &info);
	lwork = (int) work[0];
	work.resize(lwork);
	::dgetri_(&numQuad_, invVPtr, &numQuad_, iPivPtr, &(work[0]), &lwork, &info);

	TEST_FOR_EXCEPTION(info != 0, RuntimeError,
			"CellWeightsBatch::computeInvV(): inversion of V failed");

	// todo: addFlops();
	hasInvV_ = true;
}
*/
/**
 * Evaluates all basis functions of a Proriol-Koornwinder-Dubiner basis
 * of order 'order' at (x,y)
 *
 */
/*void FeketeQuadrature::evaluatePKDpolynoms(int order, double* resultPtr,
		double x, double y) const
{
	int offset = 0;

	double xy = x + y;

	double Legendre = 1.0;
	double Legendre1 = 0.0;
	for (int k = 0; k <= order; k++)
	{
		double Jacobi = 1.0;
		double Jacobi1 = 0.0;
		for (int l = 0; l <= order - k; l++)
		{
			resultPtr[offset] = Legendre * pow(xy, k) * Jacobi;

			// Update counter (implements a bijection:
			// Polynom index (k,l) -> Basis index i)
			++offset;

			// Update Jacobi polynomial
			double Jacobi2 = Jacobi1;
			Jacobi1 = Jacobi;
			int twok = 2 * k;
			int twoltwok = 2 * l + twok;
			int c1 = 2 * (l + 1) * (l + twok + 2) * (twoltwok + 1);
			int c2 = (twoltwok + 2) * (twok + 1) * (twok + 1);
			int c3 = (twoltwok + 1) * (twoltwok + 2) * (twoltwok + 3);
			int c4 = 2 * (l + twok + 1) * l * (twoltwok + 3);
			Jacobi = ((c2 + c3 * (1.0 - 2.0 * x - 2.0 * y)) * Jacobi1 - c4
					* Jacobi2) / c1;
		}

		// Update Legendre polynomial
		double Legendre2 = Legendre1;
		Legendre1 = Legendre;
		if (::fabs(xy) > 0.0)
		{
			xy = (2.0 * k + 1.0) * ((x - y) / xy) * Legendre1;
		}
		Legendre = (xy - k * Legendre2) / (k + 1);
	}
}*/
