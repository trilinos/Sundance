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
/*
 * SundanceSurfQuadrature.cpp
 *
 *  Created on: Oct 24, 2011
 *      Author: benk
 */

#include "SundanceSurfQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;


SurfQuadrature::SurfQuadrature( const QuadratureFamily& quad )
  : QuadratureFamilyBase(quad.order()) , quad_(quad)
{

}

XMLObject SurfQuadrature::toXML() const
{
  XMLObject rtn("SurfQuadrature");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}

void SurfQuadrature::getQuadRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const
{
	// IMPORTANT: this quadrature class should only be used for Surface Integrals in 3D with Brick cells
	getTriangleRule( quadPoints, quadWeights);
}


void SurfQuadrature::getTriangleRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const
{
	Array<Point> quadPoints_tmp;
	Array<double> quadWeights_tmp;
	quad_.getPoints( TriangleCell , quadPoints_tmp , quadWeights_tmp );

	// the nr. of points per line segments
	int nrPointPerLine = quadPoints_tmp.size() , ind = 0;

	// resize the point arrays and the weight arrays
	quadPoints.resize( 4 * nrPointPerLine );
	quadWeights.resize( 4 * nrPointPerLine );

	// each line segment
	for (int nrl = 0 ; nrl < 4 ; nrl++ ){
		// loop over each quadrature point
		for (int q = 0 ; q < nrPointPerLine ; q++ ){
			// copy the points and the weights several times
			quadPoints[ind] = quadPoints_tmp[q];
			quadWeights[ind] = quadWeights_tmp[q]/((double)4.0);
			ind++;
		}
	}
}
