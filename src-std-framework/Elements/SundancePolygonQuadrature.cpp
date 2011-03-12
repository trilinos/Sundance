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

#include "SundancePolygonQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;

int PolygonQuadrature::nrMaxLinePerCell_ = 5;

PolygonQuadrature::PolygonQuadrature( const QuadratureFamily& quad )
  : QuadratureFamilyBase(quad.order()) , quad_(quad)
{
  
}

XMLObject PolygonQuadrature::toXML() const
{
  XMLObject rtn("PolygonQuadrature");
  rtn.addAttribute("order", Teuchos::toString(order()));
  return rtn;
}



void PolygonQuadrature::getLineRule(Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const 
{
	Array<Point> quadPoints_tmp;
	Array<double> quadWeights_tmp;
	quad_.getPoints( LineCell , quadPoints_tmp , quadWeights_tmp );

	// the nr. of points per line segments
	int nrPointPerLine = quadPoints_tmp.size() , ind = 0;

	// resize the point arrays and the weight arrays
	quadPoints.resize( nrMaxLinePerCell_*nrPointPerLine );
	quadWeights.resize( nrMaxLinePerCell_*nrPointPerLine );

	// each line segment
	for (int nrl = 0 ; nrl < nrMaxLinePerCell_ ; nrl++ ){
		// loop over each quadrature point
		for (int q = 0 ; q < nrPointPerLine ; q++ ){
			// copy the points and the weights several times
			quadPoints[ind] = quadPoints_tmp[q];
			quadWeights[ind] = quadWeights_tmp[q]/((double)nrMaxLinePerCell_);
			ind++;
		}
	}
}
