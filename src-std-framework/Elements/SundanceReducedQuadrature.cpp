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

#include "SundanceReducedQuadrature.hpp"
#include "SundanceGaussianQuadrature.hpp"

namespace Sundance
{
using namespace Teuchos;


ReducedQuadrature::ReducedQuadrature()
  : QuadratureFamilyBase(1)
{}

XMLObject ReducedQuadrature::toXML() const 
{
  XMLObject rtn("ReducedQuadrature");
  return rtn;
}



void ReducedQuadrature::getPoints(const CellType& cellType, 
  Array<Point>& quadPoints,
  Array<double>& quadWeights) const 
{
  GaussianQuadrature basic(1);
  basic.getPoints(cellType, quadPoints, quadWeights);
}

int  ReducedQuadrature::getNumPoints( const CellType &cellType ) const 
{
  return 1;
}

void ReducedQuadrature::getAdaptedWeights(const CellType& cellType ,
  int cellDim,
  int celLID ,
  int facetIndex ,
  const Mesh& mesh ,
  const ParametrizedCurve& globalCurve ,
  Array<Point>& quadPoints ,
  Array<double>& quadWeights ,
  bool &isCut) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

}
