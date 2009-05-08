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

#include "SundanceQuadratureFamily.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace Teuchos;

int QuadratureFamily::getNumPoints( const CellType & cellType ) const
{
  const QuadratureFamilyBase* q 
    = dynamic_cast<const QuadratureFamilyBase*>(ptr().get());
  return q->getNumPoints( cellType );
}

void QuadratureFamily::getPoints(const CellType& cellType, 
                                 Array<Point>& quadPoints,
                                 Array<double>& quadWeights) const 
{
  const QuadratureFamilyBase* q 
    = dynamic_cast<const QuadratureFamilyBase*>(ptr().get());
  
  TEST_FOR_EXCEPTION(q==0, InternalError, 
                     "QuadratureFamilyStub pointer" << toXML().toString() 
                     << " could not be cast to a QuadratureFamilyBase ptr");
                     
  q->getPoints(cellType, quadPoints, quadWeights);
}

XMLObject QuadratureFamily::toXML() const 
{
  return ptr()->toXML();
}

int QuadratureFamily::order() const 
{
  return ptr()->order();
}



void QuadratureFamily::getFacetPoints(const CellType& cellType, 
                                      int facetDim,
                                      int facetIndex,
                                      Array<Point>& quadPoints,
                                      Array<double>& quadWeights) const
{
  switch(cellType)
    {
    case LineCell:
      getLineFacetQuad(facetDim, facetIndex, quadPoints, quadWeights);
      break;
    case TriangleCell:
      getTriangleFacetQuad(facetDim, facetIndex, quadPoints, quadWeights);
      break;
    case TetCell:
      getTetFacetQuad(facetDim, facetIndex, quadPoints, quadWeights);
      break;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "getFacetPoints() not implemented for cell type "
                         << cellType);
    }
}


void QuadratureFamily::getLineFacetQuad(int facetDim,
                                        int facetIndex,
                                        Array<Point>& quadPoints,
                                        Array<double>& quadWeights) const
{
  TEST_FOR_EXCEPTION(facetDim > 0, RuntimeError,
                     "Invalid facet dimension " << facetDim 
                     << " in getLineFacetQuad()");
  TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 1, RuntimeError,
                     "Invalid facet index " << facetIndex  
                     << " in getLineFacetQuad()");

  quadPoints.resize(1);
  quadWeights.resize(1);
  quadWeights[0] = 1.0;  
  quadPoints[0] = Point((double) facetIndex);
}




void QuadratureFamily::getTriangleFacetQuad(int facetDim,
                                            int facetIndex,
                                            Array<Point>& quadPoints,
                                            Array<double>& quadWeights) const
{
  TEST_FOR_EXCEPTION(facetDim > 1, RuntimeError,
                     "Invalid facet dimension " << facetDim 
                     << " in getTriangleFacetQuad()");
  TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 2, RuntimeError,
                     "Invalid facet index " << facetIndex  
                     << " in getTriangleFacetQuad()");
  if (facetDim==1)
    {
      Array<Point> facetPts;
      Array<double> facetWts;
      getPoints(LineCell, facetPts, facetWts);
      quadPoints.resize(facetPts.size());
      quadWeights.resize(facetWts.size());
      for (unsigned int i=0; i<facetPts.size(); i++)
        {
          if (facetIndex==0)
            {
              quadPoints[i] = Point(facetPts[i][0], 0.0);
              quadWeights[i] = facetWts[i];
            }
          else if (facetIndex==1)
            {
              quadPoints[i] = Point(1.0-facetPts[i][0], facetPts[i][0]);
              quadWeights[i] = facetWts[i];
            }
          else
            {
              quadPoints[i] = Point(0.0, 1.0-facetPts[i][0]);
              quadWeights[i] = facetWts[i];
            }
        }
    }
  else
    {
      quadPoints.resize(1);
      quadWeights.resize(1);
      quadWeights[0] = 1.0;  
      if (facetIndex==0) quadPoints[0] = Point(0.0, 0.0);
      else if (facetIndex==1) quadPoints[0] = Point(1.0, 0.0);
      else quadPoints[0] = Point(0.0, 1.0);
    }
}



void QuadratureFamily::getTetFacetQuad(int facetDim,
                                       int facetIndex,
                                       Array<Point>& quadPoints,
                                       Array<double>& quadWeights) const
{
  TEST_FOR_EXCEPTION(facetDim > 2, RuntimeError,
                     "Invalid facet dimension " << facetDim 
                     << " in getTetFacetQuad()");
  TEST_FOR_EXCEPTION(facetIndex < 0 || facetIndex > 4, RuntimeError,
                     "Invalid facet index " << facetIndex  
                     << " in getTetFacetQuad()");
  if (facetDim==2)
    {
      Array<Point> facetPts;
      Array<double> facetWts;
      getPoints(TriangleCell, facetPts, facetWts);
      quadPoints.resize(facetPts.size());
      quadWeights.resize(facetWts.size());
      for (unsigned int i=0; i<facetPts.size(); i++)
        {
          double s = facetPts[i][0];
          double t = facetPts[i][1];
          double x,y,z;
          if (facetIndex==0)
            {
              x = 1.0-s;
              y = 1.0-s-t;
              z = t;
            }
          else if (facetIndex==1)
            {
              x = 1.0-s;
              y = 0.0;
              z = t;
            }
          else if (facetIndex==2) 
            {
              x = 0.0;
              y = 1.0-s;
              z = t;
            }
          else
            {
              x = s;
              y = t;
              z = 0.0;
            }
          quadPoints[i] = Point(x, y, z);
          quadWeights[i] = facetWts[i];

        }
    }
  else if (facetDim==0)
    {
      quadPoints.resize(1);
      quadWeights.resize(1);
      quadWeights[0] = 1.0;  
      if (facetIndex==0) quadPoints[0] = Point(0.0, 0.0);
      else if (facetIndex==1) quadPoints[0] = Point(1.0, 0.0);
      else if (facetIndex==2) quadPoints[0] = Point(1.0, 0.0);
      else quadPoints[0] = Point(0.0, 1.0);
    }
}




