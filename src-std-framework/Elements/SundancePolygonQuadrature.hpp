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

#ifndef SUNDANCE_POLYGONQUADRATURE_H
#define SUNDANCE_POLYGONQUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceQuadratureFamily.hpp"

namespace Sundance
{
using namespace Teuchos;


/** 
 * Quadrature special for the polygon curve integration <br>.
 * There we need for each line segment inside one cell a line integral.
 * So for one cell we might need several line quadrature points, and this class
 * gives us per cell(line) a constant number of quadrature points for a line,
 * so we can integrate along the curve exactly. (each line segment exactly inside one cell) <br>
 * As an input we need a quadrature class, and the points of that quadrature rule for a line
 * will be used for each line segment of the polygon inside one cell (2D).
 */
class PolygonQuadrature : public QuadratureFamilyBase
{
public:
  /** */
  PolygonQuadrature( const QuadratureFamily& quad );

  /** */
  virtual ~PolygonQuadrature(){;}


  /** */
  virtual XMLObject toXML() const ;

  /** Describable interface */
  virtual std::string description() const 
    {return "PolygonQuadrature[order=" + Teuchos::toString(order())
        +  "]";}

  /* handleable boilerplate */
  GET_RCP(QuadratureFamilyStub);

  /** method for the user to change this value
   * too much line nr is an overhead
   * too less will cause error thrown in the code
   * @param maxNrLine */
  static void setNrMaxLinePerCell( int maxNrLine ) { nrMaxLinePerCell_ = maxNrLine; }

  /** return the maximal number of line segments inside one cell */
  static int getNrMaxLinePerCell() { return nrMaxLinePerCell_; }

protected:

  /** for polygon curve integrals only this method should be used */
  virtual void getLineRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

private:

  /** The quadrature which will be used for the polygon lines integration */
  const QuadratureFamily& quad_;

  /** maximum number of lines inside a cell, the default value is 6 but this can be changed
   * by the user */
  static int nrMaxLinePerCell_;

};
}


#endif
