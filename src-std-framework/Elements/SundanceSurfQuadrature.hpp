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
 * SundanceSurfQuadrature.hpp
 *
 *  Created on: Oct 24, 2011
 *      Author: benk
 */

#ifndef SUNDANCE_SURFQUADRATURE_H
#define SUNDANCE_SURFQUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceQuadratureFamily.hpp"

namespace Sundance
{
using namespace Teuchos;


/**
 * The surface integral class. The surface integral is build from triangles.
 * The maximum number of triangles is 4. <br>
 * IMPORTANT: this quadrature class should only be used for Surface Integrals in 3D with Brick cells
 */
class SurfQuadrature : public QuadratureFamilyBase
{
public:
  /** */
  SurfQuadrature( const QuadratureFamily& quad );

  /** */
  virtual ~SurfQuadrature(){;}


  /** */
  virtual XMLObject toXML() const ;

  /** Describable interface */
  virtual std::string description() const
    {return "SurfQuadrature[order=" + Teuchos::toString(order())
        +  "]";}

  /* handleable boilerplate */
  GET_RCP(QuadratureFamilyStub);

  /** return the maximal number of line segments inside one cell */
  static int getNrMaxTrianglePerCell() { return 4; }

protected:

  /** for surface integral integrals only this method should be used */
  virtual void getQuadRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;


  /** for surface integral integrals only this method should be used */
  virtual void getTriangleRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

private:

  /** The quadrature which will be used for the surface integration */
  const QuadratureFamily& quad_;

};
}


#endif
