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

#ifndef SUNDANCE_QUADRATUREFAMILY_H
#define SUNDANCE_QUADRATUREFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdFwk
{
  using namespace TSFExtended;
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** 
   * QuadratureFamily is a geometry-independent specification of
   * a method by which quadrature is to be carried out. For example,
   * a GaussianQuadrature family will generate Gaussian
   * quadrature points on any cell type.
   */
  class QuadratureFamily : public Handle<QuadratureFamilyStub>
  {
  public:
    /* */
    HANDLE_CTORS(QuadratureFamily, QuadratureFamilyStub);
    /** */
    XMLObject toXML() const ;

    /** */
    int order() const ;

    /** Get the quadrature points and weights for the given cell type */
    void getPoints(const CellType& cellType, 
                   Array<Point>& quadPoints,
                   Array<double>& quadWeights) const ;

    /** Get quadrature points and weights for integration on a facet of a cell */
    void getFacetPoints(const CellType& cellType, 
                        int facetDim,
                        int facetIndex,
                        Array<Point>& quadPoints,
                        Array<double>& quadWeights) const ;
  private:
    /** Get quad points for a facet of a line */
    void getLineFacetQuad(int facetDim,
                          int facetIndex,
                          Array<Point>& quadPoints,
                          Array<double>& quadWeights) const ;
    /** Get quad points for a facet of a triangle */
    void getTriangleFacetQuad(int facetDim,
                              int facetIndex,
                              Array<Point>& quadPoints,
                              Array<double>& quadWeights) const ;
  };
}

#endif
