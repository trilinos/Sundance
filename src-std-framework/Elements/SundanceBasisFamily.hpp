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

#ifndef SUNDANCE_BASISFAMILY_H
#define SUNDANCE_BASISFAMILY_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceBasisFamilyBase.hpp"
#include "SundanceOrderedHandle.hpp"
#include "Teuchos_Array.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceStdFwk::Internal;
  using namespace TSFExtended;

  /** 
   * BasisFamily is the user-level handle class for specifying the
   * basis with which a test, unknown, or discrete function is represented.
   * Basis functions can be vector-valued, as is the case with, for example,
   * the Nedelec basis in electromagnetics; the dim() method
   * returns the spatial dimension of the basis functions. Scalar-valued
   * bases naturally have dim()=1. 
   */
  class BasisFamily : public OrderedHandle<BasisFamilyBase>
    {
    public:
      /* handle ctor boilerplate */
      ORDERED_HANDLE_CTORS(BasisFamily, BasisFamilyBase);

      /** write to XML */
      XMLObject toXML() const ;

      /** return the polynomial order to which the basis is complete */
      int order() const ;

      /** return the spatial dimension of the basis elements */
      int dim() const ;

      /** return the number of nodes for this basis on the given cell type */
      int nNodes(const CellType& cellType) const ;

     //  /** */
//       void getNodalPoints(const CellType& cellType,
//                           Array<Point>& x) const ;

      /** */
      bool operator==(const BasisFamily& other) const ;

      /** Sum up the dim() values for array of bases. */
      static int size(const Array<BasisFamily>& b) ;

      /** Extract the basis from an expression */
      static BasisFamily getBasis(const Expr& expr);
    };

  /** \relates BasisFamily */
  inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b)
  {
    return tuple(a,b);
  }

  /** \relates BasisFamily */
  inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
                                 const BasisFamily& c)
  {
    return tuple(a,b,c);
  }
}

#endif
