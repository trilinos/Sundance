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
#include "SundanceBasisFamilyBase.hpp"
#include "SundanceOrderedHandle.hpp"
#include "Teuchos_Array.hpp"

namespace SundanceCore {class CommonFuncDataStub;}

namespace SundanceStdFwk
{
using Teuchos::XMLObject;
using Teuchos::tuple;
using SundanceUtils::OrderedHandle;
using SundanceCore::CommonFuncDataStub;
using TSFExtended::Handle;

/** 
 * BasisFamily is the user-level handle class for specifying the basis with
 * which a test, unknown, or discrete function is represented.  Basis
 * functions can be vector-valued, as is the case with, for example, the
 * Nedelec basis in electromagnetics; the dim() method returns the spatial
 * dimension of the basis functions. Scalar-valued bases naturally have
 * dim()=1.

*/
class BasisFamily : public OrderedHandle<BasisFamilyBase>
{
public:
  /* handle ctor boilerplate */
  ORDERED_HANDLE_CTORS(BasisFamily, BasisFamilyBase);

  /** write to XML */
  XMLObject toXML() const ;

  /** 
   * \brief Return the polynomial order of the basis functions, for use in
   * determining the quadrature rule required to integrate a product of
   * basis functions exactly. The polynomial order will be the smallest
   * integer for which all mixed partial derivatives vanish exactly.
   *  
   * Note: in H(div) and H(curl) spaces the order of accuracy is not
   * always an integer, and the relationship between the order of accuracy
   * and the return value of the order() method is not necessarily simple
   * (for instance, it can depend on things such as the convexity of the
   * boundary). Thus it is better to think of this method as specifying
   * the required order of quadrature, and not the order of accuracy of
   * approximation nor the order to which the space is complete.
   */
  int order() const ;

  /** \brief
   * Return the number of DOFs for this basis on the given 
   * reference cell type.
   */
  int nReferenceDOFs(const CellType& maximalCellType,
    const CellType& cellType) const ;

  /** 
   * \brief Return the dimension of the members of 
   * a vector-valued basis. Return 1 if the basis
   * is scalar-valued. Otherwise, return the spatial dimension.
   */
  int dim() const ;

  /** 
   * \brief Return the tensor order of the members of 
   * a basis.
   */
  int tensorOrder() const {return ptr()->tensorOrder();}

  /** */
  bool operator==(const BasisFamily& other) const ;

  

  /** \brief Inform caller as to whether I am a scalar basis. */
  bool isScalarBasis() const {return ptr()->isScalarBasis();}

  /** \brief Inform caller as to whether I am an H(div) basis. */
  bool isHDivBasis() const {return ptr()->isHDivBasis();}

  /** \brief Inform caller as to whether I am an H(curl) basis. */
  bool isHCurlBasis() const {return ptr()->isHCurlBasis();}

  /** Sum up the dim() values for array of bases. */
  static unsigned int size(const Array<BasisFamily>& b) ;

  /** Extract the basis from an expression */
  static BasisFamily getBasis(const RCP<const CommonFuncDataStub>& funcData);

  /** Extract the basis from an expression */
  static RCP<BasisDOFTopologyBase> getBasisTopology(const RCP<const CommonFuncDataStub>& funcData);

  /** */
  void refEval(
    const CellType& maximalCellType,
    const CellType& cellType,
    const Array<Point>& pts,
    const MultiIndex& deriv,
    Array<Array<Array<double> > >& result,
    int verbosity) const ;
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

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d)
{
  return tuple(a,b,c,d);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e)
{
  return tuple(a,b,c,d,e);
}


/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f)
{
  return tuple(a,b,c,d,e,f);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f, 
  const BasisFamily& g)
{
  return tuple(a,b,c,d,e,f,g);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f, 
  const BasisFamily& g, const BasisFamily& h)
{
  return tuple(a,b,c,d,e,f,g,h);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f, 
  const BasisFamily& g, const BasisFamily& h, 
  const BasisFamily& i)
{
  return tuple(a,b,c,d,e,f,g,h,i);
}

/** \relates BasisFamily */
inline Array<BasisFamily> List(const BasisFamily& a, const BasisFamily& b,
  const BasisFamily& c, const BasisFamily& d, 
  const BasisFamily& e, const BasisFamily& f, 
  const BasisFamily& g, const BasisFamily& h, 
  const BasisFamily& i, const BasisFamily& j)
{
  return tuple(a,b,c,d,e,f,g,h,i,j);
}

/** */
inline Array<BasisFamily> replicate(const BasisFamily& b, int n)
{
  Array<BasisFamily> rtn(n);
  for (int i=0; i<n; i++) rtn[i] = b;
  return rtn;
}


/** */
inline Array<BasisFamily> replicate(const Array<BasisFamily>& b, int n)
{
  Array<BasisFamily> rtn(n*b.size());
  for (unsigned int i=0; i<n*b.size(); i++) rtn[i] = b[0];
  return rtn;
}

class BasisArray : public Array<BasisFamily>
{
public:
  BasisArray() : Array<BasisFamily>() {;}

  BasisArray(int n) : Array<BasisFamily>(n) {;}

  BasisArray(const Array<BasisFamily>& a) 
    : Array<BasisFamily>(a) 
    {;}

#ifndef TRILINOS_8
  template <int n>
  BasisArray(const Teuchos::Tuple<BasisFamily, n>& a) 
    : Array<BasisFamily>(a) 
    {;}
#endif
};


/** \relates BasisFamily */
Array<std::pair<int, int> > vectorDimStructure(const Array<BasisFamily>& basis);


/** \relates BasisFamily */
Array<std::pair<int, int> > vectorDimStructure(const BasisFamily& basis);



}

#endif
